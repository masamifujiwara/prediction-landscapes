# =============================================================================
# Script: c11_gamma_state_space_with_env_covariates_parallel.R
#
# Purpose:
#   KFAS state-space forecasts for gamma (q=0,1,2) by GEAR×q using rolling-origin
#   folds from c3 (backtest_folds.rds) on the COMPLETE monthly fold grid,
#   ADDING environmental covariates as regression terms (same filtering rules as c6/c9/c10).
#
# Key fix vs earlier version:
#   - Lag warm-up months (1982) make lagged xregs NA by construction.
#   - We trim the effective training window to start at the first row where ALL xregs are complete,
#     but we still enforce that xregs are present wherever y is observed within that window.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# -----------------------------
# Required bootstrap snippet
# -----------------------------
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)

# -----------------------------
# PATHS
# -----------------------------
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Dependencies
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "KFAS", "future", "future.apply", "parallel")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(KFAS)
  library(future)
  library(future.apply)
  library(forecast)
})

# -----------------------------
# Config
# -----------------------------
script_stem <- "c11_gamma_state_space_with_env_covariates_parallel"

MODEL_NAME <- "ssm_ll_seas12_xreg"
FREQ <- 12
HORIZONS <- c(1, 12)

EVAL_TARGET_COL <- "value_deseas"
TARGET_SCALE <- "deseasonalized"

MIN_TRAIN_NON_NA <- 24
MAXIT <- 200

COVID_MISSING_MONTHS <- as.Date(c("2020-04-01", "2020-05-01"))

# -----------------------------
# Parallel plan
# -----------------------------
N_WORKERS <- min(6, max(1, parallel::detectCores(logical = FALSE) - 2))
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_WORKERS)
} else {
  future::plan(future::multicore, workers = N_WORKERS)
}
options(future.globals.maxSize = 4 * 1024^3)  # 4GB

# -----------------------------
# Directories
# -----------------------------
dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script  <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dir_logs_script, paste0("runlog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
cat("Script:", script_stem, "\nStarted:", format(Sys.time()), "\nWorkers:", N_WORKERS, "\n\n")

# -----------------------------
# Inputs (authoritative)
# -----------------------------
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "gamma_q012_deseas.rds")
IN_FOLDS <- file.path(PATHS$analysis_dir, "Data", "c3_forecast_design_backtesting", "backtest_folds.rds")
IN_ENVG  <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries", "env_gamma_index.rds")

stopifnot(file.exists(IN_GAMMA), file.exists(IN_FOLDS), file.exists(IN_ENVG))

cat("Inputs:\n")
cat("  gamma:", IN_GAMMA, "\n")
cat("  folds:", IN_FOLDS, "\n")
cat("  envg :", IN_ENVG, "\n\n")

# -----------------------------
# Outputs
# -----------------------------
OUT_FC_DERIVED <- file.path(dir_data_script,  "ssm_xreg_forecasts_gamma_by_gear_q.rds")
OUT_FC_MODEL   <- file.path(dir_model_script, "ssm_xreg_forecasts_gamma_by_gear_q.rds")
OUT_FITLOG     <- file.path(dir_diag_script,  "ssm_xreg_fit_log.csv")
OUT_XREG_MISS  <- file.path(dir_data_script,  "xreg_missingness_summary.csv")

# -----------------------------
# Helpers
# -----------------------------
expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

standardize_gamma_names <- function(gamma) {
  if (!("GEAR" %in% names(gamma))) {
    gcol <- intersect(names(gamma), c("gear", "Gear"))
    if (length(gcol) >= 1) gamma <- dplyr::rename(gamma, GEAR = all_of(gcol[1]))
  }
  if (!("q" %in% names(gamma))) {
    qcol <- intersect(names(gamma), c("Q", "hill_q", "order_q"))
    if (length(qcol) >= 1) gamma <- dplyr::rename(gamma, q = all_of(qcol[1]))
  }
  if (!("date" %in% names(gamma))) {
    if (all(c("YEAR", "MONTH") %in% names(gamma))) {
      gamma <- gamma %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(YEAR), as.integer(MONTH))))
    } else if (all(c("year", "month") %in% names(gamma))) {
      gamma <- gamma %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(.data$year), as.integer(.data$month))))
    } else {
      stop("gamma input has no date or (YEAR, MONTH) columns.")
    }
  }
  gamma
}

standardize_env_names <- function(env) {
  if (!("GEAR" %in% names(env))) {
    gcol <- intersect(names(env), c("gear", "Gear"))
    if (length(gcol) >= 1) env <- dplyr::rename(env, GEAR = all_of(gcol[1]))
  }
  if (!("date" %in% names(env))) {
    dcol <- intersect(names(env), c("DATE", "Date", "month_date"))
    if (length(dcol) >= 1) env <- dplyr::rename(env, date = all_of(dcol[1]))
  }
  env
}

collapse_dups_gamma <- function(df, value_col) {
  df %>%
    group_by(GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

collapse_dups_env <- function(env, by_cols, xreg_cols) {
  env %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by_cols))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(xreg_cols),
        ~ {
          if (all(is.na(.x))) {
            NA_real_
          } else {
            mean(.x, na.rm = TRUE)
          }
        }
      ),
      .groups = "drop"
    )
}

build_kfas_formula <- function(xreg_cols) { 
  if (length(xreg_cols) == 0) { 
    as.formula("y ~ -1 + SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type = 'dummy', Q = NA)") 
  } else { 
    rhs <- paste( 
      "-1 + SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type = 'dummy', Q = NA) +", 
      paste(xreg_cols, collapse = " + ") 
    ) 
    as.formula(paste("y ~", rhs)) 
  } 
}

fit_kfas <- function(df_train, xreg_cols) { 
  form <- build_kfas_formula(xreg_cols)
  
  # Build the initial model structure
  model0 <- SSModel( 
    formula = form, data = df_train, H = NA 
  )
  
  # Robust update function that handles the exact number of NAs safely
  upd_fn <- function(pars, model) {
    i <- 1
    if (anyNA(model$H)) {
      model$H[is.na(model$H)] <- exp(pars[i])
      i <- i + 1
    }
    if (anyNA(model$Q)) {
      nQ <- sum(is.na(model$Q))
      model$Q[is.na(model$Q)] <- exp(pars[i:(i + nQ - 1)])
    }
    model
  }
  
  # Dynamically count how many variance parameters (NAs) need to be estimated
  n_par <- sum(is.na(model0$H)) + sum(is.na(model0$Q))
  
  v0 <- var(df_train$y, na.rm = TRUE)
  if (!is.finite(v0) || v0 <= 0) v0 <- 1
  inits <- rep(log(v0), n_par)
  
  tryCatch( 
    fitSSM(model0, inits = inits, updatefn = upd_fn, method = "BFGS", control = list(maxit = MAXIT)), 
    error = function(e) e 
  ) 
}

# THE RESTORED FUNCTION
predict_kfas <- function(fit, df_new, xreg_cols, level = 0.80) { 
  if (inherits(fit, "error")) return(fit) 
  tryCatch({
    form <- build_kfas_formula(xreg_cols)
    
    # Build the future model structure, but explicitly insert the TRAINED observation variance
    mod_new <- SSModel(formula = form, data = df_new, H = fit$model$H)
    
    # Copy the TRAINED state variances over so the model has no NA placeholders
    mod_new$Q <- fit$model$Q
    
    # Predict using the fully specified SSModel
    predict(fit$model, newdata = mod_new, interval = "prediction", level = level)
  }, error = function(e) e) 
}

mk_row <- function(GEAR, q, h, t0, t1, origin_date, test_date, yobs,
                   yhat = NA_real_, l80 = NA_real_, u80 = NA_real_, l95 = NA_real_, u95 = NA_real_,
                   ok = FALSE, reason = "fail") {
  data.frame(
    GEAR = as.character(GEAR),
    q = as.character(q),
    horizon = as.integer(h),
    origin_i = as.integer(t0),
    test_i = as.integer(t1),
    origin_date = as.Date(origin_date),
    test_date = as.Date(test_date),
    model = as.character(MODEL_NAME),
    y_obs = as.numeric(yobs),
    y_hat = as.numeric(yhat),
    lower80 = as.numeric(l80),
    upper80 = as.numeric(u80),
    lower95 = as.numeric(l95),
    upper95 = as.numeric(u95),
    ok = as.logical(ok),
    reason = as.character(reason),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Load + standardize gamma
# -----------------------------
gamma_raw <- readr::read_rds(IN_GAMMA)
gamma_raw <- standardize_gamma_names(gamma_raw)
expect_cols(gamma_raw, c("GEAR", "q", "date"), "gamma")

gamma_raw <- gamma_raw %>%
  mutate(GEAR = as.character(GEAR), q = as.character(q), date = as.Date(date))

if (!(EVAL_TARGET_COL %in% names(gamma_raw))) {
  stop("Requested eval_target_col not found in gamma input: ", EVAL_TARGET_COL,
       "\nAvailable cols include: ", paste(names(gamma_raw), collapse = ", "))
}
cat("Using gamma eval_target_col:", EVAL_TARGET_COL, "\n\n")

gamma_obs <- gamma_raw %>%
  transmute(GEAR, q, date, value_use = safe_num(.data[[EVAL_TARGET_COL]])) %>%
  { collapse_dups_gamma(., value_col = "value_use") } %>%
  mutate(value = ifelse(date %in% COVID_MISSING_MONTHS, NA_real_, value)) %>%
  arrange(GEAR, q, date)

# -----------------------------
# Load + standardize env
# -----------------------------
env <- readr::read_rds(IN_ENVG)
env <- standardize_env_names(env)
expect_cols(env, c("GEAR", "date"), "env_gamma")

env <- env %>%
  mutate(GEAR = as.character(GEAR), date = as.Date(date))

key_cols <- intersect(names(env), c("GEAR", "date"))
xreg_cols <- setdiff(names(env), key_cols)
xreg_cols <- xreg_cols[vapply(env[xreg_cols], is.numeric, logical(1))]

xreg_cols <- setdiff(xreg_cols, c("YEAR", "MONTH", "n_samples"))
xreg_cols <- xreg_cols[!grepl("^TURB__", xreg_cols)]
xreg_cols <- xreg_cols[!grepl("^DO__",   xreg_cols)]
xreg_cols <- xreg_cols[!grepl("__sd$", xreg_cols)]

is_lag <- grepl("__lag[0-9]+$", xreg_cols)
lag_ok <- !is_lag | grepl("__lag(1|12)$", xreg_cols)
xreg_cols <- xreg_cols[lag_ok]

keep_pat <- "^(TEMP__mean|SAL__mean)(__(lag(1|12)))?$"
xreg_cols <- xreg_cols[grepl(keep_pat, xreg_cols)]

if (length(xreg_cols) == 0) stop("After filtering, no xreg columns remain. Check env_gamma_index.rds columns.")
cat("Final xreg cols:\n  ", paste(xreg_cols, collapse = ", "), "\n\n")

env_collapsed <- collapse_dups_env(env, by_cols = c("GEAR", "date"), xreg_cols = xreg_cols)

# -----------------------------
# Load folds + gamma groups
# -----------------------------
folds <- readRDS(IN_FOLDS)
stopifnot(is.list(folds), !is.null(folds$groups))

gamma_groups <- folds$groups[vapply(folds$groups, function(g) identical(g$series, "gamma"), logical(1))]
if (length(gamma_groups) == 0) stop("No gamma groups found in backtest_folds.rds (series == 'gamma').")

# -----------------------------
# Build aligned tables (fold grid = master)
# -----------------------------
get_gamma_group_df_aligned <- function(g, gamma_obs_g, env_collapsed_g) {
  tibble(date = as.Date(g$dates)) %>%
    left_join(gamma_obs_g, by = "date") %>%
    left_join(env_collapsed_g, by = "date") %>%
    mutate(GEAR = as.character(g$GEAR), q = as.character(g$q)) %>%
    arrange(date)
}

aligned_by_gid <- lapply(names(gamma_groups), function(gid) {
  g <- gamma_groups[[gid]]
  if (!isTRUE(g$ok)) return(NULL)
  
  df_obs_g <- gamma_obs %>% filter(GEAR == as.character(g$GEAR), q == as.character(g$q)) %>% select(date, value)
  df_env_g <- env_collapsed %>% filter(GEAR == as.character(g$GEAR)) %>% select(date, all_of(xreg_cols))
  
  df_g <- get_gamma_group_df_aligned(g, df_obs_g, df_env_g)
  
  ### ADD THIS LINE: Safely interpolate NAs on the calendar grid
  df_g <- df_g %>% mutate(across(all_of(xreg_cols), ~ as.numeric(zoo::na.approx(.x, na.rm = FALSE))))
  
  if (nrow(df_g) != length(g$dates)) {
    stop("Aligned grid length mismatch for ", gid,
         " (aligned rows=", nrow(df_g), ", folds dates=", length(g$dates), ").")
  }
  df_g
})
names(aligned_by_gid) <- names(gamma_groups)

# first index where ALL xregs are complete (fixes lag warm-up issue)
xreg_start_i_by_gid <- sapply(names(aligned_by_gid), function(gid) {
  dfg <- aligned_by_gid[[gid]]
  if (is.null(dfg)) return(NA_integer_)
  ok <- stats::complete.cases(dfg[, xreg_cols, drop = FALSE])
  i0 <- which(ok)[1]
  if (length(i0) == 0 || is.na(i0)) NA_integer_ else as.integer(i0)
})

# -----------------------------
# Job list
# -----------------------------
jobs <- list()

for (gid in names(gamma_groups)) {
  g <- gamma_groups[[gid]]
  if (!isTRUE(g$ok)) next
  df_g <- aligned_by_gid[[gid]]
  if (is.null(df_g)) next
  n_g <- nrow(df_g)
  
  for (hname in names(g$folds)) {
    for (f in g$folds[[hname]]) {
      h <- f$horizon
      if (!(h %in% HORIZONS)) next
      
      if (max(unlist(f$train_idx), unlist(f$test_idx)) > n_g) next
      
      jobs[[length(jobs) + 1]] <- list(
        gid = gid,
        GEAR = as.character(g$GEAR),
        q = as.character(g$q),
        horizon = as.integer(h),
        origin_i = as.integer(f$origin_i),
        test_i = as.integer(f$test_i),
        train_idx = as.integer(f$train_idx)
      )
    }
  }
}

cat("Gamma groups:", length(gamma_groups), "\n")
cat("Jobs to run:", length(jobs), "\n\n")

# -----------------------------
# Run one job (CLOSURE; uses globals)
# -----------------------------
run_job <- function(job) {
  gid <- job$gid
  df_g <- aligned_by_gid[[gid]]
  
  y <- df_g$value
  n <- length(y)
  
  h  <- job$horizon
  t0 <- job$origin_i
  t1 <- job$test_i
  tr_idx <- job$train_idx
  
  origin_date <- df_g$date[t0]
  test_date   <- df_g$date[t1]
  
  if (t1 > n) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = NA_real_,
                  ok = FALSE, reason = "test_i_oob"))
  }
  
  # Effective training starts when xregs are fully available (lag warm-up)
  xreg_start_i <- xreg_start_i_by_gid[[gid]]
  if (is.na(xreg_start_i)) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = y[t1],
                  ok = FALSE, reason = "no_complete_xreg_window"))
  }
  
  tr_eff <- tr_idx[tr_idx >= xreg_start_i]
  if (length(tr_eff) == 0) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = y[t1],
                  ok = FALSE, reason = "train_window_before_xreg_start"))
  }
  
  if (sum(is.finite(y[tr_eff])) < MIN_TRAIN_NON_NA) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = y[t1],
                  ok = FALSE, reason = "too_few_non_na_in_train"))
  }
  
  df_train <- data.frame(y = y[tr_eff], stringsAsFactors = FALSE)
  df_new   <- data.frame(y = rep(NA_real_, h), stringsAsFactors = FALSE)
  
  if (length(xreg_cols) > 0) { 
    X_train <- df_g[tr_eff, xreg_cols, drop = FALSE] 
    X_new   <- df_g[(t0 + 1):(t0 + h), xreg_cols, drop = FALSE] 
    
    # ADD THESE TWO LINES: Merge the covariates so KFAS can see them!
    df_train <- cbind(df_train, X_train)
    df_new   <- cbind(df_new, X_new)
  }
  
  fit <- fit_kfas(df_train, xreg_cols)
  
  if (inherits(fit, "error")) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = y[t1],
                  ok = FALSE, reason = paste0("fit_error: ", conditionMessage(fit))))
  }
  
  pr80 <- predict_kfas(fit, df_new, xreg_cols, level = 0.80)
  pr95 <- predict_kfas(fit, df_new, xreg_cols, level = 0.95)
  
  if (inherits(pr80, "error") || inherits(pr95, "error")) {
    msg <- if (inherits(pr80, "error")) conditionMessage(pr80) else conditionMessage(pr95)
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs = y[t1],
                  ok = FALSE, reason = paste0("predict_error: ", msg)))
  }
  
  yhat <- as.numeric(pr80[h, "fit"])
  l80  <- as.numeric(pr80[h, "lwr"])
  u80  <- as.numeric(pr80[h, "upr"])
  l95  <- as.numeric(pr95[h, "lwr"])
  u95  <- as.numeric(pr95[h, "upr"])
  
  yobs <- y[t1]
  
  if (!is.finite(yobs)) {
    return(mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs,
                  yhat, l80, u80, l95, u95, ok = FALSE, reason = "missing_target_month"))
  }
  
  mk_row(job$GEAR, job$q, h, t0, t1, origin_date, test_date, yobs,
         yhat, l80, u80, l95, u95, ok = TRUE, reason = NA_character_)
}

# -----------------------------
# Run jobs (parallel)
# -----------------------------
pred_list <- future.apply::future_lapply(jobs, run_job, future.seed = TRUE)

pred <- dplyr::bind_rows(pred_list) %>%
  mutate(
    GEAR = as.character(GEAR),
    q = as.character(q),
    horizon = as.integer(horizon),
    origin_i = as.integer(origin_i),
    test_i = as.integer(test_i),
    origin_date = as.Date(origin_date),
    test_date = as.Date(test_date),
    model = as.character(model),
    ok = as.logical(ok),
    reason = as.character(reason)
  ) %>%
  arrange(GEAR, q, horizon, origin_i)

# -----------------------------
# Xreg missingness summary
# -----------------------------
xreg_miss <- bind_rows(lapply(names(aligned_by_gid), function(gid) {
  df_g <- aligned_by_gid[[gid]]
  if (is.null(df_g)) return(NULL)
  tibble(
    gid = gid,
    GEAR = unique(df_g$GEAR),
    q = unique(df_g$q),
    n_months = nrow(df_g),
    first_complete_xreg = as.character(df_g$date[xreg_start_i_by_gid[[gid]]]),
    n_missing_any_xreg = sum(!stats::complete.cases(df_g[, xreg_cols, drop = FALSE])),
    pct_missing_any_xreg = 100 * mean(!stats::complete.cases(df_g[, xreg_cols, drop = FALSE])),
    n_missing_target = sum(!is.finite(df_g$value))
  )
}))
readr::write_csv(xreg_miss, OUT_XREG_MISS)

# -----------------------------
# Save (d0/d1 expectations)
# -----------------------------
out_obj <- list(
  predictions = pred,
  eval_target_col = EVAL_TARGET_COL,
  target_scale = TARGET_SCALE,
  script_stem = script_stem,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
)

readr::write_rds(out_obj, OUT_FC_DERIVED)
readr::write_rds(out_obj, OUT_FC_MODEL)

fitlog <- pred %>%
  mutate(reason2 = ifelse(ok, "ok", ifelse(is.na(reason), "fail_unknown", reason))) %>%
  count(GEAR, q, horizon, model, reason2, name = "n") %>%
  arrange(GEAR, q, horizon, model, desc(n))
readr::write_csv(fitlog, OUT_FITLOG)

cat("\nSaved:\n")
cat("  Forecasts (derived):", OUT_FC_DERIVED, "\n")
cat("  Forecasts (model)  :", OUT_FC_MODEL, "\n")
cat("  Fit log            :", OUT_FITLOG, "\n")
cat("  Xreg missingness   :", OUT_XREG_MISS, "\n\n")
cat("Done:", format(Sys.time()), "\n")
sink()