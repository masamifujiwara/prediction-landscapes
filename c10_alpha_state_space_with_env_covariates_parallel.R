# =============================================================================
# Script: c10_alpha_state_space_with_env_covariates_parallel.R  
#
# Key fixes vs prior c10:
#   1) Uses fold$train_idx (not 1:t0) consistently
#   2) Handles lag warm-up months: drop training rows before first complete xreg row
#   3) Allows xreg NA only when y is NA (structural missing target months)
#   4) Forecast horizon: xregs must be present unless the target month is structural missing
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
# Config (match c7/c8/c9 philosophy)
# -----------------------------
script_stem <- "c10_alpha_state_space_with_env_covariates_parallel"

MODEL_NAME <- "ssm_ll_seas12_xreg"
FREQ <- 12
HORIZONS <- c(1, 12)

TRAINING_RULE <- "expanding"
MIN_TRAIN_MONTHS <- 120
FIXED_TRAIN_MONTHS <- 120
ORIGIN_STEP <- 1
HOLDOUT_LAST_N <- 0
MIN_TRAIN_NON_NA <- 24

TARGET_COL_PRIMARY <- "value_deseas"
TARGET_COL_FALLBACKS <- c("residual_component", "value")
TARGET_SCALE <- "deseasonalized"

PI_LEVELS <- c(80, 95)

MAXIT <- 200
QUIET <- TRUE

# Exclude Sabine Lake (Bay 1) going forward
DROP_BAYS <- c("1")

# Structural missing target months kept on grid
COVID_MISSING_MONTHS <- as.Date(c("2020-04-01", "2020-05-01"))

# -----------------------------
# Parallel plan
# -----------------------------
N_WORKERS <- min(6, max(1, parallel::detectCores(logical = FALSE) - 2))
future::plan(future::multisession, workers = N_WORKERS)

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

# Inputs
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds")
stopifnot(file.exists(IN_ALPHA))

IN_ENV <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries", "env_alpha_index.rds")
stopifnot(file.exists(IN_ENV))

# Outputs
OUT_FC_DERIVED  <- file.path(dir_data_script,  "ssm_xreg_forecasts_alpha_by_bay.rds")
OUT_MET_DERIVED <- file.path(dir_data_script,  "ssm_xreg_metrics_alpha_by_bay.csv")
OUT_XREG_MISS   <- file.path(dir_data_script,  "xreg_missingness_summary.csv")

OUT_FC_MODEL    <- file.path(dir_model_script, "ssm_xreg_forecasts_alpha_by_bay.rds")
OUT_MET_MODEL   <- file.path(dir_model_script, "ssm_xreg_metrics_alpha_by_bay.csv")
OUT_FITLOG      <- file.path(dir_diag_script,  "ssm_xreg_fit_log.csv")

# -----------------------------
# Helpers
# -----------------------------
# -----------------------------------------------------------------------------
# Utility function: complete monthly grid (GLOBAL) for one BAY×GEAR×q series
# -----------------------------------------------------------------------------
complete_monthly_grid_global <- function(df_i, start_date, end_date,
                                         covid_missing_months = NULL) {
  stopifnot("date" %in% names(df_i))
  stopifnot(all(c("BAY", "GEAR", "q") %in% names(df_i)))
  
  # Ensure single-series input
  keys <- df_i %>% dplyr::distinct(BAY, GEAR, q)
  if (nrow(keys) != 1) {
    stop("complete_monthly_grid_global() expects df_i to contain exactly ONE BAY×GEAR×q series.")
  }
  
  # Full monthly grid
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  
  df_grid <- tibble::tibble(date = as.Date(dates)) %>%
    dplyr::mutate(
      BAY  = as.character(keys$BAY),
      GEAR = as.character(keys$GEAR),
      q    = as.character(keys$q)
    )
  
  # Join observed values (expects df_i has 'value')
  if (!("value" %in% names(df_i))) {
    stop("df_i must have a column named 'value' before building the complete grid.")
  }
  
  out <- df_grid %>%
    dplyr::left_join(
      df_i %>% dplyr::select(BAY, GEAR, q, date, value),
      by = c("BAY", "GEAR", "q", "date")
    ) %>%
    dplyr::arrange(date)
  
  # Enforce structural missing months (optional)
  if (!is.null(covid_missing_months)) {
    out <- out %>%
      dplyr::mutate(value = ifelse(date %in% as.Date(covid_missing_months), NA_real_, value))
  }
  
  out
}


expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

pick_target_col <- function(df) {
  cands <- c(TARGET_COL_PRIMARY, TARGET_COL_FALLBACKS)
  hit <- cands[cands %in% names(df)]
  if (length(hit) == 0) stop("No target column found. Tried: ", paste(cands, collapse = ", "))
  hit[1]
}

standardize_alpha_names <- function(alpha) {
  if (!("GEAR" %in% names(alpha))) {
    gcol <- intersect(names(alpha), c("gear", "Gear"))
    if (length(gcol) >= 1) alpha <- dplyr::rename(alpha, GEAR = all_of(gcol[1]))
  }
  if (!("q" %in% names(alpha))) {
    qcol <- intersect(names(alpha), c("Q", "hill_q", "order_q"))
    if (length(qcol) >= 1) alpha <- dplyr::rename(alpha, q = all_of(qcol[1]))
  }
  if (!("date" %in% names(alpha))) {
    if (all(c("YEAR", "MONTH") %in% names(alpha))) {
      alpha <- alpha %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(YEAR), as.integer(MONTH))))
    } else if (all(c("year", "month") %in% names(alpha))) {
      alpha <- alpha %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(.data$year), as.integer(.data$month))))
    } else {
      stop("alpha input has no date or (YEAR, MONTH) columns.")
    }
  }
  if (!("BAY" %in% names(alpha))) {
    bcol <- intersect(names(alpha), c("BAY_ID", "bay", "Bay", "BAYID"))
    if (length(bcol) >= 1) alpha <- dplyr::rename(alpha, BAY = all_of(bcol[1]))
  }
  alpha
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
  if (!("BAY" %in% names(env))) {
    bcol <- intersect(names(env), c("BAY_ID", "bay", "Bay", "BAYID"))
    if (length(bcol) >= 1) env <- dplyr::rename(env, BAY = all_of(bcol[1]))
  }
  env
}

collapse_dups_alpha <- function(df, value_col) {
  df %>%
    group_by(BAY, GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

collapse_dups_env <- function(env, by_cols, xreg_cols) {
  env %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      across(
        all_of(xreg_cols),
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
build_folds_one_series <- function(n, horizons,
                                   training_rule = "expanding",
                                   min_train_months = 120,
                                   fixed_train_months = 120,
                                   origin_step = 1,
                                   holdout_last_n = 0) {
  
  n_eff <- if (holdout_last_n > 0) max(0, n - holdout_last_n) else n
  if (n_eff < (min_train_months + max(horizons) + 1)) {
    return(list(ok = FALSE, reason = "Too short for min_train + max_horizon", folds = list()))
  }
  
  folds <- list()
  for (h in horizons) {
    origins <- seq(from = min_train_months, to = (n_eff - h), by = origin_step)
    folds[[paste0("h", h)]] <- lapply(origins, function(t0) {
      train_idx <- if (training_rule == "expanding") 1:t0 else max(1, t0 - fixed_train_months + 1):t0
      list(horizon = h, origin_i = t0, test_i = t0 + h, train_idx = train_idx)
    })
  }
  list(ok = TRUE, reason = NA_character_, folds = folds)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

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


# -----------------------------
# Load + standardize alpha
# -----------------------------
alpha <- readr::read_rds(IN_ALPHA)
alpha <- standardize_alpha_names(alpha)
expect_cols(alpha, c("BAY", "GEAR", "q", "date"), "alpha")

alpha <- alpha %>%
  mutate(
    BAY = as.character(BAY),
    GEAR = as.character(GEAR),
    q = as.character(q),
    date = as.Date(date)
  ) %>%
  filter(!(BAY %in% DROP_BAYS))

target_col <- pick_target_col(alpha)
cat("Using alpha target column:", target_col, "\n\n")

alpha_collapsed <- collapse_dups_alpha(alpha, value_col = target_col)

# -----------------------------
# Load + standardize env
# -----------------------------
env <- readr::read_rds(IN_ENV)
env <- standardize_env_names(env)
expect_cols(env, c("GEAR", "date"), "env")

env <- env %>%
  mutate(GEAR = as.character(GEAR), date = as.Date(date))
if ("BAY" %in% names(env)) {
  env <- env %>% mutate(BAY = as.character(BAY)) %>% filter(!(BAY %in% DROP_BAYS))
}

# Choose xreg columns (same as c9/c10 prior)
key_cols <- intersect(names(env), c("BAY", "GEAR", "date"))
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

if (length(xreg_cols) == 0) stop("After filtering, no xreg columns remain. Check env_alpha_index.rds columns.")
cat("Final xreg cols:\n  ", paste(xreg_cols, collapse = ", "), "\n\n")

env_by <- if ("BAY" %in% names(env)) c("BAY", "GEAR", "date") else c("GEAR", "date")
env_collapsed <- collapse_dups_env(env, by_cols = env_by, xreg_cols = xreg_cols)

# -----------------------------
# Global complete monthly grid
# -----------------------------
GLOBAL_START <- min(alpha_collapsed$date, na.rm = TRUE)
GLOBAL_END   <- max(alpha_collapsed$date, na.rm = TRUE)
cat("Global grid:", format(GLOBAL_START), "to", format(GLOBAL_END), "\n\n")

# -----------------------------
# Build aligned series tables + folds + xreg start index
# -----------------------------
series_keys <- alpha_collapsed %>%
  distinct(BAY, GEAR, q) %>%
  arrange(BAY, GEAR, q)

series_list <- vector("list", nrow(series_keys))
folds_list  <- vector("list", nrow(series_keys))
xreg_start_i_by_series <- rep(NA_integer_, nrow(series_keys))
xreg_miss_rows <- vector("list", nrow(series_keys))

for (i in seq_len(nrow(series_keys))) {
  key <- series_keys[i, ]
  
  df_i <- alpha_collapsed %>%
    filter(BAY == key$BAY, GEAR == key$GEAR, q == key$q) %>%
    arrange(date)
  
  df_g <- complete_monthly_grid_global(df_i, GLOBAL_START, GLOBAL_END)
  
  if ("BAY" %in% names(env_collapsed)) { 
    env_i <- env_collapsed %>% filter(BAY == key$BAY, GEAR == key$GEAR) 
    df_g <- df_g %>% left_join(env_i, by = c("BAY", "GEAR", "date")) 
  } else { 
    env_i <- env_collapsed %>% filter(GEAR == key$GEAR) 
    df_g <- df_g %>% left_join(env_i, by = c("GEAR", "date")) 
  }
  
  # # ADD THIS LINE: Interpolate NAs generated by the calendar grid merge
  # df_g <- df_g %>% mutate(across(all_of(xreg_cols), ~ as.numeric(forecast::na.interp(.x))))
  
  ### xreg start = first row with complete xregs (this handles lag warm-up months)
  # Safely interpolate NAs on the calendar grid
  df_g <- df_g %>% mutate(across(all_of(xreg_cols), ~ as.numeric(zoo::na.approx(.x, na.rm = FALSE))))
  ok_x <- complete.cases(df_g[, xreg_cols, drop = FALSE]) 
  i0   <- which(ok_x)[2] 
  if (length(i0) == 0 || is.na(i0)) i0 <- NA_integer_
  xreg_start_i_by_series[i] <- as.integer(i0)
  
  folds_i <- build_folds_one_series(
    n = nrow(df_g),
    horizons = HORIZONS,
    training_rule = TRAINING_RULE,
    min_train_months = MIN_TRAIN_MONTHS,
    fixed_train_months = FIXED_TRAIN_MONTHS,
    origin_step = ORIGIN_STEP,
    holdout_last_n = HOLDOUT_LAST_N
  )
  
  miss_frac <- sapply(xreg_cols, function(xc) mean(is.na(df_g[[xc]])))
  xreg_miss_rows[[i]] <- data.frame(
    BAY = key$BAY, GEAR = key$GEAR, q = key$q,
    xreg = names(miss_frac),
    frac_na = as.numeric(miss_frac),
    first_complete_xreg = if (is.na(i0)) NA_character_ else as.character(df_g$date[i0]),
    stringsAsFactors = FALSE
  )
  
  series_list[[i]] <- df_g
  folds_list[[i]]  <- folds_i
}

xreg_miss_df <- do.call(rbind, xreg_miss_rows)
readr::write_csv(xreg_miss_df, OUT_XREG_MISS)

cat("Aligned alpha groups:", nrow(series_keys), "\n\n")

# -----------------------------
# Job list
# -----------------------------
jobs <- list()
for (i in seq_len(nrow(series_keys))) {
  fi <- folds_list[[i]]
  if (!isTRUE(fi$ok)) next
  for (hname in names(fi$folds)) {
    for (fold in fi$folds[[hname]]) {
      jobs[[length(jobs) + 1]] <- list(series_i = i, fold = fold)
    }
  }
}
cat("Jobs to run:", length(jobs), "\n\n")

# -----------------------------
# Fit + forecast one job
# -----------------------------
run_job <- function(job,
                    series_keys, series_list, xreg_cols, xreg_start_i_by_series,
                    covid_missing_months, min_train_non_na) {
  
  i <- job$series_i
  fold <- job$fold
  
  key <- series_keys[i, ]
  df_g <- series_list[[i]]
  
  y <- safe_num(df_g$value)
  n <- length(y)
  
  h <- fold$horizon
  t0 <- fold$origin_i
  t1 <- fold$test_i
  train_idx <- fold$train_idx
  
  origin_date <- df_g$date[t0]
  test_date   <- df_g$date[t1]
  
  if (t1 > n) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = "test_i_oob",
      stringsAsFactors = FALSE
    ))
  }
  
  # Lag warm-up handling
  xreg_start_i <- xreg_start_i_by_series[i]
  if (is.na(xreg_start_i)) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = "no_complete_xreg_window",
      stringsAsFactors = FALSE
    ))
  }
  
  tr_eff <- train_idx[train_idx >= xreg_start_i]
  if (length(tr_eff) == 0) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = "train_window_before_xreg_start",
      stringsAsFactors = FALSE
    ))
  }
  
  # Enforce min observed y in effective training window
  if (sum(is.finite(y[tr_eff])) < min_train_non_na) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = "too_few_non_na_in_train",
      stringsAsFactors = FALSE
    ))
  }
  
  # Build KFAS frames
  y_train <- y[tr_eff]
  df_train <- data.frame(y = y_train, stringsAsFactors = FALSE)
  df_new   <- data.frame(y = rep(NA_real_, h), stringsAsFactors = FALSE)
  
  if (length(xreg_cols) > 0) {
    X_train <- df_g[tr_eff, xreg_cols, drop = FALSE]
    X_new   <- df_g[(t0 + 1):(t0 + h), xreg_cols, drop = FALSE]
    
    df_train <- cbind(df_train, X_train)
    df_new   <- cbind(df_new, X_new)
    
    # STRICT where y is observed: xregs must be present
    bad_train <- anyNA(X_train[is.finite(y_train), , drop = FALSE])
    if (bad_train) {
      return(data.frame(
        BAY = key$BAY, GEAR = key$GEAR, q = key$q,
        horizon = h, origin_i = t0, test_i = t1,
        origin_date = origin_date, test_date = test_date,
        model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
        lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
        ok = FALSE, reason = "missing_xreg_in_train",
        stringsAsFactors = FALSE
      ))
    }
    
    # Allow xreg NA only where y is NA (structural missing); fill to satisfy KFAS
    if (anyNA(X_train)) X_train[is.na(X_train)] <- 0
    
    y_obs <- y[t1]
    if (anyNA(X_new)) {
      # if the target month is observed, this is a hard failure
      if (is.finite(y_obs) && !(test_date %in% covid_missing_months)) {
        return(data.frame(
          BAY = key$BAY, GEAR = key$GEAR, q = key$q,
          horizon = h, origin_i = t0, test_i = t1,
          origin_date = origin_date, test_date = test_date,
          model = MODEL_NAME, y_obs = y_obs, y_hat = NA_real_,
          lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
          ok = FALSE, reason = "missing_xreg_in_forecast_horizon",
          stringsAsFactors = FALSE
        ))
      } else {
        # structural missing target month: allow and fill
        X_new[is.na(X_new)] <- 0
      }
    }
    
    for (xc in xreg_cols) {
      df_train[[xc]] <- safe_num(X_train[[xc]])
      df_new[[xc]]   <- safe_num(X_new[[xc]])
    }
  }
  
  # Fit
  fit <- fit_kfas(df_train, xreg_cols)
  if (inherits(fit, "error")) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = paste0("fit_error: ", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }
  
  pr80 <- predict_kfas(fit, df_new, xreg_cols, level = 0.80)
  pr95 <- predict_kfas(fit, df_new, xreg_cols, level = 0.95)
  
  if (inherits(pr80, "error") || inherits(pr95, "error")) {
    msg <- if (inherits(pr80, "error")) conditionMessage(pr80) else conditionMessage(pr95)
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = y[t1], y_hat = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      ok = FALSE, reason = paste0("predict_error: ", msg),
      stringsAsFactors = FALSE
    ))
  }
  
  yhat <- as.numeric(pr80[h, "fit"])
  l80  <- as.numeric(pr80[h, "lwr"])
  u80  <- as.numeric(pr80[h, "upr"])
  l95  <- as.numeric(pr95[h, "lwr"])
  u95  <- as.numeric(pr95[h, "upr"])
  
  yobs <- y[t1]
  
  if (!is.finite(yobs)) {
    return(data.frame(
      BAY = key$BAY, GEAR = key$GEAR, q = key$q,
      horizon = h, origin_i = t0, test_i = t1,
      origin_date = origin_date, test_date = test_date,
      model = MODEL_NAME, y_obs = yobs, y_hat = yhat,
      lower80 = l80, upper80 = u80, lower95 = l95, upper95 = u95,
      ok = FALSE, reason = "missing_target_month",
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    BAY = key$BAY, GEAR = key$GEAR, q = key$q,
    horizon = h, origin_i = t0, test_i = t1,
    origin_date = origin_date, test_date = test_date,
    model = MODEL_NAME, y_obs = yobs, y_hat = yhat,
    lower80 = l80, upper80 = u80, lower95 = l95, upper95 = u95,
    ok = TRUE, reason = NA_character_,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Run jobs
# -----------------------------
pred_list <- future.apply::future_lapply(
  jobs,
  run_job,
  series_keys = series_keys,
  series_list = series_list,
  xreg_cols = xreg_cols,
  xreg_start_i_by_series = xreg_start_i_by_series,
  covid_missing_months = COVID_MISSING_MONTHS,
  min_train_non_na = MIN_TRAIN_NON_NA,
  future.seed = TRUE
)

pred <- dplyr::bind_rows(pred_list) %>%
  mutate(
    BAY = as.character(BAY),
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
  arrange(BAY, GEAR, q, horizon, origin_i)

# -----------------------------
# Metrics (lightweight; full metrics computed in d1)
# -----------------------------
met <- pred %>%
  filter(ok) %>%
  group_by(GEAR, q, horizon, model) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean((y_obs - y_hat)^2, na.rm = TRUE)),
    MAE  = mean(abs(y_obs - y_hat), na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(met, OUT_MET_DERIVED)
readr::write_csv(met, OUT_MET_MODEL)

# Fit log summary
fitlog <- pred %>%
  mutate(reason2 = ifelse(ok, "ok", ifelse(is.na(reason), "fail_unknown", reason))) %>%
  count(GEAR, q, horizon, model, reason2, name = "n") %>%
  arrange(GEAR, q, horizon, model, desc(n))
readr::write_csv(fitlog, OUT_FITLOG)

# Save output (match d0/d1)
out_obj <- list(
  predictions = pred,
  eval_target_col = target_col,
  target_scale = TARGET_SCALE,
  script_stem = script_stem,
  created = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

readr::write_rds(out_obj, OUT_FC_DERIVED)
readr::write_rds(out_obj, OUT_FC_MODEL)

cat("\nSaved forecasts:\n  ", OUT_FC_DERIVED, "\n  ", OUT_FC_MODEL, "\n")
cat("Saved metrics:\n  ", OUT_MET_DERIVED, "\n  ", OUT_MET_MODEL, "\n")
cat("Saved xreg missingness:\n  ", OUT_XREG_MISS, "\n")
cat("Saved fit log:\n  ", OUT_FITLOG, "\n")
cat("\nDone:", format(Sys.time()), "\n")

sink()