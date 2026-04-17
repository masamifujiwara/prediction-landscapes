# =============================================================================
# Script: c9_alpha_arimax_from_env_index_parallel.R  (REVISED to match c7)
#
# Purpose:
#   ARIMAX/SARIMAX alpha forecasts with env regressors for BAY×GEAR×q
#   on a GLOBAL complete monthly grid for fold-index alignment with c7/c8.
#
# Output conventions (match revised c7):
#   - Save RDS as list: predictions + eval_target_col + target_scale + script_stem + created
#   - predictions columns include: BAY,GEAR,q,horizon,origin_i,test_i,origin_date,test_date,
#     model,y_obs,y_hat,lower80,upper80,lower95,upper95,ok,reason
# =============================================================================

rm(list = ls())

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
need_pkgs <- c("dplyr", "tidyr", "readr", "forecast", "future", "future.apply", "parallel", "zoo")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forecast)
  library(future)
  library(future.apply)
})

# -----------------------------
# Config
# -----------------------------
script_stem <- "c9_alpha_arimax_from_env_index_parallel"

MODEL_NAME <- "arimax_seasonal_autoarima"
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

XREG_SCALE <- FALSE

# Drop Sabine Lake (Bay 1)
DROP_BAYS <- c("1")

# -----------------------------
# Parallel plan
# -----------------------------
# Use processes (multisession), cap at 7 workers (you have 8 physical cores)
N_WORKERS <- min(7, max(1, parallel::detectCores(logical = FALSE) - 1))
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
OUT_FC_DERIVED  <- file.path(dir_data_script,  "arimax_forecasts_alpha_by_bay.rds")
OUT_MET_DERIVED <- file.path(dir_data_script,  "arimax_metrics_alpha_by_bay.csv")
OUT_XREG_MISS   <- file.path(dir_data_script,  "xreg_missingness_summary.csv")

OUT_FC_MODEL    <- file.path(dir_model_script, "arimax_forecasts_alpha_by_bay.rds")
OUT_MET_MODEL   <- file.path(dir_model_script, "arimax_metrics_alpha_by_bay.csv")

OUT_FITLOG      <- file.path(dir_diag_script,  "arimax_fit_log.csv")

# -----------------------------
# Helpers
# -----------------------------
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
        ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))
      ),
      .groups = "drop"
    )
}

complete_monthly_grid_global <- function(df_obs, global_start, global_end) {
  grid <- tibble(date = seq.Date(global_start, global_end, by = "month"))
  grid %>% left_join(df_obs, by = "date") %>% arrange(date)
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
    folds_h <- lapply(origins, function(t0) {
      train_idx <- if (training_rule == "expanding") 1:t0 else max(1, t0 - fixed_train_months + 1):t0
      test_path_idx <- (t0 + 1):(t0 + h)
      list(horizon = h, origin_i = t0, test_i = t0 + h, train_idx = train_idx, test_path_idx = test_path_idx)
    })
    folds[[paste0("h", h)]] <- folds_h
  }
  list(ok = TRUE, reason = NA_character_, folds = folds)
}

rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

coverage_rate <- function(y, lo, hi) {
  ok <- is.finite(y) & is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(NA_real_)
  mean(y[ok] >= lo[ok] & y[ok] <= hi[ok])
}

as_xreg_mat <- function(df, cols) {
  m <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(m) <- "double"
  m
}

scale_train_apply <- function(x_train, x_future) {
  mu <- apply(x_train, 2, mean, na.rm = TRUE)
  sdv <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  list(
    x_train  = sweep(sweep(x_train,  2, mu, "-"), 2, sdv, "/"),
    x_future = sweep(sweep(x_future, 2, mu, "-"), 2, sdv, "/")
  )
}

# -----------------------------
# Load alpha + env
# -----------------------------
alpha_obj <- readRDS(IN_ALPHA)
alpha <- alpha_obj
if (is.list(alpha_obj) && !is.data.frame(alpha_obj)) {
  idx <- which(vapply(alpha_obj, is.data.frame, logical(1)))
  if (length(idx) >= 1) alpha <- alpha_obj[[idx[1]]]
}
stopifnot(is.data.frame(alpha))

alpha <- standardize_alpha_names(alpha)
expect_cols(alpha, c("BAY", "GEAR", "q", "date"), "alpha")

alpha <- alpha %>%
  mutate(BAY = as.character(BAY)) %>%
  filter(!(BAY %in% DROP_BAYS))

target_col <- pick_target_col(alpha)
EVAL_TARGET_COL <- target_col
cat("Using target column:", target_col, "\n")

alpha2 <- alpha %>%
  mutate(BAY = as.character(BAY),
         GEAR = as.character(GEAR),
         q = as.character(q),
         date = as.Date(date))

alpha_collapsed <- collapse_dups_alpha(alpha2, value_col = target_col)
GLOBAL_START <- min(alpha_collapsed$date, na.rm = TRUE)
GLOBAL_END   <- max(alpha_collapsed$date, na.rm = TRUE)

env_obj <- readRDS(IN_ENV)
env <- env_obj
if (is.list(env_obj) && !is.data.frame(env_obj)) {
  idx <- which(vapply(env_obj, is.data.frame, logical(1)))
  if (length(idx) >= 1) env <- env_obj[[idx[1]]]
}
stopifnot(is.data.frame(env))

env <- standardize_env_names(env)
expect_cols(env, c("GEAR", "date"), "env")

env <- env %>%
  mutate(GEAR = as.character(GEAR), date = as.Date(date))
if ("BAY" %in% names(env)) {
  env <- env %>%
    mutate(BAY = as.character(BAY)) %>%
    filter(!(BAY %in% DROP_BAYS))
}

# -----------------------------
# Choose xreg columns (filtered set)
# -----------------------------
key_cols <- intersect(names(env), c("BAY", "GEAR", "date"))
xreg_cols <- setdiff(names(env), key_cols)
xreg_cols <- xreg_cols[vapply(env[xreg_cols], is.numeric, logical(1))]

# Remove YEAR/MONTH/n_samples
xreg_cols <- setdiff(xreg_cols, c("YEAR", "MONTH", "n_samples"))

# Remove Turbidity + Dissolved Oxygen
xreg_cols <- xreg_cols[!grepl("^TURB__", xreg_cols)]
xreg_cols <- xreg_cols[!grepl("^DO__",   xreg_cols)]

# Remove all SD terms
xreg_cols <- xreg_cols[!grepl("__sd$", xreg_cols)]

# Keep only lags 1 and 12 (and keep non-lag variables)
is_lag <- grepl("__lag[0-9]+$", xreg_cols)
lag_ok <- !is_lag | grepl("__lag(1|12)$", xreg_cols)
xreg_cols <- xreg_cols[lag_ok]

# (Optional safety) Keep only TEMP/SAL (+ their allowed lags)
# This prevents accidentally keeping stray numeric columns.
keep_pat <- "^(TEMP__mean|SAL__mean)(__(lag(1|12)))?$"
xreg_cols <- xreg_cols[grepl(keep_pat, xreg_cols)]

if (length(xreg_cols) == 0) stop("After filtering, no xreg columns remain. Check env_alpha_index.rds columns.")
cat("Final xreg cols:\n  ", paste(xreg_cols, collapse = ", "), "\n\n")

# Collapse env to the right key
env_by <- if ("BAY" %in% names(env)) c("BAY", "GEAR", "date") else c("GEAR", "date")
env_collapsed <- collapse_dups_env(env, by_cols = env_by, xreg_cols = xreg_cols)

# -----------------------------
# Series + aligned tables
# -----------------------------
series_keys <- alpha_collapsed %>%
  distinct(BAY, GEAR, q) %>%
  arrange(BAY, GEAR, q)

series_list <- vector("list", nrow(series_keys))
folds_list  <- vector("list", nrow(series_keys))
xreg_miss_rows <- list()

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
  
  # --- Fill xreg gaps (so h=12 forecast path doesn’t die on missing xregs) ---
  df_g <- df_g %>%
    arrange(date) %>%
    mutate(across(all_of(xreg_cols),
                  ~ as.numeric(zoo::na.approx(.x, na.rm = FALSE))))
  
  folds_i <- build_folds_one_series(
    n = nrow(df_g),
    horizons = HORIZONS,
    training_rule = TRAINING_RULE,
    min_train_months = MIN_TRAIN_MONTHS,
    fixed_train_months = FIXED_TRAIN_MONTHS,
    origin_step = ORIGIN_STEP,
    holdout_last_n = HOLDOUT_LAST_N
  )
  
  miss_rate <- sapply(xreg_cols, function(cc) mean(is.na(df_g[[cc]])))
  xreg_miss_rows[[length(xreg_miss_rows) + 1]] <- tibble(
    BAY = key$BAY, GEAR = key$GEAR, q = key$q,
    xreg = names(miss_rate),
    missing_rate = as.numeric(miss_rate)
  )
  
  series_list[[i]] <- df_g
  folds_list[[i]] <- folds_i
}

xreg_miss_tbl <- bind_rows(xreg_miss_rows) %>% arrange(desc(missing_rate), GEAR, q, BAY, xreg)
readr::write_csv(xreg_miss_tbl, OUT_XREG_MISS)

# -----------------------------
# Build job table
# -----------------------------
jobs <- list()
for (i in seq_len(nrow(series_keys))) {
  key <- series_keys[i, ]
  folds_i <- folds_list[[i]]
  if (!isTRUE(folds_i$ok)) next
  
  for (h_name in names(folds_i$folds)) {
    for (f in folds_i$folds[[h_name]]) {
      jobs[[length(jobs) + 1]] <- list(
        series_i = i,
        BAY = key$BAY,
        GEAR = key$GEAR,
        q = key$q,
        horizon = f$horizon,
        origin_i = f$origin_i,
        test_i = f$test_i,
        train_idx = f$train_idx,
        test_path_idx = f$test_path_idx
      )
    }
  }
}

cat("Jobs:", length(jobs), "\n")

# -----------------------------
# Run one job (standardized row)
# -----------------------------
run_one_job <- function(job) {
  df_g <- series_list[[job$series_i]]
  y <- df_g$value
  
  y_train <- y[job$train_idx]
  n_non_na <- sum(is.finite(y_train))
  y_obs <- y[job$test_i]
  
  ok_row <- TRUE
  reason <- NA_character_
  
  out <- list(ok = FALSE, mu = NA_real_,
              lo80 = NA_real_, hi80 = NA_real_,
              lo95 = NA_real_, hi95 = NA_real_,
              reason = NA_character_)
  
  if (!is.finite(y_obs)) {
    ok_row <- FALSE
    reason <- "missing_target_month"
  }
  if (n_non_na < MIN_TRAIN_NON_NA) {
    ok_row <- FALSE
    reason <- paste0("insufficient_training_nonNA<", MIN_TRAIN_NON_NA)
  }
  
  if (ok_row) {
    x_train_df <- df_g[job$train_idx, xreg_cols, drop = FALSE]
    x_fut_df   <- df_g[job$test_path_idx, xreg_cols, drop = FALSE]
    
    need_ok <- complete.cases(x_train_df) & is.finite(y_train)
    if (!any(need_ok)) {
      ok_row <- FALSE
      reason <- "no_complete_cases_in_training_y_xreg"
    } else if (any(!complete.cases(x_fut_df))) {
      ok_row <- FALSE
      reason <- "missing_xreg_in_forecast_path"
    } else {
      y_fit <- y_train[need_ok]
      x_fit <- as_xreg_mat(x_train_df[need_ok, , drop = FALSE], xreg_cols)
      x_fut <- as_xreg_mat(x_fut_df, xreg_cols)
      
      if (XREG_SCALE) {
        scaled <- scale_train_apply(x_fit, x_fut)
        x_fit <- scaled$x_train
        x_fut <- scaled$x_future
      }
      
      out <- tryCatch({
        fit <- forecast::auto.arima(
          y_fit,
          xreg = x_fit,
          seasonal = TRUE,
          stepwise = TRUE,
          approximation = FALSE
        )
        fc <- forecast::forecast(fit, h = job$horizon, xreg = x_fut, level = PI_LEVELS)
        
        mu_h <- as.numeric(fc$mean[job$horizon])
        lo80 <- hi80 <- lo95 <- hi95 <- NA_real_
        if (!is.null(fc$lower) && !is.null(fc$upper)) {
          if ("80%" %in% colnames(fc$lower)) {
            lo80 <- as.numeric(fc$lower[job$horizon, "80%"])
            hi80 <- as.numeric(fc$upper[job$horizon, "80%"])
          }
          if ("95%" %in% colnames(fc$lower)) {
            lo95 <- as.numeric(fc$lower[job$horizon, "95%"])
            hi95 <- as.numeric(fc$upper[job$horizon, "95%"])
          }
        }
        list(ok = TRUE, mu = mu_h, lo80 = lo80, hi80 = hi80, lo95 = lo95, hi95 = hi95, reason = NA_character_)
      }, error = function(e) {
        list(ok = FALSE, mu = NA_real_, lo80 = NA_real_, hi80 = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
             reason = conditionMessage(e))
      })
      
      if (!isTRUE(out$ok)) {
        ok_row <- FALSE
        reason <- out$reason %||% "model_fit_failed"
      }
    }
  }
  
  tibble(
    BAY = job$BAY, GEAR = job$GEAR, q = job$q,
    model = MODEL_NAME,
    horizon = job$horizon,
    origin_i = job$origin_i,
    test_i = job$test_i,
    origin_date = df_g$date[job$origin_i],
    test_date   = df_g$date[job$test_i],
    y_obs = y_obs,
    y_hat = out$mu,
    lower80 = out$lo80, upper80 = out$hi80,
    lower95 = out$lo95, upper95 = out$hi95,
    ok = ok_row,
    reason = reason
  )
}

res <- future.apply::future_lapply(jobs, run_one_job, future.seed = TRUE)
pred_tbl <- bind_rows(res)

fit_log <- pred_tbl %>%
  transmute(BAY, GEAR, q, horizon, origin_i, test_i, test_date, ok, reason)
readr::write_csv(fit_log, OUT_FITLOG)

metrics_tbl <- pred_tbl %>%
  group_by(BAY, GEAR, q, horizon, model) %>%
  summarise(
    n = sum(is.finite(y_obs) & is.finite(y_hat) & ok),
    rmse = rmse(y_hat - y_obs),
    mae  = mae(y_hat - y_obs),
    coverage80 = coverage_rate(y_obs, lower80, upper80),
    coverage95 = coverage_rate(y_obs, lower95, upper95),
    frac_failed = mean(!ok),
    .groups = "drop"
  )

save_obj <- list(
  predictions     = pred_tbl,
  eval_target_col = EVAL_TARGET_COL,
  target_scale    = TARGET_SCALE,
  script_stem     = script_stem,
  created         = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

saveRDS(save_obj, OUT_FC_DERIVED)
saveRDS(save_obj, OUT_FC_MODEL)
readr::write_csv(metrics_tbl, OUT_MET_DERIVED)
readr::write_csv(metrics_tbl, OUT_MET_MODEL)

cat("\nWrote:\n")
cat("  -", OUT_FC_DERIVED, "\n")
cat("  -", OUT_FC_MODEL, "\n")
cat("  -", OUT_XREG_MISS, "\n")
cat("  -", OUT_MET_DERIVED, "\n")
cat("  -", OUT_MET_MODEL, "\n")
cat("  -", OUT_FITLOG, "\n")

future::plan(future::sequential)
cat("\nDone.\n")

sink()
cat("Log written to:", log_file, "\n")