# =============================================================================
# Script: e7_arimax_env_forecasts_abundance.R  (REVISED to match c9)
#
# Purpose:
# Fit SARIMAX (ARIMAX with seasonal terms) for population abundance (CPUE)
# using environmental predictors (temperature and salinity + lags).
#
# Key fixes vs prior e7:
#   1) series_list passed explicitly to workers via series_i index (not get())
#   2) env aggregated to BAY="ALL" via group_by+summarise (not bind_cols)
#   3) bays_to_process = "ALL" to match e4 fold design (coastwide mean)
#   4) Single na.approx call (no duplicate overwrite)
#   5) xreg NA check uses complete.cases only where y is observed (c9 pattern)
#   6) xreg_start_i trims training window for lag warm-up (c9/c10 pattern)
#
# Outputs:
# - 01_Analysis/Data/e7_arimax_env_forecasts_abundance/arimax_forecasts_abundance.rds
# - 02_Results/04_Models/e7_arimax_env_forecasts_abundance/arimax_metrics_abundance.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap
# ─────────────────────────────────────────────────────────────────────────────
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
cat("Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

options(scipen = 999, warn = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Packages & Parallel Config
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(zoo)
  library(forecast)
  library(future)
  library(future.apply)
})

N_WORKERS <- min(6, max(1, parallel::detectCores(logical = FALSE) - 2))
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_WORKERS)
} else {
  future::plan(future::multicore, workers = N_WORKERS)
}
options(future.globals.maxSize = 4 * 1024^3)
cat("Parallel plan initialized with", N_WORKERS, "workers\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Directories
# ─────────────────────────────────────────────────────────────────────────────
script_stem      <- "e7_arimax_env_forecasts_abundance"
dir_data_script  <- file.path(PATHS$analysis_dir, "Data",              script_stem)
dir_model_script <- file.path(PATHS$results_dir,  "04_Models",         script_stem)
dir_diag_script  <- file.path(PATHS$results_dir,  "03_Diagnostics",    script_stem)

dir.create(dir_data_script,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script,  recursive = TRUE, showWarnings = FALSE)

### ─────────────────────────────────────────────────────────────────────────────
### 4) Parameters & Inputs
### ─────────────────────────────────────────────────────────────────────────────
### TEST_RUN: set TRUE to process 2 species x 1 gear only (quick sanity check)
TEST_RUN <- FALSE

TARGET_COL <- "value_deseas"
TARGET_SCALE <- "deseasonalized"
PI_LEVELS <- c(80, 95)
FREQ <- 12
MIN_TRAIN_NON_NA <- 24
MIN_TRAIN_NON_ZERO <- 24

### Environmental predictors (must exist in env_gamma_index after standardization)
XREG_COLS <- c(
  "TEMP__mean", "TEMP__mean__lag1", "TEMP__mean__lag12",
  "SAL__mean",  "SAL__mean__lag1",  "SAL__mean__lag12"
)

MODEL_NAME <- "sarimax"

IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance",        "alpha_abundance_deseas.rds")
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance",        "gamma_abundance_deseas.rds")
IN_FOLDS <- file.path(PATHS$analysis_dir, "Data", "e4_forecast_design_abundance",      "backtest_folds_abundance.rds")
IN_ENV_A <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries",   "env_alpha_index.rds")
IN_ENV_G <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries",   "env_gamma_index.rds")

stopifnot(file.exists(IN_ALPHA), file.exists(IN_GAMMA), file.exists(IN_FOLDS), file.exists(IN_ENV_A), file.exists(IN_ENV_G))

OUT_FCST_DERIVED <- file.path(dir_data_script,  "arimax_forecasts_abundance.rds")
OUT_FCST_MODEL   <- file.path(dir_model_script, "arimax_forecasts_abundance.rds")
OUT_METR_DERIVED <- file.path(dir_data_script,  "arimax_metrics_abundance.csv")
OUT_METR_MODEL   <- file.path(dir_model_script, "arimax_metrics_abundance.csv")
OUT_FITLOG       <- file.path(dir_diag_script,  "arimax_fit_log.csv")
OUT_XREG_MISS    <- file.path(dir_data_script,  "xreg_missingness_summary.csv")

### ─────────────────────────────────────────────────────────────────────────────
### 5) Load & Aggregate Abundance Data 
### ─────────────────────────────────────────────────────────────────────────────
cat("Loading deseasonalized abundance data...\n")
alpha <- readRDS(IN_ALPHA)
gamma <- readRDS(IN_GAMMA)

### Alpha: Retain individual bay structure
alpha_obs <- alpha %>% filter(as.character(BAY_ID) != "1") %>% transmute(
  series       = "alpha",
  BAY          = as.character(BAY_ID),
  GEAR         = as.character(GEAR),
  SPECIES_CODE = as.character(SPECIES_CODE),
  date         = as.Date(date),
  value        = as.numeric(.data[[TARGET_COL]])
)

gamma_obs <- gamma %>% transmute(
  series       = "gamma",
  BAY          = "ALL",
  GEAR         = as.character(GEAR),
  SPECIES_CODE = as.character(SPECIES_CODE),
  date         = as.Date(date),
  value        = as.numeric(.data[[TARGET_COL]])
)

ts_obs <- bind_rows(alpha_obs, gamma_obs)

### ─────────────────────────────────────────────────────────────────────────────
### 6) Load & Aggregate Environmental Data
### ─────────────────────────────────────────────────────────────────────────────
cat("Loading environmental data...\n")

### Alpha env: Retain individual bay structure
env_alpha <- readRDS(IN_ENV_A) %>% 
  filter(as.character(BAY_ID) != "1") %>% 
  mutate(date = as.Date(DATE), GEAR = as.character(GEAR), BAY = as.character(BAY_ID)) %>% 
  group_by(BAY, GEAR, date) %>% 
  summarise(across(all_of(XREG_COLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  mutate(series = "alpha")

env_gamma <- readRDS(IN_ENV_G) %>% 
  mutate(date = as.Date(DATE), GEAR = as.character(GEAR)) %>% 
  group_by(GEAR, date) %>% 
  summarise(across(all_of(XREG_COLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  mutate(series = "gamma", BAY = "ALL")

env_obs <- bind_rows(env_alpha, env_gamma)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build series_list and jobs — following c9 architecture
#    series_i indexes into series_list so workers never use get()
# ─────────────────────────────────────────────────────────────────────────────
cat("Aligning temporal grids and building job list...\n")

folds       <- readRDS(IN_FOLDS)
series_keys <- list()   # metadata for each series
series_list <- list()   # aligned dataframes (indexed by series_i)
jobs        <- list()
xreg_miss_rows <- list()

series_i <- 0L

for (gid in names(folds$groups)) {
  g <- folds$groups[[gid]]
  if (!isTRUE(g$ok)) next
  
  series_i <- series_i + 1L
  bay_id   <- as.character(g$BAY_ID)   # Pull bay_id directly from the fold group
  dates_grid <- tibble(date = as.Date(g$dates))
  
  ### Join abundance and environment — single na.approx (no duplicate)
  df_g <- dates_grid %>% left_join(
      ts_obs %>%
        filter(series == g$series, BAY == bay_id,
               GEAR == g$GEAR, SPECIES_CODE == g$SPECIES_CODE),
      by = "date"
    ) %>%
    left_join(
      env_obs %>%
        filter(series == g$series, BAY == bay_id, GEAR == g$GEAR) %>%
        select(date, all_of(XREG_COLS)),
      by = "date"
    ) %>%
    arrange(date) %>%
    # Single interpolation pass — lag warm-up NAs at start intentionally kept
    # so xreg_start_i can identify the first complete row (c9/c10 pattern)
    mutate(across(all_of(XREG_COLS),
                  ~ as.numeric(zoo::na.approx(.x, na.rm = FALSE))))

  # First row where ALL xreg columns are complete
  ok_x        <- complete.cases(df_g[, XREG_COLS, drop = FALSE])
  xreg_start_i <- if (any(ok_x)) which(ok_x)[1] else NA_integer_

  # Missingness summary
  miss_rate <- sapply(XREG_COLS, function(cc) mean(is.na(df_g[[cc]])))
  xreg_miss_rows[[series_i]] <- tibble(
    series       = g$series,
    BAY          = bay_id,
    GEAR         = g$GEAR,
    SPECIES_CODE = g$SPECIES_CODE,
    xreg         = names(miss_rate),
    missing_rate = as.numeric(miss_rate)
  )

  series_keys[[series_i]] <- list(
    series       = g$series,
    BAY          = bay_id,
    GEAR         = g$GEAR,
    SPECIES_CODE = g$SPECIES_CODE,
    xreg_start_i = xreg_start_i
  )
  series_list[[series_i]] <- df_g

  for (hname in names(g$folds)) {
    for (f in g$folds[[hname]]) {
      test_path_idx <- (f$origin_i + 1):(f$origin_i + f$horizon)
      jobs[[length(jobs) + 1]] <- list(
        series_i      = series_i,
        horizon       = f$horizon,
        origin_i      = f$origin_i,
        test_i        = f$test_i,
        train_idx     = f$train_idx,
        test_path_idx = test_path_idx,
        origin_date   = g$dates[f$origin_i],
        test_date     = g$dates[f$test_i],
        xreg_start_i  = xreg_start_i
      )
    }
  }
}

write_csv(bind_rows(xreg_miss_rows), OUT_XREG_MISS)
cat("Total jobs:", length(jobs), "\n")

# TEST_RUN: limit to first 20 jobs
if (TEST_RUN) {
  jobs <- head(jobs, 20)
  cat("TEST_RUN: limited to", length(jobs), "jobs\n")
}

### ─────────────────────────────────────────────────────────────────────────────
### 8) Worker Function
### ─────────────────────────────────────────────────────────────────────────────
run_job <- function(job, series_keys, series_list, xreg_cols, min_train_non_na, model_name, freq, pi_levels) {
  
  key  <- series_keys[[job$series_i]]
  df_g <- series_list[[job$series_i]]
  y    <- df_g$value
  
  ### Trim training window for lag warm-up
  xreg_start_i <- job$xreg_start_i
  if (is.null(xreg_start_i) || is.na(xreg_start_i)) {
    return(tibble(
      series = key$series, BAY = key$BAY, GEAR = key$GEAR, SPECIES_CODE = key$SPECIES_CODE, model = model_name,
      horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i, origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
      y_obs = y[job$test_i], y_hat = NA_real_, lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_, ok = FALSE, reason = "no_complete_xreg_window"
    ))
  }
  
  train_eff <- job$train_idx[job$train_idx >= xreg_start_i]
  if (length(train_eff) == 0) {
    return(tibble(
      series = key$series, BAY = key$BAY, GEAR = key$GEAR, SPECIES_CODE = key$SPECIES_CODE, model = model_name,
      horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
      origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
      y_obs = y[job$test_i], y_hat = NA_real_, lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_, ok = FALSE, reason = "train_window_before_xreg_start"
    ))
  }
  
  y_train <- y[train_eff]
  y_obs   <- y[job$test_i]
  n_non_na <- sum(is.finite(y_train))
  n_non_zero <- sum(y_train > 0, na.rm = TRUE)
  
  make_fail <- function(reason) {
    tibble(
      series = key$series, BAY = key$BAY, GEAR = key$GEAR,
      SPECIES_CODE = key$SPECIES_CODE, model = model_name, horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
      origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date), y_obs = y_obs, y_hat = NA_real_, lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_, ok = FALSE, reason = reason
    )
  }
  
  x_train_df <- df_g[train_eff, xreg_cols, drop = FALSE]
  x_fut_df   <- df_g[job$test_path_idx, xreg_cols, drop = FALSE]
  
  xreg_train <- as.matrix(x_train_df)
  xreg_fut   <- as.matrix(x_fut_df)
  
  if (!is.finite(y_obs)) {
    return(make_fail("missing_target_month"))
  } else if (n_non_na < min_train_non_na) {
    return(make_fail(paste0("insufficient_training_data<", min_train_non_na)))
  } else if (n_non_zero < MIN_TRAIN_NON_ZERO) {
    return(make_fail("insufficient_non_zero_data"))
  } else if (any(!is.finite(xreg_train[is.finite(y_train), , drop = FALSE]))) {
    return(make_fail("missing_xreg_in_train"))
  } else if (any(!is.finite(xreg_fut))) {
    return(make_fail("missing_xreg_in_forecast_horizon"))
  }
  
  xreg_train[is.na(xreg_train)] <- 0
  
  need_ok <- is.finite(y_train)
  if (sum(need_ok) < min_train_non_na) {
    return(make_fail("too_few_finite_y_after_xreg_filter"))
  }
  
  y_fit <- y_train[need_ok]
  x_fit <- as.matrix(xreg_train[need_ok, , drop = FALSE])
  x_fut <- as.matrix(xreg_fut)
  storage.mode(x_fit) <- "double"
  storage.mode(x_fut) <- "double"
  
  y_ts <- stats::ts(forecast::na.interp(y_fit), frequency = freq)
  
  out <- tryCatch({
    fit <- suppressWarnings(forecast::auto.arima(
      y_ts, xreg = x_fit, seasonal = TRUE, stepwise = TRUE, approximation = FALSE
    ))
    fc <- forecast::forecast(fit, h = job$horizon, xreg = x_fut, level = pi_levels)
    list(ok = TRUE, mu = as.numeric(fc$mean[job$horizon]), lo80 = as.numeric(fc$lower[job$horizon, 1]), hi80 = as.numeric(fc$upper[job$horizon, 1]), lo95 = as.numeric(fc$lower[job$horizon, 2]), hi95 = as.numeric(fc$upper[job$horizon, 2]), reason = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, mu = NA_real_, lo80 = NA_real_, hi80 = NA_real_, lo95 = NA_real_, hi95 = NA_real_, reason = conditionMessage(e))
  })
  
  tibble(
    series = key$series, BAY = key$BAY, GEAR = key$GEAR,
    SPECIES_CODE = key$SPECIES_CODE, model = model_name, horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
    origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
    y_obs   = y_obs,
    y_hat   = out$mu, lower80 = out$lo80, upper80 = out$hi80, lower95 = out$lo95, upper95 = out$hi95, ok      = isTRUE(out$ok),
    reason  = out$reason
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 9) Execute in Parallel — pass series_list explicitly to workers
# ─────────────────────────────────────────────────────────────────────────────
cat("Fitting SARIMAX models...\n")

res_list <- future.apply::future_lapply(
  jobs, run_job,
  series_keys      = series_keys,
  series_list      = series_list,
  xreg_cols        = XREG_COLS,
  min_train_non_na = MIN_TRAIN_NON_NA,
  model_name       = MODEL_NAME,
  freq             = FREQ,
  pi_levels        = PI_LEVELS,
  future.seed      = TRUE
)

pred_tbl <- dplyr::bind_rows(res_list) %>%
  mutate(
    error     = y_obs - y_hat,
    abs_error = abs(error),
    covered80 = ifelse(is.na(lower80) | is.na(upper80), NA,
                       y_obs >= lower80 & y_obs <= upper80),
    covered95 = ifelse(is.na(lower95) | is.na(upper95), NA,
                       y_obs >= lower95 & y_obs <= upper95)
  ) %>%
  arrange(series, BAY, GEAR, SPECIES_CODE, horizon, origin_i)

# ─────────────────────────────────────────────────────────────────────────────
# 10) Metrics
# ─────────────────────────────────────────────────────────────────────────────
metrics_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(series, BAY, GEAR, SPECIES_CODE, horizon, model) %>%
  summarise(
    n          = n(),
    RMSE       = sqrt(mean(error^2)),
    MAE        = mean(abs_error),
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    .groups    = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 11) Save
# ─────────────────────────────────────────────────────────────────────────────
save_obj <- list(
  created        = Sys.time(),
  eval_target_col = TARGET_COL,
  target_scale   = TARGET_SCALE,
  PI_LEVELS      = PI_LEVELS,
  xreg_cols      = XREG_COLS,
  predictions    = pred_tbl
)

saveRDS(save_obj, OUT_FCST_DERIVED)
saveRDS(save_obj, OUT_FCST_MODEL)
write_csv(metrics_tbl, OUT_METR_DERIVED)
write_csv(metrics_tbl, OUT_METR_MODEL)

fitlog <- pred_tbl %>%
  mutate(reason2 = ifelse(ok, "ok",
                          ifelse(is.na(reason), "fail_unknown", reason))) %>%
  count(series, BAY, GEAR, SPECIES_CODE, horizon, model, reason2, name = "n") %>%
  arrange(desc(n))
write_csv(fitlog, OUT_FITLOG)

future::plan(future::sequential)

cat("\nSARIMAX forecasting complete.\n")
cat("Predictions:", OUT_FCST_MODEL, "\n")
cat("Metrics:",     OUT_METR_MODEL, "\n")
cat("Fit log:",     OUT_FITLOG,     "\n")
