# =============================================================================
# Script: e8_ssm_xreg_forecasts_abundance.R
#
# Purpose:
# Fit KFAS state-space models with environmental regressors (Temperature, 
# Salinity, and lags) for population abundance (CPUE).
# 
# This script unifies the logic of c10 (Alpha) and c11 (Gamma).
# Consistent with spatial modeling protocols, Sabine Lake (Bay 1) is excluded 
# from Alpha (bay-level) predictions.
#
# Outputs:
# - 01_Analysis/Data/e8_ssm_xreg_forecasts_abundance/ssm_xreg_forecasts_abundance.rds
# - 02_Results/04_Models/e8_ssm_xreg_forecasts_abundance/ssm_xreg_metrics_abundance.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS (Windows & Linux Safe)
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
  ofile <- tryCatch(sys.frames()[]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
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
  library(KFAS)
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
cat("🚀 Initialized parallel plan with", N_WORKERS, "workers\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Directories
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e8_ssm_xreg_forecasts_abundance"
dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)

### ─────────────────────────────────────────────────────────────────────────────
### 4) Parameters & Inputs
### ─────────────────────────────────────────────────────────────────────────────
### TEST_RUN: set TRUE to process a small subset of jobs (quick sanity check)
TEST_RUN <- FALSE

MODEL_NAME <- "ssm_ll_seas12_xreg"
TARGET_COL <- "value_deseas"
TARGET_SCALE <- "deseasonalized"
PI_LEVELS <- c(80, 95)
MAXIT <- 200
MIN_TRAIN_NON_NA <- 24
MIN_TRAIN_NON_ZERO <- 24

XREG_COLS <- c(
  "TEMP__mean", "TEMP__mean__lag1", "TEMP__mean__lag12",
  "SAL__mean",  "SAL__mean__lag1",  "SAL__mean__lag12"
)

IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
IN_FOLDS <- file.path(PATHS$analysis_dir, "Data", "e4_forecast_design_abundance", "backtest_folds_abundance.rds")
IN_ENV_A <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries", "env_alpha_index.rds")
IN_ENV_G <- file.path(PATHS$analysis_dir, "Data", "b8_build_environment_timeseries", "env_gamma_index.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 5) KFAS Helper Functions (Covariate Enabled)
# ─────────────────────────────────────────────────────────────────────────────
build_ssm_xreg <- function(y, xreg_mat) {
  form <- as.formula("y ~ -1 + xreg_mat + SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type = 'dummy', Q = NA)")
  KFAS::SSModel(formula = form, H = NA)
}

fit_ssm_ml <- function(model) {
  updatefn <- function(pars, model) {
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
  
  n_par <- sum(is.na(model$H)) + sum(is.na(model$Q))
  v0 <- stats::var(model$y, na.rm = TRUE)
  if (!is.finite(v0) || v0 <= 0) v0 <- 1
  init <- rep(log(v0), n_par)
  
  fit <- tryCatch({
    KFAS::fitSSM(model, inits = init, updatefn = updatefn, method = "BFGS", control = list(maxit = MAXIT, trace = 0))
  }, error = function(e) e)
  
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), model = model))
  
  conv <- fit$optim.out$convergence
  msg  <- fit$optim.out$message
  list(
    ok = isTRUE(conv == 0),
    reason = if (isTRUE(conv == 0)) NA_character_ else paste0("convergence=", conv, " ", msg),
    model = fit$model
  )
}

forecast_ssm_xreg <- function(fitted_model, h, xreg_fut) {
  tryCatch({
    # Build future model structure — Q = NA placeholders initially
    new_mod <- KFAS::SSModel(
      rep(NA_real_, h) ~ -1 + xreg_fut +
        SSMtrend(1, Q = NA) +
        SSMseasonal(12, sea.type = 'dummy', Q = NA),
      H = fitted_model$H   # copy fitted observation variance directly
    )
    
    # Copy fitted state variances so new_mod has no NA placeholders
    new_mod$Q <- fitted_model$Q
    
    # Now predict using the fully specified model
    pred   <- predict(fitted_model, newdata = new_mod,
                      interval = "prediction", level = 0.95)
    pred80 <- predict(fitted_model, newdata = new_mod,
                      interval = "prediction", level = 0.80)
    
    list(
      ok      = TRUE,
      reason  = NA_character_,
      mean    = as.numeric(pred[h, "fit"]),
      lower80 = as.numeric(pred80[h, "lwr"]),
      upper80 = as.numeric(pred80[h, "upr"]),
      lower95 = as.numeric(pred[h, "lwr"]),
      upper95 = as.numeric(pred[h, "upr"])
    )
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e),
         mean = NA_real_, lower80 = NA_real_, upper80 = NA_real_,
         lower95 = NA_real_, upper95 = NA_real_)
  })
}

### ─────────────────────────────────────────────────────────────────────────────
### 6) Load & Standardize Abundance and Environmental Data
### ─────────────────────────────────────────────────────────────────────────────
cat("Loading deseasonalized abundance and environmental data...\n")
alpha <- readRDS(IN_ALPHA)
gamma <- readRDS(IN_GAMMA)

alpha_obs <- alpha %>% filter(as.character(BAY_ID) != "1") %>% transmute(
  series = "alpha", 
  BAY = as.character(BAY_ID), 
  GEAR = as.character(GEAR), 
  SPECIES_CODE = as.character(SPECIES_CODE), 
  date = as.Date(date), 
  value = as.numeric(.data[[TARGET_COL]])
)

gamma_obs <- gamma %>% transmute(
  series = "gamma", 
  BAY = "ALL", 
  GEAR = as.character(GEAR), 
  SPECIES_CODE = as.character(SPECIES_CODE), 
  date = as.Date(date), 
  value = as.numeric(.data[[TARGET_COL]])
)

ts_obs <- bind_rows(alpha_obs, gamma_obs)

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

### ─────────────────────────────────────────────────────────────────────────────
### 7) Prepare Fold Jobs & Align Covariates
### ─────────────────────────────────────────────────────────────────────────────
folds <- readRDS(IN_FOLDS)
jobs <- list()
series_list <- list() # Explicit list to store dataframes for parallel workers

cat("Aligning temporal grids and covariates...\n")
for (gid in names(folds$groups)) {
  g <- folds$groups[[gid]]
  if (!isTRUE(g$ok)) next
  
  bays_to_process <- as.character(g$BAY_ID)   # Pull bay_id dynamically
  dates_grid <- tibble(date = as.Date(g$dates))
  
  for (bay_id in bays_to_process) {
    df_g <- dates_grid %>%
      left_join(
        ts_obs %>% filter(series == g$series, BAY == bay_id, GEAR == g$GEAR, SPECIES_CODE == g$SPECIES_CODE),
        by = "date"
      ) %>%
      left_join(
        env_obs %>%
          filter(series == g$series, BAY == bay_id, GEAR == g$GEAR) %>%
          select(date, all_of(XREG_COLS)),
        by = "date"
      ) %>%
      arrange(date) %>%
      # Interpolate interior NAs only (rule = 1, default).
      # Lag warm-up NAs at the start of the series are intentionally left as NA
      # so that xreg_start_i can identify the first complete row and trim the
      # training window accordingly — matching the c10 approach.
      mutate(across(all_of(XREG_COLS), ~ as.numeric(zoo::na.approx(.x, na.rm = FALSE))))

    # First row where ALL xreg columns are complete (handles lag warm-up)
    ok_x <- complete.cases(df_g[, XREG_COLS, drop = FALSE])
    xreg_start_i <- if (any(ok_x)) which(ok_x)[1] else NA_integer_

    ts_key <- paste(g$series, bay_id, g$GEAR, g$SPECIES_CODE, sep = "_")
    series_list[[ts_key]] <- df_g

    for (hname in names(g$folds)) {
      for (f in g$folds[[hname]]) {
        test_path_idx <- (f$origin_i + 1):(f$origin_i + f$horizon)

        jobs[[length(jobs) + 1]] <- list(
          ts_key = ts_key, series = g$series, BAY = bay_id, GEAR = g$GEAR, SPECIES_CODE = g$SPECIES_CODE,
          horizon = f$horizon, origin_i = f$origin_i, test_i = f$test_i,
          train_idx = f$train_idx, test_path_idx = test_path_idx,
          origin_date = g$dates[f$origin_i], test_date = g$dates[f$test_i],
          xreg_start_i = xreg_start_i
        )
      }
    }
  }
}
### TEST_RUN: limit to first 20 jobs for quick validation
if (TEST_RUN) { 
  jobs <- head(jobs, 20) 
  cat("TEST_RUN: limited to", length(jobs), "jobs\n") 
}

### ─────────────────────────────────────────────────────────────────────────────
### 8) Worker Function
### ─────────────────────────────────────────────────────────────────────────────
run_job <- function(job, series_data) {
  df <- series_data[[job$ts_key]]
  
  ### Trim training window to start at xreg_start_i (handles lag warm-up NAs)
  xreg_start_i <- job$xreg_start_i
  if (is.null(xreg_start_i) || is.na(xreg_start_i)) {
    return(tibble::tibble(
      series = job$series, BAY = job$BAY, GEAR = job$GEAR, SPECIES_CODE = job$SPECIES_CODE,
      model = MODEL_NAME, horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
      origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
      y_obs = as.numeric(df$value[job$test_i]), y_hat = NA_real_, lower80 = NA_real_, upper80 = NA_real_,
      lower95 = NA_real_, upper95 = NA_real_, ok = FALSE, reason = "no_complete_xreg_window"
    ))
  }
  
  ### Effective training indices: only rows at or after xreg_start_i
  train_eff <- job$train_idx[job$train_idx >= xreg_start_i]
  if (length(train_eff) == 0) {
    return(tibble::tibble(
      series = job$series, BAY = job$BAY, GEAR = job$GEAR, SPECIES_CODE = job$SPECIES_CODE,
      model = MODEL_NAME, horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
      origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
      y_obs = as.numeric(df$value[job$test_i]), y_hat = NA_real_, lower80 = NA_real_, upper80 = NA_real_,
      lower95 = NA_real_, upper95 = NA_real_, ok = FALSE, reason = "train_window_before_xreg_start"
    ))
  }
  
  y_train <- df$value[train_eff]
  y_test  <- df$value[job$test_i]
  
  xreg_train <- as.matrix(df[train_eff,        XREG_COLS, drop = FALSE])
  xreg_fut   <- as.matrix(df[job$test_path_idx, XREG_COLS, drop = FALSE])
  
  n_non_na <- sum(is.finite(y_train))
  n_non_zero <- sum(y_train > 0, na.rm = TRUE)
  
  ok_row <- TRUE
  reason <- NA_character_
  fc <- list(ok = FALSE, mean = NA_real_, lower80 = NA_real_, upper80 = NA_real_, lower95 = NA_real_, upper95 = NA_real_, reason = NA_character_)
  
  if (!is.finite(y_test)) {
    ok_row <- FALSE; reason <- "missing_target_month"
  } else if (n_non_na < MIN_TRAIN_NON_NA) {
    ok_row <- FALSE; reason <- paste0("insufficient_training_data<", MIN_TRAIN_NON_NA)
  } else if (n_non_zero < MIN_TRAIN_NON_ZERO) {
    ok_row <- FALSE; reason <- "insufficient_non_zero_data"
  } else if (any(!is.finite(xreg_train[is.finite(y_train), , drop = FALSE]))) {
    # Only reject if xreg is missing where y is observed (strict check from c10)
    ok_row <- FALSE; reason <- "missing_xreg_in_train"
  } else if (any(!is.finite(xreg_fut))) {
    ok_row <- FALSE; reason <- "missing_xreg_in_forecast_horizon"
  } else {
    # Fill any remaining structural NAs in training (where y is NA) with 0
    xreg_train[is.na(xreg_train)] <- 0
    
    ssm0 <- build_ssm_xreg(y_train, xreg_train)
    fit  <- fit_ssm_ml(ssm0)
    
    if (!isTRUE(fit$ok)) {
      ok_row <- FALSE
      reason <- fit$reason
    } else {
      fc_ssm <- forecast_ssm_xreg(fit$model, h = job$horizon, xreg_fut = xreg_fut)
      if (!isTRUE(fc_ssm$ok)) {
        ok_row <- FALSE
        reason <- paste0("forecast_failed: ", fc_ssm$reason)
      } else {
        fc$mean <- fc_ssm$mean
        fc$lower80 <- fc_ssm$lower80
        fc$upper80 <- fc_ssm$upper80
        fc$lower95 <- fc_ssm$lower95
        fc$upper95 <- fc_ssm$upper95
        fc$ok <- TRUE
      }
    }
  }
  
  tibble::tibble(
    series = job$series, BAY = job$BAY, GEAR = job$GEAR, SPECIES_CODE = job$SPECIES_CODE,
    model = MODEL_NAME, horizon = job$horizon, origin_i = job$origin_i, test_i = job$test_i,
    origin_date = as.Date(job$origin_date), test_date = as.Date(job$test_date),
    y_obs = as.numeric(y_test), y_hat = fc$mean, lower80 = fc$lower80, upper80 = fc$upper80,
    lower95 = fc$lower95, upper95 = fc$upper95, ok = ok_row, reason = reason
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 9) Execute Jobs Parallel & Consolidate
# ─────────────────────────────────────────────────────────────────────────────
cat("Forecasting Environmental State-Space models...\n")
# Pass the series_data list explicitly to the workers

res_list <- future.apply::future_lapply(jobs, run_job, series_data = series_list, future.seed = TRUE)
pred_tbl <- dplyr::bind_rows(res_list)

pred_tbl <- pred_tbl %>%
  mutate(
    error = y_obs - y_hat,
    abs_error = abs(error),
    covered80 = ifelse(is.na(lower80) | is.na(upper80), NA, y_obs >= lower80 & y_obs <= upper80),
    covered95 = ifelse(is.na(lower95) | is.na(upper95), NA, y_obs >= lower95 & y_obs <= upper95)
  ) %>% arrange(series, BAY, GEAR, SPECIES_CODE, horizon, origin_i, model)

# ─────────────────────────────────────────────────────────────────────────────
# 10) Calculate Basic Accuracy Metrics
# ─────────────────────────────────────────────────────────────────────────────
metrics_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(series, BAY, GEAR, SPECIES_CODE, horizon, model) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(error^2)),
    MAE  = mean(abs_error),
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 11) Save Outputs & Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
save_obj <- list(
  created = Sys.time(), eval_target_col = TARGET_COL, target_scale = TARGET_SCALE,
  PI_LEVELS = PI_LEVELS, xreg_cols = XREG_COLS, predictions = pred_tbl
)

OUT_FCST_DERIVED <- file.path(dir_data_script, "ssm_xreg_forecasts_abundance.rds")
OUT_FCST_MODEL   <- file.path(dir_model_script, "ssm_xreg_forecasts_abundance.rds")
saveRDS(save_obj, OUT_FCST_DERIVED)
saveRDS(save_obj, OUT_FCST_MODEL)

OUT_METR_DERIVED <- file.path(dir_data_script, "ssm_xreg_metrics_abundance.csv")
OUT_METR_MODEL   <- file.path(dir_model_script, "ssm_xreg_metrics_abundance.csv")
write_csv(metrics_tbl, OUT_METR_DERIVED)
write_csv(metrics_tbl, OUT_METR_MODEL)

fitlog <- pred_tbl %>%
  mutate(reason2 = ifelse(ok, "ok", ifelse(is.na(reason), "fail_unknown", reason))) %>%
  count(series, BAY, GEAR, SPECIES_CODE, horizon, model, reason2, name = "n") %>%
  arrange(desc(n))
write_csv(fitlog, file.path(dir_diag_script, "ssm_xreg_fit_log.csv"))

cat("\n✅ SSM XREG forecasting complete.\n")
cat("💾 Saved predictions to:", OUT_FCST_MODEL, "\n")
cat("💾 Saved metrics to:", OUT_METR_MODEL, "\n")
cat("🧾 Diagnostics logged in:", dir_diag_script, "\n")