# =============================================================================
# Script: c5_state_space_models_forecasts.R  (STANDALONE: NO c4p dependency)
#
# Purpose:
#   Fit a KFAS state-space model (local level + seasonal(12) + obs noise [+ optional xreg])
#   for each rolling-origin fold for γ diversity, and score forecasts.
#
# Inputs (authoritative; per project conventions):
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/gamma_q012_deseas.rds
#   - 01_Analysis/Data/c3_forecast_design_backtesting/backtest_folds.rds
#
# Outputs:
#   Derived (downstream-friendly):
#     - 01_Analysis/Data/c5_state_space_models_forecasts/ssm_fits_gamma.rds
#     - 01_Analysis/Data/c5_state_space_models_forecasts/ssm_metrics_gamma.csv
#   Results (reporting):
#     - 02_Results/04_Models/c5_state_space_models_forecasts/ssm_fits_gamma.rds
#     - 02_Results/04_Models/c5_state_space_models_forecasts/ssm_metrics_gamma.csv
#   Diagnostics:
#     - 02_Results/03_Diagnostics/c5_state_space_models_forecasts/*.csv
#     - 02_Results/03_Diagnostics/c5_state_space_models_forecasts/Figures/*.(pdf|tiff)  (if ggplot2)
#   Logs:
#     - 01_Analysis/Logs/c5_state_space_models_forecasts/runlog_*.txt
#
# Key hardening:
#   - Always align gamma to fold group's COMPLETE monthly date grid (g$dates).
#   - Collapse duplicates in gamma (GEAR,q,date) via median as safety.
#   - Prefer target column "value" if present; otherwise fall back to value_deseas/residual_component.
#   - Pre-build aligned group tables once.
#   - Save fit log with skip/failure reasons.
#
# Notes:
#   - This script computes SSM-only metrics (RMSE/MAE/coverage). It does NOT compute
#     skill vs SARIMA to avoid dependency on c4p outputs.
# =============================================================================

rm(list = ls())

# -----------------------------
# Required bootstrap snippet (must include exactly)
# -----------------------------
get_script_dir <- function() {
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  
  # 2) Rscript / R CMD BATCH: --file=...
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  
  # 3) source("script.R")
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  
  # 4) fallback: assume you launched from 01_Analysis
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)

# -----------------------------
# PATHS.RData usage (required)
# -----------------------------
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file) # loads PATHS
stopifnot(exists("PATHS"))

# Validate PATHS project_dir matches computed project_dir
if (!is.null(PATHS$project_dir)) {
  paths_proj <- normalizePath(PATHS$project_dir, winslash = "/", mustWork = TRUE)
  if (!identical(paths_proj, project_dir)) {
    stop(
      "PATHS$project_dir does not match computed project_dir.\n",
      "Computed: ", project_dir, "\n",
      "PATHS:    ", paths_proj
    )
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Minimal dependencies
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "KFAS")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(KFAS)
})

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_future  <- requireNamespace("future", quietly = TRUE)
has_future_apply <- requireNamespace("future.apply", quietly = TRUE)

# -----------------------------
# Script stem + directories
# -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c5_state_space_models_forecasts.R"))
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c5_state_space_models_forecasts"

dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script  <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_diag_fig     <- file.path(dir_diag_script, "Figures")
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Inputs (authoritative)
# -----------------------------
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "gamma_q012_deseas.rds")
IN_FOLDS <- file.path(PATHS$analysis_dir, "Data", "c3_forecast_design_backtesting", "backtest_folds.rds")

stopifnot(file.exists(IN_GAMMA), file.exists(IN_FOLDS))

# -----------------------------
# Outputs
# -----------------------------
OUT_FITS_DERIVED <- file.path(dir_data_script, "ssm_fits_gamma.rds")
OUT_METR_DERIVED <- file.path(dir_data_script, "ssm_metrics_gamma.csv")

OUT_FITS_MODEL <- file.path(dir_model_script, "ssm_fits_gamma.rds")
OUT_METR_MODEL <- file.path(dir_model_script, "ssm_metrics_gamma.csv")

# Diagnostics
OUT_FITLOG <- file.path(dir_diag_script, "ssm_fit_log.csv")

# -----------------------------
# Parallel config
# -----------------------------
RUN_PARALLEL <- TRUE
N_WORKERS <- min(6, max(1, parallel::detectCores(logical = FALSE) - 2))

if (RUN_PARALLEL && has_future && has_future_apply) {
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = N_WORKERS)
  } else {
    future::plan(future::multicore, workers = N_WORKERS)
  }
  options(future.globals.maxSize = 2 * 1024^3) # 2 GB
} else {
  if (has_future) future::plan(future::sequential)
  RUN_PARALLEL <- FALSE
}

# -----------------------------
# Run log
# -----------------------------
log_file <- file.path(dir_logs_script, paste0("runlog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
cat("Script:", script_stem, "\n")
cat("Started:", format(Sys.time()), "\n")
cat("Inputs:\n  -", IN_GAMMA, "\n  -", IN_FOLDS, "\n")
cat("Parallel:", RUN_PARALLEL, " | Workers:", if (RUN_PARALLEL) N_WORKERS else 1, "\n\n")

# -----------------------------
# Model/config
# -----------------------------
TARGET_COL_EVAL_PRIMARY <- "value"
TARGET_COL_FALLBACKS <- c("value_deseas", "residual_component")

FREQ <- 12
PI_LEVELS <- c(80, 95)

LJUNGBOX_LAG <- 12
ACF_MAX_LAG <- 36

# Optional regressors aligned to fold grid (must exist in aligned df_g)
XREG_COLS <- character(0)

MAXIT <- 200
QUIET <- TRUE

MIN_TRAIN_NON_NA <- 24

# -----------------------------
# Helpers
# -----------------------------
expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

q_label <- function(q) {
  dplyr::case_when(
    is.na(q) ~ NA_character_,
    q == 0   ~ "q = 0 (richness)",
    q == 1   ~ "q = 1 (Shannon-weighted diversity)",
    q == 2   ~ "q = 2 (dominance-weighted diversity)",
    TRUE     ~ paste0("q = ", q)
  )
}

save_fig_pdf_tiff <- function(plot, base_path_no_ext, width, height) {
  if (!has_ggplot2) return(invisible(NULL))
  ggplot2::ggsave(paste0(base_path_no_ext, ".pdf"),  plot = plot, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(base_path_no_ext, ".tiff"), plot = plot, width = width, height = height, units = "in", dpi = 300,
                  compression = "lzw")
}

pick_target_col <- function(df) {
  if (TARGET_COL_EVAL_PRIMARY %in% names(df)) return(TARGET_COL_EVAL_PRIMARY)
  fb <- TARGET_COL_FALLBACKS[TARGET_COL_FALLBACKS %in% names(df)]
  if (length(fb) >= 1) return(fb[1])
  stop("No usable target column found. Looked for: ",
       paste(c(TARGET_COL_EVAL_PRIMARY, TARGET_COL_FALLBACKS), collapse = ", "))
}

collapse_dups_gamma <- function(df, value_col) {
  df %>%
    group_by(GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- KFAS model helpers ----
# ---- KFAS model helpers ----
build_ssm <- function(y, xreg = NULL) {
  
  # NOTE: SSMtrend / SSMseasonal must NOT be namespaced inside the formula.
  # Q should be NA (not list(NA)) for local level + seasonal.
  
  if (is.null(xreg)) {
    
    SSModel(
      y ~ -1 +
        SSMtrend(1, Q = NA) +
        SSMseasonal(FREQ, sea.type = "dummy", Q = NA),
      H = NA
    )
    
  } else {
    
    # If you later add xreg: KFAS wants regression via SSMregression()
    xreg <- as.matrix(xreg)
    p <- ncol(xreg)
    SSModel(
      y ~ -1 +
        SSMtrend(1, Q = NA) +
        SSMseasonal(FREQ, sea.type = "dummy", Q = NA) +
        SSMregression(xreg, Q = diag(NA_real_, p)),
      H = NA
    )
  }
}

fit_ssm_ml <- function(model) {
  
  updatefn <- function(pars, model) {
    i <- 1
    if (anyNA(model$H)) {
      model$H[is.na(model$H)] <- exp(pars[i]); i <- i + 1
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
    KFAS::fitSSM(
      model,
      inits = init,
      updatefn = updatefn,
      method = "BFGS",
      control = list(maxit = MAXIT, trace = if (QUIET) 0 else 1)
    )
  }, error = function(e) e)
  
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(fit), model = model, optim = NULL))
  }
  
  conv <- fit$optim.out$convergence
  msg  <- fit$optim.out$message %||% NA_character_
  
  list(
    ok = isTRUE(conv == 0),
    reason = if (isTRUE(conv == 0)) NA_character_ else paste0("convergence=", conv, " ", msg),
    model = fit$model,
    optim = fit$optim.out
  )
}

forecast_ssm <- function(fitted_model, h) {
  tryCatch({
    # 1. Get the mean prediction
    pred_mean <- stats::predict(fitted_model, n.ahead = h, interval = "none")
    mu_h <- as.numeric(pred_mean[h])
    
    # 2. Helper to get specific intervals by calling predict() individually (like in c8)
    get_int <- function(level_pct) {
      lev <- level_pct / 100
      pr <- stats::predict(fitted_model, n.ahead = h, interval = "prediction", level = lev)
      c(lwr = as.numeric(pr[h, "lwr"]), upr = as.numeric(pr[h, "upr"]))
    }
    
    # 3. Extract 80% and 95% separately
    int80 <- if (80 %in% PI_LEVELS) get_int(80) else c(lwr = NA_real_, upr = NA_real_)
    int95 <- if (95 %in% PI_LEVELS) get_int(95) else c(lwr = NA_real_, upr = NA_real_)
    
    list(
      ok = TRUE, 
      reason = NA_character_,
      mean = mu_h,
      lower80 = unname(int80["lwr"]),
      upper80 = unname(int80["upr"]),
      lower95 = unname(int95["lwr"]),
      upper95 = unname(int95["upr"])
    )
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e),
         mean = NA_real_, lower80 = NA_real_, upper80 = NA_real_,
         lower95 = NA_real_, upper95 = NA_real_)
  })
}

align_to_fold_grid <- function(g, gamma_obs_g) {
  df_dates <- tibble(date = as.Date(g$dates))
  df_dates %>%
    left_join(gamma_obs_g, by = "date") %>%
    mutate(
      GEAR = as.character(g$GEAR),
      q = as.numeric(g$q),
      q_label = q_label(g$q)
    ) %>%
    arrange(date)
}

# -----------------------------
# Load inputs
# -----------------------------
gamma_raw <- readRDS(IN_GAMMA) %>%
  mutate(
    GEAR = as.character(GEAR),
    q = as.numeric(q),
    date = as.Date(date)
  ) %>%
  arrange(GEAR, q, date)

expect_cols(gamma_raw, c("GEAR", "q", "date"), "gamma")

TARGET_COL_EVAL <- pick_target_col(gamma_raw)
cat("Using target column for evaluation:", TARGET_COL_EVAL, "\n")

gamma_obs <- gamma_raw %>%
  transmute(GEAR, q, date, value_use = suppressWarnings(as.numeric(.data[[TARGET_COL_EVAL]]))) %>%
  collapse_dups_gamma(value_col = "value_use") %>%
  arrange(GEAR, q, date)

folds <- readRDS(IN_FOLDS)
stopifnot(is.list(folds), !is.null(folds$groups))

# Gamma groups from folds
gamma_groups <- folds$groups[vapply(folds$groups, function(g) identical(g$series, "gamma"), logical(1))]
if (length(gamma_groups) == 0) stop("No gamma groups found in backtest folds (series=='gamma').")

# Pre-build aligned series for each gid once
aligned_by_gid <- lapply(names(gamma_groups), function(gid) {
  g <- gamma_groups[[gid]]
  if (!isTRUE(g$ok)) return(NULL)
  
  obs_g <- gamma_obs %>% filter(GEAR == g$GEAR, q == g$q) %>% select(date, value)
  df_g <- align_to_fold_grid(g, obs_g)
  
  if (nrow(df_g) != length(g$dates)) {
    stop("Aligned grid length mismatch for ", gid,
         " (aligned rows=", nrow(df_g), ", folds dates=", length(g$dates), ").")
  }
  df_g
})
names(aligned_by_gid) <- names(gamma_groups)

# -----------------------------
# Build job list (one job = one fold forecast)
# -----------------------------
jobs <- list()
job_skips <- data.frame(reason = character(0), stringsAsFactors = FALSE)

for (gid in names(gamma_groups)) {
  g <- gamma_groups[[gid]]
  if (!isTRUE(g$ok)) next
  df_g <- aligned_by_gid[[gid]]
  if (is.null(df_g)) next
  
  n_g <- nrow(df_g)
  
  for (hname in names(g$folds)) {
    folds_h <- g$folds[[hname]]
    if (length(folds_h) == 0) next
    
    for (f in folds_h) {
      h <- f$horizon
      tr_idx <- f$train_idx
      te_idx <- f$test_idx
      
      if (max(tr_idx, te_idx) > n_g) {
        job_skips <- rbind(job_skips, data.frame(reason = paste0(gid, ": index_out_of_range"), stringsAsFactors = FALSE))
        next
      }
      
      y_train <- df_g$value[tr_idx]
      y_test  <- df_g$value[te_idx]
      
      if (is.na(y_test)) next
      if (sum(!is.na(y_train)) < MIN_TRAIN_NON_NA) next
      
      jobs[[length(jobs) + 1]] <- list(
        gid = gid,
        h = h,
        tr_idx = tr_idx,
        te_idx = te_idx
      )
    }
  }
}

cat("Jobs to run:", length(jobs), "\n")
if (nrow(job_skips) > 0) {
  job_skips_summary <- dplyr::as_tibble(job_skips) %>% count(reason, sort = TRUE)
  readr::write_csv(job_skips_summary, file.path(dir_diag_script, "job_skips_summary.csv"))
}
if (length(jobs) == 0) stop("No valid jobs were created; check folds and gamma series.")

# -----------------------------
# Worker function
# -----------------------------
run_one_job <- function(job, aligned_by_gid, XREG_COLS) {
  gid <- job$gid
  df_g <- aligned_by_gid[[gid]]
  
  h <- job$h
  tr_idx <- job$tr_idx
  te_idx <- job$te_idx
  
  y_train <- df_g$value[tr_idx]
  y_test  <- df_g$value[te_idx]
  
  fit_log <- tibble(
    group_id = gid,
    GEAR = unique(df_g$GEAR),
    q = unique(df_g$q),
    horizon = h,
    origin_date = df_g$date[max(tr_idx)],
    test_date = df_g$date[te_idx],
    ok = FALSE,
    reason = NA_character_
  )
  
  if (is.na(y_test) || sum(!is.na(y_train)) < MIN_TRAIN_NON_NA) {
    fit_log$reason <- "insufficient_train_or_missing_test"
    return(list(pred = NULL, log = fit_log))
  }
  
  xreg_train <- NULL
  if (length(XREG_COLS) > 0) {
    expect_cols(df_g, XREG_COLS, "aligned gamma (xreg)")
    xreg_train <- as.matrix(df_g[tr_idx, XREG_COLS, drop = FALSE])
  }
  
  ssm0 <- build_ssm(y_train, xreg = xreg_train)
  fit  <- fit_ssm_ml(ssm0)
  
  if (!isTRUE(fit$ok)) {
    fit_log$reason <- fit$reason %||% "fit_failed"
    return(list(pred = NULL, log = fit_log))
  }
  
  fc <- forecast_ssm(fit$model, h = h)
  if (!isTRUE(fc$ok)) {
    fit_log$reason <- paste0("forecast_failed: ", fc$reason %||% NA_character_)
    return(list(pred = NULL, log = fit_log))
  }
  
  fit_log$ok <- TRUE
  
  pred <- tibble(
    group_id = gid,
    GEAR = unique(df_g$GEAR),
    q = unique(df_g$q),
    q_label = unique(df_g$q_label),
    horizon = h,
    origin_date = df_g$date[max(tr_idx)],
    test_date = df_g$date[te_idx],
    model = "ssm",
    ok = TRUE,
    reason = NA_character_,
    y_obs = as.numeric(y_test),
    y_hat = fc$mean,
    lower80 = fc$lower80,
    upper80 = fc$upper80,
    lower95 = fc$lower95,
    upper95 = fc$upper95
  )
  
  list(pred = pred, log = fit_log)
}

# -----------------------------
# Run jobs
# -----------------------------
if (RUN_PARALLEL) {
  res_list <- future.apply::future_lapply(
    jobs,
    FUN = run_one_job,
    aligned_by_gid = aligned_by_gid,
    XREG_COLS = XREG_COLS,
    future.seed = TRUE
  )
} else {
  res_list <- lapply(
    jobs,
    run_one_job,
    aligned_by_gid = aligned_by_gid,
    XREG_COLS = XREG_COLS
  )
}

pred_tbl <- bind_rows(lapply(res_list, `[[`, "pred")) %>%
  mutate(
    error = y_obs - y_hat,
    abs_error = abs(error),
    covered80 = ifelse(is.na(lower80) | is.na(upper80), NA, y_obs >= lower80 & y_obs <= upper80),
    covered95 = ifelse(is.na(lower95) | is.na(upper95), NA, y_obs >= lower95 & y_obs <= upper95)
  )

fit_log_tbl <- bind_rows(lapply(res_list, `[[`, "log"))
readr::write_csv(fit_log_tbl, OUT_FITLOG)

saveRDS(
  list(
    created = Sys.time(),
    eval_target_col = TARGET_COL_EVAL,
    model_spec = "local level + seasonal(12) + observation noise (+ optional xreg)",
    xreg_cols = XREG_COLS,
    fit_log = fit_log_tbl,
    predictions = pred_tbl
  ),
  OUT_FITS_DERIVED
)
saveRDS(readRDS(OUT_FITS_DERIVED), OUT_FITS_MODEL)

# -----------------------------
# Metrics (SSM-only)
# -----------------------------
metrics_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(GEAR, q, q_label, horizon) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(error^2)),
    MAE  = mean(abs_error),
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(metrics_tbl, OUT_METR_DERIVED)
readr::write_csv(metrics_tbl, OUT_METR_MODEL)

# -----------------------------
# Diagnostics
# -----------------------------
if (has_ggplot2) {
  acf_groups <- pred_tbl %>%
    filter(ok, is.finite(error)) %>%
    group_by(GEAR, q, q_label, horizon) %>%
    group_split()
  
  for (d in acf_groups) {
    if (nrow(d) < 20) next
    d <- d %>% arrange(test_date)
    
    ac <- stats::acf(d$error, plot = FALSE, lag.max = ACF_MAX_LAG)
    acdf <- tibble(lag = as.numeric(ac$lag), acf = as.numeric(ac$acf))
    
    p <- ggplot2::ggplot(acdf, ggplot2::aes(x = lag, y = acf)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
      ggplot2::labs(
        x = "Lag",
        y = "ACF",
        title = "SSM forecast error ACF (across folds)",
        subtitle = paste0(unique(d$GEAR), " | ", unique(d$q_label), " | h=", unique(d$horizon))
      ) +
      ggplot2::theme_bw()
    
    fname <- sprintf("acf_ssm_errors__%s__q%s__h%s",
                     unique(d$GEAR),
                     format(unique(d$q), trim = TRUE),
                     unique(d$horizon))
    
    save_fig_pdf_tiff(p, file.path(dir_diag_fig, fname), width = 6, height = 4)
  }
  
  if (nrow(pred_tbl) > 0) {
    p_err <- ggplot2::ggplot(pred_tbl, ggplot2::aes(x = error)) +
      ggplot2::geom_histogram(bins = 40) +
      ggplot2::facet_grid(horizon ~ q_label) +
      ggplot2::labs(
        x = "Forecast error (y_obs - y_hat)",
        y = "Count",
        title = "SSM forecast error distributions (gamma)"
      ) +
      ggplot2::theme_bw()
    
    save_fig_pdf_tiff(p_err, file.path(dir_diag_fig, "ssm_error_histogram"), width = 12, height = 7)
  }
}

ljung_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(GEAR, q, q_label, horizon) %>%
  summarise(
    n = n(),
    ljungbox_p = if (dplyr::n() >= (LJUNGBOX_LAG + 5))
      stats::Box.test(error, lag = LJUNGBOX_LAG, type = "Ljung-Box")$p.value else NA_real_,
    .groups = "drop"
  )
readr::write_csv(ljung_tbl, file.path(dir_diag_script, "ssm_ljungbox_across_folds.csv"))

cat("\nFinished:", format(Sys.time()), "\n")
cat("Wrote derived:\n  -", OUT_FITS_DERIVED, "\n  -", OUT_METR_DERIVED, "\n")
cat("Wrote results:\n  -", OUT_FITS_MODEL, "\n  -", OUT_METR_MODEL, "\n")
cat("Diagnostics:\n  -", dir_diag_script, "\n")
cat("Fit log:\n  -", OUT_FITLOG, "\n")
sink()
message("✅ c5 completed (standalone; no c4p dependency). Log: ", log_file)