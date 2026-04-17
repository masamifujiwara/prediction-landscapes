# =============================================================================
# Script: c4p_baseline_forecasts_gamma_parallel.R
#
# Purpose:
#   Baseline forecasts for γ diversity (per gear × q) using rolling-origin folds
#   from c3. Models: seasonal_naive + SARIMA (auto.arima).
#
# Inputs (authoritative; per project conventions):
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/gamma_q012_deseas.rds
#   - 01_Analysis/Data/c3_forecast_design_backtesting/backtest_folds.rds
#
# Outputs:
#   Derived data (downstream-friendly):
#     - 01_Analysis/Data/c4p_baseline_forecasts_gamma_parallel/baseline_forecasts_gamma.rds
#     - 01_Analysis/Data/c4p_baseline_forecasts_gamma_parallel/baseline_metrics_gamma.csv
#   Results (reporting):
#     - 02_Results/04_Models/c4p_baseline_forecasts_gamma_parallel/baseline_metrics_gamma.csv
#     - 02_Results/04_Models/c4p_baseline_forecasts_gamma_parallel/baseline_forecasts_gamma.rds
#   Diagnostics (CSV + plots):
#     - 02_Results/03_Diagnostics/c4p_baseline_forecasts_gamma_parallel/*.csv
#     - 02_Results/03_Diagnostics/c4p_baseline_forecasts_gamma_parallel/Figures/*.(pdf|tiff)  (if ggplot2)
#   Logs:
#     - 01_Analysis/Logs/c4p_baseline_forecasts_gamma_parallel/runlog_*.txt
#
# Notes:
#   1) Uses fold group's date grid (g$dates) as master grid to avoid mismatch.
#   2) Uses target column preference:
#        primary: value
#        fallbacks: value_deseas, residual_component
#   3) Collapses any duplicate (GEAR,q,date) via median as safety.
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
# Minimal dependencies (required)
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "forecast")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forecast)
})

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_future  <- requireNamespace("future", quietly = TRUE)
has_future_apply <- requireNamespace("future.apply", quietly = TRUE)

# -----------------------------
# Script stem + directories
# -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c4p_baseline_forecasts_gamma_parallel.R"))
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c4p_baseline_forecasts_gamma_parallel"

dir_data_script <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir_diag_script <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_diag_fig    <- file.path(dir_diag_script, "Figures")
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
OUT_METR_DERIVED <- file.path(dir_data_script, "baseline_metrics_gamma.csv")
OUT_FCST_DERIVED <- file.path(dir_data_script, "baseline_forecasts_gamma.rds")

OUT_METR_MODEL <- file.path(dir_model_script, "baseline_metrics_gamma.csv")
OUT_FCST_MODEL <- file.path(dir_model_script, "baseline_forecasts_gamma.rds")

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
# Model/config
# -----------------------------
TARGET_COL_EVAL_PRIMARY <- "value"
TARGET_COL_FALLBACKS <- c("value_deseas", "residual_component")
MODELS <- c("seasonal_naive", "sarima")

FREQ <- 12
PI_LEVELS <- c(80, 95)

ACF_MAX_LAG <- 36
LJUNGBOX_LAG <- 12

ARIMA_STEPWISE <- TRUE
ARIMA_APPROX   <- FALSE

MIN_TRAIN_NON_NA <- 24

# -----------------------------
# Logging
# -----------------------------
log_file <- file.path(
  dir_logs_script,
  paste0("runlog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)
sink(log_file, split = TRUE)
cat("Script:", script_stem, "\n")
cat("Started:", format(Sys.time()), "\n")
cat("Inputs:\n  -", IN_GAMMA, "\n  -", IN_FOLDS, "\n")
cat("Parallel:", RUN_PARALLEL, " | Workers:", if (RUN_PARALLEL) N_WORKERS else 1, "\n\n")

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

pick_target_col <- function(df) {
  if (TARGET_COL_EVAL_PRIMARY %in% names(df)) return(TARGET_COL_EVAL_PRIMARY)
  fb <- TARGET_COL_FALLBACKS[TARGET_COL_FALLBACKS %in% names(df)]
  if (length(fb) >= 1) return(fb[1])
  stop("No usable target column found. Looked for: ",
       paste(c(TARGET_COL_EVAL_PRIMARY, TARGET_COL_FALLBACKS), collapse = ", "))
}

pit_from_pi <- function(y, mean, lower80, upper80) {
  z <- qnorm(0.90)
  sd_hat <- (upper80 - lower80) / (2 * z)
  ifelse(is.finite(sd_hat) & sd_hat > 0, pnorm(y, mean = mean, sd = sd_hat), NA_real_)
}

forecast_one_fold <- function(y_train_ts, h, model) {
  tryCatch({
    if (model == "seasonal_naive") {
      fc <- forecast::forecast(forecast::snaive(y_train_ts, h = h), level = PI_LEVELS)
    } else if (model == "sarima") {
      fit <- suppressWarnings(
        forecast::auto.arima(
          y_train_ts,
          seasonal = TRUE,
          stepwise = ARIMA_STEPWISE,
          approximation = ARIMA_APPROX
        )
      )
      fc <- forecast::forecast(fit, h = h, level = PI_LEVELS)
    } else {
      stop("Unknown model: ", model)
    }
    
    m <- as.numeric(fc$mean)[h]
    lo <- fc$lower
    hi <- fc$upper
    col_levels <- suppressWarnings(as.numeric(gsub("%", "", colnames(lo))))
    idx80 <- which(col_levels == 80)
    idx95 <- which(col_levels == 95)
    
    list(
      ok = TRUE,
      reason = NA_character_,
      mean = m,
      lower80 = if (length(idx80) == 1) as.numeric(lo[h, idx80]) else NA_real_,
      upper80 = if (length(idx80) == 1) as.numeric(hi[h, idx80]) else NA_real_,
      lower95 = if (length(idx95) == 1) as.numeric(lo[h, idx95]) else NA_real_,
      upper95 = if (length(idx95) == 1) as.numeric(hi[h, idx95]) else NA_real_
    )
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e),
         mean = NA_real_, lower80 = NA_real_, upper80 = NA_real_,
         lower95 = NA_real_, upper95 = NA_real_)
  })
}

collapse_dups_gamma <- function(gamma_df, value_col) {
  gamma_df %>%
    group_by(GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

get_gamma_group_df_aligned <- function(g, gamma_obs_g) {
  # master grid = fold grid
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

save_fig_pdf_tiff <- function(plot, base_path_no_ext, width, height) {
  if (!has_ggplot2) return(invisible(NULL))
  ggplot2::ggsave(paste0(base_path_no_ext, ".pdf"),  plot = plot, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(base_path_no_ext, ".tiff"), plot = plot, width = width, height = height, units = "in", dpi = 300,
                  compression = "lzw")
}

# -----------------------------
# Load inputs
# -----------------------------
gamma_raw <- readRDS(IN_GAMMA)
folds <- readRDS(IN_FOLDS)

expect_cols(gamma_raw, c("GEAR", "q", "date"), "gamma")
stopifnot(is.list(folds), !is.null(folds$groups))

gamma_raw <- gamma_raw %>%
  mutate(
    GEAR = as.character(GEAR),
    q = as.numeric(q),
    date = as.Date(date)
  )

TARGET_COL_EVAL <- pick_target_col(gamma_raw)
cat("Using target column for evaluation:", TARGET_COL_EVAL, "\n")

gamma_obs <- gamma_raw %>%
  transmute(GEAR, q, date, value_use = suppressWarnings(as.numeric(.data[[TARGET_COL_EVAL]]))) %>%
  arrange(GEAR, q, date) %>%
  { collapse_dups_gamma(., value_col = "value_use") } %>%
  arrange(GEAR, q, date)

# -----------------------------
# Filter folds to gamma groups only
# -----------------------------
gamma_groups <- folds$groups[vapply(folds$groups, function(g) identical(g$series, "gamma"), logical(1))]
if (length(gamma_groups) == 0) stop("No gamma groups found in backtest_folds.rds (series == 'gamma').")

# Pre-build aligned group tables once
aligned_by_gid <- lapply(names(gamma_groups), function(gid) {
  g <- gamma_groups[[gid]]
  if (!isTRUE(g$ok)) return(NULL)
  
  df_obs_g <- gamma_obs %>% filter(GEAR == g$GEAR, q == g$q) %>% select(date, value)
  df_g <- get_gamma_group_df_aligned(g, df_obs_g)
  
  if (nrow(df_g) != length(g$dates)) {
    stop("Aligned grid length mismatch for ", gid,
         " (aligned rows=", nrow(df_g), ", folds dates=", length(g$dates), ").")
  }
  df_g
})
names(aligned_by_gid) <- names(gamma_groups)

# -----------------------------
# Build job list (fold × model)
# -----------------------------
jobs <- list()
job_skips <- data.frame(reason = character(0), stringsAsFactors = FALSE)

for (gid in names(gamma_groups)) {
  g <- gamma_groups[[gid]]
  
  if (!isTRUE(g$ok)) {
    job_skips <- rbind(job_skips, data.frame(reason = paste0(gid, ": group ok=FALSE"), stringsAsFactors = FALSE))
    next
  }
  
  df_g <- aligned_by_gid[[gid]]
  if (is.null(df_g)) {
    job_skips <- rbind(job_skips, data.frame(reason = paste0(gid, ": aligned df missing"), stringsAsFactors = FALSE))
    next
  }
  
  n_g <- nrow(df_g)
  
  for (hname in names(g$folds)) {
    folds_h <- g$folds[[hname]]
    if (length(folds_h) == 0) next
    
    for (f in folds_h) {
      h <- f$horizon
      tr_idx <- f$train_idx
      te_idx <- f$test_idx
      
      if (max(tr_idx, te_idx) > n_g) {
        job_skips <- rbind(job_skips, data.frame(reason = paste0(gid, ": index out of range"), stringsAsFactors = FALSE))
        next
      }
      
      y_train <- df_g$value[tr_idx]
      y_test  <- df_g$value[te_idx]
      
      if (is.na(y_test)) next
      if (sum(!is.na(y_train)) < MIN_TRAIN_NON_NA) next
      
      for (m in MODELS) {
        jobs[[length(jobs) + 1]] <- list(
          gid = gid,
          h = h,
          tr_idx = tr_idx,
          te_idx = te_idx,
          model = m
        )
      }
    }
  }
}

cat("Jobs to run (fold × model):", length(jobs), "\n")
if (nrow(job_skips) > 0) {
  job_skips_summary <- job_skips %>% count(reason, sort = TRUE)
  readr::write_csv(job_skips_summary, file.path(dir_diag_script, "job_skips_summary.csv"))
  cat("Wrote job skip summary CSV.\n")
}
if (length(jobs) == 0) stop("No valid jobs were created; check folds and gamma series.")

# -----------------------------
# Run jobs
# -----------------------------
run_one_job <- function(job, aligned_by_gid) {
  gid <- job$gid
  df_g <- aligned_by_gid[[gid]]
  
  h <- job$h
  tr_idx <- job$tr_idx
  te_idx <- job$te_idx
  m <- job$model
  
  y_train <- df_g$value[tr_idx]
  y_test  <- df_g$value[te_idx]
  
  # Fill training NAs for ARIMA/SNAIVE fitting
  y_train_fill <- forecast::na.interp(y_train)
  y_train_ts <- ts(y_train_fill, frequency = FREQ)
  
  fc <- forecast_one_fold(y_train_ts = y_train_ts, h = h, model = m)
  
  tibble(
    group_id = gid,
    GEAR = unique(df_g$GEAR),
    q = unique(df_g$q),
    q_label = unique(df_g$q_label),
    horizon = h,
    origin_date = df_g$date[max(tr_idx)],
    test_date = df_g$date[te_idx],
    model = m,
    ok = fc$ok,
    reason = fc$reason,
    y_obs = as.numeric(y_test),
    y_hat = fc$mean,
    lower80 = fc$lower80,
    upper80 = fc$upper80,
    lower95 = fc$lower95,
    upper95 = fc$upper95
  )
}

if (RUN_PARALLEL) {
  res_list <- future.apply::future_lapply(
    jobs,
    FUN = run_one_job,
    aligned_by_gid = aligned_by_gid,
    future.seed = TRUE
  )
} else {
  res_list <- lapply(
    jobs,
    run_one_job,
    aligned_by_gid = aligned_by_gid
  )
}

pred_tbl <- bind_rows(res_list) %>%
  mutate(
    error = y_obs - y_hat,
    abs_error = abs(error),
    covered80 = ifelse(is.na(lower80) | is.na(upper80), NA, y_obs >= lower80 & y_obs <= upper80),
    covered95 = ifelse(is.na(lower95) | is.na(upper95), NA, y_obs >= lower95 & y_obs <= upper95),
    pit = pit_from_pi(y_obs, y_hat, lower80, upper80)
  )

# Save forecasts (derived + results)
saveRDS(
  list(
    created = Sys.time(),
    eval_target_col = TARGET_COL_EVAL,
    models = MODELS,
    PI_LEVELS = PI_LEVELS,
    predictions = pred_tbl
  ),
  OUT_FCST_DERIVED
)
saveRDS(readRDS(OUT_FCST_DERIVED), OUT_FCST_MODEL)

# -----------------------------
# Metrics + skill vs seasonal-naive
# -----------------------------
metrics_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(GEAR, q, q_label, horizon, model) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(error^2)),
    MAE  = mean(abs_error),
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    pit_mean = mean(pit, na.rm = TRUE),
    pit_var  = stats::var(pit, na.rm = TRUE),
    .groups = "drop"
  )

baseline_tbl <- metrics_tbl %>%
  select(GEAR, q, horizon, model, RMSE, MAE) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(RMSE, MAE))

metrics_tbl2 <- metrics_tbl %>%
  left_join(baseline_tbl, by = c("GEAR", "q", "horizon")) %>%
  mutate(
    RMSE_skill_vs_snaive = 1 - (RMSE / RMSE_seasonal_naive),
    MAE_skill_vs_snaive  = 1 - (MAE  / MAE_seasonal_naive)
  )

readr::write_csv(metrics_tbl2, OUT_METR_DERIVED)
readr::write_csv(metrics_tbl2, OUT_METR_MODEL)

# -----------------------------
# Diagnostics
# -----------------------------
# Ljung-Box across folds
ljung_tbl <- pred_tbl %>%
  filter(ok, is.finite(error)) %>%
  group_by(GEAR, q, q_label, horizon, model) %>%
  summarise(
    n = n(),
    ljungbox_p = if (dplyr::n() >= (LJUNGBOX_LAG + 5))
      stats::Box.test(error, lag = LJUNGBOX_LAG, type = "Ljung-Box")$p.value
    else NA_real_,
    .groups = "drop"
  )
readr::write_csv(ljung_tbl, file.path(dir_diag_script, "ljungbox_across_folds.csv"))

# Error ACF plots
if (has_ggplot2) {
  acf_groups <- pred_tbl %>%
    filter(ok, is.finite(error)) %>%
    group_by(GEAR, q, q_label, horizon, model) %>%
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
        title = "Forecast error ACF (across folds)",
        subtitle = paste0(unique(d$GEAR), " | ", unique(d$q_label),
                          " | h=", unique(d$horizon), " | ", unique(d$model))
      ) +
      ggplot2::theme_bw()
    
    fname <- sprintf("acf_errors__%s__q%s__h%s__%s",
                     unique(d$GEAR),
                     format(unique(d$q), trim = TRUE),
                     unique(d$horizon),
                     unique(d$model))
    
    save_fig_pdf_tiff(p, file.path(dir_diag_fig, fname), width = 6, height = 4)
  }
  
  # PIT histogram
  pit_plot_df <- pred_tbl %>% filter(ok, is.finite(pit))
  if (nrow(pit_plot_df) > 0) {
    p_pit <- ggplot2::ggplot(pit_plot_df, ggplot2::aes(x = pit)) +
      ggplot2::geom_histogram(bins = 20) +
      ggplot2::facet_grid(horizon ~ model) +
      ggplot2::labs(
        x = "PIT (Normal approx from 80% PI)",
        y = "Count",
        title = "PIT distribution across all gamma folds"
      ) +
      ggplot2::theme_bw()
    
    save_fig_pdf_tiff(p_pit, file.path(dir_diag_fig, "pit_histogram__all_gamma"), width = 10, height = 6)
  }
  
  # Coverage plot
  cov_plot_df <- metrics_tbl %>% mutate(label = paste(GEAR, q_label, sep = " | "))
  if (nrow(cov_plot_df) > 0) {
    p_cov <- ggplot2::ggplot(cov_plot_df, ggplot2::aes(x = label, y = coverage80)) +
      ggplot2::geom_point() +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(horizon ~ model, scales = "free_y") +
      ggplot2::labs(
        x = NULL,
        y = "Empirical coverage (80%)",
        title = "Prediction interval coverage by group"
      ) +
      ggplot2::theme_bw()
    
    save_fig_pdf_tiff(p_cov, file.path(dir_diag_fig, "coverage80_by_group"), width = 11, height = 8)
  }
}

cat("\nFinished:", format(Sys.time()), "\n")
cat("Wrote derived:\n  -", OUT_FCST_DERIVED, "\n  -", OUT_METR_DERIVED, "\n")
cat("Wrote results:\n  -", OUT_FCST_MODEL, "\n  -", OUT_METR_MODEL, "\n")
cat("Diagnostics:\n  -", dir_diag_script, "\n")
sink()
message("✅ c4p completed. Log: ", log_file)