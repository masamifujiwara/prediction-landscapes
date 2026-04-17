# =============================================================================
# Script: c7_alpha_forecasting_parallel.R
#
# Purpose:
#   Baseline bay-replicated forecasting for α diversity (q=0,1,2) by BAY×GEAR×q.
#   COMPLETE monthly grid per series (BAY×GEAR×q) so h=12 = 12 calendar months.
#   Fits baseline models per fold: seasonal_naive + SARIMA (auto.arima seasonal=TRUE).
#
# Inputs:
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/alpha_q012_deseas.rds
#   - 01_Analysis/PATHS.RData
#
# Outputs:
#   Derived:
#     - 01_Analysis/Data/c7_alpha_forecasting_multilevel_or_bay_replication/alpha_forecasts_by_bay.rds
#     - 01_Analysis/Data/c7_alpha_forecasting_multilevel_or_bay_replication/alpha_skill_by_bay.csv
#   Results:
#     - 02_Results/04_Models/c7_alpha_forecasting_multilevel_or_bay_replication/alpha_forecasts_by_bay.rds
#     - 02_Results/04_Models/c7_alpha_forecasting_multilevel_or_bay_replication/alpha_skill_by_bay.csv
#   Diagnostics:
#     - 02_Results/03_Diagnostics/c7_alpha_forecasting_multilevel_or_bay_replication/Figures/*.(pdf|tiff)
#   Logs:
#     - 01_Analysis/Logs/c7_alpha_forecasting_multilevel_or_bay_replication/runlog_*.txt
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

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Dependencies
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "forecast", "future", "future.apply", "progressr", "parallel")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forecast)
  library(future)
  library(future.apply)
  library(progressr)
})
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

# -----------------------------
# Config
# -----------------------------
script_stem <- "c7_alpha_forecasting_multilevel_or_bay_replication"

HORIZONS <- c(1, 12)
TRAINING_RULE <- "expanding"
MIN_TRAIN_MONTHS <- 120
FIXED_TRAIN_MONTHS <- 120
ORIGIN_STEP <- 1
HOLDOUT_LAST_N <- 0

FREQ <- 12
ARIMA_STEPWISE <- TRUE
ARIMA_APPROX   <- FALSE
MIN_TRAIN_NONNA <- 24
MODELS <- c("seasonal_naive", "sarima")

# Enforce deseasonalized α everywhere (scientific intent)
TARGET_COL_EVAL_PRIMARY <- "value_deseas"
TARGET_COL_FALLBACKS <- c("residual_component")  # only if you really need a fallback

TARGET_SCALE <- "deseasonalized"
EVAL_TARGET_COL <- TARGET_COL_EVAL_PRIMARY

# Exclude Sabine Lake (Bay 1) going forward
DROP_BAYS <- c("1")
# -----------------------------
# Parallel plan
# -----------------------------
num_cores_phys <- parallel::detectCores()
plan_cores <- min(max(1, num_cores_phys - 2), 6)
future::plan(multisession, workers = plan_cores)
progressr::handlers(global = TRUE)

# -----------------------------
# Directories (PATHS conventions)
# -----------------------------
dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script  <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_diag_fig     <- file.path(dir_diag_script, "Figures")
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

invisible(lapply(
  list(dir_data_script, dir_logs_script, dir_diag_script, dir_diag_fig, dir_model_script),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

log_file <- file.path(dir_logs_script, paste0("runlog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
cat("Script:", script_stem, "\nStarted:", format(Sys.time()), "\nWorkers:", plan_cores, "\n\n")

# -----------------------------
# Helpers
# -----------------------------
pick_target_col <- function(df) {
  if (TARGET_COL_EVAL_PRIMARY %in% names(df)) return(TARGET_COL_EVAL_PRIMARY)
  stop(
    "Required deseasonalized target column '", TARGET_COL_EVAL_PRIMARY, "' not found in alpha input.\n",
    "Fix: run c2_deseasonalize_alpha_gamma and use its output alpha_q012_deseas.rds."
  )
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

collapse_dups_alpha <- function(df, value_col = "value") {
  df %>%
    group_by(BAY, GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      n_samples = if ("n_samples" %in% names(df)) sum(n_samples, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}

complete_monthly_grid_global <- function(df_obs, global_start, global_end) {
  tibble(date = seq.Date(global_start, global_end, by = "month")) %>%
    left_join(df_obs, by = "date") %>%
    arrange(date)
}

build_folds_one_series <- function(n, horizons,
                                   training_rule = "expanding",
                                   min_train_months = 120,
                                   fixed_train_months = 120,
                                   origin_step = 1,
                                   holdout_last_n = 0) {
  n_eff <- if (holdout_last_n > 0) max(0, n - holdout_last_n) else n
  if (n_eff < (min_train_months + max(horizons) + 1)) return(list(ok = FALSE, reason = "Too short", folds = list()))
  
  folds <- list()
  for (h in horizons) {
    origins <- seq(from = min_train_months, to = (n_eff - h), by = origin_step)
    folds[[paste0("h", h)]] <- lapply(origins, function(t0) {
      train_idx <- if (training_rule == "expanding") 1:t0 else max(1, t0 - fixed_train_months + 1):t0
      list(horizon = h, origin_i = t0, test_i = t0 + h, train_idx = train_idx, test_idx = t0 + h)
    })
  }
  list(ok = TRUE, folds = folds)
}

forecast_one_fold <- function(y_train_ts, h, model) {
  tryCatch({
    if (model == "seasonal_naive") {
      fc <- forecast::forecast(forecast::snaive(y_train_ts, h = h), level = c(80, 95))
    } else {
      fit <- suppressWarnings(
        forecast::auto.arima(
          y_train_ts,
          seasonal = TRUE,
          stepwise = ARIMA_STEPWISE,
          approximation = ARIMA_APPROX
        )
      )
      fc <- forecast::forecast(fit, h = h, level = c(80, 95))
    }
    
    lo <- fc$lower; hi <- fc$upper
    col_levels <- suppressWarnings(as.numeric(gsub("%", "", colnames(lo))))
    idx80 <- which(col_levels == 80); idx95 <- which(col_levels == 95)
    
    list(
      ok = TRUE,
      mean = as.numeric(fc$mean)[h],
      lower80 = if (length(idx80) == 1) as.numeric(lo[h, idx80]) else NA_real_,
      upper80 = if (length(idx80) == 1) as.numeric(hi[h, idx80]) else NA_real_,
      lower95 = if (length(idx95) == 1) as.numeric(lo[h, idx95]) else NA_real_,
      upper95 = if (length(idx95) == 1) as.numeric(hi[h, idx95]) else NA_real_,
      reason = NA_character_
    )
  }, error = function(e) {
    list(
      ok = FALSE,
      mean = NA_real_,
      lower80 = NA_real_, upper80 = NA_real_,
      lower95 = NA_real_, upper95 = NA_real_,
      reason = conditionMessage(e)
    )
  })
}

save_fig_pdf_tiff <- function(plot, base_path_no_ext, width, height) {
  if (!has_ggplot2) return(invisible(NULL))
  ggplot2::ggsave(paste0(base_path_no_ext, ".pdf"),  plot = plot, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(base_path_no_ext, ".tiff"), plot = plot, width = width, height = height, units = "in", dpi = 300, compression = "lzw")
}

# -----------------------------
# Load alpha (b6 output)
# -----------------------------
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds")
stopifnot(file.exists(IN_ALPHA))

alpha_obj <- readRDS(IN_ALPHA)
alpha <- if (is.data.frame(alpha_obj)) alpha_obj else alpha_obj[[which(vapply(alpha_obj, is.data.frame, logical(1)))[1]]]
stopifnot(is.data.frame(alpha))

# Standardize key columns
if (!("GEAR" %in% names(alpha))) alpha <- rename(alpha, GEAR = !!intersect(names(alpha), c("gear", "Gear"))[1])
if (!("YEAR" %in% names(alpha))) alpha <- rename(alpha, YEAR = !!intersect(names(alpha), c("year", "Year"))[1])
if (!("MONTH" %in% names(alpha))) alpha <- rename(alpha, MONTH = !!intersect(names(alpha), c("month", "Month"))[1])

# Pivot wide q0/q1/q2 -> long if needed
if (!("q" %in% names(alpha)) && all(c("q0", "q1", "q2") %in% names(alpha))) {
  alpha <- alpha %>%
    pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "value") %>%
    mutate(q = as.numeric(gsub("^q", "", q)))
}

BAY_COL <- intersect(c("BAY_ID", "BAY", "bay", "Bay", "MAJOR_AREA"), names(alpha))[1]
if (is.na(BAY_COL)) stop("Could not find a BAY column in alpha input.")

TARGET_COL_EVAL <- pick_target_col(alpha)
EVAL_TARGET_COL <- TARGET_COL_EVAL
cat("Using target column:", TARGET_COL_EVAL, "\n\n")

alpha_bay_obs <- alpha %>%
  mutate(
    BAY = as.character(.data[[BAY_COL]]),
    GEAR = as.character(GEAR),
    q = as.character(q),  # MATCH c8/c9
    date = as.Date(sprintf("%04d-%02d-01", as.integer(YEAR), as.integer(MONTH))),
    value = suppressWarnings(as.numeric(.data[[TARGET_COL_EVAL]])),
    q_label = q_label(as.numeric(q))
  ) %>%
  filter(!(BAY %in% DROP_BAYS)) %>%
  select(BAY, GEAR, q, q_label, date, value, any_of("n_samples"))

alpha_bay_obs <- collapse_dups_alpha(alpha_bay_obs, value_col = "value") %>%
  left_join(alpha_bay_obs %>% distinct(BAY, GEAR, q, q_label), by = c("BAY", "GEAR", "q"))

GLOBAL_START <- min(alpha_bay_obs$date, na.rm = TRUE)
GLOBAL_END   <- max(alpha_bay_obs$date, na.rm = TRUE)

# -----------------------------
# Split into groups and forecast
# -----------------------------
groups <- alpha_bay_obs %>% group_by(BAY, GEAR, q, q_label) %>% group_split()
cat("Processing", length(groups), "series groups on", plan_cores, "workers...\n\n")

with_progress({
  p <- progressor(steps = length(groups))
  
  results_list <- future_lapply(groups, function(df_obs) {
    p()
    
    bay <- unique(df_obs$BAY); gear <- unique(df_obs$GEAR); qq <- unique(df_obs$q); qlab <- unique(df_obs$q_label)
    
    df_g <- complete_monthly_grid_global(df_obs, GLOBAL_START, GLOBAL_END)
    folds <- build_folds_one_series(
      n = nrow(df_g),
      horizons = HORIZONS,
      training_rule = TRAINING_RULE,
      min_train_months = MIN_TRAIN_MONTHS,
      fixed_train_months = FIXED_TRAIN_MONTHS,
      origin_step = ORIGIN_STEP,
      holdout_last_n = HOLDOUT_LAST_N
    )
    
    local_preds <- list()
    
    if (isTRUE(folds$ok)) {
      for (hname in names(folds$folds)) {
        for (f in folds$folds[[hname]]) {
          
          y_test <- df_g$value[f$test_i]
          n_train_non_na <- sum(!is.na(df_g$value[f$train_idx]))
          can_fit <- (n_train_non_na >= MIN_TRAIN_NONNA)
          
          y_ts <- NULL
          if (can_fit) {
            y_ts <- ts(forecast::na.interp(df_g$value[f$train_idx]), frequency = FREQ)
          }
          
          for (m in MODELS) {
            # Always emit a row for this key
            ok_row <- TRUE
            reason_row <- NA_character_
            
            if (!can_fit) {
              ok_row <- FALSE
              reason_row <- paste0("insufficient_training_nonNA<", MIN_TRAIN_NONNA)
            }
            
            if (is.na(y_test) || !is.finite(y_test)) {
              ok_row <- FALSE
              reason_row <- "missing_target_month"
            }
            
            fc <- list(ok = FALSE, mean = NA_real_,
                       lower80 = NA_real_, upper80 = NA_real_,
                       lower95 = NA_real_, upper95 = NA_real_,
                       reason = reason_row)
            
            # Fit only if train ok AND target exists (keeps compute down; key row still written)
            if (can_fit && !(is.na(y_test) || !is.finite(y_test))) {
              fc <- forecast_one_fold(y_ts, f$horizon, m)
              if (!isTRUE(fc$ok)) {
                ok_row <- FALSE
                reason_row <- fc$reason %||% "model_fit_failed"
              }
            }
            
            local_preds[[length(local_preds) + 1]] <- tibble(
              BAY = bay,
              GEAR = gear,
              q = qq,
              q_label = qlab,
              model = m,
              fold_rule = TRAINING_RULE,
              horizon = f$horizon,
              origin_i = f$origin_i,
              test_i = f$test_i,
              origin_date = df_g$date[f$origin_i],
              test_date = df_g$date[f$test_i],
              y_obs = y_test,
              y_hat = fc$mean,
              lower80 = fc$lower80, upper80 = fc$upper80,
              lower95 = fc$lower95, upper95 = fc$upper95,
              ok = ok_row,
              reason = reason_row
            )
          }
        }
      }
    }
    
    list(
      preds = bind_rows(local_preds),
      info = tibble(BAY = bay, GEAR = gear, q = qq, q_label = qlab,
                    ok_folds = folds$ok, n_months_grid = nrow(df_g))
    )
  }, future.seed = TRUE)
})

pred_tbl <- bind_rows(lapply(results_list, `[[`, "preds")) %>%
  mutate(
    error = y_obs - y_hat,
    abs_error = abs(error),
    covered80 = is.finite(lower80) & is.finite(upper80) & is.finite(y_obs) & (y_obs >= lower80) & (y_obs <= upper80),
    covered95 = is.finite(lower95) & is.finite(upper95) & is.finite(y_obs) & (y_obs >= lower95) & (y_obs <= upper95)
  )

# Save object wrapper (metadata + predictions)
save_obj <- list(
  predictions     = pred_tbl,
  eval_target_col = EVAL_TARGET_COL,
  target_scale    = TARGET_SCALE,
  script_stem     = script_stem,
  created         = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

out_rds_data    <- file.path(dir_data_script,  "alpha_forecasts_by_bay.rds")
out_rds_results <- file.path(dir_model_script, "alpha_forecasts_by_bay.rds")
saveRDS(save_obj, out_rds_data)
saveRDS(save_obj, out_rds_results)

# Minimal per-bay metrics table (skill vs snaive computed here for diagnostics only)
metrics_by_bay <- pred_tbl %>%
  filter(ok, is.finite(y_obs), is.finite(y_hat)) %>%
  group_by(BAY, GEAR, q, q_label, horizon, model) %>%
  summarise(
    n = dplyr::n(),
    RMSE = sqrt(mean((y_obs - y_hat)^2)),
    MAE  = mean(abs(y_obs - y_hat)),
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    .groups = "drop"
  )

baseline <- metrics_by_bay %>%
  select(BAY, GEAR, q, horizon, model, RMSE) %>%
  tidyr::pivot_wider(names_from = model, values_from = RMSE)

skill_by_bay <- metrics_by_bay %>%
  left_join(baseline, by = c("BAY", "GEAR", "q", "horizon")) %>%
  mutate(RMSE_skill_vs_snaive = 1 - (RMSE / seasonal_naive))

out_csv_data    <- file.path(dir_data_script,  "alpha_skill_by_bay.csv")
out_csv_results <- file.path(dir_model_script, "alpha_skill_by_bay.csv")
readr::write_csv(skill_by_bay, out_csv_data)
readr::write_csv(skill_by_bay, out_csv_results)

# Diagnostics plot
if (has_ggplot2 && nrow(skill_by_bay) > 0) {
  p_box <- ggplot2::ggplot(dplyr::filter(skill_by_bay, model == "sarima"),
                           ggplot2::aes(x = GEAR, y = RMSE_skill_vs_snaive)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(q_label ~ horizon) +
    ggplot2::labs(title = "Alpha baseline skill distributions across bays") +
    ggplot2::theme_bw()
  
  save_fig_pdf_tiff(p_box, file.path(dir_diag_fig, "bs_vs_bt_skill_distributions__rmse"), 11, 7)
}

cat("\nWrote:\n  -", out_rds_data, "\n  -", out_rds_results, "\n  -", out_csv_data, "\n  -", out_csv_results, "\n")
cat("\nDone.\n")
sink()
future::plan(sequential)