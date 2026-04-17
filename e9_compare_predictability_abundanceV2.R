# =============================================================================
# Script: e9_compare_predictability_abundanceV2.R
#
# Purpose:
# Compiles all abundance (CPUE) forecasting results (e5, e6, e7, e8) into a 
# single master metrics dictionary.
#
# NEW PUBLICATION-GRADE UPGRADES:
# - Calculates Error Dispersion (Standard Deviation and IQR of errors).
# - Calculates Diebold-Mariano (DM) Test p-values to prove if the skill 
#   improvements over the Seasonal Naive baseline are statistically significant.
# - Replaces deprecated cur_data() with pick(everything()).
#
# Outputs:
# - 01_Analysis/Data/e9_compare_predictability_abundance/metrics_dictionary_abundance.csv
# - 02_Results/04_Models/e9_compare_predictability_abundance/metrics_dictionary_abundance.csv
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
# 2) Packages
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forecast) # Required for dm.test
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) Directories
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e9_compare_predictability_abundance"
dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)

# Input files from previous steps
IN_E5 <- file.path(PATHS$analysis_dir, "Data", "e5_baseline_forecasts_abundance", "baseline_forecasts_abundance.rds")
IN_E6 <- file.path(PATHS$analysis_dir, "Data", "e6_state_space_forecasts_abundance", "ssm_forecasts_abundance.rds")
IN_E7 <- file.path(PATHS$analysis_dir, "Data", "e7_arimax_env_forecasts_abundance", "arimax_forecasts_abundance.rds")
IN_E8 <- file.path(PATHS$analysis_dir, "Data", "e8_ssm_xreg_forecasts_abundance", "ssm_xreg_forecasts_abundance.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Load and Consolidate Raw Predictions
# ─────────────────────────────────────────────────────────────────────────────
cat("Loading raw predictions from spatial model tiers...\n")

load_preds <- function(path) {
  if (file.exists(path)) {
    data_obj <- readRDS(path)
    
    # If the file is already a direct data frame, return it
    if (is.data.frame(data_obj)) {
      return(data_obj)
    } 
    # If the file is a list containing $predictions, return that
    else if (is.list(data_obj) && "predictions" %in% names(data_obj)) {
      return(data_obj$predictions)
    } 
    # If it's something else, throw a warning so it doesn't fail silently
    else {
      warning(paste("File found, but data format not recognized:", path))
      return(NULL)
    }
    
  } else {
    warning(paste("File missing, skipping:", path))
    return(NULL)
  }
}

pred_list <- list(
  load_preds(IN_E5),
  load_preds(IN_E6),
  load_preds(IN_E7),
  load_preds(IN_E8)
)

# Combine and remove duplicates (e.g., e7 recalculates snaive/sarima to calculate its own skill)
pred_all <- dplyr::bind_rows(pred_list) %>%
  distinct(series, BAY, GEAR, SPECIES_CODE, horizon, origin_i, model, .keep_all = TRUE)

cat("Consolidated", nrow(pred_all), "total forecasting events.\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Calculate Failure Rates using pick(everything())
# ─────────────────────────────────────────────────────────────────────────────
cat("Calculating model stability and failure rates...\n")

fail_rates <- pred_all %>%
  group_by(series, BAY, GEAR, SPECIES_CODE, horizon, model) %>%
  summarise(
    total_attempts = n(),
    frac_failed = if ("ok" %in% names(pick(everything()))) mean(!ok, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Isolate Baseline and Join for Paired Testing (Diebold-Mariano)
# ─────────────────────────────────────────────────────────────────────────────
cat("Extracting baseline cross-validation pairings for Diebold-Mariano tests...\n")

base_preds <- pred_all %>%
  filter(model == "seasonal_naive", ok) %>%
  select(series, BAY, GEAR, SPECIES_CODE, horizon, test_date, 
         error_base = error, abs_error_base = abs_error)

pred_joined <- pred_all %>%
  filter(ok, is.finite(error)) %>%
  left_join(base_preds, by = c("series", "BAY", "GEAR", "SPECIES_CODE", "horizon", "test_date"))

# --- REVISION: Save fold-level paired errors for bootstrap CIs and effect sizes ---
OUT_JOINED_RDS <- file.path(dir_data_script, "pred_joined_abundance.rds")
saveRDS(pred_joined, OUT_JOINED_RDS)
cat("Saved fold-level paired errors to:", OUT_JOINED_RDS, "\n")
# ---------------------------------------------------------------------------------

# ─────────────────────────────────────────────────────────────────────────────
# 7) Calculate Final Metrics Dictionary
# ─────────────────────────────────────────────────────────────────────────────
cat("Computing RMSE, MAE, Error Dispersion, and DM p-values...\n")

metrics_dict <- pred_joined %>%
  group_by(series, BAY, GEAR, SPECIES_CODE, horizon, model) %>%
  summarise(
    n_successful = n(),
    rmse = sqrt(mean(error^2)),
    mae  = mean(abs_error),
    
    # Error Dispersion (Publication requested)
    sd_error  = sd(error),
    iqr_error = IQR(error),
    
    coverage80 = mean(covered80, na.rm = TRUE),
    coverage95 = mean(covered95, na.rm = TRUE),
    
    # Baseline comparison averages
    rmse_base = sqrt(mean(error_base^2, na.rm = TRUE)),
    mae_base  = mean(abs_error_base, na.rm = TRUE),
    
    # Diebold-Mariano Test: statistic and p-value
    dm_p_value = {
      valid <- !is.na(error) & !is.na(error_base)
      e1 <- error[valid]; e2 <- error_base[valid]
      if (length(e1) > 3 && var(e1 - e2) > 1e-8) {
        tryCatch({
          forecast::dm.test(e1, e2, h = horizon[1], power = 2, alternative = "less")$p.value
        }, error = function(e) NA_real_)
      } else { NA_real_ }
    },
    dm_statistic = {
      valid <- !is.na(error) & !is.na(error_base)
      e1 <- error[valid]; e2 <- error_base[valid]
      if (length(e1) > 3 && var(e1 - e2) > 1e-8) {
        tryCatch({
          forecast::dm.test(e1, e2, h = horizon[1], power = 2, alternative = "less")$statistic
        }, error = function(e) NA_real_)
      } else { NA_real_ }
    },
    # Cohen's d on paired differences (baseline - advanced); positive = advanced is better
    cohens_d = {
      diffs <- error_base - error
      diffs <- diffs[!is.na(diffs)]
      if (length(diffs) > 1 && sd(diffs) > 0) mean(diffs) / sd(diffs) else NA_real_
    },
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 8) Final Merging and Skill Calculations
# ─────────────────────────────────────────────────────────────────────────────
final_dict <- metrics_dict %>%
  left_join(fail_rates, by = c("series", "BAY", "GEAR", "SPECIES_CODE", "horizon", "model")) %>%
  mutate(
    skill_rmse_ratio = ifelse(is.na(rmse_base) | rmse_base == 0, NA_real_, rmse / rmse_base),
    skill_mae_ratio  = ifelse(is.na(mae_base)  | mae_base  == 0, NA_real_, mae  / mae_base),
    skill_rmse_improve_pct = ifelse(is.na(skill_rmse_ratio), NA_real_, 100 * (1 - skill_rmse_ratio)),
    skill_mae_improve_pct  = ifelse(is.na(skill_mae_ratio),  NA_real_, 100 * (1 - skill_mae_ratio)),
    
    # Significance flag (alpha = 0.05)
    is_significantly_better = ifelse(dm_p_value < 0.05, TRUE, FALSE),
    
    metrics_version = "e9_v2",
    target_definition = "deseasonalized log(CPUE)",
    baseline_definition = "seasonal_naive"
  ) %>%
  select(
    series, BAY, GEAR, SPECIES_CODE, horizon, model, total_attempts, n_successful, frac_failed,
    rmse, mae, sd_error, iqr_error, coverage80, coverage95,
    rmse_base, mae_base, skill_rmse_improve_pct, skill_mae_improve_pct,
    dm_statistic, dm_p_value, is_significantly_better, cohens_d,
    everything()
  ) %>%
  arrange(series, BAY, GEAR, SPECIES_CODE, horizon, model)

# ─────────────────────────────────────────────────────────────────────────────
# 9) Save Outputs
# ─────────────────────────────────────────────────────────────────────────────
OUT_CSV_DERIVED <- file.path(dir_data_script, "metrics_dictionary_abundance.csv")
OUT_RDS_DERIVED <- file.path(dir_data_script, "metrics_dictionary_abundance.rds")
OUT_CSV_MODEL   <- file.path(dir_model_script, "metrics_dictionary_abundance.csv")
OUT_RDS_MODEL   <- file.path(dir_model_script, "metrics_dictionary_abundance.rds")

write_csv(final_dict, OUT_CSV_DERIVED)
saveRDS(final_dict, OUT_RDS_DERIVED)
write_csv(final_dict, OUT_CSV_MODEL)
saveRDS(final_dict, OUT_RDS_MODEL)

cat("\n✅ Abundance Predictability Dictionary Compiled Successfully!\n")
cat("💾 Saved to:", OUT_CSV_MODEL, "\n")