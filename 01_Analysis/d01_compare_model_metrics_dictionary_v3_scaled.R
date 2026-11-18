# =============================================================================
# Script:  d01_compare_model_metrics_dictionary_v3_scaled.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Brings the diversity forecasts of c07 (SARIMA family, damped trend), c08
#   (state-space), c09 (ARIMAX), c10 (state-space with covariates) and c12
#   (references) into one table for scoring in f01 (lib/join_forecasts.R).
#   One row per series, forecast date, lead and model; nothing is scored here.
#   The file name keeps "metrics_dictionary" from the earlier script name.
#
# Run after: b09, c02, c03, c07, c08, c09, c10, c12
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds                       (b09)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/seasonal_coefs_by_origin_scaled.rds (c02)
#   - 01_Analysis/Data/c3_forecast_design_backtesting/forecast_origins_diversity.rds (c03)
#   - the forecast and fit-information files of the model scripts (SOURCES):
#       01_Analysis/Data/c7_alpha_forecasting_multilevel_or_bay_replication_scaled/
#         forecasts_sarima_family_diversity.rds, fit_info_sarima_family_diversity.rds (c07)
#       01_Analysis/Data/c8_alpha_state_space_models_forecasts_parallel_scaled/
#         forecasts_ssm_diversity.rds, fit_info_ssm_diversity.rds (c08)
#       01_Analysis/Data/c9_alpha_arimax_from_env_index_parallel_scaled/
#         forecasts_arimax_diversity.rds, fit_info_arimax_diversity.rds (c09)
#       01_Analysis/Data/c10_alpha_state_space_with_env_covariates_parallel_scaled/
#         forecasts_ssm_xreg_diversity.rds, fit_info_ssm_xreg_diversity.rds (c10)
#       01_Analysis/Data/c12_naive_baseline_diversity/
#         forecasts_reference_diversity.rds, fit_info_reference_diversity.rds (c12)
#
# Outputs (01_Analysis/Data/d1_compare_model_metrics_dictionary_v3_scaled/):
#   - forecasts_joined_diversity.rds   keys, origin_i, origin_date,
#       in_common, lead, target_date, target_filled, scored, model, role,
#       source, ok, mean, lo80, hi80, lo95, hi95, y_obs, season_at_target,
#       value_obs
#   - fit_info_joined_diversity.rds    one row per series, forecast date and
#       model: ok, reason, source (part: which covariate set, for combined
#       models)
#   02_Results/03_Diagnostics/<stem>/d1_join_summary.csv  rows, fits, failures by
#       model
#   01_Analysis/Logs/<stem>/run_<timestamp>.txt
#   Earlier output files in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp>.
#
# Method notes:
#   - All models and references are kept as stored, plus the combined
#     covariate models (arimax, ssm_xreg: the lag 1 and 12 set at lead 1, the
#     lag 12 set at leads 2 to 12). Roles and partners are set in model_info().
#   - Failed fits stay as rows (ok = FALSE, NA forecasts); f01 decides how
#     they enter paired comparisons.
#   - Added columns: scored (FALSE for targets in April and May 2020),
#     season_at_target and value_obs (full scale), role, source.
#   - Checks that stop the script: every model has every forecast date and
#     lead; observed anomalies are identical across sources; value_obs equals
#     the input series; no duplicates; no NA forecast from a fit marked ok.
#   - RUN_MODE "test" reads the _TEST outputs of the model scripts for their
#     test series and the forecast dates they share, and writes under a _TEST
#     stem with _TEST file names.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

# ---- Bootstrap: find the project folder and build PATHS (same in every script)
# Looks for the marker file .project_root above, in order: the file being
# source()d (RStudio Source, source()), the file given to Rscript / R -f /
# R CMD BATCH (also under nohup), the file open in RStudio, the working folder.
local({
  a <- commandArgs(FALSE)
  here <- c(unlist(lapply(rev(sys.frames()), function(f) get0("ofile", envir = f, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", a, value = TRUE)),
            a[which(a == "-f") + 1L],
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
              tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
            getwd())
  here <- gsub("~+~", " ", here[!is.na(here) & nzchar(here)], fixed = TRUE)  # Rscript on Linux writes spaces as ~+~
  for (h in here) {
    d <- normalizePath(h, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) d <- dirname(d)
    repeat {
      if (file.exists(file.path(d, ".project_root"))) {
        assign(".PROJECT_ROOT", d, envir = globalenv())
        return(invisible())
      }
      if (dirname(d) == d) break
      d <- dirname(d)
    }
  }
  stop("Project folder not found: no .project_root above ", paste(unique(here), collapse = " or "))
})
source(file.path(.PROJECT_ROOT, "01_Analysis", "lib", "paths.R"))
PATHS <- get_paths(.PROJECT_ROOT)
setwd(PATHS$analysis_dir)
# -----------------------------------------------------------------------------
source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
source(file.path(PATHS$analysis_dir, "lib", "forecast_rules.R"))
source(file.path(PATHS$analysis_dir, "lib", "seasonal.R"))
source(file.path(PATHS$analysis_dir, "lib", "join_forecasts.R"))

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
RUN_MODE    <- "full"                                        # "full" or "test" (test runs use _TEST names)
TEST_SERIES <- c("bay|2|BS|0", "coastwide|NA|BT|2")   # used only in test mode
KEYS        <- c("scale", "BAY_ID", "GEAR", "q")

script_stem <- "d1_compare_model_metrics_dictionary_v3_scaled"   # folder stem kept from the earlier script name
SFX <- ""
if (RUN_MODE == "test") { script_stem <- paste0(script_stem, "_TEST"); SFX <- "_TEST" }

SOURCES <- list(   # name = c(script stem, forecast file, fit-information file)
  c07 = c("c7_alpha_forecasting_multilevel_or_bay_replication_scaled", "forecasts_sarima_family_diversity", "fit_info_sarima_family_diversity"),
  c08 = c("c8_alpha_state_space_models_forecasts_parallel_scaled",    "forecasts_ssm_diversity",           "fit_info_ssm_diversity"),
  c09 = c("c9_alpha_arimax_from_env_index_parallel_scaled",          "forecasts_arimax_diversity",        "fit_info_arimax_diversity"),
  c10 = c("c10_alpha_state_space_with_env_covariates_parallel_scaled", "forecasts_ssm_xreg_diversity",    "fit_info_ssm_xreg_diversity"),
  c12 = c("c12_naive_baseline_diversity",                            "forecasts_reference_diversity",     "fit_info_reference_diversity"))
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_gamma    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_coefs    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "seasonal_coefs_by_origin_scaled.rds")
in_origins  <- file.path(PATHS$data_dir, "c3_forecast_design_backtesting", "forecast_origins_diversity.rds")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 2) Folders, log, inputs
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "| run mode:", RUN_MODE, "\n")
src_files <- lapply(SOURCES, function(s) file.path(PATHS$data_dir, paste0(s[1], SFX), paste0(s[2:3], SFX, ".rds")))
stop_if_missing(c(in_settings, in_alpha, in_gamma, in_coefs, in_origins, unlist(src_files)),
                c(c("b09_forecast_settings.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 3), "c03_forecast_design_backtesting.R"), rep(paste0(names(SOURCES), " (model script)"), each = 2)))
settings <- load_settings(in_settings)
sid <- function(d) do.call(paste, c(as.list(d[, KEYS, drop = FALSE]), sep = "|"))

a <- readRDS(in_alpha); g <- readRDS(in_gamma)
series <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, q = as.integer(q), date, value),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, q = as.integer(q), date, value))
coefs   <- as.data.frame(readRDS(in_coefs))
origins <- as.data.frame(readRDS(in_origins))
sources <- lapply(src_files, function(f) list(forecasts = readRDS(f[1]), info = readRDS(f[2])))
for (s in names(SOURCES)) cat("Read", s, ":", nrow(sources[[s]]$forecasts), "forecast rows,", nrow(sources[[s]]$info), "fit records\n")
if (length(unique(sid(origins))) != 48L) stop("Origin table does not have 48 series", call. = FALSE)

if (RUN_MODE == "test") {
  # Test outputs use different origin subsets per script: keep the shared ones.
  common <- Reduce(intersect, lapply(sources, function(x) paste(sid(as.data.frame(x$info)), x$info$origin_i)))
  origins <- origins[paste(sid(origins), origins$origin_i) %in% common & sid(origins) %in% TEST_SERIES, ]
  sources <- lapply(sources, function(x) {
    list(forecasts = x$forecasts[paste(sid(as.data.frame(x$forecasts)), x$forecasts$origin_i) %in% common, ],
         info      = x$info[paste(sid(as.data.frame(x$info)), x$info$origin_i) %in% common, ])
  })
  cat("Test mode:", nrow(origins), "shared origins in", length(unique(sid(origins))), "series\n")
}

# -----------------------------------------------------------------------------
# 3) Join and check (lib/join_forecasts.R)
# -----------------------------------------------------------------------------
j <- join_forecasts(sources, KEYS, origins, coefs, as.data.frame(series), settings)
fc <- j$forecasts; info <- j$info

summary_tab <- fc %>% group_by(model, role, source) %>%
  summarise(rows = n(), max_lead = max(lead), na_forecasts = sum(is.na(mean)), rows_unscored = sum(!scored),
            .groups = "drop") %>%
  left_join(info %>% group_by(model) %>% summarise(fits = n(), failed_fits = sum(!ok), .groups = "drop"), by = "model") %>%
  arrange(match(role, c("main", "trend_test", "covariate_set", "reference")), model)

# -----------------------------------------------------------------------------
# 4) Write
# -----------------------------------------------------------------------------
out_fc   <- file.path(dirs$data, paste0("forecasts_joined_diversity", SFX, ".rds"))
out_info <- file.path(dirs$data, paste0("fit_info_joined_diversity", SFX, ".rds"))
out_sum  <- file.path(dirs$diag, paste0("d1_join_summary", SFX, ".csv"))
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE)
if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|pdf|tiff|png|txt)$", full.names = TRUE)
if (length(old_d)) archive_files(old_d)
saveRDS(fc, out_fc); saveRDS(info, out_info); readr::write_csv(summary_tab, out_sum, na = "")

cat("\nRows and fits by model:\n"); print(as.data.frame(summary_tab), row.names = FALSE)
cat("\nJoined rows:", nrow(fc), "| series:", length(unique(sid(fc))), "| models:", length(unique(fc$model)), "\n")
for (f in c(out_fc, out_info, out_sum)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("d01_compare_model_metrics_dictionary_v3_scaled.R completed\n")
