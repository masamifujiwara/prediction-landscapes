# =============================================================================
# Script:  e09_compare_predictability_abundance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Brings the abundance forecasts of e05a (SARIMA family, damped trend),
#   e05b (references), e06 (state-space), e07 (ARIMAX) and e08 (state-space
#   with covariates) into one table for scoring in f01
#   (lib/join_forecasts.R). One row per series, forecast date, lead time and
#   model; nothing is scored here.
#
# Run after: e05a, e05b, e06, e07, e08, e03, e04, lib/join_forecasts.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds                        (b09)
#   01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds     (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds     (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/seasonal_coefs_by_origin_abundance.rds (e03)
#   01_Analysis/Data/e4_forecast_design_abundance/forecast_origins_abundance.rds (e04)
#   01_Analysis/Data/e5_baseline_forecasts_abundance/
#     forecasts_sarima_family_abundance.rds, fit_info_sarima_family_abundance.rds (e05a)
#   01_Analysis/Data/e5b_naive_baseline_abundance/
#     forecasts_reference_abundance.rds, fit_info_reference_abundance.rds   (e05b)
#   01_Analysis/Data/e6_state_space_forecasts_abundance/
#     forecasts_ssm_abundance.rds, fit_info_ssm_abundance.rds               (e06)
#   01_Analysis/Data/e7_arimax_env_forecasts_abundance/
#     forecasts_arimax_abundance.rds, fit_info_arimax_abundance.rds         (e07)
#   01_Analysis/Data/e8_ssm_xreg_forecasts_abundance/
#     forecasts_ssm_xreg_abundance.rds, fit_info_ssm_xreg_abundance.rds     (e08)
#   In test mode the model outputs are read from the _TEST folders and files.
#
# Outputs:
#   01_Analysis/Data/e9_compare_predictability_abundance/
#     forecasts_joined_abundance.rds   keys, Common_Name, origin_i,
#       origin_date, in_common, lead, target_date, target_filled, scored,
#       model, role, source, ok, mean, lo80, hi80, lo95, hi95, y_obs,
#       season_at_target, value_obs
#     fit_info_joined_abundance.rds    one row per series, forecast date and
#       model: ok, reason, source (part: which covariate set, for combined
#       models)
#   02_Results/03_Diagnostics/e9_compare_predictability_abundance/e9_join_summary.csv
#     rows, fits and failures by model
#   01_Analysis/Logs/e9_compare_predictability_abundance/run_<timestamp>.txt
#   In test mode the stem and file names end in _TEST.
#   Earlier output files in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> before writing.
#
# Method notes:
#   1. All models and references as stored, plus the combined covariate
#      models (arimax, ssm_xreg: the lag 1 and 12 set at lead time 1, the lag
#      12 set at lead times 2 to 12). Roles and partners are in model_info().
#   2. Failed fits stay as rows (ok = FALSE, NA forecasts); f01 decides how
#      they enter paired comparisons.
#   3. Added columns: scored (FALSE for targets in April and May 2020),
#      season_at_target and value_obs (full scale), role, source.
#   4. Checks that stop the script: every model has every forecast date and
#      lead time; observed anomalies are identical across sources; value_obs
#      equals the input series; no duplicates; no NA forecast from a fit
#      marked ok.
#   5. RUN_MODE "test" reads the _TEST outputs of the model scripts for their
#      test series and the forecast dates they share.
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
RUN_MODE    <- "full"                                        # "full" or "test" (method note 5)
TEST_SERIES <- c("bay|2|BS|601", "coastwide|NA|BT|9605")   # used only in test mode
KEYS        <- c("scale", "BAY_ID", "GEAR", "SPECIES_CODE")

script_stem <- "e9_compare_predictability_abundance"   # stem of the output folders
SFX <- ""
if (RUN_MODE == "test") { script_stem <- paste0(script_stem, "_TEST"); SFX <- "_TEST" }

SOURCES <- list(   # name = c(script stem, forecast file, fit-information file)
  e05a = c("e5_baseline_forecasts_abundance",     "forecasts_sarima_family_abundance", "fit_info_sarima_family_abundance"),
  e05b = c("e5b_naive_baseline_abundance",        "forecasts_reference_abundance",     "fit_info_reference_abundance"),
  e06  = c("e6_state_space_forecasts_abundance",  "forecasts_ssm_abundance",           "fit_info_ssm_abundance"),
  e07  = c("e7_arimax_env_forecasts_abundance",   "forecasts_arimax_abundance",        "fit_info_arimax_abundance"),
  e08  = c("e8_ssm_xreg_forecasts_abundance",     "forecasts_ssm_xreg_abundance",      "fit_info_ssm_xreg_abundance"))
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
in_gamma    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
in_coefs    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "seasonal_coefs_by_origin_abundance.rds")
in_origins  <- file.path(PATHS$data_dir, "e4_forecast_design_abundance", "forecast_origins_abundance.rds")

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
                c(c("b09_forecast_settings.R", rep("e03_deseasonalize_abundance.R", 3), "e04_forecast_design_abundance.R"), rep(paste0(names(SOURCES), " (model script)"), each = 2)))
settings <- load_settings(in_settings)
sid <- function(d) do.call(paste, c(as.list(d[, KEYS, drop = FALSE]), sep = "|"))

a <- readRDS(in_alpha); g <- readRDS(in_gamma)
series <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, SPECIES_CODE = as.character(SPECIES_CODE), date, value),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, SPECIES_CODE = as.character(SPECIES_CODE), date, value))
coefs   <- as.data.frame(readRDS(in_coefs))
origins <- as.data.frame(readRDS(in_origins))
sources <- lapply(src_files, function(f) list(forecasts = readRDS(f[1]), info = readRDS(f[2])))
for (s in names(SOURCES)) cat("Read", s, ":", nrow(sources[[s]]$forecasts), "forecast rows,", nrow(sources[[s]]$info), "fit records\n")
if (length(unique(sid(origins))) != 128L) stop("Origin table does not have 128 series", call. = FALSE)

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
out_fc   <- file.path(dirs$data, paste0("forecasts_joined_abundance", SFX, ".rds"))
out_info <- file.path(dirs$data, paste0("fit_info_joined_abundance", SFX, ".rds"))
out_sum  <- file.path(dirs$diag, paste0("e9_join_summary", SFX, ".csv"))
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
cat("e09_compare_predictability_abundance.R completed\n")
