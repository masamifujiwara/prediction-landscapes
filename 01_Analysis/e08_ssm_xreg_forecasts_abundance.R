# =============================================================================
# Script:  e08_ssm_xreg_forecasts_abundance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   State-space forecasts with water temperature and salinity (local level
#   plus covariates) for every abundance series: 112 estuary series and 16
#   coastwide series (estuary or coastwide x sampling method x species), at
#   every forecast date of e04, lead times 1 to 12 months. Models
#   (lib/models_covariates.R), as in c10:
#     ssm_xreg_lag1_12  covariates at lags 1 and 12; used at lead time 1
#                       (lead time 1 only)
#     ssm_xreg_lag12    covariates at lag 12; used at lead times 2 to 12 (lead
#                       times 1 to 12 stored)
#   Partner without covariates: ssm_level in e06.
#
# Run after: b08, e04, e03, b09, lib/models_covariates.R, lib/seasonal.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds                       (b09)
#   01_Analysis/Data/b8_build_environment_timeseries/env_alpha_index.rds      (b08)
#   01_Analysis/Data/b8_build_environment_timeseries/env_gamma_index.rds      (b08)
#   01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds    (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds    (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/seasonal_coefs_by_origin_abundance.rds (e03)
#   01_Analysis/Data/e4_forecast_design_abundance/forecast_origins_abundance.rds (e04)
#
# Outputs:
#   01_Analysis/Data/e8_ssm_xreg_forecasts_abundance/
#     forecasts_ssm_xreg_abundance.rds   one row per series, forecast date,
#       model and lead time: scale, BAY_ID, GEAR, SPECIES_CODE, Common_Name,
#       origin_i, origin_date, in_common, model, lead, mean, lo80, hi80,
#       lo95, hi95, target_date, target_filled, y_obs (anomaly scale of that
#       forecast date)
#     fit_info_ssm_xreg_abundance.rds    one row per series, forecast date
#       and model
#     checkpoints/                       one file per task (resume)
#   02_Results/03_Diagnostics/e8_ssm_xreg_forecasts_abundance/e8_fit_summary.csv
#   01_Analysis/Logs/e8_ssm_xreg_forecasts_abundance/run_<timestamp>.txt
#   In test mode the stem and file names end in _TEST.
#   Earlier output files in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> before writing.
#
# Method notes:
#   Needs KFAS.
#   1. Anomalies at each forecast date come from the seasonal curves refitted
#      at that date (e03, on log(CPUE + 1)); forecast dates, training start
#      (month 13) and longest lead time come from e04; one fit per forecast
#      date and covariate set serves every lead time. No minimum count of
#      positive training values is applied: on centered anomalies such a rule
#      would count months above the seasonal curve.
#   2. No seasonal part in the state-space model (the series are already
#      departures from the seasonal curve). Covariates from b08: TEMP and SAL
#      monthly means of the same estuary and sampling method (coastwide: the
#      sampling method's coastwide means), lagged only, standardized with the
#      training rows' mean and SD; lag columns are checked against the base
#      columns; a missing covariate in a training or target row stops the
#      script. Covariate gaps are not interpolated over the whole series,
#      which would look ahead.
#   3. Two covariate sets (settings$covariate_sets): lags 1 and 12 for lead
#      time 1, lag 12 for lead times 2 to 12; lead times 1 to 12 only.
#   4. Model code is only in lib/models.R and lib/models_covariates.R.
#   5. Restarts and collapse checks as in c08; failed fits are kept as rows
#      with ok = FALSE and a named reason.
#   6. Targets in April and May 2020 are kept with target_filled = TRUE and
#      removed at scoring; those months are never forecast dates.
#   7. Parallel over tasks of one series and a block of forecast dates, with
#      a checkpoint per task and a signature; a rerun after an interruption
#      resumes from the checkpoints.
#   8. RUN_MODE "test" runs a few series and forecast dates under the stem
#      e8_ssm_xreg_forecasts_abundance_TEST with _TEST file names.
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
source(file.path(PATHS$analysis_dir, "lib", "models.R"))
source(file.path(PATHS$analysis_dir, "lib", "models_covariates.R"))
suppressPackageStartupMessages(library(KFAS))   # ssm_xreg() needs KFAS attached

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
RUN_MODE           <- "full"           # "full" or "test" (method note 8)
N_WORKERS_OVERRIDE <- NA_integer_      # set an integer to force a count
ORIGIN_BLOCK       <- 30L              # origins per task
TEST_SERIES        <- c("bay|2|BS|601", "coastwide|NA|BT|9605")   # used only in test mode
TEST_ORIGIN_EVERY  <- 60L                                   # used only in test mode

script_stem <- "e8_ssm_xreg_forecasts_abundance"   # stem of the output folders
SFX <- ""
if (RUN_MODE == "test") { script_stem <- paste0(script_stem, "_TEST"); SFX <- "_TEST" }

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_env_a    <- file.path(PATHS$data_dir, "b8_build_environment_timeseries", "env_alpha_index.rds")
in_env_g    <- file.path(PATHS$data_dir, "b8_build_environment_timeseries", "env_gamma_index.rds")
in_alpha    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
in_gamma    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
in_coefs    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "seasonal_coefs_by_origin_abundance.rds")
in_origins  <- file.path(PATHS$data_dir, "e4_forecast_design_abundance", "forecast_origins_abundance.rds")
this_script <- file.path(PATHS$analysis_dir, "e08_ssm_xreg_forecasts_abundance.R")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr", "KFAS", "future", "future.apply")
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
stop_if_missing(c(in_settings, in_env_a, in_env_g, in_alpha, in_gamma, in_coefs, in_origins),
                c("b09_forecast_settings.R", rep("b08_build_environment_timeseries.R", 2),
                  rep("e03_deseasonalize_abundance.R", 3), "e04_forecast_design_abundance.R"))
settings <- load_settings(in_settings)
H_MAX  <- max(settings$horizons_covariate)
SET_1  <- covariate_columns(settings$covariate_sets$lead_1)       # lead 1
SET_2  <- covariate_columns(settings$covariate_sets$leads_2_12)   # leads 2 to 12
COVARS <- union(SET_1, SET_2)
cat("Leads 1 to", H_MAX, "| lead-1 covariates:", paste(SET_1, collapse = ", "),
    "| leads 2 to 12:", paste(SET_2, collapse = ", "), "\n")

MODEL_FUN_X <- function(y, H, X_train, X_future) {
  c(pad_leads(ssm_xreg(y, 1L, X_train[, SET_1, drop = FALSE], X_future[1, SET_1, drop = FALSE], "ssm_xreg_lag1_12"), H),
    ssm_xreg(y, H,  X_train[, SET_2, drop = FALSE], X_future[, SET_2, drop = FALSE], "ssm_xreg_lag12"))
}

# Covariates: keep the base and lag columns used, and check the lags (method note 2).
lag_check <- function(env, by) {
  env <- env %>% arrange(across(all_of(c(by, "DATE"))))
  for (v in c("TEMP__mean", "SAL__mean")) for (l in unique(c(settings$covariate_sets$lead_1, settings$covariate_sets$leads_2_12))) {
    col <- paste0(v, "__lag", l)
    d <- env %>% group_by(across(all_of(by))) %>%
      mutate(.ok = is.na(.data[[col]]) | abs(.data[[col]] - dplyr::lag(.data[[v]], l)) < 1e-12,
             .gap = c(NA, diff(as.integer(format(DATE, "%Y")) * 12L + as.integer(format(DATE, "%m"))))) %>% ungroup()
    if (any(d$.gap != 1L, na.rm = TRUE)) stop("Covariate file is not a complete monthly sequence", call. = FALSE)
    if (!all(d$.ok)) stop("Covariate column ", col, " is not ", v, " shifted by ", l, " months", call. = FALSE)
  }
  invisible(TRUE)
}
env_a <- readRDS(in_env_a); env_g <- readRDS(in_env_g)
lag_check(env_a, c("BAY_ID", "GEAR")); lag_check(env_g, "GEAR")
cat("Covariate lag columns checked against the base columns: OK\n")
env <- bind_rows(
  env_a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, date = DATE, across(all_of(COVARS))),
  env_g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, date = DATE, across(all_of(COVARS))))

a <- readRDS(in_alpha); g <- readRDS(in_gamma)
ts_all <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, SPECIES_CODE = as.character(SPECIES_CODE), Common_Name, date, MONTH, filled, value),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, SPECIES_CODE = as.character(SPECIES_CODE), Common_Name, date, MONTH, filled, value))
n_before <- nrow(ts_all)
ts_all <- ts_all %>% left_join(env, by = c("scale", "BAY_ID", "GEAR", "date"))
if (nrow(ts_all) != n_before) stop("Joining the covariates changed the number of rows", call. = FALSE)
coefs   <- readRDS(in_coefs)
origins <- readRDS(in_origins)
sid <- function(d) paste(d$scale, d$BAY_ID, d$GEAR, d$SPECIES_CODE, sep = "|")
ts_all$sid <- sid(ts_all); coefs$sid <- sid(coefs); origins$sid <- sid(origins)
series_ids <- sort(unique(origins$sid))
if (length(series_ids) != 128L) stop(length(series_ids), " series in the origin table; expected 128", call. = FALSE)

# Every training and target row of every origin must have all covariates (method note 2).
first_train <- origins %>% group_by(sid) %>% summarise(start = min(train_start_i), .groups = "drop")
miss_tab <- ts_all %>% group_by(sid) %>% arrange(date, .by_group = TRUE) %>% mutate(row = row_number()) %>%
  inner_join(first_train, by = "sid") %>% filter(row >= start) %>%
  summarise(rows = n(), missing = sum(!stats::complete.cases(across(all_of(COVARS)))), .groups = "drop")
if (any(miss_tab$missing > 0)) {
  print(as.data.frame(miss_tab %>% filter(missing > 0)), row.names = FALSE)
  stop("Missing covariate values after the training start (listed above)", call. = FALSE)
}
cat("Covariates complete from the training start in all", nrow(miss_tab), "series\n")

if (RUN_MODE == "test") {
  series_ids <- intersect(series_ids, TEST_SERIES)
  origins <- origins %>% filter(sid %in% series_ids) %>% group_by(sid) %>% filter(row_number() %% TEST_ORIGIN_EVERY == 1L) %>% ungroup()
}

# Tasks: one series and a block of origins.
tasks <- bind_rows(lapply(series_ids, function(s) {
  oi <- origins$origin_i[origins$sid == s]
  data.frame(sid = s, block = ceiling(seq_along(oi) / ORIGIN_BLOCK), origin_i = oi)
})) %>% group_by(sid, block) %>% summarise(origin_i = list(origin_i), .groups = "drop")
tasks$task_id <- sprintf("%s__b%02d", gsub("[|]", "_", tasks$sid), tasks$block)
cat("Series:", length(series_ids), "| origins:", nrow(origins), "| tasks:", nrow(tasks), "\n")

# -----------------------------------------------------------------------------
# 3) Checkpoints (method note 7)
# -----------------------------------------------------------------------------
ck_dir <- file.path(dirs$data, "checkpoints")
dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
lib_files <- file.path(PATHS$analysis_dir, "lib", c("models.R", "models_covariates.R", "seasonal.R", "forecast_rules.R", "helpers.R"))
signature <- paste(c(unname(tools::md5sum(c(lib_files, in_settings, in_env_a, in_env_g, in_alpha, in_gamma, in_coefs, in_origins,
                                            if (file.exists(this_script)) this_script))),
                     RUN_MODE, ORIGIN_BLOCK, paste(R.version$major, R.version$minor),
                     as.character(utils::packageVersion("KFAS"))), collapse = "|")
sig_file <- file.path(ck_dir, "signature.txt")
if (file.exists(sig_file) && !identical(readLines(sig_file, warn = FALSE), signature)) {
  old <- list.files(ck_dir, full.names = TRUE)
  cat("Code, settings or inputs changed since the checkpoints were written: archiving", length(old), "checkpoint files\n")
  archive_files(old)
}
writeLines(signature, sig_file)
done <- sub("[.]rds$", "", list.files(ck_dir, pattern = "[.]rds$"))
todo <- tasks[!tasks$task_id %in% done, ]
cat("Checkpoints found:", sum(tasks$task_id %in% done), "| tasks to run:", nrow(todo), "\n")

# -----------------------------------------------------------------------------
# 4) Forecasts (parallel)
# -----------------------------------------------------------------------------
run_task <- function(k) {
  tk <- todo[k, ]; s <- tk$sid
  ser <- ts_all[ts_all$sid == s, ]; ser <- ser[order(ser$date), ]
  og  <- origins[origins$sid == s & origins$origin_i %in% tk$origin_i[[1]], ]
  co  <- coefs[coefs$sid == s, ]
  r <- forecast_series_xreg(ser, og, co, MODEL_FUN_X, COVARS, H_max = H_MAX)
  r$forecasts <- keep_leads(r$forecasts, c(ssm_xreg_lag1_12 = 1L))   # lead-1 model: drop the padded rows
  key <- ser[1, c("scale", "BAY_ID", "GEAR", "SPECIES_CODE", "Common_Name")]
  r$forecasts <- cbind(key[rep(1, nrow(r$forecasts)), ], r$forecasts)
  r$info      <- cbind(key[rep(1, nrow(r$info)), ], r$info)
  saveRDS(r, file.path(ck_dir, paste0(tk$task_id, ".rds")))
  tk$task_id
}
if (nrow(todo) > 0) {
  w <- setup_workers(override = N_WORKERS_OVERRIDE, globals_max_gb = 4)
  on.exit(future::plan(future::sequential), add = TRUE)
  t_fit <- Sys.time()
  invisible(future.apply::future_lapply(seq_len(nrow(todo)), run_task, future.seed = settings$seed,
                                        future.scheduling = Inf, future.packages = "KFAS"))
  cat("Fitting time:", round(as.numeric(difftime(Sys.time(), t_fit, units = "mins")), 1), "min on", w$n_workers, "workers\n")
}

# -----------------------------------------------------------------------------
# 5) Combine, check, write
# -----------------------------------------------------------------------------
parts <- lapply(file.path(ck_dir, paste0(tasks$task_id, ".rds")), function(f) {
  if (!file.exists(f)) stop("Missing checkpoint ", basename(f), call. = FALSE)
  readRDS(f)
})
fc   <- bind_rows(lapply(parts, `[[`, "forecasts"))
info <- bind_rows(lapply(parts, `[[`, "info"))
n_models <- 2L
if (nrow(info) != nrow(origins) * n_models) stop("Fit records: ", nrow(info), "; expected ", nrow(origins) * n_models, call. = FALSE)
if (anyDuplicated(fc[, c("scale", "BAY_ID", "GEAR", "SPECIES_CODE", "origin_i", "model", "lead")])) stop("Duplicate forecast rows", call. = FALSE)
if (any(fc$model == "ssm_xreg_lag1_12" & fc$lead != 1L)) stop("Lag 1 and 12 model has leads other than 1", call. = FALSE)
bad <- fc %>% semi_join(info %>% filter(ok), by = c("scale", "BAY_ID", "GEAR", "SPECIES_CODE", "origin_i", "model")) %>% filter(is.na(mean))
if (nrow(bad)) stop(nrow(bad), " NA forecasts from fits marked ok", call. = FALSE)

summary_tab <- info %>% group_by(model, scale) %>%
  summarise(fits = n(), failed = sum(!ok), restarts_used = sum(order != "tries=1", na.rm = TRUE),
            median_seconds = median(seconds, na.rm = TRUE), total_hours = sum(seconds, na.rm = TRUE) / 3600,
            warnings = sum(n_warnings), .groups = "drop")
fail_tab <- info %>% filter(!ok) %>% count(model, reason, name = "n")

out_fc   <- file.path(dirs$data, paste0("forecasts_ssm_xreg_abundance", SFX, ".rds"))
out_info <- file.path(dirs$data, paste0("fit_info_ssm_xreg_abundance", SFX, ".rds"))
out_sum  <- file.path(dirs$diag, paste0("e8_fit_summary", SFX, ".csv"))
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE)
if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|pdf|tiff|png)$", full.names = TRUE)
if (length(old_d)) archive_files(old_d)
saveRDS(fc, out_fc); saveRDS(info, out_info); readr::write_csv(summary_tab, out_sum, na = "")

cat("\nFits by model and scale:\n"); print(as.data.frame(summary_tab %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\nFailures by reason:\n"); if (nrow(fail_tab)) print(as.data.frame(fail_tab), row.names = FALSE) else cat("  none\n")
cat("\nForecast rows:", nrow(fc), "\n")
for (f in c(out_fc, out_info, out_sum)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("e08_ssm_xreg_forecasts_abundance.R completed\n")
