# =============================================================================
# Script:  e05b_naive_baseline_abundance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Reference forecasts for every abundance series (112 estuary and 16
#   coastwide series), at every forecast date of e04, lead times 1 to 24
#   months (lib/models.R), as in c12:
#     ref_naive         last month's value: the last anomaly at the forecast
#                       date
#     ref_past_average  cumulative average: mean of all sampled anomalies from
#                       the series start to the forecast date (zero by
#                       construction, so it equals the seasonal curve fitted
#                       in the training window); the primary reference for
#                       skill
#     ref_12month       12-month average: mean of the sampled anomalies of the
#                       last 12 months
#   The references are written to their own file.
#
# Run after: e04, e03, b09, lib/models.R, lib/seasonal.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds                       (b09)
#   01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds    (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds    (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/seasonal_coefs_by_origin_abundance.rds (e03)
#   01_Analysis/Data/e4_forecast_design_abundance/forecast_origins_abundance.rds (e04)
#
# Outputs:
#   01_Analysis/Data/e5b_naive_baseline_abundance/
#     forecasts_reference_abundance.rds   one row per series, forecast date,
#       reference and lead time: scale, BAY_ID, GEAR, SPECIES_CODE,
#       Common_Name, origin_i, origin_date, in_common, model, lead, mean,
#       lo80, hi80, lo95, hi95, target_date, target_filled, y_obs (anomaly
#       scale of that forecast date; no intervals, so lo/hi are NA)
#     fit_info_reference_abundance.rds    one row per series, forecast date
#       and reference
#     checkpoints/                        one file per task (resume)
#   02_Results/03_Diagnostics/e5b_naive_baseline_abundance/e5b_reference_summary.csv
#   01_Analysis/Logs/e5b_naive_baseline_abundance/run_<timestamp>.txt
#   In test mode the stem and file names end in _TEST.
#   Earlier output files in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> before writing.
#
# Method notes:
#   1. Three references: last month's value, cumulative average (primary)
#      and 12-month average. Neither the full-record mean nor the mean of the
#      test period is a reference, because both use future data.
#   2. Anomalies at each forecast date come from the seasonal curves refitted
#      at that date (e03, on log(CPUE + 1)); forecast dates come from e04.
#      References use all past months, not only the model training rows from
#      month 13.
#   3. Targets in April and May 2020 are kept with target_filled = TRUE and
#      removed at scoring (drop_unscored_targets()).
#   4. Same task, checkpoint and test-mode mechanism as e05a (the script is
#      fast, so checkpoints rarely matter).
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

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
RUN_MODE           <- "full"           # "full" or "test" (method note 4)
N_WORKERS_OVERRIDE <- NA_integer_      # set an integer to force a count
ORIGIN_BLOCK       <- 30L              # origins per task
TEST_SERIES        <- c("bay|2|BS|601", "coastwide|NA|BT|9605")   # used only in test mode
TEST_ORIGIN_EVERY  <- 60L                                   # used only in test mode

script_stem <- "e5b_naive_baseline_abundance"   # stem of the output folders
SFX <- ""
if (RUN_MODE == "test") { script_stem <- paste0(script_stem, "_TEST"); SFX <- "_TEST" }

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
in_gamma    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
in_coefs    <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "seasonal_coefs_by_origin_abundance.rds")
in_origins  <- file.path(PATHS$data_dir, "e4_forecast_design_abundance", "forecast_origins_abundance.rds")
this_script <- file.path(PATHS$analysis_dir, "e05b_naive_baseline_abundance.R")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr", "future", "future.apply")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

MODEL_FUN <- NULL   # references only (forecast_series(references = TRUE))

# -----------------------------------------------------------------------------
# 2) Folders, log, inputs
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "| run mode:", RUN_MODE, "\n")
stop_if_missing(c(in_settings, in_alpha, in_gamma, in_coefs, in_origins),
                c("b09_forecast_settings.R", rep("e03_deseasonalize_abundance.R", 3), "e04_forecast_design_abundance.R"))
settings <- load_settings(in_settings)
H_MAX <- max(settings$horizons)

a <- readRDS(in_alpha); g <- readRDS(in_gamma)
ts_all <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, SPECIES_CODE = as.character(SPECIES_CODE), Common_Name, date, MONTH, filled, value),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, SPECIES_CODE = as.character(SPECIES_CODE), Common_Name, date, MONTH, filled, value))
coefs   <- readRDS(in_coefs)
origins <- readRDS(in_origins)
sid <- function(d) paste(d$scale, d$BAY_ID, d$GEAR, d$SPECIES_CODE, sep = "|")
ts_all$sid <- sid(ts_all); coefs$sid <- sid(coefs); origins$sid <- sid(origins)
series_ids <- sort(unique(origins$sid))
if (length(series_ids) != 128L) stop(length(series_ids), " series in the origin table; expected 128", call. = FALSE)
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
cat("Series:", length(series_ids), "| origins:", nrow(origins), "| tasks:", nrow(tasks), "| leads 1 to", H_MAX, "\n")

# -----------------------------------------------------------------------------
# 3) Checkpoints (method note 4)
# -----------------------------------------------------------------------------
ck_dir <- file.path(dirs$data, "checkpoints")
dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
lib_files <- file.path(PATHS$analysis_dir, "lib", c("models.R", "seasonal.R", "forecast_rules.R", "helpers.R"))
signature <- paste(c(unname(tools::md5sum(c(lib_files, in_settings, in_alpha, in_gamma, in_coefs, in_origins,
                                            if (file.exists(this_script)) this_script))),
                     RUN_MODE, ORIGIN_BLOCK, paste(R.version$major, R.version$minor),
                     "references"), collapse = "|")
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
  r <- forecast_series(ser, og, co, H_max = H_MAX, references = TRUE)
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
                                        future.scheduling = Inf))
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
n_models <- 3L
if (nrow(info) != nrow(origins) * n_models) stop("Fit records: ", nrow(info), "; expected ", nrow(origins) * n_models, call. = FALSE)
if (anyDuplicated(fc[, c("scale", "BAY_ID", "GEAR", "SPECIES_CODE", "origin_i", "model", "lead")])) stop("Duplicate forecast rows", call. = FALSE)
if (any(is.na(fc$mean) & fc$model %in% info$model[info$ok])) {
  bad <- fc %>% semi_join(info %>% filter(ok), by = c("scale", "BAY_ID", "GEAR", "SPECIES_CODE", "origin_i", "model")) %>% filter(is.na(mean))
  if (nrow(bad)) stop(nrow(bad), " NA forecasts from fits marked ok", call. = FALSE)
}

summary_tab <- info %>% group_by(model, scale) %>%
  summarise(origins = n(), failed = sum(!ok), .groups = "drop") %>%
  left_join(fc %>% group_by(model, scale) %>% summarise(mean_forecast = mean(mean), max_abs_forecast = max(abs(mean)), .groups = "drop"),
            by = c("model", "scale"))
if (max(abs(fc$mean[fc$model == "ref_past_average"])) > 1e-8) stop("The past average is not zero on the anomaly scale", call. = FALSE)
fail_tab <- info %>% filter(!ok) %>% count(model, reason, name = "n")

out_fc   <- file.path(dirs$data, paste0("forecasts_reference_abundance", SFX, ".rds"))
out_info <- file.path(dirs$data, paste0("fit_info_reference_abundance", SFX, ".rds"))
out_sum  <- file.path(dirs$diag, paste0("e5b_reference_summary", SFX, ".csv"))
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE)
if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|pdf|tiff|png)$", full.names = TRUE)
if (length(old_d)) archive_files(old_d)
saveRDS(fc, out_fc); saveRDS(info, out_info); readr::write_csv(summary_tab, out_sum, na = "")

cat("\nReferences by scale:\n"); print(as.data.frame(summary_tab %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\nFailures by reason:\n"); if (nrow(fail_tab)) print(as.data.frame(fail_tab), row.names = FALSE) else cat("  none\n")
cat("\nForecast rows:", nrow(fc), "\n")
for (f in c(out_fc, out_info, out_sum)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("e05b_naive_baseline_abundance.R completed\n")
