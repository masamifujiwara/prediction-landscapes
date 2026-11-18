# =============================================================================
# Script:  c12_naive_baseline_diversity.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Reference forecasts for every diversity series (42 estuary and 6 coastwide
#   series), at every forecast date from c03, leads 1 to 24 (lib/models.R):
#     ref_naive         last month's value (last anomaly at the forecast date)
#     ref_past_average  cumulative average: mean of all sampled anomalies from
#                       the series start to the forecast date (zero by
#                       construction, since it equals the seasonal curve
#                       fitted to the same window); the primary reference
#                       for skill
#     ref_12month       12-month average: mean of the sampled anomalies of the
#                       last 12 months
#   The references are written to their own file. The file name keeps
#   "naive_baseline" from the earlier script name.
#
# Run after: b09, c02, c03
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds                       (b09)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/seasonal_coefs_by_origin_scaled.rds (c02)
#   - 01_Analysis/Data/c3_forecast_design_backtesting/forecast_origins_diversity.rds (c03)
#
# Outputs (01_Analysis/Data/c12_naive_baseline_diversity/):
#   - forecasts_reference_diversity.rds      one row per series, forecast date,
#       reference, lead: scale, BAY_ID, GEAR, q, origin_i, origin_date,
#       in_common, model, lead, mean, lo80, hi80, lo95, hi95, target_date,
#       target_filled, y_obs (all on the anomaly scale of that forecast date;
#       no intervals: lo/hi are NA)
#   - fit_info_reference_diversity.rds       one row per series, forecast date,
#       reference
#   - checkpoints/                           one file per task (resume)
#   02_Results/03_Diagnostics/<stem>/c12_reference_summary.csv
#   01_Analysis/Logs/<stem>/run_<timestamp>.txt
#   Earlier .rds and .csv files in the data folder and earlier files in the
#   diagnostics folder are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Three references: last month's value, the cumulative average (primary)
#     and the 12-month average. The mean over the full record is not used as
#     a reference, because it uses data after the forecast date.
#   - Anomalies (departures from the seasonal curve) at each forecast date use
#     the seasonal curves fitted up to that date (c02); forecast dates come
#     from c03. References use all past months, not only the model training
#     rows from month 13.
#   - Targets in April and May 2020 are kept with target_filled = TRUE and
#     removed at scoring (drop_unscored_targets()).
#   - Same task, checkpoint and test-mode mechanism as c07 (RUN_MODE "test"
#     writes under a _TEST stem with _TEST file names).
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
RUN_MODE           <- "full"           # "full" or "test" (test runs use _TEST names)
N_WORKERS_OVERRIDE <- NA_integer_      # set an integer to force a count
ORIGIN_BLOCK       <- 30L              # origins per task
TEST_SERIES        <- c("bay|2|BS|0", "coastwide|NA|BT|2")   # used only in test mode
TEST_ORIGIN_EVERY  <- 60L                                   # used only in test mode

script_stem <- "c12_naive_baseline_diversity"   # folder stem kept from the earlier script name
SFX <- ""
if (RUN_MODE == "test") { script_stem <- paste0(script_stem, "_TEST"); SFX <- "_TEST" }

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_gamma    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_coefs    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "seasonal_coefs_by_origin_scaled.rds")
in_origins  <- file.path(PATHS$data_dir, "c3_forecast_design_backtesting", "forecast_origins_diversity.rds")
this_script <- file.path(PATHS$analysis_dir, "c12_naive_baseline_diversity.R")

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
                c("b09_forecast_settings.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 3), "c03_forecast_design_backtesting.R"))
settings <- load_settings(in_settings)
H_MAX <- max(settings$horizons)

a <- readRDS(in_alpha); g <- readRDS(in_gamma)
ts_all <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, q = as.integer(q), date, MONTH, filled, value),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, q = as.integer(q), date, MONTH, filled, value))
coefs   <- readRDS(in_coefs)
origins <- readRDS(in_origins)
sid <- function(d) paste(d$scale, d$BAY_ID, d$GEAR, d$q, sep = "|")
ts_all$sid <- sid(ts_all); coefs$sid <- sid(coefs); origins$sid <- sid(origins)
series_ids <- sort(unique(origins$sid))
if (length(series_ids) != 48L) stop(length(series_ids), " series in the origin table; expected 48", call. = FALSE)
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
# 3) Checkpoints (resume only if code, settings and inputs are unchanged)
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
  key <- ser[1, c("scale", "BAY_ID", "GEAR", "q")]
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
if (anyDuplicated(fc[, c("scale", "BAY_ID", "GEAR", "q", "origin_i", "model", "lead")])) stop("Duplicate forecast rows", call. = FALSE)
if (any(is.na(fc$mean) & fc$model %in% info$model[info$ok])) {
  bad <- fc %>% semi_join(info %>% filter(ok), by = c("scale", "BAY_ID", "GEAR", "q", "origin_i", "model")) %>% filter(is.na(mean))
  if (nrow(bad)) stop(nrow(bad), " NA forecasts from fits marked ok", call. = FALSE)
}

summary_tab <- info %>% group_by(model, scale) %>%
  summarise(origins = n(), failed = sum(!ok), .groups = "drop") %>%
  left_join(fc %>% group_by(model, scale) %>% summarise(mean_forecast = mean(mean), max_abs_forecast = max(abs(mean)), .groups = "drop"),
            by = c("model", "scale"))
if (max(abs(fc$mean[fc$model == "ref_past_average"])) > 1e-8) stop("The past average is not zero on the anomaly scale", call. = FALSE)
fail_tab <- info %>% filter(!ok) %>% count(model, reason, name = "n")

out_fc   <- file.path(dirs$data, paste0("forecasts_reference_diversity", SFX, ".rds"))
out_info <- file.path(dirs$data, paste0("fit_info_reference_diversity", SFX, ".rds"))
out_sum  <- file.path(dirs$diag, paste0("c12_reference_summary", SFX, ".csv"))
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
cat("c12_naive_baseline_diversity.R completed\n")
