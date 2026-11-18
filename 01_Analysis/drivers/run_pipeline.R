# =============================================================================
# Script:  run_pipeline.R (drivers/)
# Author:  Masami Fujiwara
#
# Purpose:
#   Runs the top-level scripts of 01_Analysis in run order, each in its own
#   Rscript process, and stops at the first failure. The run order is
#   file-name order with the exceptions listed in lib/run_order.R (f02 runs
#   before f01). The pipeline scripts are a00, b01, b03 to b09, c00 to c03, c07 to
#   c10, c12, d01, e02 to e09, f01 to f04, f07, g01 to g08, g11 to g17, g20,
#   h01, h02, h11 to h17 and h20.
#   Run from anywhere: Rscript drivers/run_pipeline.R (nohup works too), or
#   Source it in RStudio.
#
# Run after: none (it runs the pipeline).
#
# Inputs:
#   01_Analysis/*.R (the top-level scripts)
#   01_Analysis/lib/run_order.R
#   01_Analysis/drivers/check_order.R (run first when CHECK_ORDER is TRUE)
#
# Outputs:
#   01_Analysis/Logs/run_pipeline_<timestamp>/00_check_order.log
#   01_Analysis/Logs/run_pipeline_<timestamp>/<NN>_<script>.log (one per
#     script, named by its position and file name)
#   01_Analysis/Logs/run_pipeline_<timestamp>/summary.csv (machine name,
#     status, minutes and start time per script)
#   Each run writes a new timestamped folder; nothing is archived.
#
# Method notes:
#   Settings:
#   START_AT    first script to run (file name); NA = from the beginning.
#   STOP_AFTER  last script to run; NA = to the end.
#   SKIP        scripts never run by this driver. a00_getdata.R is skipped by
#               default: it rebuilds Data.RData from the raw TPWD files outside
#               the project, and Data.RData is the frozen input of the analysis.
#   SKIP_MODELS TRUE (default): skip the forecast scripts c07 to c12 and e05a
#               to e08 (models and reference forecasts, matched by their
#               number prefix), which take days to run. The later scripts
#               (d01, e09, f, g, h) then read the forecasts already saved from
#               earlier runs. Those forecasts match only as long as the
#               upstream outputs (b, c00 to c03, e02 to e04) are unchanged.
#               FALSE: run everything (several days).
#   CHECK_ORDER run drivers/check_order.R first and stop if it fails.
#   The machine name is written to the summary, because two machines can give
#   slightly different SARIMA results.
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

START_AT    <- NA_character_          # e.g. "c00_inputs_check.R"
STOP_AFTER  <- NA_character_          # e.g. "b08_build_environment_timeseries.R"
SKIP        <- c("a00_getdata.R")
SKIP_MODELS <- TRUE
CHECK_ORDER <- TRUE

an      <- PATHS$analysis_dir
source(file.path(an, "lib", "run_order.R"))
scripts <- run_order(list.files(an, pattern = "\\.R$"))
pick    <- rep(TRUE, length(scripts))
if (!is.na(START_AT)) {
  if (!START_AT %in% scripts) stop("START_AT not found: ", START_AT, call. = FALSE)
  pick <- pick & seq_along(scripts) >= match(START_AT, scripts)
}
if (!is.na(STOP_AFTER)) {
  if (!STOP_AFTER %in% scripts) stop("STOP_AFTER not found: ", STOP_AFTER, call. = FALSE)
  pick <- pick & seq_along(scripts) <= match(STOP_AFTER, scripts)
}
# Forecast scripts are matched by their number prefix, so the list does not
# depend on the rest of the name; it selects c07 to c10, c12 and e05a to e08.
model_scripts <- grep("^(c0[5-9]|c1[0-2]|e0[5-8][ab]?)_", scripts, value = TRUE)
if (SKIP_MODELS) {
  cat("SKIP_MODELS = TRUE: skipping", length(model_scripts), "model scripts:\n")
  cat(paste0("  ", model_scripts), sep = "\n")
  SKIP <- union(SKIP, model_scripts)
}
pick <- pick & !scripts %in% SKIP
todo <- scripts[pick]
if (!length(todo)) stop("Nothing to run.", call. = FALSE)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
stopifnot(file.exists(rscript))
log_dir <- file.path(PATHS$logs_dir, format(Sys.time(), "run_pipeline_%Y%m%d_%H%M%S"))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

run_one <- function(path, log) {
  system2(rscript, c("--no-save", "--no-restore", shQuote(path)), stdout = log, stderr = log)
}

if (CHECK_ORDER) {
  rc <- run_one(file.path(an, "drivers", "check_order.R"), file.path(log_dir, "00_check_order.log"))
  if (rc != 0L) stop("check_order.R failed; see ", file.path(log_dir, "00_check_order.log"), call. = FALSE)
  cat("check_order.R: OK\n")
}

cat("Running", length(todo), "scripts on", Sys.info()[["nodename"]], "\n")
summary <- data.frame(n = seq_along(todo), script = todo, status = "not run",
                      minutes = NA_real_, started = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(todo)) {
  t0  <- Sys.time()
  log <- file.path(log_dir, sprintf("%02d_%s", i, sub("\\.R$", ".log", todo[i])))
  cat(sprintf("[%2d/%d] %-62s %s ", i, length(todo), todo[i], format(t0, "%d %b %H:%M")))
  rc  <- run_one(file.path(an, todo[i]), log)
  summary$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  summary$started[i] <- format(t0, "%Y-%m-%d %H:%M:%S")
  summary$status[i]  <- if (rc == 0L) "OK" else paste0("FAILED (", rc, ")")
  cat(sprintf("%s  %.1f min\n", summary$status[i], summary$minutes[i]))
  write.csv(cbind(machine = Sys.info()[["nodename"]], summary),
            file.path(log_dir, "summary.csv"), row.names = FALSE)
  if (rc != 0L) {
    cat("\nStopped. Last lines of ", log, ":\n", sep = "")
    cat(tail(readLines(log, warn = FALSE), 20), sep = "\n")
    cat("\nFix it, then set START_AT <- \"", todo[i], "\" and run again.\n", sep = "")
    break
  }
}
cat("\nLogs in", log_dir, "\n")
