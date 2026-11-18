# =============================================================================
# File:    lib/helpers.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Shared helper functions for the pipeline scripts: input checks, output
#   folders, run logs, archiving of earlier outputs, parallel workers,
#   warning capture, monthly grids and model-fit checks. A script loads it
#   right after the standard opening block:
#     source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
#
# Used by: every top-level script (a00 to h20); also needed by
#   lib/forecast_rules.R, lib/models.R and lib/models_covariates.R
#
# Functions:
#   %||%                    NULL default (base R has it from 4.4; defined here for older R)
#   stop_if_missing()       stop with a clear message if input files are missing
#   expect_cols()           stop if a data frame lacks required columns
#   script_dirs()           this script's Data, Logs, 03_Diagnostics, 04_Models, 01_Figures folders
#   start_log(), end_log()  run log with console copy, error handler and sessionInfo()
#   archive_previous()      move existing outputs in whole folders to ARCHIVE_<timestamp>
#   archive_files()         the same for a list of individual output files
#   setup_workers()         the standard worker block (one physical core free, at most 7)
#   with_warnings()         run code, keep its value and count its warnings
#   complete_monthly_grid() full monthly date grid, left-joined to a series
#   check_kfas_fit()        convergence and variance-collapse checks for KFAS::fitSSM
#   check_arima_fit()       optimizer code and basic checks for forecast/stats ARIMA fits
#   coerce_keys()           standard key column types (not used by the pipeline scripts)
#
# Outputs:
#   start_log() writes the run log <log_dir>/<prefix><timestamp>.txt of the
#   calling script. archive_previous() and archive_files() move earlier
#   outputs to ARCHIVE_<timestamp> folders.
#
# Method notes:
#   Needs base R; dplyr and tibble for complete_monthly_grid(); future for
#   setup_workers(). Packages are called with ::, nothing is attached.
#   Where a script still defines its own copy of a helper, the local copy
#   wins because it is defined later.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

if (!exists("PATHS", inherits = TRUE)) {
  stop("lib/helpers.R needs PATHS: source it after the standard opening block.")
}

# ---- %||% --------------------------------------------------------------------
# Base R has %||% from version 4.4. Several scripts define it themselves; this
# keeps them working on older R.
if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ---- Inputs ------------------------------------------------------------------

# Stop if any required input file is missing. `producer` names the script that
# writes each file (one name for all, or one per path), so the message says
# what to run.
stop_if_missing <- function(paths, producer = NULL) {
  miss <- !file.exists(paths)
  if (!any(miss)) return(invisible(TRUE))
  prod <- if (is.null(producer)) "" else rep_len(producer, length(paths))
  msg <- paste0("  ", paths[miss],
                if (is.null(producer)) "" else paste0("   (written by ", prod[miss], ")"),
                collapse = "\n")
  stop("Required input file(s) not found:\n", msg, call. = FALSE)
}

# Stops with a message naming the missing columns.
expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

# ---- Folders -----------------------------------------------------------------

# Folders for one script, named by its fixed stem (the script_stem the script
# sets, for example "f1_forecast_scoring").
script_dirs <- function(stem, which = c("data", "logs", "diag", "models"), create = TRUE) {
  roots <- list(
    data   = PATHS$data_dir,
    logs   = PATHS$logs_dir,
    diag   = file.path(PATHS$results_dir, "03_Diagnostics"),
    models = file.path(PATHS$results_dir, "04_Models"),
    figs   = file.path(PATHS$results_dir, "01_Figures")
  )
  bad <- setdiff(which, names(roots))
  if (length(bad)) stop("script_dirs(): unknown folder type(s): ", paste(bad, collapse = ", "))
  out <- lapply(roots[which], function(r) file.path(r, stem))
  if (create) for (d in out) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  out
}

# ---- Run log -----------------------------------------------------------------
# on.exit() does not work at the top level of a script (under source() it fires
# at once, under Rscript never), so the log is closed by end_log() at the end
# and by an error handler if the script stops early. This is the pattern used
# in c02, c07 to c10 and the f, g and h scripts.

.helpers_log <- new.env(parent = emptyenv())

start_log <- function(log_dir, stem, prefix = "run_", extra = NULL) {
  # readr prints progress lines as messages, which the log captures; switch
  # them off.
  options(readr.show_progress = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(log_dir, paste0(prefix, stamp, ".txt"))
  con <- file(log_file, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  assign("con", con, envir = .helpers_log)
  assign("file", log_file, envir = .helpers_log)

  options(error = function() {
    cat("\nERROR: ", geterrmessage(), sep = "")
    if (requireNamespace("future", quietly = TRUE))
      try(future::plan(future::sequential), silent = TRUE)
    end_log("FAILED")
    if (!interactive()) quit(save = "no", status = 1)
  })

  cat("Script:      ", stem, "\n", sep = "")
  cat("Run started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat("Machine:     ", Sys.info()[["nodename"]], "\n", sep = "")
  cat("R:           ", R.version.string, "\n", sep = "")
  cat("Log file:    ", log_file, "\n", sep = "")
  if (length(extra)) cat(paste0(names(extra), ": ", unlist(extra), collapse = "\n"), "\n", sep = "")
  cat("\n")
  invisible(log_file)
}

# Prints the status and sessionInfo(), closes the sinks and the file, and puts
# the default error handler back. Safe to call twice.
end_log <- function(status = "finished") {
  con <- get0("con", envir = .helpers_log, inherits = FALSE)
  if (is.null(con) || !isTRUE(tryCatch(isOpen(con), error = function(e) FALSE))) return(invisible(NULL))
  cat("\nRun ", status, ": ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat("\n---- sessionInfo ----\n")
  print(sessionInfo())
  while (sink.number(type = "message") != 2L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(con)
  rm(list = "con", envir = .helpers_log)
  options(error = NULL)
  invisible(NULL)
}

# ---- Archive before writing ----------------------------------------------------

# Moves every file (not subfolder) directly inside each folder into
# <folder>/ARCHIVE_<tag>/. Pass subfolders such as Figures/ as their own
# entries in `dirs` if they hold outputs too. Stops if any move fails.
archive_previous <- function(dirs, tag = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  n_total <- 0L
  for (d in dirs) {
    if (!dir.exists(d)) next
    fs <- list.files(d, full.names = TRUE, recursive = FALSE)
    fs <- fs[!dir.exists(fs)]
    if (length(fs) == 0) next
    adir <- file.path(d, paste0("ARCHIVE_", tag))
    dir.create(adir, recursive = TRUE, showWarnings = FALSE)
    ok <- file.rename(fs, file.path(adir, basename(fs)))
    if (any(!ok)) {
      stop("Could not archive existing output(s): ",
           paste(basename(fs)[!ok], collapse = ", "),
           ". Refusing to run and risk overwriting them.")
    }
    cat("Archived", sum(ok), "file(s) to", adir, "\n")
    n_total <- n_total + sum(ok)
  }
  invisible(n_total)
}

# Moves only the listed files that exist, each into ARCHIVE_<tag>/ beside it.
# Use this when a folder also holds files the script must not move (for
# example inputs or checkpoints).
archive_files <- function(paths, tag = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  paths <- paths[file.exists(paths)]
  for (p in paths) {
    arch <- file.path(dirname(p), paste0("ARCHIVE_", tag))
    dir.create(arch, showWarnings = FALSE)
    ok <- file.rename(p, file.path(arch, basename(p)))
    if (!isTRUE(ok)) stop("Could not archive ", p, ". Stopping rather than overwrite.")
    cat("Archived: ", p, " -> ", arch, "\n", sep = "")
  }
  invisible(length(paths))
}

# ---- Parallel workers ----------------------------------------------------------

# The standard worker block: one physical core free, capped at 7, an override
# to force a count, and the counts printed so they reach the run log. Windows
# uses multisession, other systems multicore.
setup_workers <- function(override = NA_integer_, cap = 7L, set_plan = TRUE,
                          globals_max_gb = NULL) {
  n_phys <- parallel::detectCores(logical = FALSE)
  if (!is.finite(n_phys) || n_phys < 1) n_phys <- parallel::detectCores(logical = TRUE)
  if (!is.finite(n_phys) || n_phys < 1) n_phys <- 2L
  n_workers <- if (is.finite(override) && override >= 1) {
    as.integer(override)
  } else {
    max(1L, min(as.integer(cap), as.integer(n_phys) - 1L))
  }
  if (set_plan) {
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = n_workers)
    } else {
      future::plan(future::multicore, workers = n_workers)
    }
  }
  if (!is.null(globals_max_gb)) options(future.globals.maxSize = globals_max_gb * 1024^3)
  cat("Physical cores detected:", n_phys, " Workers:", n_workers,
      if (is.finite(override)) "(override)" else "", "\n")
  invisible(list(n_phys = as.integer(n_phys), n_workers = n_workers))
}

# ---- Warnings ------------------------------------------------------------------

# Runs `expr`, returns list(value, warnings). Use instead of suppressWarnings()
# around model fits, so the number and text of warnings can be logged.
with_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# ---- Monthly grid --------------------------------------------------------------

# Full monthly grid from `start` to `end`, with the series left-joined by `date`
# and sorted by date. With start and end left NULL it spans the series itself;
# with both given it is a common grid for all series. Months with no row in
# `df` get NA in every other column, keys included.
complete_monthly_grid <- function(df, start = NULL, end = NULL) {
  if (!"date" %in% names(df)) stop("complete_monthly_grid(): df has no 'date' column")
  if (is.null(start) || is.null(end)) {
    df <- dplyr::arrange(df, date)
    if (is.null(start)) start <- min(df$date, na.rm = TRUE)
    if (is.null(end))   end   <- max(df$date, na.rm = TRUE)
  }
  grid <- tibble::tibble(date = seq.Date(as.Date(start), as.Date(end), by = "month"))
  dplyr::arrange(dplyr::left_join(grid, df, by = "date"), date)
}

# ---- Fit checks ----------------------------------------------------------------

# For a KFAS::fitSSM() result: optimizer code, H, largest Q, variance collapse
# and log-likelihood, plus ok and a named reason, so failures can be counted
# by cause.
# v0: training variance; var_tol: VAR_COLLAPSE_TOL used by the script.
check_kfas_fit <- function(fit, v0, var_tol) {
  if (inherits(fit, "error") || is.null(fit)) {
    return(list(conv = NA_integer_, H = NA_real_, Qmax = NA_real_,
                var_collapsed = NA, loglik = NA_real_, ok = FALSE,
                reason = if (inherits(fit, "error")) paste("fit_error:", conditionMessage(fit)) else "fit_error"))
  }
  Hv <- suppressWarnings(as.numeric(fit$model$H)[1])
  Qv <- tryCatch({
    Qa <- fit$model$Q
    d  <- dim(Qa)[1]
    diag(matrix(as.numeric(Qa[, , 1]), nrow = d, ncol = d))
  }, error = function(e) NA_real_)
  vmax <- suppressWarnings(max(c(Hv, Qv), na.rm = TRUE))
  conv <- tryCatch(as.integer(fit$optim.out$convergence), error = function(e) NA_integer_)
  var_collapsed <- is.finite(vmax) && (vmax < var_tol * v0)
  reason <- if (is.na(conv)) "no_convergence_code" else
            if (conv != 0L) paste0("convergence=", conv) else
            if (var_collapsed) "variance_collapsed" else NA_character_
  list(
    conv = conv,
    H = Hv,
    Qmax = suppressWarnings(max(Qv, na.rm = TRUE)),
    var_collapsed = var_collapsed,
    loglik = tryCatch(as.numeric(-fit$optim.out$value), error = function(e) NA_real_),
    ok = is.na(reason),
    reason = reason
  )
}

# For forecast::Arima / auto.arima or stats::arima results: optimizer code
# (0 = converged), the chosen order, and finite coefficients and sigma2.
check_arima_fit <- function(fit) {
  if (inherits(fit, "error") || is.null(fit)) {
    return(list(code = NA_integer_, order = NA_character_, ok = FALSE,
                reason = if (inherits(fit, "error")) paste("fit_error:", conditionMessage(fit)) else "fit_error"))
  }
  code <- if (!is.null(fit$code)) as.integer(fit$code) else NA_integer_
  ord  <- tryCatch(paste(forecast::arimaorder(fit), collapse = ","), error = function(e) NA_character_)
  finite_ok <- all(is.finite(stats::coef(fit))) && is.finite(fit$sigma2)
  reason <- if (!is.na(code) && code != 0L) paste0("optim_code=", code) else
            if (!finite_ok) "nonfinite_coef_or_sigma2" else NA_character_
  list(code = code, order = ord, ok = is.na(reason), reason = reason)
}

# ---- Key column types ------------------------------------------------------------

# Sets the standard types on whichever key columns are present: estuary IDs,
# YEAR, MONTH, q, horizon, origin/test indices, fold ID and SPECIES_CODE to
# integer; GEAR to character "BS"/"BT"; date columns to Date. Bay value
# "ALL" (the coastwide/gamma marker) becomes NA_integer_, and a `scale`
# column ("gamma" / "alpha") is added if missing. Stops, rather than
# introducing NA, if a value cannot be converted exactly. Not used by the
# current pipeline scripts.
coerce_keys <- function(df, name = "data") {
  int_cols  <- c("BAY_ID", "BAY", "MAJOR_AREA", "YEAR", "MONTH", "q", "horizon",
                 "origin_i", "test_i", "fold", "fold_id", "SPECIES_CODE")
  date_cols <- c("date", "DATE", "test_date", "origin_date")

  to_int <- function(x, col) {
    if (is.integer(x)) return(x)
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) x <- trimws(x)
    num <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & (is.na(num) | abs(num - round(num)) > 1e-9)
    if (any(bad)) stop(name, ": column ", col, " has values that are not whole numbers: ",
                       paste(utils::head(unique(x[bad]), 5), collapse = ", "), call. = FALSE)
    as.integer(round(num))
  }

  for (bc in intersect(c("BAY_ID", "BAY"), names(df))) {
    if (is.character(df[[bc]]) || is.factor(df[[bc]])) {
      is_all <- !is.na(df[[bc]]) & as.character(df[[bc]]) == "ALL"
      if (any(is_all)) {
        if (!"scale" %in% names(df)) df$scale <- ifelse(is_all, "gamma", "alpha")
        df[[bc]] <- ifelse(is_all, NA, as.character(df[[bc]]))
      }
    }
  }
  for (col in intersect(int_cols, names(df))) df[[col]] <- to_int(df[[col]], col)

  if ("GEAR" %in% names(df)) {
    df$GEAR <- as.character(df$GEAR)
    bad <- setdiff(unique(df$GEAR[!is.na(df$GEAR)]), c("BS", "BT"))
    if (length(bad)) stop(name, ": GEAR has values other than BS and BT: ",
                          paste(bad, collapse = ", "), call. = FALSE)
  }
  for (col in intersect(date_cols, names(df))) {
    if (!inherits(df[[col]], "Date")) {
      d <- as.Date(df[[col]])
      if (any(is.na(d) & !is.na(df[[col]]))) stop(name, ": column ", col, " could not be read as dates", call. = FALSE)
      df[[col]] <- d
    }
  }
  df
}
