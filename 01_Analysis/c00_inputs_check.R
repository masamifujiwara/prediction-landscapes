# =============================================================================
# Script:  c00_inputs_check.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Checks the diversity series that c01, c02 and the model scripts read (the
#   filled series from b09) and stops if anything is wrong. Then writes the
#   missing-month table per series (months observed, filled, missing before
#   the first or after the last month, and gaps inside the series before
#   filling), a one-row summary per input, and a report of real but unusual
#   cells (report only; nothing is changed).
#
# Run after: b09
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   01_Analysis/Data/b9_forecast_settings/alpha_q012_filled.rds (b09, from b06)
#   01_Analysis/Data/b9_forecast_settings/gamma_q012_filled.rds (b09, from b07)
#
# Outputs:
#   02_Results/03_Diagnostics/c1_inputs_check/
#     c1_inputcheck__target_missingness.csv   one row per series
#     c1_inputcheck__alpha_q012__summary.csv  one row
#     c1_inputcheck__gamma_q012__summary.csv  one row
#     c1_inputcheck__qc_cells.csv             cells with q0 below Q0_LOW, and
#                                             cells flagged extrap_over_2n
#     c1_inputcheck__period_summary.csv       2010 to 2019 against 2020 to
#                                             2024, by scale and sampling method
#   01_Analysis/Logs/c1_inputs_check/run_<timestamp>.txt
#   Earlier outputs are moved to an ARCHIVE_<timestamp> subfolder first.
#
# Method notes:
#   - c00 fills nothing itself; the inputs are the b09 filled files.
#   - One series is one estuary and sampling method (bay level) or one
#     sampling method (coastwide); every check keys on the estuary as well.
#   - Checks (all failures are listed, then the script stops): sampling
#     methods and estuaries as in settings.rds; q0, q1, q2 present, finite and
#     above 0; no duplicate cell; months 1 to 12; no gap inside a series;
#     filled = TRUE exactly on settings$gap_months in every series (n_samples
#     0 there, above 0 elsewhere); q0 >= q1 >= q2.
#   - Missing months are counted per series, not per q (every q has the same
#     months). Months missing at the start and end are counted against the
#     common monthly grid, which runs from the first to the last month of any
#     series.
#   - Cells dominated by a single taxon (q0 below 2), cells extrapolated
#     beyond twice their catch, and the higher bag seine diversity in 2020 to
#     2024 are real values; c00 lists and summarizes them only.
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

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
script_stem <- "c1_inputs_check"   # output folder stem (short form c1_)

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "b9_forecast_settings", "alpha_q012_filled.rds")
in_gamma    <- file.path(PATHS$data_dir, "b9_forecast_settings", "gamma_q012_filled.rds")

Q0_LOW  <- 2                                     # QC report: cells with q0 below this
PERIODS <- list("2010-2019" = 2010:2019,         # QC report: the bag seine shift
                "2020-2024" = 2020:2024)         # (real values, report only)
Q_COLS  <- c("q0", "q1", "q2")

options(scipen = 999, warn = 1, width = 200)   # wide tables print on one line in the log

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("logs", "diag"))
archive_previous(dirs$diag)
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

# -----------------------------------------------------------------------------
# 4) Inputs
# -----------------------------------------------------------------------------
stop_if_missing(c(in_settings, in_alpha, in_gamma), "b09_forecast_settings.R")
settings <- load_settings(in_settings)
alpha <- readRDS(in_alpha)
gamma <- readRDS(in_gamma)
key_cols <- c("GEAR", "YEAR", "MONTH", "n_samples", Q_COLS, "filled")
expect_cols(alpha, c("BAY_ID", key_cols), "alpha_q012_filled")
expect_cols(gamma, key_cols, "gamma_q012_filled")
cat("\nalpha rows:", nrow(alpha), "| gamma rows:", nrow(gamma), "\n")

# Both scales in one table: one series = scale x BAY_ID x GEAR.
targets <- bind_rows(
  alpha %>% mutate(scale = "bay"),
  gamma %>% mutate(scale = "coastwide", BAY_ID = NA_integer_)
) %>%
  mutate(t = month_index(YEAR, MONTH)) %>%
  arrange(scale, BAY_ID, GEAR, t)
series_keys <- c("scale", "BAY_ID", "GEAR")

# -----------------------------------------------------------------------------
# 5) Checks. All failures are listed, then the script stops.
# -----------------------------------------------------------------------------
problems <- character(0)
fail <- function(...) problems <<- c(problems, paste0(...))

if (!setequal(unique(alpha$GEAR), settings$gears)) fail("alpha gears: ", paste(sort(unique(alpha$GEAR)), collapse = ", "))
if (!setequal(unique(gamma$GEAR), settings$gears)) fail("gamma gears: ", paste(sort(unique(gamma$GEAR)), collapse = ", "))
if (!setequal(unique(alpha$BAY_ID), settings$bays)) fail("alpha bays: ", paste(sort(unique(alpha$BAY_ID)), collapse = ", "),
                                                        "; settings: ", paste(settings$bays, collapse = ", "))
n_bay_gear <- nrow(distinct(alpha, BAY_ID, GEAR))
if (n_bay_gear != length(settings$bays) * length(settings$gears)) {
  fail("alpha has ", n_bay_gear, " bay x gear series; expected ", length(settings$bays) * length(settings$gears))
}
if (any(!targets$MONTH %in% 1:12)) fail("months outside 1 to 12: ", sum(!targets$MONTH %in% 1:12), " rows")

q_mat <- as.matrix(targets[, Q_COLS])
if (anyNA(q_mat)) fail("NA in q0, q1 or q2: ", sum(is.na(q_mat)), " values")
if (any(!is.finite(q_mat[!is.na(q_mat)]) | q_mat[!is.na(q_mat)] <= 0)) fail("q0, q1 or q2 not finite or not above 0")
bad_order <- with(targets, q1 > q0 + 1e-9 | q2 > q1 + 1e-9)
if (any(bad_order, na.rm = TRUE)) fail("q0 >= q1 >= q2 fails in ", sum(bad_order, na.rm = TRUE), " rows")

n_dup <- targets %>% count(across(all_of(series_keys)), YEAR, MONTH) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) fail("duplicate cells (scale, bay, gear, year, month): ", n_dup)

per_series <- targets %>%
  group_by(across(all_of(series_keys))) %>%
  summarise(
    complete       = all(diff(t) == 1L),
    filled_ok      = setequal(month_date(YEAR[filled], MONTH[filled]), settings$gap_months) &&
                     sum(filled) == length(settings$gap_months),
    n_samples_ok   = all(n_samples[filled] == 0L) && all(n_samples[!filled] > 0L),
    .groups = "drop"
  )
lab <- function(d) ifelse(is.na(d$BAY_ID), paste("coastwide", d$GEAR), paste("bay", d$BAY_ID, d$GEAR))
for (chk in c("complete", "filled_ok", "n_samples_ok")) {
  bad <- per_series[!per_series[[chk]], ]
  if (nrow(bad)) fail(chk, " fails in: ", paste(lab(bad), collapse = "; "))
}

if (length(problems)) {
  cat("\nInput check FAILED:\n", paste0("  - ", problems, collapse = "\n"), "\n", sep = "")
  stop("c00: ", length(problems), " input check(s) failed; see the list above.", call. = FALSE)
}
cat("\nAll input checks passed:", nrow(per_series), "series, both gears, bays",
    paste(settings$bays, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 6) Missing months per series
# -----------------------------------------------------------------------------
# Interior gaps are counted on the observed rows, that is before the b09 fill.
grid_first <- min(targets$t)
grid_last  <- max(targets$t)
fmt_t <- function(t) sprintf("%04d-%02d", t %/% 12L, t %% 12L + 1L)

missingness <- targets %>%
  group_by(across(all_of(series_keys))) %>%
  summarise(
    first_month          = fmt_t(min(t)),
    last_month           = fmt_t(max(t)),
    n_months             = n(),
    n_observed           = sum(!filled),
    n_filled             = sum(filled),
    filled_months        = paste(fmt_t(t[filled]), collapse = " "),
    n_leading_missing    = min(t) - grid_first,
    n_trailing_missing   = grid_last - max(t),
    n_interior_missing   = sum(diff(t[!filled]) - 1L),
    n_interior_gaps      = sum(diff(t[!filled]) > 1L),
    longest_interior_gap = max(c(0L, diff(t[!filled]) - 1L)),
    .groups = "drop"
  ) %>%
  mutate(grid = paste(fmt_t(grid_first), "to", fmt_t(grid_last)))

cat("\n=== Target missingness (common grid ", fmt_t(grid_first), " to ", fmt_t(grid_last), ") ===\n", sep = "")
print(as.data.frame(missingness %>% select(-grid)), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7) One-row summary per input
# -----------------------------------------------------------------------------
summarise_input <- function(d) {
  tibble::tibble(
    n_rows      = nrow(d),
    n_series    = nrow(distinct(d, across(all_of(series_keys)))),
    gears       = paste(sort(unique(d$GEAR)), collapse = " "),
    bays        = paste(sort(unique(d$BAY_ID)), collapse = " "),
    first_month = fmt_t(min(d$t)),
    last_month  = fmt_t(max(d$t)),
    n_observed  = sum(!d$filled),
    n_filled    = sum(d$filled),
    q0_min = min(d$q0), q0_median = stats::median(d$q0), q0_max = max(d$q0),
    q1_min = min(d$q1), q1_median = stats::median(d$q1), q1_max = max(d$q1),
    q2_min = min(d$q2), q2_median = stats::median(d$q2), q2_max = max(d$q2)
  )
}
sum_alpha <- summarise_input(filter(targets, scale == "bay"))
sum_gamma <- summarise_input(filter(targets, scale == "coastwide"))

# -----------------------------------------------------------------------------
# 8) QC report, report only
# -----------------------------------------------------------------------------
obs <- filter(targets, !filled)
qc_cols <- c("scale", "BAY_ID", "GEAR", "YEAR", "MONTH", "n_samples", "N_total", "S_obs",
             "m_at_C_target", Q_COLS)
qc_cells <- bind_rows(
  obs %>% filter(q0 < Q0_LOW) %>% mutate(check = paste0("q0 below ", Q0_LOW)),
  obs %>% filter(extrap_over_2n %in% TRUE) %>% mutate(check = "extrapolated beyond twice the catch")
) %>%
  select(check, any_of(qc_cols)) %>%
  arrange(check, scale, BAY_ID, GEAR, YEAR, MONTH)

period_summary <- bind_rows(lapply(names(PERIODS), function(p) {
  obs %>%
    filter(YEAR %in% PERIODS[[p]]) %>%
    group_by(scale, GEAR) %>%
    summarise(period = p,
              n_cells = n(),
              median_individuals_per_sample = stats::median(N_total / n_samples),
              median_S_obs = stats::median(S_obs),
              median_q0 = stats::median(q0), median_q1 = stats::median(q1), median_q2 = stats::median(q2),
              .groups = "drop")
})) %>%
  select(scale, GEAR, period, everything()) %>%
  arrange(scale, GEAR, period)

cat("\n=== QC report (real values, nothing changed) ===\n")
print(as.data.frame(qc_cells %>% count(check, scale, GEAR)), row.names = FALSE)
cat("\nPeriods (observed cells; individuals per sample = N_total / n_samples per cell):\n")
print(as.data.frame(period_summary %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)

# -----------------------------------------------------------------------------
# 9) Outputs
# -----------------------------------------------------------------------------
f <- function(x) file.path(dirs$diag, paste0("c1_inputcheck__", x))
out <- c(missingness = f("target_missingness.csv"),
         alpha       = f("alpha_q012__summary.csv"),
         gamma       = f("gamma_q012__summary.csv"),
         qc          = f("qc_cells.csv"),
         periods     = f("period_summary.csv"))
readr::write_csv(missingness,    out[["missingness"]], na = "")
readr::write_csv(sum_alpha,      out[["alpha"]],       na = "")
readr::write_csv(sum_gamma,      out[["gamma"]],       na = "")
readr::write_csv(qc_cells,       out[["qc"]],          na = "")
readr::write_csv(period_summary, out[["periods"]],     na = "")
for (p in out) cat("Saved:", p, "\n")

cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("c00_inputs_check.R completed\n")
