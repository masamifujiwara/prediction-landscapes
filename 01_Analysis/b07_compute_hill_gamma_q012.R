# =============================================================================
# Script:  b07_compute_hill_gamma_q012.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Computes coastwide (gamma) diversity as coverage-standardized Hill numbers
#   (q = 0, 1, 2) with iNEXT for every GEAR x YEAR x MONTH cell, pooling the
#   catch of all estuaries within each sampling method, and writes the cell
#   coverage diagnostics and a plot of the series. The method is the same as
#   in b06 (estuary scale), with the Hill-number functions shared in
#   lib/hill.R.
#
# Run after: b03, b04, b05, lib/hill.R.
#
# Inputs:
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds, bt_events.rds (b04)
#   01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds (b05):
#     keys$valid_bays
#   01_Analysis/lib/hill.R (shared Hill-number functions)
#
# Outputs:
#   01_Analysis/Data/01_TimeSeries/gamma_q012.rds, gamma_q012.csv
#     (the csv without the abund_vec list column); read by b09, f07. The folder
#     name does not follow the script stem.
#   02_Results/03_Diagnostics/b7_compute_hill_gamma_q012/
#     b7_gamma_coverage_diagnostics__cells.csv
#     b7_gamma_cells_that_fail_estimateD.csv (only if a cell fails)
#     b7_gamma_q012_timeseries__by_gear.png
#   Run log: 01_Analysis/Logs/b7_compute_hill_gamma_q012/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> first (in 01_TimeSeries
#   only gamma_q012.rds and .csv are moved).
#
# Method notes:
#   - Same Hill-number method as b06: one common coverage target, C_target =
#     10th percentile of the cells' observed coverage, capped at C_cap (0.90);
#     Chao and Jost (2012) coverage correction when no species is caught
#     twice; no bootstrap by default (NBOOT = 0), because the diversity values
#     do not depend on it and nothing downstream reads its confidence limits.
#   - Placeholder rows (SPECIES_CODE = NA, catch = 0) are left out of the
#     abundance vectors. A cell with effort but no catch gets q0 = q1 = q2 = 0
#     and Chat_obs = 0. Pooled catches are rounded to whole numbers for iNEXT
#     (catches are already whole numbers).
#   - A coastwide cell is kept only when every estuary was sampled with that
#     sampling method in that month (n_bays); otherwise it describes only part
#     of the coast. This drops the first months of the bay trawl record, before
#     every estuary was trawled. The script stops if a month inside the record
#     lacks an estuary.
#   - Estuaries are pooled by sampling effort, not with equal weight: every
#     haul or tow of the sampling method in the month counts once. Bag seine
#     effort is equal across estuaries. In the bay trawl, the two Laguna Madre
#     estuaries (bays 7 and 8, 10 tows per month) count half as much as the
#     others (20); TPWD sets effort roughly by estuary size. The shares are
#     constant over the record, so the weighting cannot create changes over
#     time. An equal-weight version (the same number of random tows per
#     estuary, averaged over many draws) is a possible sensitivity check and
#     is not used here.
#   - extrap_over_2n marks cells whose target coverage needs more than twice
#     the observed catch (round(m_at_C_target) > 2 x N_total, iNEXT's warning
#     rule); such cells are kept and reported, as in b06.
#   - A cell that iNEXT cannot estimate is listed, and the script stops at the
#     end, after trying every cell.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

set.seed(102)

# -----------------------------------------------------------------------------
# 0) Script directory (RStudio / Rscript / source / nohup) + setwd()
# -----------------------------------------------------------------------------
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
source(file.path(PATHS$analysis_dir, "lib", "hill.R"))


# -----------------------------------------------------------------------------
# 2) Options + logging helpers (nohup-safe; no emoji by default)
# -----------------------------------------------------------------------------
options(scipen = 999, warn = 1)

msg <- function(...) cat(paste0(..., collapse = ""), "\n")

# -----------------------------------------------------------------------------
# 3) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "ggplot2", "iNEXT")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 4) Parameters
# -----------------------------------------------------------------------------
# Stem used for the log and diagnostics folders.
script_stem <- "b7_compute_hill_gamma_q012"

# Bootstrap resamples per cell for the confidence limits. 0 = none (default).
NBOOT <- 0
nboot <- NBOOT
conf  <- 0.95

C_cap      <- 0.90
C_quantile <- 0.10
C_min      <- 0.50  # diagnostics flag only

strict_noninteger_check <- FALSE

# -----------------------------------------------------------------------------
# 5) Paths + load inputs (from 01_Analysis/Data/<script_stem>/...)
# -----------------------------------------------------------------------------
in_station_core <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_bs_events    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bs_events.rds")
in_bt_events    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bt_events.rds")
in_index_keys   <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.rds")

# Derived output folder for the coastwide time series (does not follow the stem; see header)
ts_dir <- file.path(PATHS$data_dir, "01_TimeSeries")
dir.create(ts_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(ts_dir, "gamma_q012.rds")
out_csv <- file.path(ts_dir, "gamma_q012.csv")

dirs <- script_dirs(script_stem, c("logs", "diag"))
archive_files(c(out_rds, out_csv))
archive_previous(dirs$diag)
start_log(dirs$logs, script_stem)
msg("Project directory: ", PATHS$project_dir)

stop_if_missing(c(in_station_core, in_bs_events, in_bt_events, in_index_keys),
                c("b03_harmonize_station_metadata.R", "b04_join_organisms_to_station.R",
                  "b04_join_organisms_to_station.R", "b05_define_spatiotemporal_units.R"))

station_core <- readRDS(in_station_core)
bs_events    <- readRDS(in_bs_events)
bt_events    <- readRDS(in_bt_events)
index_keys   <- readRDS(in_index_keys)

station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays
if (is.null(valid_bays) || length(valid_bays) == 0) {
  stop("index_keys.rds has no keys$valid_bays (rerun b05_define_spatiotemporal_units.R)", call. = FALSE)
}
msg("Bays pooled: ", paste(valid_bays, collapse = " "))

diag_dir <- dirs$diag

# -----------------------------------------------------------------------------
# 6) Build gamma "effort grid" (gear x year x month) from station_core (keeps effort-only cells)
# -----------------------------------------------------------------------------
station_all <- bind_rows(station_bs_core, station_bt_core) %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  ) %>%
  filter(BAY_ID %in% valid_bays) %>%
  select(SAMPLE_ID, GEAR, BAY_ID, YEAR, MONTH)

gamma_grid <- station_all %>%
  group_by(GEAR, YEAR, MONTH) %>%
  summarise(n_samples = n(), n_bays = n_distinct(BAY_ID), .groups = "drop") %>%
  arrange(GEAR, YEAR, MONTH) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

# Keep a coastwide cell only when every bay was sampled with that gear in that
# month; otherwise it describes part of the coast. Only months before a gear's
# first complete month may be dropped; a gap inside the record stops the script.
incomplete <- gamma_grid %>% filter(n_bays < length(valid_bays))
if (nrow(incomplete) > 0) {
  first_full <- gamma_grid %>% filter(n_bays == length(valid_bays)) %>%
    group_by(GEAR) %>% summarise(first_ym = min(YEAR_MONTH), .groups = "drop")
  inside <- incomplete %>% left_join(first_full, by = "GEAR") %>% filter(YEAR_MONTH > first_ym)
  if (nrow(inside) > 0) {
    stop("Coastwide months inside the record with fewer than all bays sampled: ",
         paste(inside$GEAR, inside$YEAR_MONTH, collapse = ", "), call. = FALSE)
  }
  msg("Coastwide cells dropped (fewer than ", length(valid_bays), " bays sampled): ",
      paste0(incomplete$GEAR, " ", incomplete$YEAR_MONTH, " (", incomplete$n_bays, " bays)", collapse = ", "))
  gamma_grid <- gamma_grid %>% filter(n_bays == length(valid_bays))
}

# -----------------------------------------------------------------------------
# 7) Utilities: safe_chat(), sanitize_abundance() and
#    estimate_q012_at_coverage() are in lib/hill.R (sourced above).
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 8) Prepare events and pool abundances coastwide by (gear x year x month)
#     IMPORTANT: remove sentinel NA species rows before pooling
# -----------------------------------------------------------------------------
prep_events <- function(events_df, gear_label, valid_bays) {
  events_df %>%
    mutate(
      GEAR = gear_label,
      BAY_ID = as.integer(MAJOR_AREA),
      YEAR  = as.integer(YEAR_i),
      MONTH = as.integer(MONTH_i),
      SPECIES_CODE = as.character(SPECIES_CODE),
      catch = as.numeric(catch)
    ) %>%
    filter(BAY_ID %in% valid_bays)
}

bs_events2 <- prep_events(bs_events, "BS", valid_bays)
bt_events2 <- prep_events(bt_events, "BT", valid_bays)
events_all <- bind_rows(bs_events2, bt_events2)

gamma_species <- events_all %>%
  filter(!is.na(SPECIES_CODE)) %>%  # sentinel removed
  group_by(GEAR, YEAR, MONTH, SPECIES_CODE) %>%
  summarise(abund = sum(catch, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(abund) & abund > 0)

gamma_vectors <- gamma_species %>%
  group_by(GEAR, YEAR, MONTH) %>%
  summarise(
    abund_vec = list(setNames(abund, SPECIES_CODE)),
    N_total   = sum(abund),
    S_obs     = n(),
    .groups   = "drop"
  )

gamma_cells <- gamma_grid %>%
  left_join(gamma_vectors, by = c("GEAR", "YEAR", "MONTH")) %>%
  mutate(
    N_total = coalesce(N_total, 0),
    S_obs   = coalesce(S_obs, 0)
  ) %>%
  mutate(
    abund_vec = map2(abund_vec, N_total, ~{
      if (is.null(.x) || is.na(.y) || .y == 0) setNames(numeric(0), character(0)) else .x
    })
  )

gamma_cells <- gamma_cells %>%
  mutate(Chat_obs = map_dbl(abund_vec, safe_chat))

chat_nonempty <- gamma_cells$Chat_obs[gamma_cells$N_total > 0]
if (length(chat_nonempty) == 0) stop("No non-empty gamma cells found. Check input events.", call. = FALSE)

C_target <- min(C_cap, as.numeric(stats::quantile(chat_nonempty, probs = C_quantile, na.rm = TRUE)))
C_target <- max(min(C_target, 0.999), 0.001)

msg("Coverage target (C_target): ", sprintf("%.3f", C_target),
    " (cap=", C_cap, ", quantile=", C_quantile, ")")

gamma_cells <- gamma_cells %>%
  mutate(
    extrapolated  = (N_total > 0) & (C_target > Chat_obs),
    coverage_fail = (N_total == 0) | (Chat_obs < C_min),
    C_target      = C_target
  )

# -----------------------------------------------------------------------------
# 9) Diagnostics: non-integer check
# -----------------------------------------------------------------------------
nonint_summary <- gamma_vectors %>%
  mutate(any_nonint = map_lgl(abund_vec, ~{
    v <- as.numeric(.x)
    v <- v[is.finite(v) & v > 0]
    if (length(v) == 0) return(FALSE)
    any(abs(v - round(v)) > 1e-8)
  })) %>%
  summarise(
    n_cells = n(),
    n_any_nonint = sum(any_nonint),
    prop_any_nonint = mean(any_nonint)
  )

msg("\n=== Abundance integrality check (before rounding) ===")
print(nonint_summary)

if (isTRUE(strict_noninteger_check) && nonint_summary$n_any_nonint > 0) {
  stop("Non-integer pooled abundances detected. Set strict_noninteger_check=FALSE to allow rounding, ",
       "or revisit what 'catch' represents in events.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 10) Estimate gamma q0/q1/q2 per cell
# -----------------------------------------------------------------------------
msg("\n=== Estimating gamma diversity per cell (nboot=", nboot, ") ===")

est_list <- vector("list", nrow(gamma_cells))
fail_msg <- rep(NA_character_, nrow(gamma_cells))
for (i in seq_len(nrow(gamma_cells))) {
  est_list[[i]] <- tryCatch(
    estimate_q012_at_coverage(
      x = gamma_cells$abund_vec[[i]],
      C_target = C_target,
      nboot = nboot,
      conf = conf,
      strict_noninteger_check = strict_noninteger_check
    ),
    error = function(e) { fail_msg[i] <<- conditionMessage(e); NULL }
  )
  if (i %% 200 == 0) msg("... estimate progress: ", i, "/", nrow(gamma_cells))
}

# Every cell has been tried; stop now if any failed, with the list saved.
if (any(!is.na(fail_msg))) {
  bad <- gamma_cells[!is.na(fail_msg), ] %>%
    select(GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs, C_target) %>%
    mutate(error = fail_msg[!is.na(fail_msg)])
  bad_path <- file.path(diag_dir, "b7_gamma_cells_that_fail_estimateD.csv")
  readr::write_csv(bad, bad_path)
  stop(nrow(bad), " cell(s) failed in iNEXT::estimateD(); list saved to ", bad_path, call. = FALSE)
}

gamma_est <- bind_cols(gamma_cells, bind_rows(est_list)) %>%
  # Same rule as iNEXT's warning: round(m) > 2 x observed catch.
  mutate(extrap_over_2n = !is.na(m_at_C_target) & round(m_at_C_target) > 2 * N_total)

over2n <- gamma_est %>% filter(extrap_over_2n)
msg("\nCells extrapolated to more than twice the observed catch: ", nrow(over2n),
    " of ", nrow(gamma_est), " (kept; q0 less reliable)")
if (nrow(over2n) > 0) {
  print(as.data.frame(over2n %>%
    transmute(GEAR, YEAR, MONTH, N_total, S_obs, Chat_obs = round(Chat_obs, 3),
              m_at_C_target = round(m_at_C_target), times_catch = round(m_at_C_target / N_total, 2),
              q0 = round(q0, 1))), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 11) Diagnostics outputs (ALL under 02_Results/03_Diagnostics/<script_stem>/)
# -----------------------------------------------------------------------------
diag_cells_path <- file.path(diag_dir, "b7_gamma_coverage_diagnostics__cells.csv")
readr::write_csv(
  gamma_est %>%
    select(GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs,
           C_target, extrapolated, coverage_fail, m_at_C_target, extrap_over_2n,
           q0, q1, q2, any_of(c("lcl_q0", "lcl_q1", "lcl_q2", "ucl_q0", "ucl_q1", "ucl_q2"))),
  diag_cells_path
)
msg("Saved: ", diag_cells_path)

diag_summary <- gamma_est %>%
  summarise(
    n_cells       = n(),
    n_empty       = sum(N_total == 0),
    prop_empty    = mean(N_total == 0),
    n_fail_cov    = sum(coverage_fail),
    prop_fail_cov = mean(coverage_fail),
    n_extrap_over_2n = sum(extrap_over_2n),
    median_chat   = median(Chat_obs),
    p10_chat      = as.numeric(quantile(Chat_obs, 0.10)),
    p90_chat      = as.numeric(quantile(Chat_obs, 0.90))
  )

msg("\n=== Coverage diagnostics summary (gamma cells) ===")
print(diag_summary)

plot_df <- gamma_est %>%
  select(GEAR, YEAR_MONTH, q0, q1, q2) %>%
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "gamma") %>%
  mutate(date = as.Date(paste0(YEAR_MONTH, "-01")))

p <- ggplot(plot_df, aes(x = date, y = gamma)) +
  geom_line() +
  facet_grid(q ~ GEAR, scales = "free_y") +
  labs(
    x = "Year-month",
    y = "Coverage-standardized gamma diversity",
    title = paste0("Gamma diversity time series (C_target=", sprintf("%.3f", C_target), ")")
  )

plot_path <- file.path(diag_dir, "b7_gamma_q012_timeseries__by_gear.png")
ggsave(plot_path, p, width = 12, height = 7, dpi = 300)
msg("Saved: ", plot_path)

# -----------------------------------------------------------------------------
# 12) Save outputs (derived data)
# -----------------------------------------------------------------------------
# ---- Checks before writing ----
if (!setequal(unique(gamma_est$GEAR), c("BS", "BT"))) stop("Expected gears BS and BT", call. = FALSE)
nonempty <- gamma_est$N_total > 0
if (anyNA(gamma_est$q0[nonempty]) || anyNA(gamma_est$q1[nonempty]) || anyNA(gamma_est$q2[nonempty])) {
  stop("Missing q0, q1 or q2 in non-empty cells", call. = FALSE)
}
msg("Cells: ", nrow(gamma_est), " (", sum(!nonempty), " with effort but no catch)")

saveRDS(gamma_est, out_rds)
readr::write_csv(gamma_est %>% select(-abund_vec), out_csv)  # a list column cannot go in a csv

msg("Saved: ", out_rds)
msg("Saved: ", out_csv)

msg("\nNOTE:")
msg("  - Sentinel rows (SPECIES_CODE = NA, catch = 0) are excluded from abundance vectors.")
msg("  - Time cells with effort but N_total = 0 are kept and assigned q0=q1=q2=0, Chat_obs=0.")
msg("  - Abundances are rounded to integers before iNEXT::estimateD().")
msg("  - Bootstrap resamples per cell (NBOOT): ", NBOOT, if (NBOOT == 0) " (no confidence limits)" else "")
end_log("finished")
msg("b07_compute_hill_gamma_q012.R completed")
