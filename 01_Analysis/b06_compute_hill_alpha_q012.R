# =============================================================================
# Script:  b06_compute_hill_alpha_q012.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Computes estuary-scale (alpha) diversity as coverage-standardized Hill
#   numbers (q = 0 richness, q = 1 Shannon type, q = 2 Simpson type) with iNEXT
#   for every BAY_ID x GEAR x YEAR x MONTH cell with sampling effort, and
#   writes the cell coverage diagnostics and an example plot.
#
# Run after: b03, b04, b05, lib/hill.R.
#
# Inputs:
#   01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds, bt_events.rds (b04)
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds (b05):
#     keys$valid_bays
#   01_Analysis/lib/hill.R (shared Hill-number functions)
#
# Outputs:
#   01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds, alpha_q012.csv
#     (the csv without the abund_vec list column); read by b09, f07, h12
#   02_Results/03_Diagnostics/b6_compute_hill_alpha_q012/
#     b6_alpha_coverage_diagnostics__cells.csv
#     b6_alpha_cells_that_fail_estimateD.csv (only if a cell fails)
#     b6_alpha_q012_timeseries_examples.png
#   Run log: 01_Analysis/Logs/b6_compute_hill_alpha_q012/run_<timestamp>.txt
#   Earlier outputs in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> first. Output folders use the stem
#   b6_compute_hill_alpha_q012.
#
# Method notes:
#   - All cells are estimated at one common sample coverage, C_target: the
#     10th percentile of the cells' observed coverage, capped at C_cap (0.90).
#     A lower target that would let every cell stay within twice its catch was
#     considered and not used: it would change every value to accommodate a
#     few cells, while most cells support a much higher coverage.
#   - Placeholder rows (SPECIES_CODE = NA, catch = 0) are left out of the
#     abundance vectors. A cell with effort but no catch gets q0 = q1 = q2 = 0
#     and Chat_obs = 0.
#   - iNEXT::estimateD(datatype = "abundance") needs whole-number counts, so
#     pooled catches are rounded (catches are already whole numbers).
#   - Coverage (Chat_obs) uses the Chao and Jost (2012) correction when no
#     species is caught twice.
#   - No bootstrap by default (NBOOT = 0): the diversity values do not depend
#     on it, and nothing downstream reads its confidence limits. Set NBOOT to
#     a positive number to add lcl_q0 to ucl_q2 (serial and slow).
#   - extrap_over_2n marks cells whose target coverage needs a sample of more
#     than twice the observed catch (round(m_at_C_target) > 2 x N_total, the
#     rule iNEXT uses for its warning), where q0 may be biased (Chao et al.
#     2014). These cells are kept, listed in the log and reported in Online
#     Resource 1; they are not set to missing.
#   - All taxa that pass a00 are kept; there is no screen by how many samples a
#     taxon occurs in, because such a screen would bias the coverage estimate.
#   - A cell that iNEXT cannot estimate is listed, and the script stops at the
#     end, after trying every cell.
#   - Effort weighting does not arise at estuary scale; see b07 for the coast.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

set.seed(101)

msg <- function(...) cat(paste0(..., collapse = ""), "\n")

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

options(scipen = 999, warn = 1)

# -----------------------------------------------------------------------------
# 2) Packages
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
# 3) Parameters
# -----------------------------------------------------------------------------
# Stem used for the output folders.
script_stem <- "b6_compute_hill_alpha_q012"

# Bootstrap resamples per cell for the confidence limits. 0 = none (default):
# the diversity values do not depend on it.
NBOOT <- 0
nboot <- NBOOT
conf  <- 0.95

C_cap      <- 0.90
C_quantile <- 0.10
C_min      <- 0.50

example_bays <- c(2, 5, 8)  # Galveston, Aransas, Lower Laguna Madre (plot only)

# If TRUE, stop if any non-integer abundances are detected before rounding
strict_noninteger_check <- FALSE

# -----------------------------------------------------------------------------
# 4) Paths + load derived inputs
# -----------------------------------------------------------------------------
in_bs_events <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bs_events.rds")
in_bt_events <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bt_events.rds")
in_station   <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_keys      <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.rds")

dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
msg("Project directory: ", PATHS$project_dir)

stop_if_missing(c(in_bs_events, in_bt_events, in_station, in_keys),
                c("b04_join_organisms_to_station.R", "b04_join_organisms_to_station.R",
                  "b03_harmonize_station_metadata.R", "b05_define_spatiotemporal_units.R"))

bs_events <- readRDS(in_bs_events)
bt_events <- readRDS(in_bt_events)
station_core <- readRDS(in_station)
index_keys <- readRDS(in_keys)

station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays
if (is.null(valid_bays) || length(valid_bays) == 0) {
  stop("index_keys.rds has no keys$valid_bays (rerun b05_define_spatiotemporal_units.R)", call. = FALSE)
}
msg("Bays: ", paste(valid_bays, collapse = " "))

# Station effort table for alpha grid (includes zero-catch stations, excludes no-effort months)
station_all <- bind_rows(station_bs_core, station_bt_core) %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  ) %>%
  filter(BAY_ID %in% valid_bays) %>%
  select(SAMPLE_ID, GEAR, BAY_ID, YEAR, MONTH)

if (anyDuplicated(station_bs_core$SAMPLE_ID) || anyDuplicated(station_bt_core$SAMPLE_ID)) {
  stop("Duplicated SAMPLE_ID in station_core", call. = FALSE)
}

out_dir  <- dirs$data
diag_dir <- dirs$diag

# -----------------------------------------------------------------------------
# 5) Build pooled abundance vectors per alpha cell (BAY_ID x GEAR x YEAR x MONTH)
#    IMPORTANT: remove sentinel rows (SPECIES_CODE = NA) before pooling
# -----------------------------------------------------------------------------
prep_events <- function(events_df, gear_label, valid_bays) {
  events_df %>%
    mutate(
      GEAR  = gear_label,
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

cell_species <- events_all %>%
  filter(!is.na(SPECIES_CODE)) %>%   # sentinel removed here
  group_by(BAY_ID, GEAR, YEAR, MONTH, SPECIES_CODE) %>%
  summarise(abund = sum(catch, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(abund) & abund > 0)

cell_n_samples <- station_all %>%
  count(BAY_ID, GEAR, YEAR, MONTH, name = "n_samples")

alpha_grid <- cell_n_samples %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH)) %>%
  arrange(BAY_ID, GEAR, YEAR, MONTH)

# -----------------------------------------------------------------------------
# 6) Utilities: safe_chat(), sanitize_abundance() and
#    estimate_q012_at_coverage() are in lib/hill.R (sourced above).
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 7) Compute observed coverage per cell, define global C_target
# -----------------------------------------------------------------------------
cell_vectors <- cell_species %>%
  group_by(BAY_ID, GEAR, YEAR, MONTH) %>%
  summarise(
    abund_vec = list(setNames(abund, SPECIES_CODE)),
    N_total   = sum(abund),
    S_obs     = n(),
    .groups   = "drop"
  )

alpha_cells <- alpha_grid %>%
  left_join(cell_vectors, by = c("BAY_ID", "GEAR", "YEAR", "MONTH")) %>%
  mutate(
    N_total = coalesce(N_total, 0),
    S_obs   = coalesce(S_obs, 0)
  ) %>%
  mutate(
    abund_vec = map2(abund_vec, N_total, ~{
      if (is.null(.x) || is.na(.y) || .y == 0) setNames(numeric(0), character(0)) else .x
    })
  ) %>%
  mutate(Chat_obs = map_dbl(abund_vec, safe_chat))

chat_nonempty <- alpha_cells$Chat_obs[alpha_cells$N_total > 0]
if (length(chat_nonempty) == 0) stop("No non-empty cells found. Check input events.", call. = FALSE)

C_target <- min(C_cap, as.numeric(stats::quantile(chat_nonempty, probs = C_quantile, na.rm = TRUE)))
C_target <- max(min(C_target, 0.999), 0.001)

msg("Coverage target (C_target): ", sprintf("%.3f", C_target),
    " (cap=", C_cap, ", quantile=", C_quantile, ")")

alpha_cells <- alpha_cells %>%
  mutate(
    extrapolated  = (N_total > 0) & (C_target > Chat_obs),
    coverage_fail = (N_total == 0) | (Chat_obs < C_min),
    C_target      = C_target
  )

# -----------------------------------------------------------------------------
# 7b) Diagnostics: non-integer check (informational)
# -----------------------------------------------------------------------------
nonint_summary <- cell_vectors %>%
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
       "or revisit what 'catch' represents.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 8) Estimate alpha q0/q1/q2 at C_target per cell
# -----------------------------------------------------------------------------
msg("\n=== Estimating alpha diversity per cell (nboot=", nboot, ") ===")

est_list <- vector("list", nrow(alpha_cells))
fail_msg <- rep(NA_character_, nrow(alpha_cells))
for (i in seq_len(nrow(alpha_cells))) {
  est_list[[i]] <- tryCatch(
    estimate_q012_at_coverage(
      x = alpha_cells$abund_vec[[i]],
      C_target = C_target,
      nboot = nboot,
      conf = conf,
      strict_noninteger_check = strict_noninteger_check
    ),
    error = function(e) { fail_msg[i] <<- conditionMessage(e); NULL }
  )
  if (i %% 500 == 0) msg("... estimate progress: ", i, "/", nrow(alpha_cells))
}

# Every cell has been tried; stop now if any failed, with the list saved.
if (any(!is.na(fail_msg))) {
  bad <- alpha_cells[!is.na(fail_msg), ] %>%
    select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs, C_target) %>%
    mutate(error = fail_msg[!is.na(fail_msg)])
  bad_path <- file.path(diag_dir, "b6_alpha_cells_that_fail_estimateD.csv")
  readr::write_csv(bad, bad_path)
  stop(nrow(bad), " cell(s) failed in iNEXT::estimateD(); list saved to ", bad_path, call. = FALSE)
}

alpha_est <- bind_cols(alpha_cells, bind_rows(est_list)) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH),
         # Same rule as iNEXT's warning: round(m) > 2 x observed catch.
         extrap_over_2n = !is.na(m_at_C_target) & round(m_at_C_target) > 2 * N_total)

# Cells where reaching C_target needs more than twice the observed catch.
# iNEXT warns that q0 may be biased there; they are kept and reported.
over2n <- alpha_est %>% filter(extrap_over_2n)
msg("\nCells extrapolated to more than twice the observed catch: ", nrow(over2n),
    " of ", nrow(alpha_est), " (kept; q0 less reliable)")
if (nrow(over2n) > 0) {
  print(as.data.frame(over2n %>%
    transmute(BAY_ID, GEAR, YEAR, MONTH, N_total, S_obs, Chat_obs = round(Chat_obs, 3),
              m_at_C_target = round(m_at_C_target), times_catch = round(m_at_C_target / N_total, 2),
              q0 = round(q0, 1))), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 9) Diagnostics summary + cell CSV
# -----------------------------------------------------------------------------
diag_summary <- alpha_est %>%
  summarise(
    n_cells       = n(),
    n_empty       = sum(N_total == 0),
    prop_empty    = mean(N_total == 0),
    n_fail_cov    = sum(coverage_fail),
    prop_fail_cov = mean(coverage_fail),
    n_extrap_over_2n = sum(extrap_over_2n),
    median_chat   = stats::median(Chat_obs),
    p10_chat      = as.numeric(stats::quantile(Chat_obs, 0.10)),
    p90_chat      = as.numeric(stats::quantile(Chat_obs, 0.90))
  )

msg("\n=== Coverage diagnostics summary (alpha cells) ===")
print(diag_summary)

diag_cells_path <- file.path(diag_dir, "b6_alpha_coverage_diagnostics__cells.csv")
readr::write_csv(
  alpha_est %>%
    select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs,
           C_target, extrapolated, coverage_fail, m_at_C_target, extrap_over_2n,
           q0, q1, q2, any_of(c("lcl_q0", "lcl_q1", "lcl_q2", "ucl_q0", "ucl_q1", "ucl_q2"))),
  diag_cells_path
)
msg("Saved: ", diag_cells_path)

# -----------------------------------------------------------------------------
# 10) Quick plot examples for a few bays (q0/q1/q2)
# -----------------------------------------------------------------------------
plot_df <- alpha_est %>%
  filter(BAY_ID %in% example_bays) %>%
  select(BAY_ID, GEAR, YEAR_MONTH, q0, q1, q2) %>%
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "alpha") %>%
  mutate(date = as.Date(paste0(YEAR_MONTH, "-01")))

p <- ggplot(plot_df, aes(x = date, y = alpha)) +
  geom_line() +
  facet_grid(q ~ BAY_ID + GEAR, scales = "free_y") +
  labs(
    x = "Year-month",
    y = "Coverage-standardized alpha diversity",
    title = paste0("Alpha diversity time series examples (C_target=", sprintf("%.3f", C_target), ")")
  )

plot_path <- file.path(diag_dir, "b6_alpha_q012_timeseries_examples.png")
ggsave(plot_path, p, width = 12, height = 8, dpi = 300)
msg("Saved: ", plot_path)

# -----------------------------------------------------------------------------
# 11) Save outputs -> 01_Analysis/Data/<script_stem>/
# -----------------------------------------------------------------------------
out_rds <- file.path(out_dir, "alpha_q012.rds")
out_csv <- file.path(out_dir, "alpha_q012.csv")

# ---- Checks before writing ----
for (g in c("BS", "BT")) {
  bays_g <- sort(unique(alpha_est$BAY_ID[alpha_est$GEAR == g]))
  if (!identical(as.integer(bays_g), as.integer(sort(valid_bays)))) {
    stop("Gear ", g, " has bays ", paste(bays_g, collapse = " "), "; expected ",
         paste(sort(valid_bays), collapse = " "), call. = FALSE)
  }
}
nonempty <- alpha_est$N_total > 0
if (anyNA(alpha_est$q0[nonempty]) || anyNA(alpha_est$q1[nonempty]) || anyNA(alpha_est$q2[nonempty])) {
  stop("Missing q0, q1 or q2 in non-empty cells", call. = FALSE)
}
msg("Cells: ", nrow(alpha_est), " (", sum(!nonempty), " with effort but no catch)")

saveRDS(alpha_est, out_rds)
readr::write_csv(alpha_est %>% select(-abund_vec), out_csv)  # a list column cannot go in a csv

msg("Saved: ", out_rds)
msg("Saved: ", out_csv)

msg("\nNOTE:")
msg("  - Sentinel rows (SPECIES_CODE = NA, catch = 0) are excluded from abundance vectors.")
msg("  - Cells with sampling effort but N_total = 0 are kept and assigned q0=q1=q2=0, Chat_obs=0.")
msg("  - Abundances are rounded to integers before iNEXT::estimateD().")
msg("  - Bootstrap resamples per cell (NBOOT): ", NBOOT, if (NBOOT == 0) " (no confidence limits)" else "")
end_log("finished")
msg("b06_compute_hill_alpha_q012.R completed")
