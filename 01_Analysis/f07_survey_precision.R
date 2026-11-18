# =============================================================================
# Script:  f07_survey_precision.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Sampling uncertainty of every monthly value (biodiversity and abundance,
#   each estuary and coastwide, each sampling method): how much a month's
#   value depends on which hauls or tows happened to be taken. Gives the
#   numbers behind the tolerance (Methods; Online Resource 1, Sect. 1.10):
#   the tolerance (settings$tolerance, a share of the series mean for
#   biodiversity and log units for abundance) was set at about twice the
#   typical sampling uncertainty of one estuary and month, and this script
#   computes that uncertainty within the pipeline on the current data.
#
# Run after: a00, b03, b04, b06, b07, b09, c02, e02, lib/helpers.R,
#   lib/forecast_rules.R, lib/hill.R (independent of the forecasts)
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   - 01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   - 01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds, bt_events.rds (b04)
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds (b06; check)
#   - 01_Analysis/Data/01_TimeSeries/gamma_q012.rds (b07; check)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds,
#     gamma_q012_deseas_scaled.rds (c02; series means)
#   - 01_Analysis/Data/e2_aggregate_abundance/alpha_abundance_targets.rds,
#     gamma_abundance_targets.rds (e02; check)
#   - 01_Analysis/Data.RData (a00; catches, aggregated as in e02)
#
# Outputs:
#   01_Analysis/Data/f7_survey_precision/, each table as .rds and .csv:
#   f7_cell_precision     one row per series and month: n_samples, value,
#                         sd_sampling, sd_pct_value and sd_pct_mean (biodiversity)
#   f7_series_precision   one row per series: median sampling SD (and
#                         percentages), the tolerance, tolerance / median SD
#   02_Results/03_Diagnostics/f7_survey_precision/f7_summary.csv  medians over
#     the series by class, scale and sampling method (and q or taxon)
#   01_Analysis/Logs/f7_survey_precision/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   1. Bootstrap of the hauls or tows within each estuary and month: the
#      samples of the month are drawn with replacement (as many as were
#      taken), the monthly value is calculated again from the drawn samples,
#      and the standard deviation over the resamples is the sampling
#      uncertainty of that value.
#   2. Coastwide: the hauls or tows are drawn within each estuary (the
#      estuaries are the strata of the survey), so each estuary keeps its own
#      number of samples, and the drawn samples are then pooled as in b07 and
#      e02.
#   3. Biodiversity is calculated exactly as in b06 and b07: catches pooled by
#      taxon, Hill numbers q = 0, 1, 2 at the common coverage C_target with
#      estimate_q012_at_coverage() (lib/hill.R; iNEXT::estimateD, nboot = 0).
#      Abundance as in e02: mean catch per haul or tow including zeros
#      (menhaden: codes 604, 658 and 161 summed), then log(CPUE + 1).
#   4. The value from the original samples must reproduce b06, b07 and e02
#      (difference below CHECK_TOL); otherwise the script stops.
#   5. Only the target months that were evaluated (the month after the first
#      forecast date to December 2024) are resampled and summarized, so the
#      numbers describe the months on which the tolerance was used. Earlier
#      months had fewer bag seine hauls (Online Resource 1, Sect. 1.1).
#   6. Reported per series: the median sampling SD, for biodiversity also as a
#      percentage of the month's value and of the series mean (the mean of all
#      sampled months over the whole record, as in f01, from c02). The
#      tolerance is compared with the series-mean version because the
#      tolerance is a percentage of the series mean.
#   7. Resamples: N_RESAMPLE_DIV for biodiversity (iNEXT is slow) and
#      N_RESAMPLE_ABUND for abundance. Seed from settings$seed.
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
source(file.path(PATHS$analysis_dir, "lib", "hill.R"))

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
N_RESAMPLE_DIV     <- 100L       # bootstrap resamples per month, biodiversity
N_RESAMPLE_ABUND   <- 1000L      # bootstrap resamples per month, abundance
CHECK_TOL          <- 1e-8       # reproduction of b06, b07 and e02
N_WORKERS_OVERRIDE <- NA_integer_  # set an integer to force a count

TARGET_CODES <- c("602" = "602", "608" = "608", "633" = "633", "9605" = "9605", "9618" = "9618",
                  "601" = "601", "9600" = "9600",
                  "604" = "161", "658" = "161", "161" = "161")   # member code -> target (as e02)
TARGET_NAMES <- c("602" = "Atlantic croaker", "608" = "Spot", "633" = "Pinfish", "9605" = "Blue crab",
                  "9618" = "Brown shrimp", "601" = "Bay anchovy", "9600" = "White shrimp", "161" = "Menhaden")

script_stem <- "f7_survey_precision"          # a test copy uses its own stem, e.g. "f7_survey_precision_TEST"
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_station  <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_bs_ev    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bs_events.rds")
in_bt_ev    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bt_events.rds")
in_b06      <- file.path(PATHS$data_dir, "b6_compute_hill_alpha_q012", "alpha_q012.rds")
in_b07      <- file.path(PATHS$data_dir, "01_TimeSeries", "gamma_q012.rds")
in_c02_a    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_c02_g    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_e02_a    <- file.path(PATHS$data_dir, "e2_aggregate_abundance", "alpha_abundance_targets.rds")
in_e02_g    <- file.path(PATHS$data_dir, "e2_aggregate_abundance", "gamma_abundance_targets.rds")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "tidyr", "readr", "iNEXT", "future", "future.apply")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 2) Folders, log, inputs
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")
stop_if_missing(c(in_settings, in_station, in_bs_ev, in_bt_ev, in_b06, in_b07, in_c02_a, in_c02_g,
                  in_e02_a, in_e02_g, PATHS$clean_data_rdata),
                c("b09_forecast_settings.R", "b03_harmonize_station_metadata.R",
                  rep("b04_join_organisms_to_station.R", 2), "b06_compute_hill_alpha_q012.R",
                  "b07_compute_hill_gamma_q012.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 2),
                  rep("e02_aggregate_abundance.R", 2), "a00_getdata.R"))
settings <- load_settings(in_settings)
bays <- as.integer(settings$bays)

# Evaluated target months: from the month after the first forecast date to December 2024.
eval_start <- seq(as.Date(settings$first_origin_date), by = "month", length.out = 2)[2]
eval_end   <- as.Date("2024-12-01")
in_eval <- function(y, m) { d <- as.Date(sprintf("%04d-%02d-01", y, m)); d >= eval_start & d <= eval_end }
cat("Evaluated target months:", format(eval_start, "%Y-%m"), "to", format(eval_end, "%Y-%m"), "\n")
cat("Resamples: biodiversity", N_RESAMPLE_DIV, "| abundance", N_RESAMPLE_ABUND, "| seed", settings$seed, "\n")

station_core <- readRDS(in_station)
samples <- bind_rows(
  station_core$station_bs_core %>% mutate(GEAR = "BS"),
  station_core$station_bt_core %>% mutate(GEAR = "BT")) %>%
  transmute(SAMPLE_ID, GEAR, BAY_ID = as.integer(MAJOR_AREA), YEAR = as.integer(YEAR_i), MONTH = as.integer(MONTH_i)) %>%
  filter(BAY_ID %in% bays) %>%
  mutate(key = paste(GEAR, sprintf("%.0f", SAMPLE_ID)))   # text key, the same on every worker
if (anyDuplicated(samples[c("SAMPLE_ID", "GEAR")])) stop("Duplicated SAMPLE_ID in station_core", call. = FALSE)

# Catches by sample and taxon for biodiversity (as b06: placeholder rows left out).
events <- bind_rows(readRDS(in_bs_ev) %>% mutate(GEAR = "BS"), readRDS(in_bt_ev) %>% mutate(GEAR = "BT"))
expect_cols(events, c("SAMPLE_ID", "MAJOR_AREA", "YEAR_i", "MONTH_i", "SPECIES_CODE", "catch"), "events")
div_catch <- events %>%
  filter(!is.na(SPECIES_CODE), as.integer(MAJOR_AREA) %in% bays) %>%
  transmute(key = paste(GEAR, sprintf("%.0f", SAMPLE_ID)), SPECIES_CODE = as.character(SPECIES_CODE), catch = as.numeric(catch))

b06 <- readRDS(in_b06); b07 <- readRDS(in_b07)
C_target <- unique(c(b06$C_target, b07$C_target))
if (length(C_target) != 1L) stop("b06 and b07 do not share one C_target", call. = FALSE)
cat("C_target:", C_target, "\n")

# Catches by sample and target for abundance (as e02).
load(PATHS$clean_data_rdata)
ab_catch <- bind_rows(bagseine_organisms %>% transmute(SAMPLE_ID, GEAR = "BS", code = as.character(SPECIES_CODE), catch = as.numeric(catch)),
                      baytrawl_organisms %>% transmute(SAMPLE_ID, GEAR = "BT", code = as.character(SPECIES_CODE), catch = as.numeric(catch))) %>%
  filter(code %in% names(TARGET_CODES)) %>%
  mutate(target = unname(TARGET_CODES[code])) %>%
  group_by(SAMPLE_ID, GEAR, target) %>% summarise(catch = sum(catch), .groups = "drop") %>%
  mutate(key = paste(GEAR, sprintf("%.0f", SAMPLE_ID)))
if (anyNA(ab_catch$catch)) stop("Missing catch in a target record", call. = FALSE)

# Series means of biodiversity (as f01: all sampled months of the whole record).
lev <- bind_rows(
  readRDS(in_c02_a) %>% filter(!filled) %>% group_by(BAY_ID, GEAR, q) %>% summarise(level = mean(value), .groups = "drop") %>%
    transmute(scale = "estuary", BAY_ID = as.integer(BAY_ID), GEAR, q = as.integer(q), level),
  readRDS(in_c02_g) %>% filter(!filled) %>% group_by(GEAR, q) %>% summarise(level = mean(value), .groups = "drop") %>%
    transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, q = as.integer(q), level))

# -----------------------------------------------------------------------------
# 3) Cells to resample
# -----------------------------------------------------------------------------
# A cell is one estuary (or the coast), sampling method and month; strata are
# the estuaries (one stratum for an estuary cell, every estuary for a coastwide cell).
cells_bay <- samples %>% filter(in_eval(YEAR, MONTH)) %>% distinct(GEAR, BAY_ID, YEAR, MONTH) %>%
  mutate(scale = "estuary")
cells_coast <- samples %>% filter(in_eval(YEAR, MONTH)) %>% group_by(GEAR, YEAR, MONTH) %>%
  summarise(n_bays = n_distinct(BAY_ID), .groups = "drop")
if (any(cells_coast$n_bays != length(bays))) stop("A coastwide month in the evaluated period lacks an estuary", call. = FALSE)
cells_coast <- cells_coast %>% transmute(GEAR, BAY_ID = NA_integer_, YEAR, MONTH, scale = "coastwide")
cells <- bind_rows(cells_bay, cells_coast) %>% arrange(scale, GEAR, BAY_ID, YEAR, MONTH)
cells$cell <- seq_len(nrow(cells))
cat("Cells: estuary", nrow(cells_bay), "| coastwide", nrow(cells_coast), "\n")

# Sample lists per cell, by stratum.
samp_split <- split(samples, list(samples$GEAR, samples$YEAR, samples$MONTH), drop = TRUE)
cell_strata <- lapply(seq_len(nrow(cells)), function(i) {
  s <- samp_split[[paste(cells$GEAR[i], cells$YEAR[i], cells$MONTH[i], sep = ".")]]
  if (!is.na(cells$BAY_ID[i])) s <- s[s$BAY_ID == cells$BAY_ID[i], ]
  split(s$key, s$BAY_ID)
})
cells$n_samples <- vapply(cell_strata, function(z) sum(lengths(z)), integer(1))

# One draw: the samples of every stratum drawn with replacement (draw = FALSE: the samples as taken).
draw_ids <- function(strata, draw) {
  unlist(lapply(strata, function(ids) if (draw) ids[sample.int(length(ids), length(ids), replace = TRUE)] else ids),
         use.names = FALSE)
}

# -----------------------------------------------------------------------------
# 4) Biodiversity
# -----------------------------------------------------------------------------
div_by_sample <- split(div_catch[c("SPECIES_CODE", "catch")], div_catch$key)
hill_of <- function(keys) {
  z <- div_by_sample[keys]
  z <- z[!vapply(z, is.null, logical(1))]
  if (!length(z)) return(c(0, 0, 0))
  z <- do.call(rbind, z)
  r <- estimate_q012_at_coverage(tapply(z$catch, z$SPECIES_CODE, sum), C_target = C_target, nboot = 0, conf = 0.95)
  c(r$q0, r$q1, r$q2)
}
div_cell <- function(i) {
  obs <- hill_of(draw_ids(cell_strata[[i]], FALSE))
  bs <- vapply(seq_len(N_RESAMPLE_DIV), function(b) hill_of(draw_ids(cell_strata[[i]], TRUE)), numeric(3))
  data.frame(cell = i, q = 0:2, value = obs, sd_sampling = apply(bs, 1, stats::sd), n_na = rowSums(is.na(bs)))
}

wk <- setup_workers(N_WORKERS_OVERRIDE, globals_max_gb = 4)
cat("\nBiodiversity: resampling", nrow(cells), "cells ...\n")
t0 <- Sys.time()
div_res <- future.apply::future_lapply(seq_len(nrow(cells)), div_cell, future.seed = as.integer(settings$seed),
                                       future.scheduling = 4)
div_res <- do.call(rbind, div_res)
cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
future::plan(future::sequential)
if (any(div_res$n_na > 0)) stop("estimateD() failed on ", sum(div_res$n_na), " resamples", call. = FALSE)

div <- cells %>% inner_join(div_res, by = "cell")
# Check against b06 and b07.
ref_div <- bind_rows(
  b06 %>% transmute(scale = "estuary", GEAR, BAY_ID = as.integer(BAY_ID), YEAR = as.integer(YEAR), MONTH = as.integer(MONTH), q0, q1, q2),
  b07 %>% transmute(scale = "coastwide", GEAR, BAY_ID = NA_integer_, YEAR = as.integer(YEAR), MONTH = as.integer(MONTH), q0, q1, q2)) %>%
  tidyr::pivot_longer(c(q0, q1, q2), names_to = "q", values_to = "value_ref") %>% mutate(q = as.integer(sub("q", "", q)))
div <- div %>% left_join(ref_div, by = c("scale", "GEAR", "BAY_ID", "YEAR", "MONTH", "q"))
d_div <- max(abs(div$value - div$value_ref))
cat("Biodiversity: largest difference from b06/b07:", signif(d_div, 3), "\n")
if (!is.finite(d_div) || d_div > CHECK_TOL) stop("Biodiversity values do not reproduce b06/b07", call. = FALSE)
div <- div %>% left_join(lev, by = c("scale", "BAY_ID", "GEAR", "q"))
if (anyNA(div$level)) stop("Missing series mean for a biodiversity series", call. = FALSE)
div <- div %>% mutate(class = "biodiversity", target = paste0("q", q),
                      sd_pct_value = ifelse(value > 0, 100 * sd_sampling / value, NA_real_),
                      sd_pct_mean = 100 * sd_sampling / level)

# -----------------------------------------------------------------------------
# 5) Abundance
# -----------------------------------------------------------------------------
set.seed(as.integer(settings$seed))
tg <- names(TARGET_NAMES)
ab_by_sample <- split(ab_catch[c("target", "catch")], ab_catch$key)
ab_mat <- function(keys) {                  # samples x targets, zeros included
  m <- matrix(0, length(keys), length(tg), dimnames = list(NULL, tg))
  z <- ab_by_sample[keys]; n <- vapply(z, function(x) if (is.null(x)) 0L else nrow(x), integer(1))
  if (sum(n)) { k <- do.call(rbind, z[n > 0]); m[cbind(rep(seq_along(keys), n), match(k$target, tg))] <- k$catch }
  m
}
cat("\nAbundance: resampling", nrow(cells), "cells ...\n")
ab_rows <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  st <- cell_strata[[i]]; ids <- unlist(st, use.names = FALSE)
  M <- ab_mat(ids)
  pos <- split(seq_along(ids), rep(seq_along(st), lengths(st)))
  # Draw within each stratum: B x n matrix of row indices.
  idx <- do.call(cbind, lapply(pos, function(p) matrix(p[sample.int(length(p), length(p) * N_RESAMPLE_ABUND, replace = TRUE)],
                                                       nrow = N_RESAMPLE_ABUND)))
  bs <- vapply(tg, function(t) { v <- M[, t]; log1p(rowMeans(matrix(v[idx], nrow = N_RESAMPLE_ABUND))) }, numeric(N_RESAMPLE_ABUND))
  ab_rows[[i]] <- data.frame(cell = i, target = tg, value = log1p(colMeans(M)), sd_sampling = apply(bs, 2, stats::sd))
}
ab <- cells %>% inner_join(do.call(rbind, ab_rows), by = "cell")
ref_ab <- bind_rows(
  readRDS(in_e02_a) %>% filter(!filled) %>% transmute(scale = "estuary", GEAR, BAY_ID = as.integer(BAY_ID), YEAR = as.integer(YEAR),
                                                     MONTH = as.integer(MONTH), target = as.character(SPECIES_CODE), value_ref = log_cpue),
  readRDS(in_e02_g) %>% filter(!filled) %>% transmute(scale = "coastwide", GEAR, BAY_ID = NA_integer_, YEAR = as.integer(YEAR),
                                                     MONTH = as.integer(MONTH), target = as.character(SPECIES_CODE), value_ref = log_cpue))
ab <- ab %>% left_join(ref_ab, by = c("scale", "GEAR", "BAY_ID", "YEAR", "MONTH", "target"))
d_ab <- max(abs(ab$value - ab$value_ref))
cat("Abundance: largest difference from e02:", signif(d_ab, 3), "\n")
if (!is.finite(d_ab) || d_ab > CHECK_TOL) stop("Abundance values do not reproduce e02", call. = FALSE)
ab <- ab %>% mutate(class = "abundance", q = NA_integer_, level = NA_real_, sd_pct_value = NA_real_, sd_pct_mean = NA_real_,
                    target = unname(TARGET_NAMES[target]))

# -----------------------------------------------------------------------------
# 6) Tables
# -----------------------------------------------------------------------------
cell_tab <- bind_rows(div, ab) %>%
  transmute(class, scale, GEAR, BAY_ID, target, YEAR, MONTH, n_samples, value, sd_sampling, sd_pct_value, sd_pct_mean) %>%
  arrange(class, scale, GEAR, BAY_ID, target, YEAR, MONTH)

series_tab <- cell_tab %>% group_by(class, scale, GEAR, BAY_ID, target) %>%
  summarise(months = n(), n_samples_median = stats::median(n_samples),
            sd_median = stats::median(sd_sampling), sd_pct_value_median = stats::median(sd_pct_value, na.rm = TRUE),
            sd_pct_mean_median = stats::median(sd_pct_mean), .groups = "drop") %>%
  mutate(tolerance = ifelse(class == "biodiversity", 100 * settings$tolerance$diversity_rel, settings$tolerance$abundance_log),
         sd_ref = ifelse(class == "biodiversity", sd_pct_mean_median, sd_median),
         tol_over_sd = ifelse(sd_ref > 0, tolerance / sd_ref, NA_real_)) %>%   # NA: no sampling variation (almost no catch)
  select(-sd_ref)

med <- function(x) round(stats::median(x, na.rm = TRUE), 3)
summ_one <- function(df) df %>% summarise(series = n(), n_samples = med(n_samples_median), sd = med(sd_median),
                                          sd_pct_value = med(sd_pct_value_median), sd_pct_mean = med(sd_pct_mean_median),
                                          tol_over_sd = med(tol_over_sd), .groups = "drop")
summ <- bind_rows(
  summ_one(series_tab %>% group_by(class, scale, GEAR)) %>% mutate(target = "all"),
  summ_one(series_tab %>% group_by(class, scale, GEAR, target))) %>%
  arrange(class, scale, GEAR, target != "all", target)
cat("\nMedians over the series (sd: log units for abundance; sd_pct_*: % for biodiversity):\n")
print(as.data.frame(summ), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7) Save
# -----------------------------------------------------------------------------
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE); if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv)$", full.names = TRUE); if (length(old_d)) archive_files(old_d)
saveRDS(cell_tab, file.path(dirs$data, "f7_cell_precision.rds"))
readr::write_csv(cell_tab, file.path(dirs$data, "f7_cell_precision.csv"), na = "")
saveRDS(series_tab, file.path(dirs$data, "f7_series_precision.rds"))
readr::write_csv(series_tab, file.path(dirs$data, "f7_series_precision.csv"), na = "")
readr::write_csv(summ, file.path(dirs$diag, "f7_summary.csv"), na = "")
cat("\nSaved to", dirs$data, "and", dirs$diag, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
