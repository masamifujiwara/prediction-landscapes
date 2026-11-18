# =============================================================================
# Script:  f01_forecast_scoring.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Scores every out-of-sample forecast of every model and reference for the
#   48 diversity and 128 abundance series at leads 1 to 24 (lib/scoring.R):
#   absolute accuracy, skill, uncertainty ranges, horizons, quadrants, the
#   variance partition, the three periods, the trend and covariate tests, the
#   check against the published accuracy measure, the distribution of the
#   monthly misses, the largest misses, summaries by group, and rank
#   correlations with the slow-pattern diagnostics of f02.
#
# Run after: b09, c02, d01, e09, f02 (runs before f01, see lib/run_order.R),
#   lib/helpers.R, lib/forecast_rules.R, lib/join_forecasts.R, lib/scoring.R
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   - 01_Analysis/Data/d1_compare_model_metrics_dictionary_v3_scaled/forecasts_joined_diversity.rds (d01)
#   - 01_Analysis/Data/e9_compare_predictability_abundance/forecasts_joined_abundance.rds (e09)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds and
#     gamma_q012_deseas_scaled.rds (c02; mean level of each diversity series)
#   - 01_Analysis/Data/f2_slow_pattern_diagnostics/f2_series_diagnostics.rds (f02)
#
# Outputs:
#   01_Analysis/Data/f1_forecast_scoring/, each table as .rds and .csv:
#   f1_by_lead               series x model x lead: n, scaled RMSE with range,
#                            bias, miss quantiles, hit rate, interval coverage
#                            and width, skill against the cumulative average,
#                            the 12-month average, last month's value and the
#                            partner model (with ranges), full-scale variance
#                            explained, skill retained
#   f1_horizons              series x model: absolute and relative horizons
#                            (range and value rules), return of skill
#   f1_horizon_by_tolerance  series x model x tolerance
#   f1_quadrants             series x model x lead (1, 12, 24)
#   f1_variance_budget       series x model
#   f1_periods               series x model x period x lead
#   f1_period_horizons       series x model x period
#   f1_published_check       series x model x lead (1, 12), all forecast dates
#   f1_failed_cases          series x model
#   f1_miss_distribution     class x scale x sampling method x model x lead
#                            (pooled)
#   f1_miss_distribution_by_target  the same by q or taxon
#   f1_largest_misses        the LARGEST_N largest standardized misses per class
#   f1_coherent_miss_months  target months ranked by the share of series with
#                            a large miss, with counts and shares of large
#                            misses where the observed value was higher or
#                            lower than forecast
#   f1_group_summary         medians by class, scale, sampling method, target,
#                            model and lead
#   f1_bay_vs_coast          coastwide against the estuaries, by cell, model
#                            and lead
#   f1_diag_correlations     Spearman correlations of f02 diagnostics with
#                            SARIMA results
#   01_Analysis/Data/f1_forecast_scoring/partial/f1_partial_<class>.rds  saved
#     result per class (reused on a rerun only if code, settings and input match)
#   02_Results/03_Diagnostics/f1_forecast_scoring/f1_headline.csv  short
#     headline table
#   01_Analysis/Logs/f1_forecast_scoring/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   1. The main analysis uses the common forecast dates (in_common); the
#      periods use all forecast dates. Targets in April and May 2020 are never
#      scored.
#   2. Tolerance from settings.rds: a share of the series' mean level for
#      diversity (mean over all sampled months) and a fixed number of log
#      units for abundance, the same at every lead. Horizon-against-tolerance
#      curves run over settings$tolerance_range in steps of TOL_STEP.
#   3. Failed fits: the case is left out of every comparison involving the
#      model, for both members of the pair; counts are written to
#      f1_failed_cases.
#   4. Horizons by the uncertainty-range rule (primary) and by value; return
#      of skill after the relative horizon; quadrants at leads 1, 12 and 24.
#   5. Moving-block bootstrap: blocks of 12 forecast dates, 2,000 resamples,
#      seed settings$seed plus the series number, 95% ranges.
#   6. Variance partition for every model: seasonal, slow (skill at leads 12
#      to the longest lead), short-lived (skill at lead 1 minus slow),
#      unpredictable.
#   7. Monthly-miss distribution pooled by class, scale, sampling method,
#      model and lead (quantiles and share within the tolerance), also by
#      target (q or taxon).
#   8. Largest misses: SARIMA (automatic) at lead 1, all forecast dates. A
#      "large" miss is more than LARGE_SD times the series' own RMSE; months
#      in which many series miss at once are listed (coherent months).
#   9. Group summaries: by Hill order and taxon, estuary against coastwide,
#      skill retained over leads.
#  10. Series are scored in parallel (setup_workers(); N_WORKERS_OVERRIDE).
#      The result of each class is saved and reused only if the code,
#      settings and input are unchanged, so a rerun skips a finished class.
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
source(file.path(PATHS$analysis_dir, "lib", "join_forecasts.R"))   # model_info()
source(file.path(PATHS$analysis_dir, "lib", "scoring.R"))

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
N_WORKERS_OVERRIDE <- NA_integer_   # set an integer to force a count
B_BOOT      <- 2000L
BLOCK_LEN   <- 12L
LARGE_SD    <- 2      # a miss larger than this many times the series' RMSE is "large"
LARGEST_N   <- 200L
TOL_STEP    <- 0.1
RUN_CLASSES <- c("diversity", "abundance")

script_stem <- "f1_forecast_scoring"
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_div      <- file.path(PATHS$data_dir, "d1_compare_model_metrics_dictionary_v3_scaled", "forecasts_joined_diversity.rds")
in_abd      <- file.path(PATHS$data_dir, "e9_compare_predictability_abundance", "forecasts_joined_abundance.rds")
in_div_a    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_div_g    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_f2       <- file.path(PATHS$data_dir, "f2_slow_pattern_diagnostics", "f2_series_diagnostics.rds")

options(scipen = 999, warn = 1, width = 220)
required_pkgs <- c("dplyr", "readr", "tidyr", "future", "future.apply")
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
stop_if_missing(c(in_settings, in_div, in_abd, in_div_a, in_div_g, in_f2),
                c("b09_forecast_settings.R", "d01_compare_model_metrics_dictionary_v3_scaled.R", "e09_compare_predictability_abundance.R",
                  rep("c02_deseasonalize_alpha_gamma_scaled.R", 2), "f02_slow_pattern_diagnostics.R"))
settings <- load_settings(in_settings)
mi <- model_info(); partner_of <- setNames(mi$partner, mi$model)
f2 <- readRDS(in_f2)

# Mean level of each diversity series over all sampled months (method note 2).
lev_div <- bind_rows(
  readRDS(in_div_a) %>% filter(!filled) %>% group_by(BAY_ID, GEAR, q) %>% summarise(level = mean(value), .groups = "drop") %>%
    transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, q = as.integer(q), level),
  readRDS(in_div_g) %>% filter(!filled) %>% group_by(GEAR, q) %>% summarise(level = mean(value), .groups = "drop") %>%
    transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, q = as.integer(q), level))

cls <- list(
  diversity = list(file = in_div, keys = c("scale", "BAY_ID", "GEAR", "q"), tol = settings$tolerance$diversity_rel,
                   grid = seq(settings$tolerance_range$diversity_rel[1], settings$tolerance_range$diversity_rel[2], TOL_STEP), n = 48L),
  abundance = list(file = in_abd, keys = c("scale", "BAY_ID", "GEAR", "SPECIES_CODE"), tol = settings$tolerance$abundance_log,
                   grid = seq(settings$tolerance_range$abundance_log[1], settings$tolerance_range$abundance_log[2], TOL_STEP), n = 128L))
cat("Tolerances: diversity", cls$diversity$tol, "of the mean level; abundance", cls$abundance$tol, "log units\n")

w <- setup_workers(override = N_WORKERS_OVERRIDE, globals_max_gb = 8)
on.exit(future::plan(future::sequential), add = TRUE)

# -----------------------------------------------------------------------------
# 3) Score every series, one class at a time
# -----------------------------------------------------------------------------
res <- list()
this_script <- file.path(PATHS$analysis_dir, "f01_forecast_scoring.R")
for (cl in RUN_CLASSES) {
  cc <- cls[[cl]]; t1 <- Sys.time()
  # Saved result per class (method note 10): reused only if code, settings and
  # input are unchanged, so a rerun after a reboot skips a finished class.
  sig <- paste(c(unname(tools::md5sum(c(cc$file, in_settings, in_div_a, in_div_g, file.path(PATHS$analysis_dir, "lib", c("scoring.R", "join_forecasts.R")),
                                      if (file.exists(this_script)) this_script))), B_BOOT, BLOCK_LEN, LARGE_SD), collapse = "|")
  cache <- file.path(dirs$data, "partial", paste0("f1_partial_", cl, ".rds"))
  if (file.exists(cache)) {
    pc <- readRDS(cache)
    if (identical(pc$sig, sig)) { res[[cl]] <- pc$res; cat("\n", cl, ": saved result reused (", cache, ")\n", sep = ""); next }
    cat("\n", cl, ": saved result out of date, scoring again\n", sep = "")
  }
  fc <- readRDS(cc$file)
  keep <- c(cc$keys, intersect("Common_Name", names(fc)), "origin_i", "in_common", "lead", "target_date", "scored", "model", "ok",
            "mean", "lo80", "hi80", "lo95", "hi95", "y_obs", "value_obs")
  fc <- fc[, keep]
  sid <- do.call(paste, c(as.list(fc[, cc$keys]), sep = "|"))
  r <- rle(sid)
  if (length(r$values) != cc$n || anyDuplicated(r$values)) stop(cl, ": series are not stored in contiguous blocks of ", cc$n, " series", call. = FALSE)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  pieces <- lapply(seq_along(starts), function(i) fc[starts[i]:ends[i], ])
  rm(fc, sid); invisible(gc())
  keys_tab <- do.call(rbind, lapply(pieces, function(p) p[1, c(cc$keys, intersect("Common_Name", names(p)))])); rownames(keys_tab) <- NULL
  levels_v <- if (cl == "diversity") {
    k <- keys_tab %>% left_join(lev_div, by = cc$keys); if (anyNA(k$level)) stop("Missing mean level for a diversity series", call. = FALSE); k$level
  } else rep(1, nrow(keys_tab))
  cat("\n", cl, ":", length(pieces), "series read and split in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min\n")
  t1 <- Sys.time()
  # Each worker receives only its own series (not the whole table).
  out <- future.apply::future_mapply(function(p, lv, i) {
    score_series(p, lv, cc$tol, cc$grid, B = B_BOOT, L = BLOCK_LEN, seed = settings$seed + i, partner_of = partner_of)
  }, pieces, levels_v, seq_along(pieces), SIMPLIFY = FALSE, future.seed = TRUE, future.scheduling = Inf,
     future.globals = list(cc = cc, B_BOOT = B_BOOT, BLOCK_LEN = BLOCK_LEN, settings = settings, partner_of = partner_of,
                           score_series = score_series, block_counts = block_counts, q_lo = q_lo, q_hi = q_hi, count_horizon = count_horizon))
  cat(cl, ": scored in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min on", w$n_workers, "workers\n")
  rm(pieces); invisible(gc())
  add_keys <- function(name) do.call(rbind, lapply(seq_along(out), function(i) {
    x <- out[[i]][[name]]; if (is.null(x) || !nrow(x)) return(NULL); cbind(class = cl, keys_tab[rep(i, nrow(x)), , drop = FALSE], level = levels_v[i], x, row.names = NULL) }))
  res[[cl]] <- lapply(setNames(nm = c("by_lead", "horizons", "horizon_by_tolerance", "quadrants", "budget", "periods", "period_horizons",
                                     "published", "failed", "big")), add_keys)
  # Pooled monthly-miss distribution (method note 7), then drop the raw misses.
  ms <- add_keys("miss")
  tgt <- if (cl == "diversity") paste0("q", ms$q) else ms$Common_Name
  qs <- function(g) ms %>% mutate(target = tgt) %>% group_by(across(all_of(g))) %>%
    summarise(n = n(), q10 = stats::quantile(miss_scaled, 0.10), q25 = stats::quantile(miss_scaled, 0.25), median = stats::median(miss_scaled),
              q75 = stats::quantile(miss_scaled, 0.75), q90 = stats::quantile(miss_scaled, 0.90), hit_rate = mean(miss_scaled <= cc$tol), .groups = "drop")
  res[[cl]]$miss_distribution <- qs(c("class", "scale", "GEAR", "model", "lead"))
  res[[cl]]$miss_distribution_by_target <- qs(c("class", "scale", "GEAR", "target", "model", "lead"))
  rm(ms, out); invisible(gc())
  dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(sig = sig, res = res[[cl]]), cache)
}
future::plan(future::sequential)

bind <- function(name) {
  x <- lapply(res, `[[`, name)
  x <- lapply(x, function(d) { if (is.null(d)) return(NULL); if (!"q" %in% names(d)) d$q <- NA_integer_
    if (!"SPECIES_CODE" %in% names(d)) d$SPECIES_CODE <- NA_character_; if (!"Common_Name" %in% names(d)) d$Common_Name <- NA_character_; d })
  d <- bind_rows(x)
  d$target <- ifelse(d$class == "diversity", paste0("q", d$q), d$Common_Name)
  d %>% relocate(class, scale, BAY_ID, GEAR, target, q, SPECIES_CODE, Common_Name)
}
by_lead <- bind("by_lead"); horizons <- bind("horizons"); hz_tol <- bind("horizon_by_tolerance"); quadrants <- bind("quadrants")
budget <- bind("budget"); periods <- bind("periods"); period_hz <- bind("period_horizons"); published <- bind("published")
failed <- bind("failed"); big <- bind("big")
miss_dist <- bind_rows(lapply(res, `[[`, "miss_distribution")); miss_dist_t <- bind_rows(lapply(res, `[[`, "miss_distribution_by_target"))

# Checks.
n_expected <- sum(vapply(cls[RUN_CLASSES], `[[`, integer(1), "n"))
n_series <- nrow(distinct(by_lead, class, scale, BAY_ID, GEAR, target))
if (n_series != n_expected) stop("Scored ", n_series, " series; expected ", n_expected, call. = FALSE)
if (anyDuplicated(by_lead[, c("class", "scale", "BAY_ID", "GEAR", "target", "model", "lead")])) stop("Duplicate rows in f1_by_lead", call. = FALSE)
cat("\nScored", n_series, "series;", nrow(by_lead), "series x model x lead rows\n")

# -----------------------------------------------------------------------------
# 4) Largest misses and coherent months (method note 8)
# -----------------------------------------------------------------------------
big$large <- abs(big$miss_sd) > LARGE_SD
largest <- big %>% group_by(class) %>% arrange(desc(abs(miss_sd)), .by_group = TRUE) %>% slice_head(n = LARGEST_N) %>% ungroup()
# Errors are forecast minus observed, so miss_sd < 0 means the observed value
# was higher than forecast. Shares are of all series scored in that month;
# n_large_* count the large misses by direction. The flag `large` must stay
# logical until the shares are computed.
coherent <- big %>% group_by(class, target_date) %>% summarise(series = n(),
  n_large = sum(large), n_large_obs_higher = sum(large & miss_sd < 0), n_large_obs_lower = sum(large & miss_sd > 0),
  .groups = "drop") %>%
  mutate(share_large = n_large / series, share_large_obs_higher = n_large_obs_higher / series,
         share_large_obs_lower = n_large_obs_lower / series) %>%
  arrange(class, desc(share_large))
cat("\nShare of lead-1 SARIMA misses larger than", LARGE_SD, "times the series' RMSE:", round(mean(big$large), 3), "\n")
cat("Target months with the largest share of series missing at once (top 5 per class):\n")
print(as.data.frame(coherent %>% group_by(class) %>% slice_head(n = 5) %>% ungroup() %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5) Group summaries (method note 9) and correlations with f02
# -----------------------------------------------------------------------------
group_summary <- by_lead %>% group_by(class, scale, GEAR, target, model, lead) %>%
  summarise(series = n(), rmse_scaled = stats::median(rmse_scaled), hit_rate = stats::median(hit_rate), skill_past = stats::median(skill_past),
            skill_m12 = stats::median(skill_m12), skill_naive = stats::median(skill_naive), skill_partner = stats::median(skill_partner),
            share_skill_past_above0 = mean(skill_past_lo > 0), skill_retained = stats::median(skill_retained), .groups = "drop")
coast <- by_lead %>% filter(scale == "coastwide") %>% select(class, GEAR, target, model, lead, coast_rmse = rmse_scaled, coast_skill = skill_past)
bay_vs_coast <- by_lead %>% filter(scale == "bay") %>% group_by(class, GEAR, target, model, lead) %>%
  summarise(bays = n(), bay_median_rmse = stats::median(rmse_scaled), bay_median_skill = stats::median(skill_past), .groups = "keep") %>%
  left_join(coast, by = c("class", "GEAR", "target", "model", "lead")) %>%
  left_join(by_lead %>% filter(scale == "bay") %>% select(class, GEAR, target, model, lead, rmse_scaled) %>% group_by(class, GEAR, target, model, lead) %>%
              summarise(rm = list(rmse_scaled), .groups = "drop"), by = c("class", "GEAR", "target", "model", "lead")) %>%
  mutate(bays_more_accurate_than_coast = vapply(seq_along(rm), function(i) sum(rm[[i]] < coast_rmse[i]), integer(1))) %>% select(-rm) %>% ungroup()

sar <- by_lead %>% filter(model == "sarima_auto") %>% select(class, scale, BAY_ID, GEAR, target, lead, rmse_scaled, skill_past, skill_m12)
res_series <- sar %>% filter(lead %in% c(1, 12, 24)) %>% mutate(gap = skill_past - skill_m12) %>%
  tidyr::pivot_wider(id_cols = c(class, scale, BAY_ID, GEAR, target), names_from = lead, values_from = c(rmse_scaled, skill_past, gap)) %>%
  left_join(horizons %>% filter(model == "sarima_auto") %>% select(class, scale, BAY_ID, GEAR, target, rel_horizon_past, abs_horizon),
            by = c("class", "scale", "BAY_ID", "GEAR", "target")) %>%
  left_join(periods %>% filter(model == "sarima_auto", lead %in% c(1, 12)) %>% select(class, scale, BAY_ID, GEAR, target, period, lead, skill_past) %>%
              tidyr::pivot_wider(names_from = c(period, lead), values_from = skill_past, names_prefix = "p_") %>%
              transmute(class, scale, BAY_ID, GEAR, target, period_diff_1 = `p_2020-2024_1` - `p_1992-2007_1`,
                        period_diff_12 = `p_2020-2024_12` - `p_1992-2007_12`), by = c("class", "scale", "BAY_ID", "GEAR", "target"))
dg <- f2 %>% mutate(abs_sen_sd = abs(sen_per_decade_sd)) %>%
  select(class, scale, BAY_ID, GEAR, target, seasonal_share, slow_share, slow_ratio, abs_sen_sd, acf1, acf12, step_2020_sd, swing_ratio)
both <- inner_join(res_series, dg, by = c("class", "scale", "BAY_ID", "GEAR", "target"))
if (nrow(both) != n_expected) stop("f02 diagnostics do not match the ", n_expected, " scored series", call. = FALSE)
dvars <- c("seasonal_share", "slow_share", "slow_ratio", "abs_sen_sd", "acf1", "acf12", "step_2020_sd", "swing_ratio")
rvars <- setdiff(names(res_series), c("class", "scale", "BAY_ID", "GEAR", "target"))
diag_cor <- bind_rows(lapply(split(both, list(both$class, both$scale), drop = TRUE), function(x) {
  expand.grid(diagnostic = dvars, result = rvars, stringsAsFactors = FALSE) %>% rowwise() %>%
    mutate(class = x$class[1], scale = x$scale[1], n = sum(is.finite(x[[diagnostic]]) & is.finite(x[[result]])),
           spearman = suppressWarnings(stats::cor(x[[diagnostic]], x[[result]], method = "spearman", use = "complete.obs"))) %>% ungroup()
}))

headline <- by_lead %>% filter(model %in% c("sarima_auto", "ssm_level", "arimax", "ssm_xreg", "damped_trend"), lead %in% c(1, 6, 12, 24)) %>%
  group_by(class, scale, model, lead) %>% summarise(series = n(), rmse_scaled = stats::median(rmse_scaled), hit_rate = stats::median(hit_rate),
  skill_past = stats::median(skill_past), skill_m12 = stats::median(skill_m12), share_skilful = mean(skill_past_lo > 0), .groups = "drop")

# -----------------------------------------------------------------------------
# 6) Write
# -----------------------------------------------------------------------------
tabs <- list(f1_by_lead = by_lead, f1_horizons = horizons, f1_horizon_by_tolerance = hz_tol, f1_quadrants = quadrants,
             f1_variance_budget = budget, f1_periods = periods, f1_period_horizons = period_hz, f1_published_check = published,
             f1_failed_cases = failed, f1_miss_distribution = miss_dist, f1_miss_distribution_by_target = miss_dist_t,
             f1_largest_misses = largest, f1_coherent_miss_months = coherent, f1_group_summary = group_summary,
             f1_bay_vs_coast = bay_vs_coast, f1_diag_correlations = diag_cor)
old <- list.files(dirs$data, pattern = "^f1_.*[.](rds|csv)$", full.names = TRUE); if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|png|pdf)$", full.names = TRUE); if (length(old_d)) archive_files(old_d)
for (nm in names(tabs)) { saveRDS(tabs[[nm]], file.path(dirs$data, paste0(nm, ".rds"))); readr::write_csv(tabs[[nm]], file.path(dirs$data, paste0(nm, ".csv")), na = "") }
readr::write_csv(headline, file.path(dirs$diag, "f1_headline.csv"), na = "")

cat("\nHeadline (medians over series; skill against the past average; share skilful = lower end of the range above zero):\n")
print(as.data.frame(headline %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)
cat("\nCases removed by failed fits (common forecast dates):\n")
print(as.data.frame(failed %>% group_by(class, model) %>% summarise(cases = sum(cases_failed), .groups = "drop") %>% filter(cases > 0)), row.names = FALSE)
cat("\nSaved", length(tabs), "tables (rds and csv) in", dirs$data, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("f01_forecast_scoring.R completed\n")
