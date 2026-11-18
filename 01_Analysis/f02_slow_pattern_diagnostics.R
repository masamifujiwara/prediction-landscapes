# =============================================================================
# Script:  f02_slow_pattern_diagnostics.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Describes the slow, non-seasonal patterns of every target series from the
#   data alone (no forecasts): 48 diversity series (q0, q1, q2; seven
#   estuaries and coastwide; bag seine and bay trawl) and 128 abundance series
#   (8 taxa). Gives the slow-pattern diagnostics and the correlation-by-lag
#   curves for the Online Resources. f01 relates them to the forecast results.
#
# Run after: b09, c02, e03, lib/helpers.R, lib/forecast_rules.R,
#   lib/seasonal.R (runs before f01, see lib/run_order.R)
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds (e03)
#   - 01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds (e03)
#
# Outputs:
#   01_Analysis/Data/f2_slow_pattern_diagnostics/:
#   - f2_series_diagnostics.rds / .csv  one row per series (keys: class,
#                                       scale, BAY_ID, GEAR, target)
#   - f2_acf_by_lag.rds / .csv          one row per series and lag
#   - f2_yearly_means.rds               yearly mean departures, for figures
#   02_Results/03_Diagnostics/f2_slow_pattern_diagnostics/f2_summary_by_group.csv
#     medians and quantiles by class, scale and sampling method
#   01_Analysis/Logs/f2_slow_pattern_diagnostics/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   1. Computed on the departures from the seasonal curve fitted to the whole
#      record (value_deseas of c02 and e03). This is for description only;
#      the forecasts use curves fitted at each forecast date. April and May
#      2020 (filled, never sampled) are set to missing everywhere.
#   2. Measures per series:
#      - seasonal share: share of the total variance explained by the seasonal
#        curve (1 - variance of the departures / variance of the values);
#      - trend: Mann-Kendall test and Sen slope on the yearly means of the
#        departures (years with at least MIN_MONTHS_YEAR sampled months), the
#        slope also in SD units per decade; straight-line trend share of the
#        monthly departures and its slope; KPSS differencing order
#        (forecast::ndiffs, for reference only);
#      - slow share: share of the variance of the departures in a 25-month
#        running mean, against N_SIM simulated AR(1) series with the same
#        lag-1 correlation (95th percentile, ratio to the median, and the share
#        of simulations at or above the observed value);
#      - breaks: shifts in the mean of the yearly means, at most two, chosen by
#        BIC with segments of at least MIN_SEGMENT_YEARS years (exhaustive
#        search, same criterion as strucchange::breakpoints);
#      - the 2020 to 2024 step: mean departure 2020 to 2024 minus 2010 to 2019,
#        in SD units of the departures;
#      - seasonal swing and peak month in the first and second half of the
#        record (three-harmonic curve fitted to each half);
#      - correlation of the departures at lags 1 to MAX_LAG months.
#   3. In the summary by group, slow_above_ar1_q95 is counted before
#      slow_share is replaced by its text summary.
#   4. Rank correlations with the forecast results (the gap between skill
#      against the cumulative average and against the 12-month average,
#      period differences, horizons) are computed in f01, which reads this
#      output.
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

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
MIN_MONTHS_YEAR   <- 6L      # sampled months for a yearly mean
RUNMEAN_MONTHS    <- 25L     # running-mean window for the slow share
N_SIM             <- 1000L   # AR(1) simulations per series
MIN_SEGMENT_YEARS <- 5L      # shortest segment between breaks
MAX_BREAKS        <- 2L
MAX_LAG           <- 36L
STEP_LATE         <- c(2020L, 2024L)
STEP_EARLY        <- c(2010L, 2019L)
KPSS_ALPHA        <- 0.05

script_stem <- "f2_slow_pattern_diagnostics"
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_div_a <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_div_g <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_abd_a <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
in_abd_g <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr", "forecast", "zoo")
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
stop_if_missing(c(in_settings, in_div_a, in_div_g, in_abd_a, in_abd_g),
                c("b09_forecast_settings.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 2), rep("e03_deseasonalize_abundance.R", 2)))
settings <- load_settings(in_settings)

series <- bind_rows(
  readRDS(in_div_a) %>% transmute(class = "diversity", scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, target = paste0("q", q),
                                  date, MONTH, filled, value, value_deseas),
  readRDS(in_div_g) %>% transmute(class = "diversity", scale = "coastwide", BAY_ID = NA_integer_, GEAR, target = paste0("q", q),
                                  date, MONTH, filled, value, value_deseas),
  readRDS(in_abd_a) %>% transmute(class = "abundance", scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, target = Common_Name,
                                  date, MONTH, filled, value, value_deseas),
  readRDS(in_abd_g) %>% transmute(class = "abundance", scale = "coastwide", BAY_ID = NA_integer_, GEAR, target = Common_Name,
                                  date, MONTH, filled, value, value_deseas)) %>%
  mutate(YEAR = as.integer(format(date, "%Y")),
         sampled = !filled & !(format(date, "%Y-%m") %in% format(settings$unscored_target_months, "%Y-%m")),
         dep = ifelse(sampled, value_deseas, NA_real_))
KEYS <- c("class", "scale", "BAY_ID", "GEAR", "target")
series$sid <- do.call(paste, c(as.list(series[, KEYS]), sep = "|"))
n_series <- length(unique(series$sid))
if (n_series != 176L) stop(n_series, " series; expected 176 (48 diversity, 128 abundance)", call. = FALSE)
if (any(series$filled & !(format(series$date, "%Y-%m") %in% format(settings$gap_months, "%Y-%m")))) {
  stop("Filled months other than the gap months", call. = FALSE)
}
cat("Series:", n_series, "| rows:", nrow(series), "\n")

# -----------------------------------------------------------------------------
# 3) Helpers
# -----------------------------------------------------------------------------
# Mann-Kendall test (via Kendall's tau against time) and Sen slope per year.
mann_kendall <- function(y, t) {
  ok <- is.finite(y); y <- y[ok]; t <- t[ok]
  ct <- suppressWarnings(stats::cor.test(t, y, method = "kendall", exact = FALSE))
  ij <- utils::combn(length(y), 2)
  sen <- stats::median((y[ij[2, ]] - y[ij[1, ]]) / (t[ij[2, ]] - t[ij[1, ]]))
  c(tau = unname(ct$estimate), p = ct$p.value, sen_per_year = sen)
}
# Mean-shift breakpoints by BIC, at most k_max breaks, segments >= h values.
# The BIC counts each segment mean, each break position and the variance, as
# strucchange::breakpoints does for a constant-only model.
breaks_bic <- function(y, h, k_max) {
  n <- length(y); cs <- c(0, cumsum(y)); cs2 <- c(0, cumsum(y^2))
  rss <- function(a, b) { m <- b - a + 1; s <- cs[b + 1] - cs[a]; cs2[b + 1] - cs2[a] - s^2 / m }
  bic <- function(r, k) n * log(r / n) + ((k + 1) + k + 1) * log(n)
  best <- list(k = 0L, at = integer(0), bic = bic(rss(1, n), 0L))
  if (k_max >= 1 && n >= 2 * h) for (b1 in h:(n - h)) {
    v <- bic(rss(1, b1) + rss(b1 + 1, n), 1L)
    if (v < best$bic) best <- list(k = 1L, at = b1, bic = v)
  }
  if (k_max >= 2 && n >= 3 * h) for (b1 in h:(n - 2 * h)) for (b2 in (b1 + h):(n - h)) {
    v <- bic(rss(1, b1) + rss(b1 + 1, b2) + rss(b2 + 1, n), 2L)
    if (v < best$bic) best <- list(k = 2L, at = c(b1, b2), bic = v)
  }
  best
}
# Share of variance in a centered running mean (series interpolated over gaps).
runmean_share <- function(x, w) {
  xi <- as.numeric(zoo::na.approx(x, na.rm = FALSE)); xi <- xi - mean(xi, na.rm = TRUE)
  rm <- as.numeric(stats::filter(xi, rep(1 / w, w), sides = 2))
  stats::var(rm, na.rm = TRUE) / stats::var(xi, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 4) Per series
# -----------------------------------------------------------------------------
set.seed(settings$seed)
ids <- sort(unique(series$sid))
res <- vector("list", length(ids)); acf_rows <- vector("list", length(ids)); yr_rows <- vector("list", length(ids))
for (i in seq_along(ids)) {
  d <- series[series$sid == ids[i], ]; d <- d[order(d$date), ]
  if (any(diff(d$YEAR * 12L + as.integer(format(d$date, "%m"))) != 1L)) stop(ids[i], ": not a complete monthly series", call. = FALSE)
  dep <- d$dep; sdd <- stats::sd(dep, na.rm = TRUE)
  # seasonal share
  seas_share <- 1 - stats::var(dep, na.rm = TRUE) / stats::var(d$value[d$sampled])
  # yearly means and trend
  yr <- d %>% filter(sampled) %>% group_by(YEAR) %>% summarise(n = n(), mean_dep = mean(dep), .groups = "drop") %>%
    filter(n >= MIN_MONTHS_YEAR)
  mk <- mann_kendall(yr$mean_dep, yr$YEAR)
  t_mo <- seq_along(dep) / 120
  lmf <- stats::lm(dep ~ t_mo)
  kp <- tryCatch(forecast::ndiffs(as.numeric(zoo::na.approx(dep, na.rm = FALSE))[is.finite(zoo::na.approx(dep, na.rm = FALSE))],
                                  alpha = KPSS_ALPHA, test = "kpss"), error = function(e) NA_integer_)
  # slow share against AR(1)
  ac <- stats::acf(dep, lag.max = MAX_LAG, plot = FALSE, na.action = stats::na.pass)$acf[-1]
  phi <- ac[1]
  slow <- runmean_share(dep, RUNMEAN_MONTHS)
  nn <- length(dep); miss <- is.na(dep)
  sims <- vapply(seq_len(N_SIM), function(s) {
    z <- as.numeric(stats::arima.sim(list(ar = phi), n = nn)); z[miss] <- NA; runmean_share(z, RUNMEAN_MONTHS)
  }, numeric(1))
  # breaks in the yearly means
  bk <- breaks_bic(yr$mean_dep, MIN_SEGMENT_YEARS, MAX_BREAKS)
  bk_years <- if (bk$k > 0) yr$YEAR[bk$at] + 1L else integer(0)   # first year of each new segment
  seg_id <- findInterval(seq_len(nrow(yr)), c(1, bk$at + 1L))
  seg_means <- tapply(yr$mean_dep, seg_id, mean)
  # 2020 to 2024 step
  late  <- mean(dep[d$YEAR >= STEP_LATE[1] & d$YEAR <= STEP_LATE[2]], na.rm = TRUE)
  early <- mean(dep[d$YEAR >= STEP_EARLY[1] & d$YEAR <= STEP_EARLY[2]], na.rm = TRUE)
  # seasonal swing by half of the record
  mid <- d$date[d$sampled][ceiling(sum(d$sampled) / 2)]
  sw <- lapply(list(first = d$sampled & d$date <= mid, second = d$sampled & d$date > mid), function(k) {
    cv <- curve_from_fit(fit_harmonic(d$MONTH[k], d$value[k]))
    c(swing = max(cv) - min(cv), peak = which.max(cv))
  })
  res[[i]] <- data.frame(d[1, KEYS], n_months = nn, n_sampled = sum(d$sampled), first_month = min(d$date), last_month = max(d$date),
    sd_departure = sdd, seasonal_share = seas_share,
    mk_tau = mk[["tau"]], mk_p = mk[["p"]], sen_per_decade = 10 * mk[["sen_per_year"]], sen_per_decade_sd = 10 * mk[["sen_per_year"]] / sdd,
    n_years = nrow(yr), trend_share = summary(lmf)$r.squared, trend_per_decade_sd = unname(stats::coef(lmf)[2]) / sdd,
    kpss_ndiffs = kp, acf1 = ac[1], acf12 = ac[12],
    slow_share = slow, slow_null_median = stats::median(sims), slow_null_q95 = unname(stats::quantile(sims, 0.95)),
    slow_ratio = slow / stats::median(sims), slow_p = mean(sims >= slow),
    n_breaks = bk$k, break_years = paste(bk_years, collapse = ";"),
    largest_shift_sd = if (bk$k > 0) max(abs(diff(as.numeric(seg_means)))) / sdd else 0,
    step_2020_sd = (late - early) / sdd,
    swing_first = sw$first[["swing"]], swing_second = sw$second[["swing"]], swing_ratio = sw$second[["swing"]] / sw$first[["swing"]],
    peak_first = sw$first[["peak"]], peak_second = sw$second[["peak"]], half_split = mid,
    stringsAsFactors = FALSE)
  acf_rows[[i]] <- data.frame(d[1, KEYS], lag = seq_len(MAX_LAG), acf = ac, row.names = NULL)
  yr_rows[[i]] <- data.frame(d[1, KEYS], yr, row.names = NULL)
}
diag_tab <- do.call(rbind, res); rownames(diag_tab) <- NULL
acf_tab <- do.call(rbind, acf_rows); yr_tab <- do.call(rbind, yr_rows)

# Checks: every series once, finite core measures.
if (nrow(diag_tab) != 176L || anyDuplicated(diag_tab[, KEYS])) stop("Diagnostics table does not have 176 unique series", call. = FALSE)
core <- c("seasonal_share", "mk_tau", "sen_per_decade_sd", "trend_share", "slow_share", "slow_null_q95", "acf1", "step_2020_sd", "swing_ratio")
bad <- sapply(core, function(v) sum(!is.finite(diag_tab[[v]])))
if (any(bad > 0)) { print(bad[bad > 0]); stop("Non-finite diagnostics (listed above)", call. = FALSE) }

# -----------------------------------------------------------------------------
# 5) Summary by group, write
# -----------------------------------------------------------------------------
q3 <- function(x) sprintf("%.2f [%.2f, %.2f]", stats::median(x, na.rm = TRUE), stats::quantile(x, 0.1, na.rm = TRUE), stats::quantile(x, 0.9, na.rm = TRUE))
summ <- diag_tab %>% group_by(class, scale, GEAR) %>% summarise(
  series = n(), seasonal_share = q3(seasonal_share),
  trend_mk_p05 = sum(mk_p < 0.05), trend_up = sum(mk_p < 0.05 & mk_tau > 0), trend_down = sum(mk_p < 0.05 & mk_tau < 0),
  sen_per_decade_sd = q3(sen_per_decade_sd), slow_above_ar1_q95 = sum(slow_share > slow_null_q95), slow_share = q3(slow_share),
  with_breaks = sum(n_breaks > 0), step_2020_sd = q3(step_2020_sd), swing_ratio = q3(swing_ratio),
  acf1 = q3(acf1), acf12 = q3(acf12), .groups = "drop")
acf_summ <- acf_tab %>% group_by(class, scale, lag) %>% summarise(median = stats::median(acf), q10 = stats::quantile(acf, 0.1),
                                                                   q90 = stats::quantile(acf, 0.9), .groups = "drop")

out_rds <- file.path(dirs$data, "f2_series_diagnostics.rds"); out_csv <- file.path(dirs$data, "f2_series_diagnostics.csv")
out_acf <- file.path(dirs$data, "f2_acf_by_lag.rds"); out_acf_csv <- file.path(dirs$data, "f2_acf_by_lag.csv")
out_yr  <- file.path(dirs$data, "f2_yearly_means.rds"); out_sum <- file.path(dirs$diag, "f2_summary_by_group.csv")
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE); if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|png|pdf)$", full.names = TRUE); if (length(old_d)) archive_files(old_d)
saveRDS(diag_tab, out_rds); readr::write_csv(diag_tab, out_csv, na = "")
saveRDS(acf_tab, out_acf); readr::write_csv(acf_tab, out_acf_csv, na = "")
saveRDS(yr_tab, out_yr); readr::write_csv(summ, out_sum, na = "")

cat("\nSummary by class, scale and gear (median [10th, 90th percentile]):\n"); print(as.data.frame(summ), row.names = FALSE)
cat("\nCorrelation of the departures by lag (median over series):\n")
print(as.data.frame(acf_summ %>% filter(lag %in% c(1, 2, 3, 6, 12, 24, 36)) %>% select(class, scale, lag, median) %>%
  mutate(median = round(median, 2)) %>% tidyr::pivot_wider(names_from = lag, values_from = median)), row.names = FALSE)
for (f in c(out_rds, out_csv, out_acf, out_acf_csv, out_yr, out_sum)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("f02_slow_pattern_diagnostics.R completed\n")
