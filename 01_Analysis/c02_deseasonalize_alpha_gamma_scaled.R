# =============================================================================
# Script:  c02_deseasonalize_alpha_gamma_scaled.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Removes the seasonal cycle from every diversity series and writes the
#   departures from the seasonal curve (value_deseas).
#   - Coastwide (gamma, 6 series): the curve fitted in c01 is subtracted.
#   - Estuaries (alpha, 42 series): each estuary, sampling method and q gets
#     its own harmonic curve, fitted here, and that curve is subtracted.
#   In both cases the full fitted monthly mean (intercept plus curve) is
#   subtracted, so each series of departures is centered on zero. ("_scaled"
#   in the names means estuary-specific curves.)
#   These full-record departures describe the series. The forecasts use
#   curves refitted at every possible forecast date on the months up to that
#   date: their coefficients are written to seasonal_coefs_by_origin_scaled.rds,
#   and the model scripts rebuild the departures for a forecast date with
#   anomalies_at_origin() (lib/seasonal.R).
#
# Run after: c01, b09, lib/seasonal.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   01_Analysis/Data/b9_forecast_settings/alpha_q012_filled.rds (b09, from b06)
#   01_Analysis/Data/b9_forecast_settings/gamma_q012_filled.rds (b09, from b07)
#   01_Analysis/Data/c1_seasonality_characterization/seasonality_fits.rds (c01)
#
# Outputs:
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/
#     alpha_q012_deseas_scaled.rds          BAY_ID, GEAR, q, q_label, YEAR,
#                                           MONTH, date, filled, value,
#                                           seasonal_component,
#                                           residual_component, value_deseas
#     gamma_q012_deseas_scaled.rds          the same without BAY_ID
#     seasonal_map_scaled.rds               coastwide curves (GEAR, q, MONTH,
#                                           seasonal)
#     seasonal_map_alpha_bay_scaled.rds     estuary curves (BAY_ID, GEAR, q,
#                                           MONTH, seasonal)
#     seasonal_coefs_by_origin_scaled.rds   curves at each forecast date:
#                                           scale, BAY_ID (NA coastwide),
#                                           GEAR, q, origin_i, origin_date,
#                                           allowed_origin, n_train, b0, sin1,
#                                           cos1, ..., cos3
#     The series of departures are read by c03, c07 to c10, c12, d01, f01 to
#     f03, f07, g11, h13 and h14; the per-forecast-date curves by c03, c07 to
#     c10, c12 and d01.
#   02_Results/03_Diagnostics/c2_deseasonalize_alpha_gamma_scaled/
#     alpha_bay_harmonic_fits_scaled.csv        one row per estuary series (h14)
#     alpha_residual_seasonality_scaled.csv     own curve, coastwide curve,
#                                               2-harmonic comparison (h14)
#     alpha_residual_seasonality_summary_scaled.csv
#     alpha_residual_acf_scaled.csv             lags 1 and 12, own and
#                                               coastwide curve (h14)
#     seasonal_by_origin_checks_scaled.csv      per series: forecast dates,
#                                               last curve against the full
#                                               record, largest change between
#                                               consecutive forecast dates (h14)
#   02_Results/01_Figures/c2_deseasonalize_alpha_gamma_scaled/
#     alpha_bay_seasonal_curves_scaled.pdf and .tiff  estuary curves and the
#                                                     coastwide curve
#   01_Analysis/Logs/c2_deseasonalize_alpha_gamma_scaled/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   - Model per series, as in c01: value ~ sin1 + cos1 + sin2 + cos2 + sin3 +
#     cos3 on calendar month (lib/seasonal.R), natural scale, no year or
#     trend term, lm() on the observed (sampled) months: of the whole record
#     for the description, of the months up to each forecast date for the
#     forecasts.
#   - Curves are fitted on observed months only (filled == FALSE). April and
#     May 2020 keep their rows, with departures from the filled values and a
#     `filled` column, so every model gets the complete series and never
#     fills.
#   - The number of harmonics comes from lib/seasonal.R; the script stops if
#     the c01 fits use another number. Two harmonics are kept as a comparison
#     (month_r2_sens2). The chance level of leftover seasonality after a
#     K-harmonic fit is (11 - 2K) / (n - 1).
#   - Leftover seasonality and residual autocorrelation are computed with the
#     filled months set to NA on the complete monthly grid, so lags are
#     calendar months. The leftover seasonality when the coastwide curve is
#     subtracted from an estuary series is also reported (columns
#     *_coastwide_curve).
#   - Curves for the forecasts are fitted inside the training window, at
#     every month that could be a forecast date (from the first month with
#     settings$min_train_months sampled months to the month before the last;
#     see lib/seasonal.R), for estuary and coastwide series alike. Months in
#     settings$skip_origin_months are flagged allowed_origin = FALSE. Which
#     forecast dates are used, and the training start, are set in b09 and c03.
#   - Estuaries and sampling methods come from settings.rds.
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
script_stem <- "c2_deseasonalize_alpha_gamma_scaled"   # output folder stem (short form c2_)
SUFFIX      <- "_scaled"

N_HARMONICS    <- SEASONAL_N_HARMONICS   # 3, from lib/seasonal.R
SENS_HARMONICS <- 2L    # comparison only; not used for the output
MIN_OBS        <- 36L   # minimum observed months per bay series
MIN_PER_MONTH  <- 2L    # every calendar month observed at least this often
Q_LEVELS       <- 0:2

in_settings  <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha     <- file.path(PATHS$data_dir, "b9_forecast_settings", "alpha_q012_filled.rds")
in_gamma     <- file.path(PATHS$data_dir, "b9_forecast_settings", "gamma_q012_filled.rds")
in_seas_fits <- file.path(PATHS$data_dir, "c1_seasonality_characterization", "seasonality_fits.rds")

options(scipen = 999, warn = 1, width = 200)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "tibble", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag", "figs"))
archive_previous(c(dirs$data, dirs$diag, dirs$figs))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")
cat("Harmonics:", N_HARMONICS, "(sensitivity", SENS_HARMONICS, ")\n")

out_alpha       <- file.path(dirs$data, paste0("alpha_q012_deseas", SUFFIX, ".rds"))
out_gamma       <- file.path(dirs$data, paste0("gamma_q012_deseas", SUFFIX, ".rds"))
out_seasmap     <- file.path(dirs$data, paste0("seasonal_map", SUFFIX, ".rds"))
out_seasmap_bay <- file.path(dirs$data, paste0("seasonal_map_alpha_bay", SUFFIX, ".rds"))
out_coefs       <- file.path(dirs$data, paste0("seasonal_coefs_by_origin", SUFFIX, ".rds"))
out_origin_csv  <- file.path(dirs$diag, paste0("seasonal_by_origin_checks", SUFFIX, ".csv"))
out_fits_csv    <- file.path(dirs$diag, paste0("alpha_bay_harmonic_fits", SUFFIX, ".csv"))
out_resid_csv   <- file.path(dirs$diag, paste0("alpha_residual_seasonality", SUFFIX, ".csv"))
out_ressum_csv  <- file.path(dirs$diag, paste0("alpha_residual_seasonality_summary", SUFFIX, ".csv"))
out_acf_csv     <- file.path(dirs$diag, paste0("alpha_residual_acf", SUFFIX, ".csv"))
out_fig_base    <- file.path(dirs$figs, paste0("alpha_bay_seasonal_curves", SUFFIX))

# -----------------------------------------------------------------------------
# 4) Utilities
# -----------------------------------------------------------------------------
q_label <- function(q) {
  dplyr::case_when(
    q == 0 ~ "q = 0 (richness)",
    q == 1 ~ "q = 1 (Shannon-weighted diversity)",
    q == 2 ~ "q = 2 (dominance-weighted diversity)",
    TRUE   ~ paste0("q = ", q)
  )
}

# Wide b09 file (q0, q1, q2) to long: [BAY_ID], GEAR, q, YEAR, MONTH, filled, value
to_long <- function(df, name) {
  expect_cols(df, c("GEAR", "YEAR", "MONTH", "filled", paste0("q", Q_LEVELS)), name)
  df %>%
    select(any_of("BAY_ID"), GEAR, YEAR, MONTH, filled, all_of(paste0("q", Q_LEVELS))) %>%
    mutate(GEAR = as.character(GEAR), YEAR = as.integer(YEAR), MONTH = as.integer(MONTH)) %>%
    pivot_longer(cols = all_of(paste0("q", Q_LEVELS)), names_to = "q_raw", values_to = "value") %>%
    mutate(q = as.integer(sub("^q", "", q_raw)), value = as.numeric(value)) %>%
    select(-q_raw)
}

# harmonic_X(), fit_harmonic(), curve_from_fit(), seasonal_coefs_by_origin()
# and seasonal_level() come from lib/seasonal.R.

# Share of variance explained by calendar-month means (leftover seasonality).
# Chance level after a K-harmonic fit: (11 - 2K) / (n - 1). NA values (the
# filled months) are left out.
month_r2 <- function(x, month) {
  ok <- is.finite(x)
  x <- x[ok]; month <- month[ok]
  if (length(x) < 13) return(NA_real_)
  mm  <- ave(x, month, FUN = mean)
  tot <- sum((x - mean(x))^2)
  if (tot <= 0) return(NA_real_)
  1 - sum((x - mm)^2) / tot
}

# ACF on the complete monthly grid with the filled months as NA.
acf_grid <- function(x, lags = c(1, 12)) {
  r <- stats::acf(x, lag.max = max(lags), plot = FALSE, na.action = stats::na.pass)
  setNames(as.numeric(r$acf[lags + 1]), paste0("acf", lags))
}

# Coastwide seasonal map from the c01 fits (fitted monthly mean, value_hat).
build_seasonal_map <- function(seas_fits) {
  nm <- names(seas_fits)
  if (is.null(nm) || any(!nzchar(nm))) stop("seasonality_fits must be a named list like 'BS__q0'.", call. = FALSE)
  rows <- lapply(seq_along(seas_fits), function(i) {
    fit <- seas_fits[[i]]
    if (!isTRUE(fit$ok)) stop("c01 fit ", nm[[i]], " is not ok", call. = FALSE)
    k_c01 <- (length(stats::coef(fit$model)) - 1L) / 2L
    if (k_c01 != N_HARMONICS) stop("c01 fit ", nm[[i]], " has ", k_c01, " harmonics; c02 uses ", N_HARMONICS, call. = FALSE)
    parts <- strsplit(nm[[i]], "__q", fixed = TRUE)[[1]]
    sc <- fit$seasonal_curve
    expect_cols(sc, c("MONTH", "value_hat"), paste("seasonal_curve of", nm[[i]]))
    tibble::tibble(GEAR = parts[1], q = as.integer(parts[2]),
                   MONTH = as.integer(sc$MONTH), seasonal = as.numeric(sc$value_hat))
  })
  out <- bind_rows(rows)
  if (anyDuplicated(out[, c("GEAR", "q", "MONTH")])) stop("Coastwide seasonal map has duplicated GEAR x q x MONTH keys.", call. = FALSE)
  out
}

finish_ts <- function(df) {
  df %>%
    mutate(date = month_date(YEAR, MONTH),
           q_label = q_label(q),
           seasonal_component = seasonal,
           residual_component = value - seasonal_component,
           value_deseas = residual_component) %>%
    arrange(GEAR, q, YEAR, MONTH) %>%
    select(any_of("BAY_ID"), GEAR, q, q_label, YEAR, MONTH, date, filled,
           value, seasonal_component, residual_component, value_deseas)
}

# -----------------------------------------------------------------------------
# 5) Inputs and checks
# -----------------------------------------------------------------------------
stop_if_missing(c(in_settings, in_alpha, in_gamma, in_seas_fits),
                c("b09_forecast_settings.R", "b09_forecast_settings.R", "b09_forecast_settings.R",
                  "c01_seasonality_characterization.R"))
settings  <- load_settings(in_settings)
alpha_ts  <- to_long(readRDS(in_alpha), "alpha_q012_filled")
gamma_ts  <- to_long(readRDS(in_gamma), "gamma_q012_filled")
seas_fits <- readRDS(in_seas_fits)

if (!"BAY_ID" %in% names(alpha_ts)) stop("alpha input has no BAY_ID column.", call. = FALSE)
alpha_ts$BAY_ID <- as.integer(alpha_ts$BAY_ID)
if (anyDuplicated(alpha_ts[, c("BAY_ID", "GEAR", "q", "YEAR", "MONTH")])) stop("alpha input has duplicated keys.", call. = FALSE)
if (anyDuplicated(gamma_ts[, c("GEAR", "q", "YEAR", "MONTH")])) stop("gamma input has duplicated keys.", call. = FALSE)
if (!setequal(unique(alpha_ts$BAY_ID), settings$bays)) stop("Bays in alpha differ from settings$bays", call. = FALSE)
if (!setequal(unique(alpha_ts$GEAR), settings$gears) || !setequal(unique(gamma_ts$GEAR), settings$gears)) {
  stop("Gears differ from settings$gears", call. = FALSE)
}
n_alpha_series <- nrow(distinct(alpha_ts, BAY_ID, GEAR, q))
n_gamma_series <- nrow(distinct(gamma_ts, GEAR, q))
if (n_alpha_series != length(settings$bays) * length(settings$gears) * length(Q_LEVELS)) {
  stop("alpha has ", n_alpha_series, " series; expected ", length(settings$bays) * length(settings$gears) * length(Q_LEVELS), call. = FALSE)
}
if (n_gamma_series != length(settings$gears) * length(Q_LEVELS)) {
  stop("gamma has ", n_gamma_series, " series; expected ", length(settings$gears) * length(Q_LEVELS), call. = FALSE)
}
cat("\nAlpha rows:", nrow(alpha_ts), "| series:", n_alpha_series, "| filled rows:", sum(alpha_ts$filled), "\n")
cat("Gamma rows:", nrow(gamma_ts), "| series:", n_gamma_series, "| filled rows:", sum(gamma_ts$filled), "\n\n")

# -----------------------------------------------------------------------------
# 6) Coastwide: the c01 curve
# -----------------------------------------------------------------------------
seasonal_map <- build_seasonal_map(seas_fits)
if (nrow(seasonal_map) != n_gamma_series * 12L) stop("Coastwide seasonal map has ", nrow(seasonal_map), " rows; expected ", n_gamma_series * 12L, call. = FALSE)
gamma_deseas <- gamma_ts %>%
  left_join(seasonal_map, by = c("GEAR", "q", "MONTH")) %>%
  finish_ts()

# -----------------------------------------------------------------------------
# 7) Bays: one harmonic fit per BAY_ID x GEAR x q, on observed months
# -----------------------------------------------------------------------------
keys <- alpha_ts %>% distinct(BAY_ID, GEAR, q) %>% arrange(BAY_ID, GEAR, q)
fit_rows <- vector("list", nrow(keys))
map_rows <- vector("list", nrow(keys))

for (i in seq_len(nrow(keys))) {
  k   <- keys[i, ]
  obs <- alpha_ts %>% filter(BAY_ID == k$BAY_ID, GEAR == k$GEAR, q == k$q, !filled, is.finite(value))
  per_month <- tabulate(obs$MONTH, nbins = 12)

  tag <- sprintf("BAY %s, %s, q = %s", k$BAY_ID, k$GEAR, k$q)
  if (nrow(obs) < MIN_OBS) stop(tag, ": only ", nrow(obs), " observed months (minimum ", MIN_OBS, ").", call. = FALSE)
  if (any(per_month < MIN_PER_MONTH)) {
    stop(tag, ": calendar month(s) ", paste(which(per_month < MIN_PER_MONTH), collapse = ", "),
         " observed fewer than ", MIN_PER_MONTH, " times.", call. = FALSE)
  }

  fit <- fit_harmonic(obs$MONTH, obs$value, N_HARMONICS)
  co <- stats::coef(fit)
  if (any(!is.finite(co))) stop(tag, ": non-finite harmonic coefficients.", call. = FALSE)
  curve <- curve_from_fit(fit, N_HARMONICS)

  fit_rows[[i]] <- tibble::tibble(
    BAY_ID = k$BAY_ID, GEAR = k$GEAR, q = k$q,
    n_obs = nrow(obs), min_obs_per_month = min(per_month),
    intercept = unname(co["(Intercept)"]),
    amp_annual = sqrt(unname(co["sin1"])^2 + unname(co["cos1"])^2),
    amp_semiannual = sqrt(unname(co["sin2"])^2 + unname(co["cos2"])^2),
    peak_to_trough = max(curve) - min(curve),
    peak_month = which.max(curve), trough_month = which.min(curve),
    r2_seasonal = summary(fit)$r.squared,
    rel_amp_annual = sqrt(unname(co["sin1"])^2 + unname(co["cos1"])^2) / unname(co["(Intercept)"])
  )
  map_rows[[i]] <- tibble::tibble(BAY_ID = k$BAY_ID, GEAR = k$GEAR, q = k$q, MONTH = 1:12, seasonal = curve)
}

bay_fits         <- bind_rows(fit_rows)
seasonal_map_bay <- bind_rows(map_rows)

alpha_deseas <- alpha_ts %>%
  left_join(seasonal_map_bay, by = c("BAY_ID", "GEAR", "q", "MONTH")) %>%
  finish_ts()

# -----------------------------------------------------------------------------
# 8) Output checks
# -----------------------------------------------------------------------------
if (nrow(alpha_deseas) != nrow(alpha_ts) || nrow(gamma_deseas) != nrow(gamma_ts)) stop("Row count changed during deseasonalization.", call. = FALSE)
for (d in list(alpha = alpha_deseas, gamma = gamma_deseas)) {
  if (anyNA(d$seasonal_component)) stop("Missing seasonal component in some rows.", call. = FALSE)
  if (anyNA(d$value_deseas)) stop("Missing anomaly in some rows (the b09 series are complete).", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 9) Diagnostics: leftover seasonality (own curve, coastwide curve, 2 harmonics)
# -----------------------------------------------------------------------------
# All on observed months; the filled months are NA on the monthly grid.
coast_on_bay <- alpha_ts %>%
  left_join(seasonal_map, by = c("GEAR", "q", "MONTH")) %>%
  mutate(anom_coast = value - seasonal) %>%
  select(BAY_ID, GEAR, q, YEAR, MONTH, anom_coast)
if (anyNA(coast_on_bay$anom_coast)) stop("Coastwide curve missing for some bay rows.", call. = FALSE)

resid_rows <- vector("list", nrow(keys))
acf_rows   <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  g <- alpha_deseas %>% filter(BAY_ID == k$BAY_ID, GEAR == k$GEAR, q == k$q) %>% arrange(YEAR, MONTH) %>%
    left_join(coast_on_bay, by = c("BAY_ID", "GEAR", "q", "YEAR", "MONTH"))
  own   <- ifelse(g$filled, NA_real_, g$value_deseas)
  coast <- ifelse(g$filled, NA_real_, g$anom_coast)

  obs  <- g %>% filter(!filled)
  fit2 <- fit_harmonic(obs$MONTH, obs$value, SENS_HARMONICS)
  res2 <- obs$value - curve_from_fit(fit2, SENS_HARMONICS)[obs$MONTH]

  n_ok <- sum(is.finite(own))
  resid_rows[[i]] <- tibble::tibble(
    BAY_ID = k$BAY_ID, GEAR = k$GEAR, q = k$q, n_obs = n_ok,
    chance_level = (11 - 2 * N_HARMONICS) / (n_ok - 1),
    month_r2_new = month_r2(own, g$MONTH),
    month_r2_sens2 = month_r2(res2, obs$MONTH),
    chance_level_sens2 = (11 - 2 * SENS_HARMONICS) / (n_ok - 1),
    mean_anom_new = mean(own, na.rm = TRUE),
    sd_anom_new = stats::sd(own, na.rm = TRUE),
    month_r2_coastwide_curve = month_r2(coast, g$MONTH),
    mean_anom_coastwide_curve = mean(coast, na.rm = TRUE),
    sd_anom_coastwide_curve = stats::sd(coast, na.rm = TRUE)
  )
  a_own <- acf_grid(own)
  a_coast <- acf_grid(coast)
  acf_rows[[i]] <- tibble::tibble(BAY_ID = k$BAY_ID, GEAR = k$GEAR, q = k$q,
                                  acf1_new = a_own[["acf1"]], acf12_new = a_own[["acf12"]],
                                  acf1_coastwide_curve = a_coast[["acf1"]], acf12_coastwide_curve = a_coast[["acf12"]])
}
resid_diag <- bind_rows(resid_rows)
acf_diag   <- bind_rows(acf_rows)

resid_summary <- resid_diag %>%
  group_by(GEAR, q) %>%
  summarise(n_series = n(),
            across(c(month_r2_coastwide_curve, month_r2_new, month_r2_sens2,
                     mean_anom_coastwide_curve, mean_anom_new),
                   list(mean = ~ mean(.x), max = ~ max(.x))),
            chance_level = mean(chance_level),
            .groups = "drop")

# -----------------------------------------------------------------------------
# 10) Curves at every possible forecast date
# -----------------------------------------------------------------------------
# One coefficient row per series and origin month. Checks per series: the
# curve at the last origin against the full-record curve, and the largest
# change of the curve between consecutive origins, both in units of the SD of
# the full-record anomalies.
by_origin <- function(df, key, full_curve) {
  df <- df %>% arrange(YEAR, MONTH)
  co <- seasonal_coefs_by_origin(df$date, df$MONTH, df$value, sampled = !df$filled,
                                 min_train = settings$min_train_months, K = N_HARMONICS)
  co$allowed_origin <- is_allowed_origin(co$origin_date, settings)
  B <- as.matrix(co[, seasonal_coef_names(N_HARMONICS)])
  curves <- B %*% t(as.matrix(harmonic_X(1:12, N_HARMONICS)))       # origins x 12 months
  sd_anom <- stats::sd(df$value_deseas[!df$filled])
  # Training anomalies average zero at every origin (see lib/seasonal.R)
  tr_mean <- vapply(c(1L, nrow(co)), function(j) {
    rows <- seq_len(co$origin_i[j]); rows <- rows[!df$filled[rows]]
    mean(df$value[rows] - seasonal_level(B[j, ], df$MONTH[rows]))
  }, numeric(1))
  check <- tibble::tibble(
    n_origins = nrow(co), first_origin = min(co$origin_date), last_origin = max(co$origin_date),
    last_vs_full_sd = max(abs(curves[nrow(curves), ] - full_curve)) / sd_anom,
    max_step_sd = if (nrow(curves) > 1) max(abs(diff(curves))) / sd_anom else NA_real_,
    median_step_sd = if (nrow(curves) > 1) stats::median(apply(abs(diff(curves)), 1, max)) / sd_anom else NA_real_,
    train_mean_abs_max = max(abs(tr_mean))
  )
  list(coefs = bind_cols(key, co), check = bind_cols(key, check))
}

origin_parts <- list()
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  df <- alpha_deseas %>% filter(BAY_ID == k$BAY_ID, GEAR == k$GEAR, q == k$q)
  full <- seasonal_map_bay %>% filter(BAY_ID == k$BAY_ID, GEAR == k$GEAR, q == k$q) %>% arrange(MONTH)
  origin_parts[[length(origin_parts) + 1]] <- by_origin(df, tibble::tibble(scale = "bay", BAY_ID = k$BAY_ID, GEAR = k$GEAR, q = k$q), full$seasonal)
}
gkeys <- gamma_deseas %>% distinct(GEAR, q) %>% arrange(GEAR, q)
for (i in seq_len(nrow(gkeys))) {
  k <- gkeys[i, ]
  df <- gamma_deseas %>% filter(GEAR == k$GEAR, q == k$q)
  full <- seasonal_map %>% filter(GEAR == k$GEAR, q == k$q) %>% arrange(MONTH)
  origin_parts[[length(origin_parts) + 1]] <- by_origin(df, tibble::tibble(scale = "coastwide", BAY_ID = NA_integer_, GEAR = k$GEAR, q = k$q), full$seasonal)
}
coefs_by_origin <- bind_rows(lapply(origin_parts, `[[`, "coefs")) %>%
  select(scale, BAY_ID, GEAR, q, origin_i, origin_date, allowed_origin, n_train, all_of(seasonal_coef_names(N_HARMONICS)))
origin_checks <- bind_rows(lapply(origin_parts, `[[`, "check"))
if (nrow(origin_checks) != n_alpha_series + n_gamma_series) stop("Per-origin curves missing for some series.", call. = FALSE)
if (any(origin_checks$train_mean_abs_max > 1e-8)) stop("Training anomalies do not average zero at some origin.", call. = FALSE)
if (anyNA(select(coefs_by_origin, -BAY_ID))) stop("NA in the per-origin coefficients.", call. = FALSE)

# -----------------------------------------------------------------------------
# 11) Write
# -----------------------------------------------------------------------------
saveRDS(alpha_deseas, out_alpha)
saveRDS(gamma_deseas, out_gamma)
saveRDS(seasonal_map, out_seasmap)
saveRDS(seasonal_map_bay, out_seasmap_bay)
saveRDS(coefs_by_origin, out_coefs)
readr::write_csv(origin_checks, out_origin_csv)
readr::write_csv(bay_fits, out_fits_csv)
readr::write_csv(resid_diag, out_resid_csv)
readr::write_csv(resid_summary, out_ressum_csv)
readr::write_csv(acf_diag, out_acf_csv)

# Figure: bay curves against the coastwide curve (diagnostic, not a manuscript figure)
plot_df <- seasonal_map_bay %>%
  mutate(series = paste("Bay", BAY_ID)) %>%
  bind_rows(seasonal_map %>% mutate(BAY_ID = NA_integer_, series = "Coastwide")) %>%
  mutate(q_lab = paste0("q = ", q), is_coast = series == "Coastwide")
p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MONTH, y = seasonal, group = series)) +
  ggplot2::geom_line(data = dplyr::filter(plot_df, !is_coast), colour = "grey55", linewidth = 0.4) +
  ggplot2::geom_line(data = dplyr::filter(plot_df, is_coast), colour = "black", linewidth = 1.0) +
  ggplot2::facet_grid(q_lab ~ GEAR, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = 1:12) +
  ggplot2::labs(x = "Month", y = "Fitted seasonal level (effective species)",
                title = "Bay seasonal curves (grey) against the coastwide curve (black)") +
  ggplot2::theme_bw(base_size = 9)
W_MM <- 170; H_MM <- 180
stopifnot(H_MM / W_MM <= 8.0 / 6.5)
ggplot2::ggsave(paste0(out_fig_base, ".pdf"), p, width = W_MM, height = H_MM, units = "mm")
ggplot2::ggsave(paste0(out_fig_base, ".tiff"), p, width = W_MM, height = H_MM, units = "mm",
                dpi = 300, compression = "lzw")

# -----------------------------------------------------------------------------
# 12) Console summary
# -----------------------------------------------------------------------------
cat("\nLeftover seasonality in bay diversity anomalies (share of variance explained by\n")
cat("calendar-month means, observed months; chance is about ", sprintf("%.3f", mean(resid_diag$chance_level)), "):\n", sep = "")
print(as.data.frame(resid_summary %>% select(GEAR, q, month_r2_coastwide_curve_mean, month_r2_new_mean,
                                             month_r2_new_max, month_r2_sens2_mean, mean_anom_coastwide_curve_mean) %>%
                      mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)

cat("\nCurves at every possible forecast origin:", nrow(coefs_by_origin), "rows,",
    nrow(origin_checks), "series; origins per series", min(origin_checks$n_origins), "to", max(origin_checks$n_origins), "\n")
cat("Last origin against the full-record curve (max |difference| / anomaly SD): median",
    round(stats::median(origin_checks$last_vs_full_sd), 4), ", max", round(max(origin_checks$last_vs_full_sd), 4), "\n")
cat("Largest change between consecutive origins (/ anomaly SD): median of series maxima",
    round(stats::median(origin_checks$max_step_sd), 4), ", max", round(max(origin_checks$max_step_sd), 4), "\n")

cat("\nSaved:\n")
for (f in c(out_alpha, out_gamma, out_seasmap, out_seasmap_bay, out_coefs, out_fits_csv, out_resid_csv,
            out_ressum_csv, out_acf_csv, out_origin_csv, paste0(out_fig_base, c(".pdf", ".tiff")))) cat("  ", f, "\n", sep = "")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("c02_deseasonalize_alpha_gamma_scaled.R completed\n")
