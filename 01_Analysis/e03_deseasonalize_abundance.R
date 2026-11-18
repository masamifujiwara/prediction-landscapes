# =============================================================================
# Script:  e03_deseasonalize_abundance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Removes the seasonal cycle from every abundance series (8 targets, per
#   estuary and coastwide, both sampling methods) and writes the departures
#   from the seasonal curve (value_deseas) that the abundance forecasting
#   models use as targets. Each series gets its own harmonic curve, fitted on
#   its sampled months; the full fitted monthly mean (intercept plus curve) is
#   subtracted, so each anomaly series is centered on zero, as for diversity
#   in c02. These full-record anomalies describe the series; the forecasts use
#   curves refitted at every possible forecast date on the months up to that
#   date, whose coefficients are written to
#   seasonal_coefs_by_origin_abundance.rds (rebuild the anomalies for a
#   forecast date with anomalies_at_origin() in lib/seasonal.R).
#
# Run after: e02, b09, lib/seasonal.R
#
# Inputs:
#   01_Analysis/Data/e2_aggregate_abundance/alpha_abundance_targets.rds (e02)
#   01_Analysis/Data/e2_aggregate_abundance/gamma_abundance_targets.rds (e02)
#   01_Analysis/Data/b9_forecast_settings/settings.rds                  (b09)
#
# Outputs:
#   01_Analysis/Data/e3_deseasonalize_abundance/
#     alpha_abundance_deseas.rds   BAY_ID, GEAR, SPECIES_CODE, Common_Name,
#                                  YEAR, MONTH, date, filled, n_samples,
#                                  mean_cpue, value, seasonal_component,
#                                  residual_component, value_deseas
#     gamma_abundance_deseas.rds   the same without BAY_ID, plus n_bays
#     seasonal_map_abundance.rds   fitted monthly level per series
#                                  (scale, BAY_ID, GEAR, SPECIES_CODE, MONTH,
#                                  seasonal)
#     seasonal_coefs_by_origin_abundance.rds  curves at each forecast date:
#                                  scale, BAY_ID, GEAR, SPECIES_CODE,
#                                  origin_i, origin_date, allowed_origin,
#                                  n_train, b0, sin1, ..., cos3
#   02_Results/03_Diagnostics/e3_deseasonalize_abundance/
#     seasonality_summary.csv          per scale, sampling method and target:
#                                      amplitudes, R2, leftover seasonality
#                                      (mean, max)
#     e3_harmonic_fits.csv             one row per series
#     e3_residual_seasonality.csv      one row per series: leftover
#                                      seasonality (3 and 2 harmonics),
#                                      chance level, lag-1 and lag-12 ACF
#     e3_seasonal_by_origin_checks.csv one row per series: forecast dates,
#                                      last date against the full record,
#                                      largest change between dates
#   01_Analysis/Logs/e3_deseasonalize_abundance/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   Model per series: value ~ harmonics 1 to 3 of calendar month
#   (lib/seasonal.R), where value = log(CPUE + 1) from e02; lm() on the
#   sampled months (filled == FALSE), of the whole record for the description
#   and of the months up to each forecast date for the forecasts.
#   1. The curve is fitted on log(CPUE + 1), the scale of the forecast
#      target; raw CPUE fits a harmonic curve poorly.
#   2. Three harmonics (periods of 12, 6 and 4 months), as for diversity; the
#      number is set in lib/seasonal.R. With two harmonics, calendar month
#      still explained more of the anomaly variance than expected by chance;
#      with three it did not. Three harmonics also forecast as well as two at
#      the earliest forecast dates, so there is no special rule for short
#      training windows. The two-harmonic result is written as a comparison.
#   3. The full fitted monthly mean is subtracted, as in c02, so the anomalies
#      do not keep the series level. seasonal_component holds the fitted
#      monthly level.
#   4. Curves are fitted on sampled months only (filled == FALSE). The filled
#      months (April and May 2020, from e02) keep their rows, with anomalies
#      and the `filled` column, so every model gets the complete series and
#      never fills. Sampled months include months with zero catch.
#   5. Diagnostics over months (leftover seasonality, ACF) use the complete
#      monthly grid with the filled months as NA, so lags are calendar months.
#      The chance level of leftover seasonality after a K-harmonic fit is
#      (11 - 2K) / (n - 1), not 11 / (n - 1).
#   6. No year or trend term, as in c01 and c02. Curves for the forecasts are
#      fitted inside the training window at every month that could be a
#      forecast date (lib/seasonal.R); dates in settings$skip_origin_months
#      are flagged allowed_origin = FALSE. The full-record fit is kept for the
#      description.
#   7. The script stops if a series has fewer than 36 sampled months, a
#      calendar month sampled fewer than twice, or non-finite coefficients.
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
script_stem <- "e3_deseasonalize_abundance"   # stem of the output folders

N_HARMONICS    <- SEASONAL_N_HARMONICS   # 3, from lib/seasonal.R (method note 2)
SENS_HARMONICS <- 2L    # comparison only; not used for the output
MIN_OBS        <- 36L   # minimum sampled months per series (as c02)
MIN_PER_MONTH  <- 2L    # every calendar month sampled at least this often (as c02)

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "e2_aggregate_abundance", "alpha_abundance_targets.rds")
in_gamma    <- file.path(PATHS$data_dir, "e2_aggregate_abundance", "gamma_abundance_targets.rds")

options(scipen = 999, warn = 1, width = 200)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")
cat("Harmonics:", N_HARMONICS, "(comparison", SENS_HARMONICS, ")\n")

out_alpha  <- file.path(dirs$data, "alpha_abundance_deseas.rds")
out_gamma  <- file.path(dirs$data, "gamma_abundance_deseas.rds")
out_map    <- file.path(dirs$data, "seasonal_map_abundance.rds")
out_summ   <- file.path(dirs$diag, "seasonality_summary.csv")
out_fits   <- file.path(dirs$diag, "e3_harmonic_fits.csv")
out_resid  <- file.path(dirs$diag, "e3_residual_seasonality.csv")
out_coefs  <- file.path(dirs$data, "seasonal_coefs_by_origin_abundance.rds")
out_origin <- file.path(dirs$diag, "e3_seasonal_by_origin_checks.csv")

# -----------------------------------------------------------------------------
# 4) Utilities
# -----------------------------------------------------------------------------
# harmonic_X(), fit_harmonic(), curve_from_fit(), seasonal_coefs_by_origin()
# and seasonal_level() come from lib/seasonal.R.

# Share of variance explained by calendar-month means (leftover seasonality).
# NA values (the filled months) are left out. Chance level: (11 - 2K) / (n - 1).
month_r2 <- function(x, month) {
  ok <- is.finite(x)
  x <- x[ok]; month <- month[ok]
  if (length(x) < 13) return(NA_real_)
  mm  <- ave(x, month, FUN = mean)
  tot <- sum((x - mean(x))^2)
  if (tot <= 0) return(NA_real_)
  1 - sum((x - mm)^2) / tot
}

# ACF on the complete monthly grid with the filled months as NA (method note 5).
acf_grid <- function(x, lags = c(1, 12)) {
  r <- stats::acf(x, lag.max = max(lags), plot = FALSE, na.action = stats::na.pass)
  setNames(as.numeric(r$acf[lags + 1]), paste0("acf", lags))
}

amp <- function(co, k) sqrt(unname(co[paste0("sin", k)])^2 + unname(co[paste0("cos", k)])^2)

# -----------------------------------------------------------------------------
# 5) Inputs and checks
# -----------------------------------------------------------------------------
stop_if_missing(c(in_settings, in_alpha, in_gamma),
                c("b09_forecast_settings.R", "e02_aggregate_abundance.R", "e02_aggregate_abundance.R"))
settings <- load_settings(in_settings)
key_cols <- c("GEAR", "YEAR", "MONTH", "date", "SPECIES_CODE", "Common_Name", "filled", "n_samples", "mean_cpue", "log_cpue")
alpha_in <- readRDS(in_alpha); expect_cols(alpha_in, c("BAY_ID", key_cols), "alpha_abundance_targets")
gamma_in <- readRDS(in_gamma); expect_cols(gamma_in, c(key_cols, "n_bays"), "gamma_abundance_targets")

ts_all <- bind_rows(
  alpha_in %>% mutate(scale = "bay", BAY_ID = as.integer(BAY_ID)),
  gamma_in %>% mutate(scale = "coastwide", BAY_ID = NA_integer_)
) %>%
  mutate(value = log_cpue)

if (!setequal(unique(alpha_in$BAY_ID), settings$bays)) stop("Bays in e02 output differ from settings$bays", call. = FALSE)
if (!setequal(unique(ts_all$GEAR), settings$gears)) stop("Gears differ from settings$gears", call. = FALSE)
if (anyNA(ts_all$value)) stop("NA in log_cpue (e02 output should be complete)", call. = FALSE)
filled_dates <- unique(ts_all$date[ts_all$filled])
if (!setequal(filled_dates, as.Date(settings$gap_months))) stop("Filled months in e02 output differ from settings$gap_months", call. = FALSE)

keys <- ts_all %>% distinct(scale, BAY_ID, GEAR, SPECIES_CODE, Common_Name) %>% arrange(scale, BAY_ID, GEAR, SPECIES_CODE)
cat("Series:", nrow(keys), "(bay", sum(keys$scale == "bay"), ", coastwide", sum(keys$scale == "coastwide"), ")",
    "| rows:", nrow(ts_all), "| filled rows:", sum(ts_all$filled), "\n\n")

# -----------------------------------------------------------------------------
# 6) One harmonic fit per series, on sampled months
# -----------------------------------------------------------------------------
fit_rows <- vector("list", nrow(keys))
map_rows <- vector("list", nrow(keys))
res_rows <- vector("list", nrow(keys))
out_rows <- vector("list", nrow(keys))
coef_rows  <- vector("list", nrow(keys))
check_rows <- vector("list", nrow(keys))

for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  g <- ts_all %>%
    filter(scale == k$scale, GEAR == k$GEAR, SPECIES_CODE == k$SPECIES_CODE,
           (is.na(BAY_ID) & is.na(k$BAY_ID)) | (BAY_ID %in% k$BAY_ID)) %>%
    arrange(YEAR, MONTH)
  tag <- sprintf("%s %s, %s, %s", k$scale, ifelse(is.na(k$BAY_ID), "", k$BAY_ID), k$GEAR, k$Common_Name)
  if (any(diff(month_index(g$YEAR, g$MONTH)) != 1L)) stop(tag, ": not a complete monthly series", call. = FALSE)

  smp <- g %>% filter(!filled)
  per_month <- tabulate(smp$MONTH, nbins = 12)
  if (nrow(smp) < MIN_OBS) stop(tag, ": only ", nrow(smp), " sampled months (minimum ", MIN_OBS, ").", call. = FALSE)
  if (any(per_month < MIN_PER_MONTH)) {
    stop(tag, ": calendar month(s) ", paste(which(per_month < MIN_PER_MONTH), collapse = ", "),
         " sampled fewer than ", MIN_PER_MONTH, " times.", call. = FALSE)
  }

  fit <- fit_harmonic(smp$MONTH, smp$value, N_HARMONICS)
  co <- stats::coef(fit)
  if (any(!is.finite(co))) stop(tag, ": non-finite harmonic coefficients.", call. = FALSE)
  curve <- curve_from_fit(fit, N_HARMONICS)

  g <- g %>% mutate(seasonal_component = curve[MONTH],
                    residual_component = value - seasonal_component,
                    value_deseas = residual_component)
  out_rows[[i]] <- g

  fit_rows[[i]] <- tibble::tibble(
    scale = k$scale, BAY_ID = k$BAY_ID, GEAR = k$GEAR, SPECIES_CODE = k$SPECIES_CODE, Common_Name = k$Common_Name,
    n_sampled = nrow(smp), min_sampled_per_month = min(per_month),
    intercept = unname(co["(Intercept)"]),
    amp_annual = amp(co, 1), amp_semiannual = amp(co, 2), amp_4month = amp(co, 3),
    peak_to_trough = max(curve) - min(curve),
    peak_month = which.max(curve), trough_month = which.min(curve),
    r2_seasonal = summary(fit)$r.squared
  )
  map_rows[[i]] <- tibble::tibble(scale = k$scale, BAY_ID = k$BAY_ID, GEAR = k$GEAR,
                                  SPECIES_CODE = k$SPECIES_CODE, MONTH = 1:12, seasonal = curve)

  # Diagnostics on calendar months; filled months as NA (method note 5).
  own <- ifelse(g$filled, NA_real_, g$value_deseas)
  fit2 <- fit_harmonic(smp$MONTH, smp$value, SENS_HARMONICS)
  res2 <- ifelse(g$filled, NA_real_, g$value - curve_from_fit(fit2, SENS_HARMONICS)[g$MONTH])
  a_own <- acf_grid(own)
  n_ok <- sum(is.finite(own))
  res_rows[[i]] <- tibble::tibble(
    scale = k$scale, BAY_ID = k$BAY_ID, GEAR = k$GEAR, SPECIES_CODE = k$SPECIES_CODE, Common_Name = k$Common_Name,
    n_sampled = n_ok,
    month_r2 = month_r2(own, g$MONTH),
    chance_level = (11 - 2 * N_HARMONICS) / (n_ok - 1),
    month_r2_2harm = month_r2(res2, g$MONTH),
    chance_level_2harm = (11 - 2 * SENS_HARMONICS) / (n_ok - 1),
    mean_anom = mean(own, na.rm = TRUE), sd_anom = stats::sd(own, na.rm = TRUE),
    acf1 = a_own[["acf1"]], acf12 = a_own[["acf12"]]
  )

  # Curves at every possible forecast date (method note 6). Checks: last date
  # against the full-record curve, largest change between consecutive dates
  # (both / SD of the full-record anomalies), training anomalies average zero.
  co <- seasonal_coefs_by_origin(g$date, g$MONTH, g$value, sampled = !g$filled,
                                 min_train = settings$min_train_months, K = N_HARMONICS)
  co$allowed_origin <- is_allowed_origin(co$origin_date, settings)
  B <- as.matrix(co[, seasonal_coef_names(N_HARMONICS)])
  curves <- B %*% t(as.matrix(harmonic_X(1:12, N_HARMONICS)))
  sd_anom <- stats::sd(own, na.rm = TRUE)
  tr_mean <- vapply(c(1L, nrow(co)), function(j) {
    rows <- seq_len(co$origin_i[j]); rows <- rows[!g$filled[rows]]
    mean(g$value[rows] - seasonal_level(B[j, ], g$MONTH[rows]))
  }, numeric(1))
  key <- tibble::tibble(scale = k$scale, BAY_ID = k$BAY_ID, GEAR = k$GEAR, SPECIES_CODE = k$SPECIES_CODE)
  coef_rows[[i]] <- bind_cols(key, co)
  check_rows[[i]] <- bind_cols(key, tibble::tibble(
    Common_Name = k$Common_Name, n_origins = nrow(co),
    first_origin = min(co$origin_date), last_origin = max(co$origin_date),
    last_vs_full_sd = max(abs(curves[nrow(curves), ] - curve)) / sd_anom,
    max_step_sd = max(abs(diff(curves))) / sd_anom,
    median_step_sd = stats::median(apply(abs(diff(curves)), 1, max)) / sd_anom,
    train_mean_abs_max = max(abs(tr_mean))))
}

fits      <- bind_rows(fit_rows)
seas_map  <- bind_rows(map_rows)
resid_tab <- bind_rows(res_rows)
deseas    <- bind_rows(out_rows)
coefs_by_origin <- bind_rows(coef_rows) %>%
  select(scale, BAY_ID, GEAR, SPECIES_CODE, origin_i, origin_date, allowed_origin, n_train,
         all_of(seasonal_coef_names(N_HARMONICS)))
origin_checks <- bind_rows(check_rows)

out_cols <- c("GEAR", "SPECIES_CODE", "Common_Name", "YEAR", "MONTH", "date", "filled", "n_samples",
              "mean_cpue", "value", "seasonal_component", "residual_component", "value_deseas")
alpha_deseas <- deseas %>% filter(scale == "bay") %>%
  select(BAY_ID, all_of(out_cols)) %>% arrange(BAY_ID, GEAR, SPECIES_CODE, date)
gamma_deseas <- deseas %>% filter(scale == "coastwide") %>%
  select(all_of(out_cols[1:8]), n_bays, all_of(out_cols[9:13])) %>% arrange(GEAR, SPECIES_CODE, date)

# -----------------------------------------------------------------------------
# 7) Output checks
# -----------------------------------------------------------------------------
if (nrow(alpha_deseas) != nrow(alpha_in) || nrow(gamma_deseas) != nrow(gamma_in)) stop("Row count changed.", call. = FALSE)
for (d in list(alpha_deseas, gamma_deseas)) {
  if (anyNA(d)) stop("NA in the output (the e02 series are complete).", call. = FALSE)
}
if (nrow(fits) != nrow(keys)) stop("Not every series was fitted.", call. = FALSE)
if (nrow(origin_checks) != nrow(keys)) stop("Per-origin curves missing for some series.", call. = FALSE)
if (any(origin_checks$train_mean_abs_max > 1e-8)) stop("Training anomalies do not average zero at some origin.", call. = FALSE)
if (anyNA(select(coefs_by_origin, -BAY_ID))) stop("NA in the per-origin coefficients.", call. = FALSE)

summary_tab <- fits %>%
  left_join(resid_tab %>% select(scale, BAY_ID, GEAR, SPECIES_CODE, month_r2, month_r2_2harm, chance_level, acf1),
            by = c("scale", "BAY_ID", "GEAR", "SPECIES_CODE")) %>%
  group_by(scale, GEAR, SPECIES_CODE, Common_Name) %>%
  summarise(n_series = n(),
            amp_annual = mean(amp_annual), peak_to_trough = mean(peak_to_trough),
            r2_seasonal = mean(r2_seasonal),
            leftover_mean = mean(month_r2), leftover_max = max(month_r2),
            leftover_2harm_mean = mean(month_r2_2harm), chance_level = mean(chance_level),
            acf1_mean = mean(acf1), .groups = "drop") %>%
  arrange(scale, GEAR, desc(peak_to_trough))

# -----------------------------------------------------------------------------
# 8) Write
# -----------------------------------------------------------------------------
saveRDS(alpha_deseas, out_alpha)
saveRDS(gamma_deseas, out_gamma)
saveRDS(seas_map, out_map)
saveRDS(coefs_by_origin, out_coefs)
readr::write_csv(origin_checks, out_origin, na = "")
readr::write_csv(summary_tab, out_summ, na = "")
readr::write_csv(fits, out_fits, na = "")
readr::write_csv(resid_tab, out_resid, na = "")

cat("Leftover seasonality in the anomalies (share of variance explained by calendar\n")
cat("month, sampled months; chance about (11 - 2K) / (n - 1)):\n")
print(as.data.frame(resid_tab %>% group_by(scale) %>%
  summarise(series = n(), harm3_mean = mean(month_r2), harm3_max = max(month_r2), chance3 = mean(chance_level),
            harm2_mean = mean(month_r2_2harm), harm2_max = max(month_r2_2harm), chance2 = mean(chance_level_2harm),
            acf1_median = median(acf1), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)

cat("\nCurves at every possible forecast origin:", nrow(coefs_by_origin), "rows,", nrow(origin_checks),
    "series; origins per series", min(origin_checks$n_origins), "to", max(origin_checks$n_origins), "\n")
cat("Last origin against the full-record curve (max |difference| / anomaly SD): median",
    round(stats::median(origin_checks$last_vs_full_sd), 4), ", max", round(max(origin_checks$last_vs_full_sd), 4), "\n")
cat("Largest change between consecutive origins (/ anomaly SD): median of series maxima",
    round(stats::median(origin_checks$max_step_sd), 4), ", max", round(max(origin_checks$max_step_sd), 4), "\n")

cat("\nSaved:\n")
for (f in c(out_alpha, out_gamma, out_map, out_coefs, out_summ, out_fits, out_resid, out_origin)) cat("  ", f, "\n", sep = "")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("e03_deseasonalize_abundance.R completed\n")
