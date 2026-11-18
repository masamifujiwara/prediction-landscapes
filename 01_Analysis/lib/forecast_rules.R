# =============================================================================
# File:    lib/forecast_rules.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Shared forecast rules: month indexing, reading the forecast settings,
#   filling the April and May 2020 gap in a target series, deciding which
#   target months are scored, and choosing the forecast dates (origins) of a
#   series. There is one copy only: never copy these functions into a script.
#
# Used by: b09, c00 to c03, c07 to c10, c12, d01, e02 to e09, f01 to f04,
#   f07, g05, g08, g14, g16, h11 to h13; also lib/seasonal.R
#
# Functions:
#   month_index()                 consecutive integers for consecutive months
#   month_date()                  first day of the month, as Date
#   load_settings()               reads the b09 settings and prints them
#   print_settings()              prints a settings list (into the run log)
#   fill_gap_seasonal_anomaly()   fills missing months of one target series
#   is_scored_target()            TRUE for target months that are scored
#   is_allowed_origin()           TRUE for months that may be forecast dates
#   drop_unscored_targets()       removes forecasts of unscored target months
#   forecast_origins()            forecast dates of one complete series
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09), read by
#   load_settings(). Needs stop_if_missing() from lib/helpers.R.
#
# Method notes:
#   The COVID rule:
#   - April and May 2020: no sampling anywhere, so they are the only months
#     missing inside any series.
#   - Covariates: filled by linear interpolation in b08 (filled = TRUE).
#   - Targets: filled once, by fill_gap_seasonal_anomaly(), and marked
#     filled = TRUE. Diversity is filled in b09 and abundance in e02. Model
#     scripts never fill by themselves (no na.interp, no dropping rows and
#     rebuilding ts()).
#   - Scoring: a forecast whose target month is in
#     settings$unscored_target_months is removed from every skill score,
#     error measure and summary (every model, every lead time, the reference
#     forecasts included). Use drop_unscored_targets(); never write the months
#     into a script.
#   - Forecast dates: no forecasts are made from months in
#     settings$skip_origin_months (is_allowed_origin()).
#   Forecast dates: forecast_origins() is the one place that decides the
#   forecast dates of a series. Every model and reference script reads the
#   origin tables written by c03 (diversity) and e04 (abundance), which call
#   it; no script builds its own.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

# Month index: consecutive integers for consecutive months.
month_index <- function(YEAR, MONTH) as.integer(YEAR) * 12L + as.integer(MONTH) - 1L

# First day of the month, as Date.
month_date <- function(YEAR, MONTH) as.Date(sprintf("%04d-%02d-01", as.integer(YEAR), as.integer(MONTH)))

# Reads settings.rds from b09. With show = TRUE, prints the settings (into the
# run log when one is open).
load_settings <- function(path = file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds"),
                          show = TRUE) {
  stop_if_missing(path, "b09_forecast_settings.R")
  s <- readRDS(path)
  if (show) {
    cat("Settings loaded from:", path, "\n")
    print_settings(s)
  }
  s
}

print_settings <- function(s) {
  fmt <- function(x) {
    if (inherits(x, "Date")) x <- format(x, "%Y-%m")
    if (is.list(x)) return(paste(names(x), vapply(x, function(v) paste(fmt(v), collapse = ", "), ""), sep = " = ", collapse = "; "))
    paste(x, collapse = ", ")
  }
  for (nm in names(s)) cat(sprintf("  %-24s %s\n", nm, fmt(s[[nm]])))
  invisible(s)
}

# Fills the rows `fill` of one monthly series by seasonal-anomaly interpolation:
#   z = transform(value); for each calendar month, the mean of z over the
#   climatology rows; the anomaly z - mean is interpolated linearly across
#   the gap; then the mean is added back and the result back-transformed.
# Climatology rows: observed (not filled, not NA) rows up to and including
# clim_end (a Date, the first observed month after the gap), so the fill uses
# no data from later than clim_end.
# YEAR and MONTH must form one complete monthly series (no missing rows).
fill_gap_seasonal_anomaly <- function(value, YEAR, MONTH, fill, clim_end,
                                      transform = log, inverse = exp) {
  t <- month_index(YEAR, MONTH)
  if (any(diff(t) != 1L)) stop("fill_gap_seasonal_anomaly(): the series is not a complete monthly sequence")
  if (any(!is.na(value[fill]))) stop("fill_gap_seasonal_anomaly(): rows to fill must be NA")
  obs <- !fill & !is.na(value)
  z <- rep(NA_real_, length(value))
  z[obs] <- transform(value[obs])
  if (any(!is.finite(z[obs]))) stop("fill_gap_seasonal_anomaly(): transform gives non-finite values")
  clim_rows <- obs & month_date(YEAR, MONTH) <= clim_end
  clim <- tapply(z[clim_rows], factor(MONTH[clim_rows], levels = 1:12), mean)
  if (anyNA(clim[MONTH[fill]])) stop("fill_gap_seasonal_anomaly(): no climatology for a month to fill")
  anom <- z - clim[MONTH]
  anom_f <- stats::approx(t[obs], anom[obs], xout = t[fill], rule = 1)$y
  if (anyNA(anom_f)) stop("fill_gap_seasonal_anomaly(): a gap is not bounded by observed months")
  out <- value
  out[fill] <- inverse(anom_f + clim[MONTH[fill]])
  out
}

# TRUE for target months that count in skill scores.
is_scored_target <- function(date, settings) !(as.Date(date) %in% settings$unscored_target_months)

# TRUE for months that may be used as forecast origins.
is_allowed_origin <- function(date, settings) !(as.Date(date) %in% settings$skip_origin_months)

# Removes forecasts whose target month is not scored. `date_col` names the
# column holding the target month (Date, first day of month). Prints how many
# rows were removed, so the log shows it.
drop_unscored_targets <- function(df, settings, date_col = "test_date") {
  if (!date_col %in% names(df)) stop("drop_unscored_targets(): no column ", date_col)
  keep <- is_scored_target(df[[date_col]], settings)
  cat("Forecasts removed from scoring (target month in ",
      paste(format(settings$unscored_target_months, "%Y-%m"), collapse = ", "), "): ",
      sum(!keep), " of ", length(keep), "\n", sep = "")
  df[keep, , drop = FALSE]
}

# Forecast origins of ONE complete monthly series (rows in date order).
# date: first day of each month; filled: TRUE for filled months (never
# training-only gaps, see the COVID rule). Returns one row per origin:
#   origin_i        row of the origin in the series
#   origin_date     the origin month (the last month the forecast may use)
#   train_start_i   first training row (settings$training_start_index)
#   n_train_sampled sampled months from train_start_i to the origin
#   max_lead        longest lead with a target inside the record (<= max(settings$horizons))
#   in_common       TRUE for the common origin set used at every lead time
#                   (origin at least max(settings$horizons) months before the end)
# Rules (settings from b09): origins from settings$first_origin_date
# (the same for every series) to the month before the series end; at least
# settings$min_train_months sampled months from the series start; never in
# settings$skip_origin_months; every model trains from row
# settings$training_start_index of its series.
forecast_origins <- function(date, filled, settings) {
  date <- as.Date(date); n <- length(date)
  t <- as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m"))
  if (n < 2L || any(diff(t) != 1L)) stop("forecast_origins(): not a complete monthly series", call. = FALSE)
  need <- c("first_origin_date", "min_train_months", "training_start_index", "horizons", "skip_origin_months", "origin_step")
  miss <- setdiff(need, names(settings))
  if (length(miss)) stop("forecast_origins(): settings lack ", paste(miss, collapse = ", "), " (rerun b09)", call. = FALSE)
  if (settings$origin_step != 1L) stop("forecast_origins(): only origin_step = 1 is implemented", call. = FALSE)
  sampled_cum <- cumsum(!filled)
  i <- seq_len(n - 1L)
  keep <- date[i] >= as.Date(settings$first_origin_date) &
          sampled_cum[i] >= settings$min_train_months &
          is_allowed_origin(date[i], settings)
  i <- i[keep]
  if (!length(i)) stop("forecast_origins(): no origin satisfies the rules", call. = FALSE)
  s0 <- as.integer(settings$training_start_index)
  if (any(i < s0)) stop("forecast_origins(): an origin precedes the training start", call. = FALSE)
  hmax <- max(settings$horizons)
  data.frame(origin_i = i, origin_date = date[i], train_start_i = s0,
             n_train_sampled = sampled_cum[i] - c(0L, sampled_cum)[s0],
             max_lead = pmin(hmax, n - i), in_common = (n - i) >= hmax)
}
