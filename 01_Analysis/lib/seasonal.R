# =============================================================================
# File:    lib/seasonal.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Shared seasonal-curve code: the harmonic regression on calendar month,
#   its coefficients at every possible forecast date (origin), and the
#   departures from the seasonal curve (anomalies) as seen from one forecast
#   date. There is one copy only: never copy these functions into a script.
#
# Used by: c01, c02, c07 to c10, c12, d01, e03, e05a to e09, f02, f04; also
#   lib/models.R (anomalies_at_origin())
#
# Functions:
#   SEASONAL_N_HARMONICS        number of harmonics (the one place it is set)
#   harmonic_X()                harmonic design matrix with intercept
#   seasonal_coef_names()       names of the coefficient columns
#   fit_harmonic()              lm() fit of the curve on sampled months
#   curve_from_fit()            fitted monthly level, January to December
#   seasonal_level()            fitted level for any months from coefficients
#   seasonal_coefs_by_origin()  coefficients at every possible forecast date
#   anomalies_at_origin()       departures from the curve fitted at one date
#
# Method notes:
#   The seasonal curve of a series is a harmonic regression on calendar
#   month, fitted with an intercept on the sampled months (filled == FALSE):
#     value = b0 + sum over k = 1..K of [ s_k sin(2 pi k m / 12) + c_k cos(2 pi k m / 12) ]
#   with m = MONTH (1 to 12). The anomaly is value minus the fitted monthly
#   level (b0 plus the harmonics), so the anomalies are centered on zero over
#   the months the curve was fitted on.
#   1. Three harmonics (periods of 12, 6 and 4 months) for every target,
#      diversity and abundance: for diversity, three forecast slightly better
#      than two and leave only chance-level seasonality.
#   2. The curve is refitted at every possible forecast date on the sampled
#      months up to and including that date (in-window fit), so a forecast
#      never uses a curve that saw the future. c02 and e03 store the
#      coefficients per series and forecast date; model scripts rebuild the
#      anomalies of a forecast date with anomalies_at_origin().
#   3. No time or trend term in the curve, and no seasonal swing that changes
#      over time (tested and dropped: it made seasonal-only forecasts worse).
#   4. Possible forecast dates: every month of a series, from the first month
#      at which at least settings$min_train_months sampled months are
#      available, to the month before the last. Which of them are used, and
#      from which month the training starts, are set by forecast_origins()
#      in lib/forecast_rules.R. April and May 2020 are never forecast dates
#      (the COVID rule, lib/forecast_rules.R).
#   5. Property used by the reference forecasts: with an intercept in the
#      curve, the training anomalies at a forecast date average exactly zero,
#      so the cumulative-average reference equals the in-window seasonal
#      climatology.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

SEASONAL_N_HARMONICS <- 3L

# Harmonic design matrix with intercept; columns (Intercept), sin1, cos1, ...
harmonic_X <- function(month, K = SEASONAL_N_HARMONICS) {
  m <- as.numeric(month)
  cols <- list(`(Intercept)` = rep(1, length(m)))
  for (k in seq_len(K)) {
    cols[[paste0("sin", k)]] <- sin(2 * pi * k * m / 12)
    cols[[paste0("cos", k)]] <- cos(2 * pi * k * m / 12)
  }
  as.data.frame(cols, check.names = FALSE)
}

# Names of the coefficient columns for K harmonics.
seasonal_coef_names <- function(K = SEASONAL_N_HARMONICS) {
  c("b0", as.vector(rbind(paste0("sin", seq_len(K)), paste0("cos", seq_len(K)))))
}

# lm() fit of the curve; month and value are the sampled months only.
fit_harmonic <- function(month, value, K = SEASONAL_N_HARMONICS) {
  ok <- is.finite(value)
  dat <- cbind(data.frame(y = value[ok]), harmonic_X(month[ok], K)[, -1, drop = FALSE])
  fml <- stats::as.formula(paste("y ~", paste(setdiff(names(dat), "y"), collapse = " + ")))
  stats::lm(fml, data = dat)
}

# Fitted monthly level (12 values, January to December) of an lm() fit.
curve_from_fit <- function(fit, K = SEASONAL_N_HARMONICS) {
  as.numeric(stats::predict(fit, newdata = harmonic_X(1:12, K)[, -1, drop = FALSE]))
}

# Fitted level for any months from a coefficient vector in seasonal_coef_names() order.
seasonal_level <- function(coefs, month) {
  coefs <- as.numeric(coefs)
  K <- (length(coefs) - 1L) %/% 2L
  as.numeric(as.matrix(harmonic_X(month, K)) %*% coefs)
}

# Coefficients of the curve at every possible origin of ONE complete monthly
# series (rows in date order). sampled: TRUE for sampled months (filled ==
# FALSE). Returns one row per origin: origin_i (row of the origin in the
# series), origin_date, n_train (sampled months used), and the coefficients.
seasonal_coefs_by_origin <- function(date, month, value, sampled, min_train,
                                     K = SEASONAL_N_HARMONICS) {
  n <- length(value)
  if (any(diff(as.numeric(format(date, "%Y")) * 12 + as.numeric(format(date, "%m"))) != 1)) {
    stop("seasonal_coefs_by_origin(): not a complete monthly series", call. = FALSE)
  }
  use <- sampled & is.finite(value)
  cum <- cumsum(use)
  first <- which(cum >= min_train)[1]
  if (is.na(first) || first >= n) stop("seasonal_coefs_by_origin(): fewer than ", min_train, " sampled months before the last month", call. = FALSE)
  X <- as.matrix(harmonic_X(month, K))
  origins <- first:(n - 1L)
  B <- matrix(NA_real_, length(origins), 2L * K + 1L, dimnames = list(NULL, seasonal_coef_names(K)))
  for (j in seq_along(origins)) {
    rows <- which(use[seq_len(origins[j])])
    per_month <- tabulate(month[rows], nbins = 12)
    if (any(per_month < 2L)) stop("seasonal_coefs_by_origin(): a calendar month has fewer than 2 sampled values before origin ", format(date[origins[j]]), call. = FALSE)
    fit <- stats::lm.fit(X[rows, , drop = FALSE], value[rows])
    if (any(!is.finite(fit$coefficients))) stop("seasonal_coefs_by_origin(): non-finite coefficients at origin ", format(date[origins[j]]), call. = FALSE)
    B[j, ] <- fit$coefficients
  }
  data.frame(origin_i = origins, origin_date = date[origins], n_train = cum[origins], B, check.names = FALSE)
}

# Anomalies of ONE series as seen from one origin: every row (training and
# test months alike) minus the curve fitted at that origin. coef_row: one row
# of the coefficient table (columns b0, sin1, cos1, ...). Returns the series
# with seasonal_component and value_deseas replaced.
anomalies_at_origin <- function(series, coef_row, value_col = "value") {
  cn <- grep("^(b0|sin[0-9]+|cos[0-9]+)$", names(coef_row), value = TRUE)
  K <- (length(cn) - 1L) %/% 2L
  coefs <- unlist(coef_row[1, seasonal_coef_names(K)], use.names = FALSE)
  series$seasonal_component <- seasonal_level(coefs, series$MONTH)
  series$value_deseas <- series[[value_col]] - series$seasonal_component
  series
}
