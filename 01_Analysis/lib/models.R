# =============================================================================
# File:    lib/models.R
# Author:  Masami Fujiwara
#
# Purpose:
#   The forecasting models and reference forecasts, written once and used by
#   every model and reference-forecast script. There is one copy only: never
#   copy these functions into a script.
#   Every model works on the anomalies of ONE series as seen from ONE
#   forecast date (origin): value minus the seasonal curve fitted on the
#   months up to that date (anomalies_at_origin(), lib/seasonal.R). It is
#   fitted on the training rows (settings$training_start_index to the
#   forecast date) and forecasts leads 1 to H.
#
# Used by: c07 to c10, c12 (diversity); e05a to e08 (abundance); also
#   lib/models_covariates.R
#
# Functions:
#   MODEL_PI_LEVELS, SSM_*   prediction-interval levels and state-space fit settings
#   na_forecast()            empty forecast frame (failed fits)
#   fc_frame()               forecast frame from a forecast::forecast() object
#   model_result()           forecasts plus a one-row info record
#   sarima_family()          SARIMA automatic, no trend and linear trend
#   damped_trend()           damped-trend exponential smoothing
#   ssm_trend()              state-space local level or local linear trend
#   reference_forecasts()    last month's value, cumulative average, 12-month average
#   forecast_series()        loop over the forecast dates of one series
#   Needs lib/helpers.R (with_warnings(), check_arima_fit(), check_kfas_fit())
#   and lib/seasonal.R; KFAS must be attached for ssm_trend().
#
# Method notes:
#   1. Eight forecasting models. Headline models: SARIMA automatic
#      (auto.arima), state-space local level, ARIMAX, state-space with
#      covariates (the last two in lib/models_covariates.R). Trend test:
#      SARIMA no trend and SARIMA linear trend (both with the automatic
#      model's orders), damped trend (exponential smoothing), state-space
#      local linear trend. The word "drift" is not used in the paper: a
#      linear trend is a constant change per month (with differencing) or a
#      time term (without).
#   2. SARIMA automatic: forecast::auto.arima(seasonal = TRUE, stepwise = TRUE,
#      approximation = FALSE) on the anomalies. Whether it chose a trend term
#      is recorded at every forecast date (has_trend).
#   3. SARIMA no trend / linear trend: the automatic model's orders refitted
#      with the trend term forbidden / always included. When the automatic
#      model already has (lacks) the trend, its forecasts are reused for
#      linear trend (no trend) rather than refitted. With two or more
#      differences no linear trend can be fitted; that is recorded.
#   4. State-space models have NO seasonal part: the seasonal cycle is
#      removed beforehand, and a seasonal part fitted noise and made forecasts
#      worse. Local level (SSMtrend degree 1) and local linear trend (degree
#      2); variances estimated by fitSSM (BFGS) with restarts from several
#      starting values; a fit is rejected if it did not converge, if all
#      variances collapsed (below 1e-6 of the training variance) or if the
#      forecast SD collapsed (below 1% of the training SD).
#   5. One fit per forecast date gives every lead time.
#   6. Reference forecasts on the anomalies: last month's value (naive: the
#      last anomaly at the forecast date), cumulative average (mean of all
#      sampled anomalies from the series start to the forecast date; zero by
#      construction, so it equals the in-window seasonal climatology),
#      12-month average (mean of the sampled anomalies of the last 12
#      months). They use all past months, not only the model training rows.
#   7. Training rows include the filled months (April and May 2020) with
#      their filled values (the COVID rule, lib/forecast_rules.R); targets in
#      those months are kept in the output with target_filled = TRUE and
#      removed at scoring (drop_unscored_targets()).
#   8. Failed fits produce NA forecasts with a named reason; nothing is
#      silently replaced.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

MODEL_PI_LEVELS     <- c(80, 95)
SSM_VAR_COLLAPSE    <- 1e-6
SSM_PI_COLLAPSE     <- 0.01
SSM_MAXIT           <- 200L
SSM_RESTART_OFFSETS <- list(c(0, 0, 0), c(-3, 0, 0), c(0, -3, -3), c(2, -2, -2), c(-6, 1, 1))

# Empty forecast frame for leads 1..H (used for failed fits).
na_forecast <- function(H) {
  data.frame(lead = seq_len(H), mean = NA_real_, lo80 = NA_real_, hi80 = NA_real_,
             lo95 = NA_real_, hi95 = NA_real_)
}

# Forecast frame from a forecast::forecast() object.
fc_frame <- function(fc, H) {
  lo <- as.matrix(fc$lower); hi <- as.matrix(fc$upper)
  data.frame(lead = seq_len(H), mean = as.numeric(fc$mean)[seq_len(H)],
             lo80 = lo[seq_len(H), 1], hi80 = hi[seq_len(H), 1],
             lo95 = lo[seq_len(H), 2], hi95 = hi[seq_len(H), 2])
}

# One model result: forecasts plus a one-row info record.
model_result <- function(fc, ok, reason = NA_character_, order = NA_character_,
                         has_trend = NA, n_warnings = 0L, seconds = NA_real_) {
  list(fc = fc, info = data.frame(ok = ok, reason = reason, order = order, has_trend = has_trend,
                                  n_warnings = as.integer(n_warnings), seconds = seconds))
}

# ---- SARIMA family -------------------------------------------------------------

# Returns a named list of model results: sarima_auto, sarima_notrend, sarima_lintrend.
sarima_family <- function(y, H, freq = 12L) {
  yt <- stats::ts(as.numeric(y), frequency = freq)
  t0 <- proc.time()[[3]]
  r <- with_warnings(tryCatch(forecast::auto.arima(yt, seasonal = TRUE, stepwise = TRUE, approximation = FALSE),
                              error = function(e) e))
  auto <- r$value; sec_auto <- proc.time()[[3]] - t0
  chk <- check_arima_fit(auto)
  if (!chk$ok) {
    bad <- model_result(na_forecast(H), FALSE, chk$reason, chk$order, NA, length(r$warnings), sec_auto)
    return(list(sarima_auto = bad, sarima_notrend = bad, sarima_lintrend = bad))
  }
  has_trend <- "drift" %in% names(stats::coef(auto))
  fc_auto <- fc_frame(forecast::forecast(auto, h = H, level = MODEL_PI_LEVELS), H)
  out <- list(sarima_auto = model_result(fc_auto, TRUE, NA, chk$order, has_trend, length(r$warnings), sec_auto))
  ord <- forecast::arimaorder(auto)
  nonseas <- ord[1:3]
  seas <- if (length(ord) >= 6) ord[4:6] else c(0, 0, 0)
  has_mean <- "intercept" %in% names(stats::coef(auto))
  refit <- function(drift) {
    t1 <- proc.time()[[3]]
    rr <- with_warnings(tryCatch(
      forecast::Arima(yt, order = nonseas, seasonal = list(order = seas, period = freq),
                      include.mean = has_mean, include.drift = drift),
      error = function(e) e))
    f <- rr$value; s <- proc.time()[[3]] - t1; ck <- check_arima_fit(f)
    if (!ck$ok) return(model_result(na_forecast(H), FALSE, ck$reason, ck$order, drift, length(rr$warnings), s))
    got <- "drift" %in% names(stats::coef(f))
    if (drift && !got) return(model_result(na_forecast(H), FALSE, "linear_trend_not_possible_two_differences",
                                            ck$order, FALSE, length(rr$warnings), s))
    model_result(fc_frame(forecast::forecast(f, h = H, level = MODEL_PI_LEVELS), H), TRUE, NA, ck$order, got,
                 length(rr$warnings), s)
  }
  if (has_trend) {
    out$sarima_lintrend <- model_result(fc_auto, TRUE, "same_as_automatic", chk$order, TRUE, 0L, 0)
    out$sarima_notrend  <- refit(FALSE)
  } else {
    out$sarima_notrend  <- model_result(fc_auto, TRUE, "same_as_automatic", chk$order, FALSE, 0L, 0)
    out$sarima_lintrend <- refit(TRUE)
  }
  out[c("sarima_auto", "sarima_notrend", "sarima_lintrend")]
}

# ---- Damped trend ------------------------------------------------------------------

damped_trend <- function(y, H, freq = 12L) {
  t0 <- proc.time()[[3]]
  r <- with_warnings(tryCatch(forecast::ets(stats::ts(as.numeric(y), frequency = freq), model = "AAN", damped = TRUE),
                              error = function(e) e))
  s <- proc.time()[[3]] - t0
  if (inherits(r$value, "error")) return(list(damped_trend = model_result(na_forecast(H), FALSE,
                                                paste("fit_error:", conditionMessage(r$value)), NA, NA, length(r$warnings), s)))
  f <- r$value
  ok <- all(is.finite(f$par)) && is.finite(f$sigma2)
  if (!ok) return(list(damped_trend = model_result(na_forecast(H), FALSE, "nonfinite_parameters", NA, TRUE, length(r$warnings), s)))
  list(damped_trend = model_result(fc_frame(forecast::forecast(f, h = H, level = MODEL_PI_LEVELS), H), TRUE, NA,
                                   paste0("phi=", signif(f$par[["phi"]], 3)), TRUE, length(r$warnings), s))
}

# ---- State-space (KFAS), no seasonal part ---------------------------------------------

# degree 1 = local level, 2 = local linear trend.
ssm_trend <- function(y, H, degree = 1L, freq = 12L) {
  if (!"package:KFAS" %in% search()) stop("ssm_trend(): attach KFAS first (library(KFAS))", call. = FALSE)
  y <- as.numeric(y); v0 <- stats::var(y)
  if (!is.finite(v0) || v0 <= 0) v0 <- 1
  t0 <- proc.time()[[3]]
  build <- function() {
    Q <- if (degree == 1L) list(matrix(NA)) else list(matrix(NA), matrix(NA))
    # KFAS recognises its model terms by name inside the formula, so SSMtrend is
    # written without the KFAS:: prefix (the package must be attached).
    yts <- stats::ts(y, frequency = freq)
    KFAS::SSModel(yts ~ SSMtrend(degree, Q = Q), H = matrix(NA))
  }
  n_par <- degree + 1L
  best <- NULL; best_ll <- -Inf; reason <- "fit_error"; n_tries <- 0L
  for (off in SSM_RESTART_OFFSETS) {
    n_tries <- n_tries + 1L
    o <- rep(0, n_par); o[seq_len(min(n_par, length(off)))] <- off[seq_len(min(n_par, length(off)))]
    fit <- tryCatch(KFAS::fitSSM(build(), inits = log(v0) + o, method = "BFGS",
                                 control = list(maxit = SSM_MAXIT)), error = function(e) e)
    ck <- check_kfas_fit(fit, v0, SSM_VAR_COLLAPSE)
    if (!ck$ok) { reason <- ck$reason; next }
    pr95 <- tryCatch(stats::predict(fit$model, n.ahead = H, interval = "prediction", level = 0.95), error = function(e) e)
    pr80 <- tryCatch(stats::predict(fit$model, n.ahead = H, interval = "prediction", level = 0.80), error = function(e) e)
    if (inherits(pr95, "error") || inherits(pr80, "error")) { reason <- "predict_error"; next }
    sdv <- (as.numeric(pr95[, "upr"]) - as.numeric(pr95[, "lwr"])) / (2 * stats::qnorm(0.975))
    if (any(!is.finite(sdv)) || min(sdv) < SSM_PI_COLLAPSE * sqrt(v0)) { reason <- "forecast_sd_collapsed"; next }
    if (ck$loglik > best_ll) { best <- list(pr95 = pr95, pr80 = pr80); best_ll <- ck$loglik }
    break
  }
  s <- proc.time()[[3]] - t0
  name <- if (degree == 1L) "ssm_level" else "ssm_lintrend"
  if (is.null(best)) return(setNames(list(model_result(na_forecast(H), FALSE, reason, NA, degree == 2L, 0L, s)), name))
  fc <- data.frame(lead = seq_len(H), mean = as.numeric(best$pr95[, "fit"]),
                   lo80 = as.numeric(best$pr80[, "lwr"]), hi80 = as.numeric(best$pr80[, "upr"]),
                   lo95 = as.numeric(best$pr95[, "lwr"]), hi95 = as.numeric(best$pr95[, "upr"]))
  setNames(list(model_result(fc, TRUE, NA, paste0("tries=", n_tries), degree == 2L, 0L, s)), name)
}

# ---- Reference forecasts ---------------------------------------------------------------

# a: anomalies of the whole series at this origin; sampled: !filled; o: origin row.
reference_forecasts <- function(a, sampled, o, H) {
  past <- which(sampled[seq_len(o)])
  last12 <- past[past > o - 12L]
  mk <- function(v) model_result(data.frame(lead = seq_len(H), mean = rep(v, H), lo80 = NA_real_, hi80 = NA_real_,
                                            lo95 = NA_real_, hi95 = NA_real_), is.finite(v))
  list(ref_naive = mk(a[o]), ref_past_average = mk(mean(a[past])), ref_12month = mk(mean(a[last12])))
}

# ---- Shared loop over the origins of one series -------------------------------------------

# series: one complete monthly series in date order (columns date, MONTH,
#   filled, value). origins: its rows of the c03/e04 origin table. coefs: its
#   rows of the c02/e03 per-origin coefficient table. model_fun(y_train, H)
#   returns a named list of model results. H_max: longest lead for these models
#   (settings$horizons or settings$horizons_covariate). With references = TRUE,
#   model_fun is ignored and the reference forecasts are made.
# Returns list(forecasts, info): forecasts one row per origin, model and lead
#   (origin_i, origin_date, in_common, model, lead, target_date, target_filled,
#   y_obs, mean, lo80, hi80, lo95, hi95); info one row per origin and model.
forecast_series <- function(series, origins, coefs, model_fun = NULL, H_max, references = FALSE) {
  series <- series[order(series$date), ]
  sampled <- !series$filled
  fc_rows <- vector("list", nrow(origins)); info_rows <- vector("list", nrow(origins))
  for (j in seq_len(nrow(origins))) {
    o <- origins$origin_i[j]; H <- min(H_max, origins$max_lead[j])
    cr <- coefs[coefs$origin_i == o, , drop = FALSE]
    if (nrow(cr) != 1L) stop("forecast_series(): no unique seasonal curve for origin row ", o, call. = FALSE)
    a <- anomalies_at_origin(series, cr)$value_deseas
    res <- if (references) reference_forecasts(a, sampled, o, H) else
      model_fun(a[origins$train_start_i[j]:o], H)
    tgt <- o + seq_len(H)
    fc_rows[[j]] <- do.call(rbind, lapply(names(res), function(m) {
      data.frame(origin_i = o, origin_date = origins$origin_date[j], in_common = origins$in_common[j], model = m,
                 res[[m]]$fc, target_date = series$date[tgt], target_filled = series$filled[tgt], y_obs = a[tgt])
    }))
    info_rows[[j]] <- do.call(rbind, lapply(names(res), function(m) {
      data.frame(origin_i = o, origin_date = origins$origin_date[j], model = m,
                 n_train = o - origins$train_start_i[j] + 1L, res[[m]]$info)
    }))
  }
  list(forecasts = do.call(rbind, fc_rows), info = do.call(rbind, info_rows))
}
