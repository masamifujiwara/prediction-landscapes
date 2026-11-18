# =============================================================================
# File:    lib/models_covariates.R
# Author:  Masami Fujiwara
#
# Purpose:
#   The covariate models (ARIMAX and state-space with covariates), written
#   once and used by the covariate model scripts. Kept in a file of its own
#   so that lib/models.R, and the checkpoints of the scripts that use it,
#   stay unchanged.
#
# Used by: c09, c10 (diversity); e07, e08 (abundance)
#
# Functions:
#   covariate_columns()      column names of the b08 covariates for a set of lags
#   arimax()                 ARIMAX for one forecast date
#   ssm_xreg()               state-space local level plus covariates (note 6)
#   pad_leads()              pads a result made for fewer leads with NA rows
#   keep_leads()             drops rows beyond each model's longest lead
#   forecast_series_xreg()   loop over the forecast dates of one series
#   Needs lib/models.R (model_result(), fc_frame(), na_forecast(),
#   MODEL_PI_LEVELS, SSM_* settings, forecast_series()) and lib/helpers.R
#   (with_warnings(), check_arima_fit(), check_kfas_fit()); KFAS must be
#   attached for ssm_xreg().
#
# Method notes:
#   1. Covariates: monthly mean water temperature and salinity of the same
#      estuary and sampling method (coastwide: the sampling method's
#      coastwide means, b08), on their own scale (not deseasonalized, not
#      standardized), as lagged values only, so every value a forecast uses
#      is known at its forecast date.
#   2. Two covariate sets (settings$covariate_sets): lags 1 and 12 for lead 1;
#      lag 12 for leads 2 to 12. Each set is one fit per forecast date. The
#      lag-12 fit also gives a lead-1 forecast; it is stored too (no extra
#      cost), and the scoring scripts use the lag 1 and 12 fit at lead 1.
#   3. Leads 1 to 12 only (settings$horizons_covariate); no lag-24 set.
#   4. ARIMAX: forecast::auto.arima(xreg = ..., seasonal = TRUE, stepwise = TRUE,
#      approximation = FALSE) on the anomalies, so it is the SARIMA automatic
#      model plus the covariates (regression with SARIMA errors). Whether it
#      chose a trend term is recorded.
#   5. A missing covariate value in a training or forecast row stops the
#      script (b08 leaves none after the training start); nothing is filled
#      here.
#   6. State-space with covariates (c10, e08): the local level of
#      lib/models.R (no seasonal part, as in every state-space model here)
#      plus the covariates as regression terms with fixed coefficients
#      (diffuse prior, estimated by the Kalman filter). Covariates are
#      standardized with the mean and SD of the training rows (an exact
#      reparameterization that only helps the optimizer). Variances by fitSSM
#      (BFGS) with the restarts and collapse checks of ssm_trend() in
#      lib/models.R.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

# Column names in the b08 covariate files for a set of lags.
covariate_columns <- function(lags, vars = c("TEMP__mean", "SAL__mean")) {
  as.vector(t(outer(vars, lags, function(v, l) paste0(v, "__lag", l))))
}

# ARIMAX for one origin. y: training anomalies; X_train: covariates of the same
# rows; X_future: covariates of the H target rows. Returns a named list of one
# model result.
arimax <- function(y, H, X_train, X_future, name, freq = 12L) {
  X_train <- as.matrix(X_train); X_future <- as.matrix(X_future)[seq_len(H), , drop = FALSE]
  if (nrow(X_train) != length(y)) stop("arimax(): X_train has ", nrow(X_train), " rows for ", length(y), " values", call. = FALSE)
  if (anyNA(X_train) || anyNA(X_future)) stop("arimax(): missing covariate value (", name, ")", call. = FALSE)
  t0 <- proc.time()[[3]]
  r <- with_warnings(tryCatch(
    forecast::auto.arima(stats::ts(as.numeric(y), frequency = freq), xreg = X_train,
                         seasonal = TRUE, stepwise = TRUE, approximation = FALSE),
    error = function(e) e))
  s <- proc.time()[[3]] - t0
  chk <- check_arima_fit(r$value)
  if (!chk$ok) return(setNames(list(model_result(na_forecast(H), FALSE, chk$reason, chk$order, NA,
                                                 length(r$warnings), s)), name))
  fit <- r$value
  fc <- tryCatch(forecast::forecast(fit, xreg = X_future, level = MODEL_PI_LEVELS), error = function(e) e)
  if (inherits(fc, "error")) return(setNames(list(model_result(na_forecast(H), FALSE,
                                     paste("forecast_error:", conditionMessage(fc)), chk$order, NA, length(r$warnings), s)), name))
  setNames(list(model_result(fc_frame(fc, H), TRUE, NA, chk$order, "drift" %in% names(stats::coef(fit)),
                             length(r$warnings), s)), name)
}

# State-space local level plus covariates for one origin (decision 6).
ssm_xreg <- function(y, H, X_train, X_future, name, freq = 12L) {
  if (!"package:KFAS" %in% search()) stop("ssm_xreg(): attach KFAS first (library(KFAS))", call. = FALSE)
  X_train <- as.matrix(X_train); X_future <- as.matrix(X_future)[seq_len(H), , drop = FALSE]
  y <- as.numeric(y)
  if (nrow(X_train) != length(y)) stop("ssm_xreg(): X_train has ", nrow(X_train), " rows for ", length(y), " values", call. = FALSE)
  if (anyNA(X_train) || anyNA(X_future)) stop("ssm_xreg(): missing covariate value (", name, ")", call. = FALSE)
  mu <- colMeans(X_train); sdv <- apply(X_train, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv <= 0] <- 1
  Z_train  <- sweep(sweep(X_train, 2, mu, "-"), 2, sdv, "/")
  Z_future <- sweep(sweep(X_future, 2, mu, "-"), 2, sdv, "/")
  cn <- paste0("x", seq_len(ncol(Z_train)))
  build <- function(yv, Z, Qm, Hm) {
    d <- as.data.frame(Z); names(d) <- cn
    d$yy <- stats::ts(yv, frequency = freq)
    # KFAS recognises SSMtrend by name inside the formula (package attached).
    form <- stats::as.formula(paste("yy ~ SSMtrend(1, Q = list(Qm)) +", paste(cn, collapse = " + ")), env = environment())
    KFAS::SSModel(form, data = d, H = Hm)
  }
  v0 <- stats::var(y); if (!is.finite(v0) || v0 <= 0) v0 <- 1
  t0 <- proc.time()[[3]]
  best <- NULL; reason <- "fit_error"; n_tries <- 0L
  for (off in SSM_RESTART_OFFSETS) {
    n_tries <- n_tries + 1L
    fit <- tryCatch(KFAS::fitSSM(build(y, Z_train, matrix(NA), matrix(NA)), inits = log(v0) + off[1:2],
                                 method = "BFGS", control = list(maxit = SSM_MAXIT)), error = function(e) e)
    ck <- check_kfas_fit(fit, v0, SSM_VAR_COLLAPSE)
    if (!ck$ok) { reason <- ck$reason; next }
    new <- tryCatch(build(rep(NA_real_, H), Z_future, matrix(fit$model$Q[1, 1, 1]), fit$model$H), error = function(e) e)
    if (inherits(new, "error")) { reason <- "predict_error"; next }
    pr95 <- tryCatch(stats::predict(fit$model, newdata = new, interval = "prediction", level = 0.95), error = function(e) e)
    pr80 <- tryCatch(stats::predict(fit$model, newdata = new, interval = "prediction", level = 0.80), error = function(e) e)
    if (inherits(pr95, "error") || inherits(pr80, "error")) { reason <- "predict_error"; next }
    sdf <- (as.numeric(pr95[, "upr"]) - as.numeric(pr95[, "lwr"])) / (2 * stats::qnorm(0.975))
    if (length(sdf) != H || any(!is.finite(sdf)) || min(sdf) < SSM_PI_COLLAPSE * sqrt(v0)) { reason <- "forecast_sd_collapsed"; next }
    best <- list(pr95 = pr95, pr80 = pr80)
    break
  }
  s <- proc.time()[[3]] - t0
  if (is.null(best)) return(setNames(list(model_result(na_forecast(H), FALSE, reason, NA, FALSE, 0L, s)), name))
  fc <- data.frame(lead = seq_len(H), mean = as.numeric(best$pr95[, "fit"]),
                   lo80 = as.numeric(best$pr80[, "lwr"]), hi80 = as.numeric(best$pr80[, "upr"]),
                   lo95 = as.numeric(best$pr95[, "lwr"]), hi95 = as.numeric(best$pr95[, "upr"]))
  setNames(list(model_result(fc, TRUE, NA, paste0("tries=", n_tries), FALSE, 0L, s)), name)
}

# Pads a model result made for fewer leads (the lead-1 covariate set) to H
# rows of NA, because forecast_series() expects H rows from every model. The
# calling script removes the padded rows (keep_leads()).
pad_leads <- function(res, H) {
  lapply(res, function(m) { k <- nrow(m$fc); if (k < H) m$fc <- rbind(m$fc, na_forecast(H)[(k + 1L):H, ]); m })
}

# Drops forecast rows beyond each model's own longest lead (named vector).
keep_leads <- function(fc, max_lead_by_model) {
  lim <- max_lead_by_model[fc$model]
  fc[is.na(lim) | fc$lead <= lim, , drop = FALSE]
}

# Loop over the origins of one series for models with covariates.
# series: one complete monthly series in date order with the covariate
#   columns joined by date. covar_cols: every covariate column the models use.
# model_fun_x(y, H, X_train, X_future) returns a named list of model results.
# Each origin is passed to forecast_series() (lib/models.R) on its own, so the
# anomalies, training rows and output layout are exactly those of the models
# without covariates.
forecast_series_xreg <- function(series, origins, coefs, model_fun_x, covar_cols, H_max) {
  series <- series[order(series$date), ]
  miss <- setdiff(covar_cols, names(series))
  if (length(miss)) stop("forecast_series_xreg(): covariate columns missing: ", paste(miss, collapse = ", "), call. = FALSE)
  X <- as.matrix(series[, covar_cols, drop = FALSE])
  parts <- lapply(seq_len(nrow(origins)), function(j) {
    o <- origins$origin_i[j]
    rows <- origins$train_start_i[j]:o
    # H is the lead count forecast_series() passes for this origin.
    f <- function(y, H) model_fun_x(y, H, X[rows, , drop = FALSE], X[o + seq_len(H), , drop = FALSE])
    forecast_series(series, origins[j, , drop = FALSE], coefs, f, H_max = H_max)
  })
  list(forecasts = do.call(rbind, lapply(parts, `[[`, "forecasts")),
       info      = do.call(rbind, lapply(parts, `[[`, "info")))
}
