# =============================================================================
# File:    lib/join_forecasts.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Brings the forecasts of every model script and the reference forecasts of
#   one target class (diversity: c07 to c10 and c12; abundance: e05a to e08)
#   into one table for scoring in f01. One copy only.
#
# Used by: d01 and e09 (join_forecasts()), f01 (model_info())
#
# Functions:
#   model_info()               role and comparison partner of every model name
#   combine_covariate_sets()   adds the combined covariate models (note 2)
#   rows_of_series()           rows of a table that belong to one series
#   join_forecasts()           joins and checks all sources, one series at a time
#
# Method notes:
#   1. One row per series, forecast date (origin), lead time and model, on
#      the anomaly scale of that forecast date, for every model and reference
#      as stored by the model scripts. Nothing is scored here; f01 computes
#      every measure.
#   2. The covariate models are also stored as one combined model each:
#      "arimax" = arimax_lag1_12 at lead 1 and arimax_lag12 at leads 2 to 12;
#      "ssm_xreg" likewise. The two covariate-set versions stay in the table
#      under their own names.
#   3. Roles (model_info()): the four headline forecasting models (role
#      "main": sarima_auto, ssm_level, arimax, ssm_xreg); trend test
#      (sarima_notrend, sarima_lintrend, damped_trend against sarima_auto;
#      ssm_lintrend against ssm_level); covariate test (arimax against
#      sarima_auto, ssm_xreg against ssm_level); references (ref_naive, last
#      month's value; ref_past_average, the cumulative average, primary;
#      ref_12month, the 12-month average).
#   4. Failed fits stay as rows with ok = FALSE and NA forecasts (their reason
#      is in the fit-information files). How a failed fit enters a paired
#      comparison is decided in f01.
#   5. Added columns: scored (FALSE when the target month is April or May
#      2020, is_scored_target()), season_at_target (the seasonal curve fitted
#      at the forecast date, evaluated at the target month) and value_obs
#      (y_obs plus season_at_target, the observed value on the target's own
#      scale), so f01 can also score on the full scale (seasonal share).
#   6. Checks that stop the script: every source has every forecast date of
#      the origin table and the expected leads; the observed anomaly of a
#      series, forecast date and lead is the same in every source; value_obs
#      equals the observed value in the input series (proves the seasonal
#      curve is rebuilt correctly); no duplicate rows; no NA forecast from a
#      fit marked ok.
#   7. Works one series at a time (rows selected by their key values), so
#      memory stays small.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

# Roles and partners of every model name that can appear.
model_info <- function() {
  data.frame(
    model   = c("sarima_auto", "ssm_level", "arimax", "ssm_xreg",
                "sarima_notrend", "sarima_lintrend", "damped_trend", "ssm_lintrend",
                "arimax_lag1_12", "arimax_lag12", "ssm_xreg_lag1_12", "ssm_xreg_lag12",
                "ref_naive", "ref_past_average", "ref_12month"),
    role    = c("main", "main", "main", "main",
                "trend_test", "trend_test", "trend_test", "trend_test",
                "covariate_set", "covariate_set", "covariate_set", "covariate_set",
                "reference", "reference", "reference"),
    partner = c(NA, NA, "sarima_auto", "ssm_level",
                "sarima_auto", "sarima_auto", "sarima_auto", "ssm_level",
                "sarima_auto", "sarima_auto", "ssm_level", "ssm_level",
                NA, NA, NA),
    stringsAsFactors = FALSE)
}

# Adds the combined covariate models (decision 2) to a forecast table and to
# its fit-information table.
combine_covariate_sets <- function(fc, info, prefix) {
  a <- fc[fc$model == paste0(prefix, "_lag1_12") & fc$lead == 1L, , drop = FALSE]
  b <- fc[fc$model == paste0(prefix, "_lag12") & fc$lead >= 2L, , drop = FALSE]
  if (!nrow(a) || !nrow(b)) stop("combine_covariate_sets(): no ", prefix, " rows", call. = FALSE)
  comb <- rbind(a, b); comb$model <- prefix
  ia <- info[info$model == paste0(prefix, "_lag1_12"), , drop = FALSE]; ia$part <- "lead_1"
  ib <- info[info$model == paste0(prefix, "_lag12"), , drop = FALSE];   ib$part <- "leads_2_12"
  ic <- rbind(ia, ib); ic$model <- prefix
  if (!"part" %in% names(info)) info$part <- NA_character_
  list(fc = rbind(fc, comb), info = rbind(info, ic))
}

# Rows of a data frame belonging to one series (key values in a one-row frame;
# NA matches NA, for BAY_ID of the coastwide series).
rows_of_series <- function(df, keys, key_row) {
  sel <- rep(TRUE, nrow(df))
  for (k in keys) {
    v <- key_row[[k]][1]
    sel <- sel & (if (is.na(v)) is.na(df[[k]]) else (!is.na(df[[k]]) & df[[k]] == v))
  }
  which(sel)
}

# keys: series key columns (diversity: scale, BAY_ID, GEAR, q; abundance:
# scale, BAY_ID, GEAR, SPECIES_CODE). sources: named list of list(forecasts,
# info) read from the model scripts. origins: the c03/e04 origin table.
# coefs: the c02/e03 per-origin coefficient table. series: the input series
# (keys, date, value) used for the value_obs check. Returns list(forecasts,
# info) for all series.
join_forecasts <- function(sources, keys, origins, coefs, series, settings, tol = 1e-8) {
  for (s in names(sources)) {
    sources[[s]]$forecasts <- as.data.frame(sources[[s]]$forecasts); sources[[s]]$forecasts$source <- s
    sources[[s]]$info <- as.data.frame(sources[[s]]$info); sources[[s]]$info$source <- s
  }
  cat("Forecast rows read:", sum(sapply(sources, function(x) nrow(x$forecasts))),
      "| fit records:", sum(sapply(sources, function(x) nrow(x$info))), "\n")
  # Longest lead of each model over the whole run.
  H_by_model <- unlist(lapply(sources, function(x) tapply(x$forecasts$lead, x$forecasts$model, max)))
  names(H_by_model) <- sub("^[^.]+[.]", "", names(H_by_model))
  print(H_by_model)
  mi <- model_info()
  unknown <- setdiff(names(H_by_model), mi$model)
  if (length(unknown)) stop("Unknown model names: ", paste(unknown, collapse = ", "), call. = FALSE)
  cn <- grep("^(b0|sin[0-9]+|cos[0-9]+)$", names(coefs), value = TRUE)
  K <- (length(cn) - 1L) %/% 2L

  skeys <- unique(origins[, keys, drop = FALSE])
  skeys <- skeys[do.call(order, skeys), , drop = FALSE]
  out_fc <- vector("list", nrow(skeys)); out_info <- vector("list", nrow(skeys))
  max_dy <- 0; max_dv <- 0
  for (r in seq_len(nrow(skeys))) {
    kr <- skeys[r, , drop = FALSE]
    lab <- paste(unlist(kr), collapse = "|")
    og <- origins[rows_of_series(origins, keys, kr), , drop = FALSE]
    fc <- do.call(rbind, lapply(sources, function(x) x$forecasts[rows_of_series(x$forecasts, keys, kr), , drop = FALSE]))
    info <- do.call(rbind, lapply(sources, function(x) x$info[rows_of_series(x$info, keys, kr), c(keys, "origin_i", "model", "ok", "reason", "source")]))
    rownames(fc) <- NULL; rownames(info) <- NULL

    # Completeness: every model has every origin and leads 1..min(H, max_lead).
    for (m in names(H_by_model)) {
      exp_rows <- sum(pmin(H_by_model[[m]], og$max_lead))
      got <- sum(fc$model == m); n_o <- sum(info$model == m)
      if (got != exp_rows || n_o != nrow(og) || !setequal(info$origin_i[info$model == m], og$origin_i)) {
        stop("Series ", lab, ", model ", m, ": ", got, " forecast rows (expected ", exp_rows, "), ", n_o,
             " fit records (expected ", nrow(og), ")", call. = FALSE)
      }
    }
    # Observed anomaly identical across sources (decision 6).
    yk <- fc$origin_i * 1000L + fc$lead
    dy <- abs(fc$y_obs - fc$y_obs[match(yk, yk)])
    if (any(dy > tol)) stop("Series ", lab, ": observed anomalies differ between sources", call. = FALSE)
    max_dy <- max(max_dy, dy)
    # ok flag from the fit records.
    ik <- paste(info$origin_i, info$model)
    if (anyDuplicated(ik)) stop("Series ", lab, ": duplicate fit records", call. = FALSE)
    fc$ok <- info$ok[match(paste(fc$origin_i, fc$model), ik)]
    if (anyNA(fc$ok)) stop("Series ", lab, ": forecast rows without a fit record", call. = FALSE)
    if (any(fc$ok & is.na(fc$mean))) stop("Series ", lab, ": NA forecasts from fits marked ok", call. = FALSE)
    # Combined covariate models (decision 2).
    for (p in c("arimax", "ssm_xreg")) {
      if (any(fc$model == paste0(p, "_lag12"))) { cc <- combine_covariate_sets(fc, info, p); fc <- cc$fc; info <- cc$info }
    }
    if (!"part" %in% names(info)) info$part <- NA_character_
    if (anyDuplicated(fc[, c("origin_i", "lead", "model")])) stop("Series ", lab, ": duplicate forecast rows", call. = FALSE)
    # Seasonal curve at the target month from the origin's fit (decision 5).
    co <- coefs[rows_of_series(coefs, keys, kr), , drop = FALSE]
    ci <- match(fc$origin_i, co$origin_i)
    if (anyNA(ci)) stop("Series ", lab, ": forecast rows without a seasonal curve for their origin", call. = FALSE)
    mo <- as.integer(format(fc$target_date, "%m"))
    sv <- co$b0[ci]
    for (k in seq_len(K)) sv <- sv + co[[paste0("sin", k)]][ci] * sin(2 * pi * k * mo / 12) +
                                     co[[paste0("cos", k)]][ci] * cos(2 * pi * k * mo / 12)
    fc$season_at_target <- sv
    fc$value_obs <- fc$y_obs + sv
    se <- series[rows_of_series(series, keys, kr), , drop = FALSE]
    v <- se$value[match(fc$target_date, se$date)]
    if (anyNA(v)) stop("Series ", lab, ": target months missing from the input series", call. = FALSE)
    dv <- max(abs(v - fc$value_obs))
    if (dv > tol) stop("Series ", lab, ": value_obs differs from the input series by up to ", signif(dv, 3), call. = FALSE)
    max_dv <- max(max_dv, dv)
    out_fc[[r]] <- fc; out_info[[r]] <- info
  }
  cat("Every model has every origin and its full set of leads in all", nrow(skeys), "series\n")
  cat("Observed anomalies identical across all sources (max difference ", signif(max_dy, 3), ")\n", sep = "")
  cat("Seasonal curve rebuilt: observed values equal the input series (max difference ", signif(max_dv, 3), ")\n", sep = "")
  fc <- do.call(rbind, out_fc); rm(out_fc)
  info <- do.call(rbind, out_info)
  fc$scored <- is_scored_target(fc$target_date, settings)
  fc$role <- mi$role[match(fc$model, mi$model)]
  out_cols <- c(keys, intersect("Common_Name", names(fc)), "origin_i", "origin_date", "in_common", "lead", "target_date",
                "target_filled", "scored", "model", "role", "source", "ok", "mean", "lo80", "hi80", "lo95", "hi95",
                "y_obs", "season_at_target", "value_obs")
  fc <- fc[, out_cols]
  rownames(fc) <- NULL; rownames(info) <- NULL
  list(forecasts = fc, info = info)
}
