# =============================================================================
# File:    lib/scoring.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Scoring of the out-of-sample forecasts of ONE series (all models,
#   forecast dates and lead times): accuracy, skill, uncertainty ranges,
#   horizons, management situations, variance partition, accuracy by period
#   and the second accuracy measure. One copy only.
#
# Used by: f01
#
# Functions:
#   block_counts()    moving-block bootstrap counts of each forecast date
#   q_lo(), q_hi()    2.5% and 97.5% quantiles by column
#   count_horizon()   leads from 1 up to the first failing lead
#   score_series()    all scores of one series; returns a list of tables
#
# Method notes:
#   1. Prediction error = forecast minus observed value at the target month
#      (anomaly scale of the forecast date; identical on the full scale).
#      Absolute accuracy: RMSE; for diversity divided by the series' mean
#      level over all sampled months, for abundance in log(CPUE + 1) units
#      ("scaled" RMSE below means either). Supporting: bias, median and
#      quantiles of the size of the monthly miss, share of monthly misses
#      within the tolerance (hit rate), coverage and width of the 80% and 95%
#      intervals.
#   2. Skill = 1 - SSE(model) / SSE(reference) over the same cases, against
#      the cumulative average (primary), the 12-month average and last
#      month's value, and against each model's partner (trend and covariate
#      tests).
#   3. A failed fit removes that forecast date and lead from every comparison
#      that involves the model, for both members of the pair; targets in
#      April and May 2020 are never scored.
#   4. Uncertainty: moving-block bootstrap over forecast dates (blocks of
#      BLOCK_LEN consecutive forecast dates, B resamples, seed fixed per
#      series), the same resamples for every model and lead of a series; 95%
#      ranges.
#   5. Horizons count leads from 1 to the first failing lead (0 if lead 1
#      fails). Absolute: scaled RMSE within the tolerance; by the
#      uncertainty-range rule (upper end of the range within the tolerance,
#      primary) and by value. Relative: skill above zero against the
#      cumulative average and against the 12-month average; range rule
#      (lower end above zero, primary) and value. Return of skill: first lead
#      after the relative horizon at which the lower end is above zero again.
#   6. Management situations at leads 1, 12, 24: accurate (absolute rule) x
#      skill above zero against the cumulative average: forecast-led (both),
#      mean-led (accurate only), forecast with buffer (skill only),
#      precaution (neither); by range (primary) and by value.
#   7. Variance partition: seasonal share at the target months (1 - variance
#      of the observed departures / variance of the observed values, common
#      forecast dates, lead 1), and of the rest: slow = skill against the
#      cumulative average averaged over the long leads (12 to the longest
#      lead), at least 0; short-lived = skill at lead 1 minus slow, at least
#      0; unpredictable = the remainder.
#   8. Periods of the target month: 1992-2007, 2008-2019 and 2020-2024 (all
#      forecast dates).
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

# Moving-block bootstrap: a B x n matrix of how often each origin (in time
# order) is drawn.
block_counts <- function(n, B, L, seed) {
  set.seed(seed)
  L <- min(L, n); k <- ceiling(n / L)
  C <- matrix(0, B, n)
  for (b in seq_len(B)) {
    st <- sample.int(n - L + 1L, k, replace = TRUE)
    idx <- as.vector(outer(0:(L - 1L), st, "+"))[seq_len(n)]
    C[b, ] <- tabulate(idx, nbins = n)
  }
  C
}

q_lo <- function(M) apply(M, 2, stats::quantile, 0.025, na.rm = TRUE, names = FALSE)
q_hi <- function(M) apply(M, 2, stats::quantile, 0.975, na.rm = TRUE, names = FALSE)

# Leads from 1 up to the first lead where ok is FALSE (NA counts as failing).
count_horizon <- function(ok) { ok[is.na(ok)] <- FALSE; f <- which(!ok); if (length(f)) f[1] - 1L else length(ok) }

# d: the joined rows of one series (one class). level: the series' mean level
# (1 for abundance). tol: tolerance on the scaled RMSE. tol_grid: tolerances
# for the horizon-against-tolerance curves. Returns a list of tables.
score_series <- function(d, level, tol, tol_grid, B, L, seed, partner_of, long_from = 12L) {
  origins <- sort(unique(d$origin_i)); n_o <- length(origins); H <- max(d$lead)
  models <- sort(unique(d$model))
  oi <- match(d$origin_i, origins)
  mk <- function(col, m) { M <- matrix(NA_real_, n_o, H); s <- d$model == m; M[cbind(oi[s], d$lead[s])] <- d[[col]][s]; M }
  # Case information shared by all models (from any model with all leads).
  ref_m <- if ("ref_past_average" %in% models) "ref_past_average" else models[1]
  s0 <- d$model == ref_m
  Y <- mk("y_obs", ref_m); V <- mk("value_obs", ref_m)
  SC <- matrix(FALSE, n_o, H); SC[cbind(oi[s0], d$lead[s0])] <- d$scored[s0]
  TY <- matrix(NA_integer_, n_o, H); TY[cbind(oi[s0], d$lead[s0])] <- as.integer(format(d$target_date[s0], "%Y"))
  TD <- matrix(NA_real_, n_o, H); TD[cbind(oi[s0], d$lead[s0])] <- as.numeric(d$target_date[s0])
  common <- tapply(d$in_common[s0], oi[s0], function(x) x[1])[as.character(seq_len(n_o))]
  common <- as.logical(common)

  E <- list(); OK <- list(); LO80 <- list(); HI80 <- list(); LO95 <- list(); HI95 <- list()
  for (m in models) {
    E[[m]] <- mk("mean", m) - Y
    okm <- matrix(FALSE, n_o, H); s <- d$model == m; okm[cbind(oi[s], d$lead[s])] <- d$ok[s] & !is.na(d$mean[s])
    OK[[m]] <- okm & SC
    LO80[[m]] <- mk("lo80", m); HI80[[m]] <- mk("hi80", m); LO95[[m]] <- mk("lo95", m); HI95[[m]] <- mk("hi95", m)
  }
  Hm <- vapply(models, function(m) max(d$lead[d$model == m]), integer(1))

  # Bootstrap counts: common origins (main) and all origins (periods).
  ic <- which(common)
  Cc <- block_counts(length(ic), B, L, seed)
  Ca <- block_counts(n_o, B, L, seed + 1L)

  # Sums of squared errors over valid cases: point value and bootstrap.
  sse <- function(Em, Vm, rows, C) {
    E2 <- Em[rows, , drop = FALSE]^2; E2[!Vm[rows, , drop = FALSE]] <- 0; E2[is.na(E2)] <- 0
    N <- Vm[rows, , drop = FALSE] * 1
    list(p = colSums(E2), n = colSums(N), b = C %*% E2, bn = C %*% N)
  }
  by_lead <- list(); budget <- list()
  refs <- c(past = "ref_past_average", m12 = "ref_12month", naive = "ref_naive")
  for (m in models) {
    h <- seq_len(Hm[[m]])
    Vm <- OK[[m]]
    a <- sse(E[[m]], Vm, ic, Cc)
    rmse <- sqrt(a$p / a$n) / level
    rb <- sqrt(a$b / a$bn) / level
    em <- E[[m]][ic, , drop = FALSE]; vm <- Vm[ic, , drop = FALSE]
    am <- abs(em) / level; am[!vm] <- NA
    lv <- function(f) apply(am, 2, f)
    out <- data.frame(model = m, lead = h, n = a$n[h], rmse_scaled = rmse[h], rmse_lo = q_lo(rb)[h], rmse_hi = q_hi(rb)[h],
      bias_scaled = (colSums(ifelse(vm, em, 0)) / a$n / level)[h],
      miss_q10 = lv(function(x) stats::quantile(x, 0.10, na.rm = TRUE))[h], miss_q25 = lv(function(x) stats::quantile(x, 0.25, na.rm = TRUE))[h],
      miss_median = lv(function(x) stats::median(x, na.rm = TRUE))[h], miss_q75 = lv(function(x) stats::quantile(x, 0.75, na.rm = TRUE))[h],
      miss_q90 = lv(function(x) stats::quantile(x, 0.90, na.rm = TRUE))[h], hit_rate = lv(function(x) mean(x <= tol, na.rm = TRUE))[h],
      sd_obs_scaled = (apply(ifelse(vm, Y[ic, , drop = FALSE], NA), 2, stats::sd, na.rm = TRUE) / level)[h])
    if (!grepl("^ref_", m) || !all(is.na(LO80[[m]]))) {
      cov <- function(lo, hi) { inside <- Y[ic, , drop = FALSE] >= lo[ic, , drop = FALSE] & Y[ic, , drop = FALSE] <= hi[ic, , drop = FALSE]
        inside[!vm] <- NA; colMeans(inside, na.rm = TRUE) }
      wid <- function(lo, hi) { w <- (hi - lo)[ic, , drop = FALSE]; w[!vm] <- NA; colMeans(w, na.rm = TRUE) / level }
      out$coverage80 <- cov(LO80[[m]], HI80[[m]])[h]; out$coverage95 <- cov(LO95[[m]], HI95[[m]])[h]
      out$width80_scaled <- wid(LO80[[m]], HI80[[m]])[h]; out$width95_scaled <- wid(LO95[[m]], HI95[[m]])[h]
    } else { out$coverage80 <- NA_real_; out$coverage95 <- NA_real_; out$width80_scaled <- NA_real_; out$width95_scaled <- NA_real_ }
    # Skill against the references and the partner (paired cases).
    comp <- c(refs, partner = if (!is.na(partner_of[m])) unname(partner_of[m]) else NA_character_)
    for (nm in names(comp)) {
      r <- comp[[nm]]
      if (is.na(r) || !(r %in% models) || r == m) { out[[paste0("skill_", nm)]] <- NA_real_; out[[paste0("skill_", nm, "_lo")]] <- NA_real_
        out[[paste0("skill_", nm, "_hi")]] <- NA_real_; out[[paste0("n_", nm)]] <- NA_integer_; next }
      Vp <- Vm & OK[[r]]
      sm <- sse(E[[m]], Vp, ic, Cc); sr <- sse(E[[r]], Vp, ic, Cc)
      sk <- 1 - sm$p / sr$p; skb <- 1 - sm$b / sr$b
      out[[paste0("skill_", nm)]] <- sk[h]; out[[paste0("skill_", nm, "_lo")]] <- q_lo(skb)[h]; out[[paste0("skill_", nm, "_hi")]] <- q_hi(skb)[h]
      out[[paste0("n_", nm)]] <- as.integer(sm$n[h])
    }
    # Full-scale variance explained (seasonal share), lead by lead.
    vv <- V[ic, , drop = FALSE]; vv[!vm] <- NA
    var_full <- apply(vv, 2, stats::var, na.rm = TRUE)
    out$r2_full <- (1 - (a$p / a$n) / var_full)[h]
    by_lead[[m]] <- out
  }
  bl <- do.call(rbind, by_lead); rownames(bl) <- NULL
  s1_by_model <- bl$skill_past[bl$lead == 1][match(bl$model, bl$model[bl$lead == 1])]
  bl$skill_retained <- ifelse(is.finite(s1_by_model) & s1_by_model > 0, bl$skill_past / s1_by_model, NA_real_)

  # Horizons and quadrants.
  hz <- list(); quad <- list(); tolc <- list()
  for (m in models) {
    x <- bl[bl$model == m, ]
    rr <- function(v) count_horizon(v)
    rel_h <- rr(x$skill_past_lo > 0); rel_v <- rr(x$skill_past > 0)
    ret <- if (rel_h < nrow(x)) { w <- which(x$skill_past_lo > 0 & x$lead > rel_h + 1L); if (length(w)) x$lead[w[1]] else NA_integer_ } else NA_integer_
    hz[[m]] <- data.frame(model = m, max_lead = nrow(x),
      abs_horizon = rr(x$rmse_hi <= tol), abs_horizon_value = rr(x$rmse_scaled <= tol),
      rel_horizon_past = rel_h, rel_horizon_past_value = rel_v,
      rel_horizon_12m = rr(x$skill_m12_lo > 0), rel_horizon_12m_value = rr(x$skill_m12 > 0),
      return_of_skill_lead = ret)
    tolc[[m]] <- do.call(rbind, lapply(tol_grid, function(t) data.frame(model = m, tolerance = t,
      abs_horizon = rr(x$rmse_hi <= t), abs_horizon_value = rr(x$rmse_scaled <= t))))
    for (l in intersect(c(1L, 12L, 24L), x$lead)) {
      y <- x[x$lead == l, ]
      lab <- function(acc, sk) ifelse(acc & sk, "forecast-led", ifelse(acc, "mean-led", ifelse(sk, "forecast with buffer", "precaution")))
      quad[[paste(m, l)]] <- data.frame(model = m, lead = l, accurate = y$rmse_hi <= tol, skilful = isTRUE(y$skill_past_lo > 0),
        accurate_value = y$rmse_scaled <= tol, skilful_value = isTRUE(y$skill_past > 0),
        quadrant = lab(y$rmse_hi <= tol, isTRUE(y$skill_past_lo > 0)), quadrant_value = lab(y$rmse_scaled <= tol, isTRUE(y$skill_past > 0)))
    }
    # Variance budget.
    s1 <- x$skill_past[x$lead == 1]
    sl <- mean(x$skill_past[x$lead >= min(long_from, nrow(x))], na.rm = TRUE)
    yy <- Y[ic, 1]; vv <- V[ic, 1]; ok1 <- OK[[m]][ic, 1]
    seas <- 1 - stats::var(yy[ok1]) / stats::var(vv[ok1])
    slow <- max(0, sl); short <- max(0, s1 - slow); unp <- max(0, 1 - slow - short)
    budget[[m]] <- data.frame(model = m, seasonal_share = seas, skill_lead1 = s1, skill_long = sl,
      seasonal = seas, slow = (1 - seas) * slow, short_lived = (1 - seas) * short, unpredictable = (1 - seas) * unp)
  }

  # Periods of the target month (all origins), value and bootstrap range.
  per <- list()
  pcut <- function(y) cut(y, c(-Inf, 2007, 2019, Inf), labels = c("1992-2007", "2008-2019", "2020-2024"))
  allr <- seq_len(n_o)
  for (m in models) {
    h <- seq_len(Hm[[m]])
    for (p in levels(pcut(2000))) {
      inP <- matrix(pcut(TY) == p, n_o, H); inP[is.na(inP)] <- FALSE
      Vm <- OK[[m]] & inP
      a <- sse(E[[m]], Vm, allr, Ca)
      out <- data.frame(model = m, period = p, lead = h, n = a$n[h], rmse_scaled = (sqrt(a$p / a$n) / level)[h],
                        rmse_lo = q_lo(sqrt(a$b / a$bn) / level)[h], rmse_hi = q_hi(sqrt(a$b / a$bn) / level)[h])
      em <- abs(E[[m]]) / level; em[!Vm] <- NA
      out$hit_rate <- colMeans(em <= tol, na.rm = TRUE)[h]
      for (nm in c("past", "m12")) {
        r <- refs[[nm]]
        if (!(r %in% models) || r == m) { out[[paste0("skill_", nm)]] <- NA_real_; out[[paste0("skill_", nm, "_lo")]] <- NA_real_; out[[paste0("skill_", nm, "_hi")]] <- NA_real_; next }
        Vp <- Vm & OK[[r]]
        sm <- sse(E[[m]], Vp, allr, Ca); sr <- sse(E[[r]], Vp, allr, Ca)
        out[[paste0("skill_", nm)]] <- (1 - sm$p / sr$p)[h]
        skb <- 1 - sm$b / sr$b
        out[[paste0("skill_", nm, "_lo")]] <- q_lo(skb)[h]; out[[paste0("skill_", nm, "_hi")]] <- q_hi(skb)[h]
      }
      per[[paste(m, p)]] <- out
    }
  }
  pr <- do.call(rbind, per); rownames(pr) <- NULL
  ph <- do.call(rbind, lapply(split(pr, list(pr$model, pr$period), drop = TRUE), function(x) {
    x <- x[order(x$lead), ]
    data.frame(model = x$model[1], period = x$period[1], abs_horizon = count_horizon(x$rmse_hi <= tol),
               abs_horizon_value = count_horizon(x$rmse_scaled <= tol), rel_horizon_past = count_horizon(x$skill_past_lo > 0),
               rel_horizon_past_value = count_horizon(x$skill_past > 0))
  }))

  # Second accuracy measure: leads 1 and 12, all origins; RMSE / SD of the
  # observed values over the same cases, and skill against last month's value.
  pub <- do.call(rbind, lapply(models, function(m) {
    l <- intersect(c(1L, 12L), seq_len(Hm[[m]]))
    do.call(rbind, lapply(l, function(ld) {
      v <- OK[[m]][, ld]; e <- E[[m]][v, ld]; y <- Y[v, ld]
      vn <- v & OK[["ref_naive"]][, ld]
      data.frame(model = m, lead = ld, n = sum(v), rmse_over_sd = sqrt(mean(e^2)) / stats::sd(y),
                 skill_naive = 1 - sum(E[[m]][vn, ld]^2) / sum(E[["ref_naive"]][vn, ld]^2))
    }))
  }))

  # Cases removed by failed fits, per model (common origins, scored targets).
  failed <- data.frame(model = models, cases_failed = vapply(models, function(m) {
    sum(SC[ic, seq_len(Hm[[m]])] & !OK[[m]][ic, seq_len(Hm[[m]])]) }, numeric(1)))

  # Monthly misses of every model (common origins), for pooled distributions;
  # SARIMA lead-1 misses on all origins for the largest-miss list.
  miss <- do.call(rbind, lapply(models, function(m) {
    h <- seq_len(Hm[[m]]); em <- E[[m]][ic, h, drop = FALSE]; vm <- OK[[m]][ic, h, drop = FALSE]
    data.frame(model = m, lead = rep(h, each = length(ic))[as.vector(vm)], miss_scaled = (abs(em) / level)[vm])
  }))
  big <- NULL
  if ("sarima_auto" %in% models) {
    v <- OK[["sarima_auto"]][, 1]
    big <- data.frame(origin_i = origins[v], target_date = as.Date(TD[v, 1], origin = "1970-01-01"), miss_scaled = E[["sarima_auto"]][v, 1] / level,
                      miss_sd = E[["sarima_auto"]][v, 1] / sqrt(mean(E[["sarima_auto"]][v, 1]^2)))
  }
  list(by_lead = bl, horizons = do.call(rbind, hz), horizon_by_tolerance = do.call(rbind, tolc), quadrants = do.call(rbind, quad),
       budget = do.call(rbind, budget), periods = pr, period_horizons = ph, published = pub, failed = failed, miss = miss, big = big)
}
