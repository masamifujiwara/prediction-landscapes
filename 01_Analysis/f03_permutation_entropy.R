# =============================================================================
# Script:  f03_permutation_entropy.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Model-free measure of predictability for every target series: 48
#   diversity series (q0, q1, q2; seven estuaries and coastwide; bag seine and
#   bay trawl) and 128 abundance series (8 taxa). Weighted permutation entropy
#   (WPE) says how much of a series' ups and downs follow recognizable
#   ordering patterns: 0 fully patterned, 1 no pattern (like random noise).
#   Two versions: a step of one month (short-term ups and downs) and a step of
#   twelve months (the same calendar month across years, which shows slow
#   change). Compared with shuffled copies and with simulated simple
#   carry-over (AR(1)), and related to the forecast results of f01 and the
#   slow-pattern diagnostics of f02. Reported in the Online Resources (final
#   figure by g17); in the main text only as support for the reading that
#   predictability comes from the series' own past, not from the covariates.
#   Never presented as evidence of chaos or density dependence.
#
# Run after: b09, c02, e03, f01, f02, lib/helpers.R, lib/forecast_rules.R
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   - 01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds (e03)
#   - 01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds (e03)
#   - 01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds, f1_horizons.rds (f01)
#   - 01_Analysis/Data/f2_slow_pattern_diagnostics/f2_series_diagnostics.rds (f02)
#
# Outputs:
#   01_Analysis/Data/f3_permutation_entropy/:
#   - f3_entropy_by_series.rds / .csv   one row per series, step and m
#   - f3_series_checks.rds / .csv       one row per series: zero catches,
#                                       ties, sensitivity versions
#   - f3_correlations.rds / .csv        rank correlations by class and scale
#   02_Results/03_Diagnostics/f3_permutation_entropy/:
#   - f3_summary_by_group.csv           medians and counts by class, scale,
#                                       sampling method, step and m
#   - f3_checks_by_group.csv            zero catches, ties and sensitivity
#                                       versions by class, scale and
#                                       sampling method
#   - f3_entropy_vs_skill_draft.png     draft figure (final figure by g17)
#   01_Analysis/Logs/f3_permutation_entropy/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   1. Computed on the departures from the seasonal curve fitted to the whole
#      record (value_deseas of c02 and e03), as in f02; without removing the
#      seasonal cycle the entropy would mostly measure seasonality. April and
#      May 2020 (filled, never sampled) are set to missing; any window that
#      contains a missing month is left out (no interpolation, following
#      Pennekamp et al. 2019, Box 3).
#   2. Weighted permutation entropy (Fadlallah et al. 2013, as used by
#      Pennekamp et al. 2019), normalized by log(m!), pattern lengths m = 3,
#      4, 5; the unweighted version alongside. Headline m = 4 (24 possible
#      patterns; m = 5 has 120 and is borderline).
#      Window weight: the variance of the m values in the window.
#   3. Two steps between the values of a pattern: 1 month (primary) and 12
#      months (the same calendar month across years; it sees slow change in
#      level, which the one-month step cannot).
#   4. Comparison with noise, per series, step and m: N_SURR shuffled copies
#      (same values, same missing months) and N_SURR simulated AR(1) series
#      with the series' own lag-1 correlation (simple carry-over; same missing
#      months). Reported: median and 2.5th percentile of each, and whether the
#      observed entropy lies below the 2.5th percentile.
#   5. Zero catches (abundance). In the one-month version, consecutive
#      zero-catch months (for example a taxon absent in winter) are equal in
#      the data, but their departures differ only by the seasonal curve, so
#      the curve, not the taxon, decides their order ("hidden ties").
#      Reported per series: share of zero-catch months and share of one-month
#      windows with two or more zero months; series above ZERO_FLAG are
#      flagged and correlations are given with and without them; sensitivity:
#      entropy with those windows left out, only where at least MIN_WINDOWS
#      windows remain. True ties (equal departures) can only occur in the
#      same calendar month of different years, so they never fall in a
#      one-month window (checked in the run).
#   6. True ties in the twelve-month version (same calendar month, same catch
#      rate) are broken at random: a jitter far smaller than the smallest
#      difference between distinct values. Reported: share of twelve-month
#      windows with a tie, and the spread of the entropy over N_TIEBREAK
#      different random tie-breaks; series above TIE12_FLAG are flagged.
#   7. Sensitivity of the one-month version (m = 4), observed series only:
#      (a) departures from the twelve monthly means instead of the
#      three-harmonic curve (removes every fixed seasonal feature);
#      (b) departures divided by each calendar month's own spread (so erratic
#      months do not dominate).
#   8. Rank (Spearman) correlations of 1 - WPE (m = 4) by class and scale:
#      one-month version with SARIMA automatic skill against the cumulative
#      average at leads 1 and 12 and against the 12-month average at lead 1,
#      the horizon against the cumulative average (range rule), RMSE at lead
#      1, lag-1 correlation, slow share; twelve-month version with skill at
#      leads 12 and 1, slow share, the 2020 to 2024 step and the lag-12
#      correlation.
#
# References: Bandt and Pompe (2002) Phys. Rev. Lett. 88: 174102; Fadlallah
# et al. (2013) Phys. Rev. E 87: 022911; Pennekamp et al. (2019) Ecol.
# Monogr. 89: e01359.
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

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
M_ORDERS    <- 3:5     # pattern lengths
M_HEADLINE  <- 4L      # pattern length for the correlations, checks and figure
STEPS       <- c(1L, 12L)  # months between the values of a pattern
N_SURR      <- 1000L   # shuffled and AR(1) copies per series
ZERO_FLAG   <- 0.25    # share of one-month windows with 2+ zero months above which a series is flagged
MIN_WINDOWS <- 100L    # fewest windows for the entropy without zero windows
TIE12_FLAG  <- 0.50    # share of twelve-month windows with a tie above which a series is flagged
N_TIEBREAK  <- 20L     # random tie-breaks for the spread of the twelve-month entropy
MAIN_MODEL  <- "sarima_auto"

script_stem <- "f3_permutation_entropy"
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_div_a <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_div_g <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_abd_a <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
in_abd_g <- file.path(PATHS$data_dir, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
in_f1_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
in_f1_hor  <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_horizons.rds")
in_f2      <- file.path(PATHS$data_dir, "f2_slow_pattern_diagnostics", "f2_series_diagnostics.rds")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr", "ggplot2")
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
stop_if_missing(c(in_settings, in_div_a, in_div_g, in_abd_a, in_abd_g, in_f1_lead, in_f1_hor, in_f2),
                c("b09_forecast_settings.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 2),
                  rep("e03_deseasonalize_abundance.R", 2), rep("f01_forecast_scoring.R", 2), "f02_slow_pattern_diagnostics.R"))
settings <- load_settings(in_settings)

KEYS <- c("class", "scale", "BAY_ID", "GEAR", "target")
series <- bind_rows(
  readRDS(in_div_a) %>% transmute(class = "diversity", scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, target = paste0("q", q),
                                  date, MONTH, filled, value, value_deseas, zero = FALSE),
  readRDS(in_div_g) %>% transmute(class = "diversity", scale = "coastwide", BAY_ID = NA_integer_, GEAR, target = paste0("q", q),
                                  date, MONTH, filled, value, value_deseas, zero = FALSE),
  readRDS(in_abd_a) %>% transmute(class = "abundance", scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, target = Common_Name,
                                  date, MONTH, filled, value, value_deseas, zero = !filled & mean_cpue == 0),
  readRDS(in_abd_g) %>% transmute(class = "abundance", scale = "coastwide", BAY_ID = NA_integer_, GEAR, target = Common_Name,
                                  date, MONTH, filled, value, value_deseas, zero = !filled & mean_cpue == 0)) %>%
  mutate(sampled = !filled & !(format(date, "%Y-%m") %in% format(settings$unscored_target_months, "%Y-%m")),
         dep = ifelse(sampled, value_deseas, NA_real_), val = ifelse(sampled, value, NA_real_), zero = zero & sampled)
series$sid <- do.call(paste, c(as.list(series[, KEYS]), sep = "|"))
ids <- sort(unique(series$sid))
if (length(ids) != 176L) stop(length(ids), " series; expected 176 (48 diversity, 128 abundance)", call. = FALSE)
if (anyNA(series$zero)) stop("Zero-catch flag missing for some rows", call. = FALSE)
cat("Series:", length(ids), "| rows:", nrow(series), "\n")

# -----------------------------------------------------------------------------
# 3) Permutation entropy
# -----------------------------------------------------------------------------
# Row indices of the windows of m values with step d: one list element per
# position in the window.
window_index <- function(n, m, d) {
  nw <- n - (m - 1L) * d
  lapply(seq_len(m), function(j) (1L + (j - 1L) * d):(nw + (j - 1L) * d))
}
# Weighted and unweighted permutation entropy, normalized by log(m!), for
# every column of X (one series per column; all columns share the same
# missing months). Windows that contain a missing month, or that `drop`
# marks, are left out. The ordinal pattern of a window is coded by the
# results of all pairwise comparisons, one code per ordering of m distinct
# values, so values must have no ties (see break_ties()). Window weight:
# variance of its values.
perm_entropy_mat <- function(X, m, d = 1L, drop = NULL) {
  X <- as.matrix(X)
  idx <- window_index(nrow(X), m, d)
  keep <- Reduce(`&`, lapply(idx, function(r) !is.na(X[r, 1L])))
  if (!is.null(drop)) keep <- keep & !drop
  idx <- lapply(idx, function(r) r[keep]); nw <- sum(keep)
  if (nw < 2L) return(list(wpe = rep(NA_real_, ncol(X)), pe = rep(NA_real_, ncol(X)), n_windows = nw))
  W <- lapply(idx, function(r) X[r, , drop = FALSE])
  code <- matrix(0, nw, ncol(X)); b <- 0L
  for (i in 1:(m - 1L)) for (j in (i + 1L):m) { code <- code + (W[[i]] < W[[j]]) * 2^b; b <- b + 1L }
  mu <- Reduce(`+`, W) / m
  w <- Reduce(`+`, lapply(W, function(v) (v - mu)^2)) / m
  sw <- colSums(w); hw <- numeric(ncol(X)); hu <- numeric(ncol(X))
  for (cc in unique(as.vector(code))) {
    k <- code == cc
    pw <- colSums(w * k) / sw; pu <- colSums(k) / nw
    hw <- hw - ifelse(pw > 0, pw * log(pw), 0); hu <- hu - ifelse(pu > 0, pu * log(pu), 0)
  }
  list(wpe = hw / log(factorial(m)), pe = hu / log(factorial(m)), n_windows = nw)
}
perm_entropy <- function(x, m, d = 1L, drop = NULL) unlist(perm_entropy_mat(matrix(x, ncol = 1L), m, d, drop))
# Random tie-breaking: a jitter far below the smallest gap between distinct
# values, so only the order of tied values is decided at random.
break_ties <- function(x) {
  v <- sort(unique(x[is.finite(x)]))
  eps <- if (length(v) > 1L) min(diff(v)) * 1e-3 else 1e-9
  x + stats::runif(length(x), -eps, eps)
}
# Share of complete windows (m values, step d) in which a condition holds:
# `count` = number of flagged values in the window (zero months), or `tie` =
# two values in the window are exactly equal.
window_share <- function(x, flag, m, d, what = c("count2", "tie")) {
  what <- match.arg(what)
  idx <- window_index(length(x), m, d)
  ok <- Reduce(`&`, lapply(idx, function(r) !is.na(x[r])))
  hit <- if (what == "count2") Reduce(`+`, lapply(idx, function(r) flag[r])) >= 2L else {
    X <- sapply(idx, function(r) x[r]); apply(X, 1, function(v) anyDuplicated(v) > 0L)
  }
  list(share = if (any(ok)) mean(hit[ok]) else NA_real_, hit = hit & ok, n = sum(ok))
}

# Checks of the entropy code on known cases.
chk_up <- perm_entropy(as.numeric(1:200), 4L)
set.seed(1L); chk_noise <- perm_entropy(stats::rnorm(1e5), 4L)
if (abs(chk_up[["pe"]]) > 1e-12 || abs(chk_up[["wpe"]]) > 1e-12) stop("Entropy check failed: a monotone series must give 0", call. = FALSE)
if (abs(chk_noise[["pe"]] - 1) > 0.005) stop("Entropy check failed: white noise must give about 1 (", round(chk_noise[["pe"]], 4), ")", call. = FALSE)
x_ex <- c(4, 7, 9, 10, 6, 11, 3)    # Bandt and Pompe (2002) example, m = 3
if (abs(perm_entropy(x_ex, 3L)[["pe"]] * log(6) - (-(2 / 5) * log(2 / 5) * 2 - (1 / 5) * log(1 / 5))) > 1e-12) {
  stop("Entropy check failed: Bandt and Pompe example", call. = FALSE)
}
x12 <- rep(x_ex, each = 12)[1:72]   # the same example read with a twelve-month step
if (abs(perm_entropy(x12, 3L, 12L)[["pe"]] - perm_entropy(x_ex[1:6], 3L)[["pe"]]) > 1e-12) stop("Entropy check failed: twelve-month step", call. = FALSE)
cat("Entropy code checks passed (monotone 0, white noise", round(chk_noise[["pe"]], 4), ", Bandt and Pompe example, twelve-month step)\n")

# -----------------------------------------------------------------------------
# 4) Per series: observed entropy, surrogates, zero catches, ties, sensitivity
# -----------------------------------------------------------------------------
set.seed(settings$seed)
res <- vector("list", length(ids)); chk <- vector("list", length(ids))
for (i in seq_along(ids)) {
  d <- series[series$sid == ids[i], ]; d <- d[order(d$date), ]
  x <- d$dep; miss <- is.na(x); nn <- length(x); obs <- x[!miss]
  phi <- stats::acf(x, lag.max = 1L, plot = FALSE, na.action = stats::na.pass)$acf[2]
  # zero catches and ties
  z1 <- window_share(x, d$zero, M_HEADLINE, 1L, "count2")
  t1 <- window_share(x, NULL, M_HEADLINE, 1L, "tie")
  t12 <- window_share(x, NULL, M_HEADLINE, 12L, "tie")
  if (t1$share > 0) stop(ids[i], ": a one-month window contains tied values", call. = FALSE)
  xb <- break_ties(x)
  tb12 <- replicate(N_TIEBREAK, perm_entropy(break_ties(x), M_HEADLINE, 12L)[["wpe"]])
  # surrogates (generated once, used for every step and m)
  shuf <- replicate(N_SURR, { s <- x; s[!miss] <- sample(obs); break_ties(s) })
  ar1  <- replicate(N_SURR, { s <- as.numeric(stats::arima.sim(list(ar = phi), n = nn)); s[miss] <- NA; s })
  out <- list()
  for (st in STEPS) for (m in M_ORDERS) {
    o  <- perm_entropy(xb, m, st)
    sw <- perm_entropy_mat(shuf, m, st)$wpe
    aw <- perm_entropy_mat(ar1, m, st)$wpe
    out[[length(out) + 1L]] <- data.frame(d[1, KEYS], step = st, m = m, n_months = nn, n_sampled = sum(!miss),
      n_windows = o[["n_windows"]], acf1 = phi, wpe = o[["wpe"]], pe = o[["pe"]], predictability = 1 - o[["wpe"]],
      wpe_shuffle_median = stats::median(sw), wpe_shuffle_q025 = unname(stats::quantile(sw, 0.025)),
      wpe_ar1_median = stats::median(aw), wpe_ar1_q025 = unname(stats::quantile(aw, 0.025)),
      below_shuffle = o[["wpe"]] < stats::quantile(sw, 0.025), below_ar1 = o[["wpe"]] < stats::quantile(aw, 0.025),
      beyond_ar1 = stats::median(aw) - o[["wpe"]], row.names = NULL, stringsAsFactors = FALSE)
  }
  res[[i]] <- do.call(rbind, out)
  # sensitivity versions (one-month step, headline m)
  mm  <- ave(d$val, d$MONTH, FUN = function(v) mean(v, na.rm = TRUE)); x_mm <- d$val - mm
  msd <- ave(x, d$MONTH, FUN = function(v) stats::sd(v, na.rm = TRUE))
  x_sd <- ifelse(is.finite(msd) & msd > 0, x / msd, ifelse(is.na(x), NA_real_, 0))
  nz <- perm_entropy(xb, M_HEADLINE, 1L, drop = z1$hit)
  chk[[i]] <- data.frame(d[1, KEYS], zero_month_share = mean(d$zero[!miss]),
    zero_window_share = z1$share, zero_flag = z1$share > ZERO_FLAG,
    wpe_without_zero_windows = if (nz[["n_windows"]] >= MIN_WINDOWS) nz[["wpe"]] else NA_real_,
    n_windows_without_zero = nz[["n_windows"]],
    tie12_window_share = t12$share, tie12_flag = t12$share > TIE12_FLAG,
    wpe12_tiebreak_sd = stats::sd(tb12), wpe12_tiebreak_range = diff(range(tb12)),
    wpe_monthly_means = perm_entropy(break_ties(x_mm), M_HEADLINE, 1L)[["wpe"]],
    wpe_month_spread = perm_entropy(break_ties(x_sd), M_HEADLINE, 1L)[["wpe"]],
    row.names = NULL, stringsAsFactors = FALSE)
  if (i %% 20 == 0) cat("  series", i, "of", length(ids), "\n")
}
ent <- do.call(rbind, res); rownames(ent) <- NULL
chk_tab <- do.call(rbind, chk); rownames(chk_tab) <- NULL

# Checks: every series, step and m once; entropies within 0 and 1.
if (nrow(ent) != 176L * length(STEPS) * length(M_ORDERS) || anyDuplicated(ent[, c(KEYS, "step", "m")])) {
  stop("Entropy table does not have one row per series, step and m", call. = FALSE)
}
if (any(!is.finite(ent$wpe)) || any(ent$wpe < 0 | ent$wpe > 1)) stop("Entropy outside 0 to 1 or missing", call. = FALSE)
if (nrow(chk_tab) != 176L || any(!is.finite(chk_tab$wpe_monthly_means)) || any(!is.finite(chk_tab$wpe_month_spread))) {
  stop("Sensitivity versions missing for some series", call. = FALSE)
}
cat("No one-month window contains tied values (checked in every series)\n")

# -----------------------------------------------------------------------------
# 5) Link to the forecast results (f01) and the data diagnostics (f02)
# -----------------------------------------------------------------------------
bl <- readRDS(in_f1_lead) %>% filter(model == MAIN_MODEL)
hz <- readRDS(in_f1_hor) %>% filter(model == MAIN_MODEL)
f2 <- readRDS(in_f2)
pick <- function(df, h, col, new) df %>% filter(lead == h) %>% select(all_of(KEYS), all_of(col)) %>% rename(!!new := all_of(col))
res_tab <- pick(bl, 1L, "skill_past", "skill_past_1") %>%
  left_join(pick(bl, 12L, "skill_past", "skill_past_12"), by = KEYS) %>%
  left_join(pick(bl, 1L, "skill_m12", "skill_m12_1"), by = KEYS) %>%
  left_join(pick(bl, 1L, "rmse_scaled", "rmse_1"), by = KEYS) %>%
  left_join(hz %>% select(all_of(KEYS), rel_horizon_past), by = KEYS) %>%
  left_join(f2 %>% select(all_of(KEYS), slow_share, step_2020_sd, acf12), by = KEYS)
link <- ent %>% filter(m == M_HEADLINE) %>% select(all_of(KEYS), step, predictability, acf1) %>%
  left_join(res_tab, by = KEYS) %>% left_join(chk_tab %>% select(all_of(KEYS), zero_flag, tie12_flag), by = KEYS)
need <- c("skill_past_1", "skill_past_12", "skill_m12_1", "rmse_1", "rel_horizon_past", "slow_share", "step_2020_sd", "acf12")
if (nrow(link) != 176L * length(STEPS) || anyNA(link[, need])) stop("Not every series matched its f01 and f02 results", call. = FALSE)
results_by_step <- list(`1` = c("skill_past_1", "skill_past_12", "skill_m12_1", "rel_horizon_past", "rmse_1", "acf1", "slow_share"),
                        `12` = c("skill_past_12", "skill_past_1", "slow_share", "step_2020_sd", "acf12"))
rho <- function(a, b) if (length(a) >= 5L) suppressWarnings(stats::cor(a, b, method = "spearman")) else NA_real_
cors <- do.call(rbind, lapply(split(link, list(link$class, link$scale, link$step), drop = TRUE), function(g) {
  st <- g$step[1]; fl <- if (st == 1L) g$zero_flag else g$tie12_flag
  do.call(rbind, lapply(results_by_step[[as.character(st)]], function(r) {
    ct <- suppressWarnings(stats::cor.test(g$predictability, g[[r]], method = "spearman", exact = FALSE))
    data.frame(class = g$class[1], scale = g$scale[1], step = st, m = M_HEADLINE, result = r, n = nrow(g),
               spearman = unname(ct$estimate), p = ct$p.value, n_unflagged = sum(!fl),
               spearman_unflagged = rho(g$predictability[!fl], g[[r]][!fl]), row.names = NULL)
  }))
}))
rownames(cors) <- NULL

# -----------------------------------------------------------------------------
# 6) Summary, draft figure, write
# -----------------------------------------------------------------------------
q3 <- function(x) sprintf("%.3f [%.3f, %.3f]", stats::median(x, na.rm = TRUE), stats::quantile(x, 0.1, na.rm = TRUE), stats::quantile(x, 0.9, na.rm = TRUE))
summ <- ent %>% group_by(class, scale, GEAR, step, m) %>% summarise(
  series = n(), wpe = q3(wpe), pe = q3(pe), wpe_shuffle_median = q3(wpe_shuffle_median), wpe_ar1_median = q3(wpe_ar1_median),
  below_shuffle = sum(below_shuffle), below_ar1 = sum(below_ar1), .groups = "drop")
summ_chk <- chk_tab %>% left_join(ent %>% filter(step == 1L, m == M_HEADLINE) %>% select(all_of(KEYS), wpe), by = KEYS) %>%
  group_by(class, scale, GEAR) %>% summarise(series = n(), zero_month_share = q3(zero_month_share), zero_window_share = q3(zero_window_share),
  zero_flagged = sum(zero_flag), wpe = q3(wpe), wpe_without_zero_windows = q3(wpe_without_zero_windows),
  tie12_window_share = q3(tie12_window_share), tie12_flagged = sum(tie12_flag), wpe12_tiebreak_sd = q3(wpe12_tiebreak_sd),
  wpe_monthly_means = q3(wpe_monthly_means), wpe_month_spread = q3(wpe_month_spread), .groups = "drop")

fig_dat <- bind_rows(
  link %>% filter(step == 1L) %>% transmute(class, scale, GEAR, x = predictability, y = skill_past_1, flagged = zero_flag,
                                            version = "One-month step against skill 1 month ahead"),
  link %>% filter(step == 12L) %>% transmute(class, scale, GEAR, x = predictability, y = skill_past_12, flagged = tie12_flag,
                                             version = "Twelve-month step against skill 12 months ahead"))
p <- ggplot2::ggplot(fig_dat, ggplot2::aes(x, y, colour = interaction(scale, GEAR, sep = ", "), shape = flagged)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey50") + ggplot2::geom_point(alpha = 0.8) +
  ggplot2::facet_grid(version ~ class, labeller = ggplot2::label_wrap_gen(28)) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("no", "yes")) +
  ggplot2::labs(x = paste0("Intrinsic predictability (1 - weighted permutation entropy, m = ", M_HEADLINE, ")"),
                y = "Skill against the cumulative average", colour = "Scale, gear", shape = "Flagged",
                title = "Draft: permutation entropy against forecast skill (SARIMA automatic)") +
  ggplot2::theme_bw(base_size = 10)

out_rds <- file.path(dirs$data, "f3_entropy_by_series.rds"); out_csv <- file.path(dirs$data, "f3_entropy_by_series.csv")
out_chk <- file.path(dirs$data, "f3_series_checks.rds");   out_chk_csv <- file.path(dirs$data, "f3_series_checks.csv")
out_cor <- file.path(dirs$data, "f3_correlations.rds");    out_cor_csv <- file.path(dirs$data, "f3_correlations.csv")
out_sum <- file.path(dirs$diag, "f3_summary_by_group.csv"); out_sum2 <- file.path(dirs$diag, "f3_checks_by_group.csv")
out_png <- file.path(dirs$diag, "f3_entropy_vs_skill_draft.png")
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE); if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv|png|pdf)$", full.names = TRUE); if (length(old_d)) archive_files(old_d)
saveRDS(ent, out_rds); readr::write_csv(ent, out_csv, na = "")
saveRDS(chk_tab, out_chk); readr::write_csv(chk_tab, out_chk_csv, na = "")
saveRDS(cors, out_cor); readr::write_csv(cors, out_cor_csv, na = "")
readr::write_csv(summ, out_sum, na = ""); readr::write_csv(summ_chk, out_sum2, na = "")
ggplot2::ggsave(out_png, p, width = 170, height = 150, units = "mm", dpi = 200)

cat("\nEntropy by class, scale, gear, step and m (median [10th, 90th percentile]; counts of series):\n")
print(as.data.frame(summ), row.names = FALSE)
cat("\nZero catches, ties and sensitivity versions (one-month step, m = ", M_HEADLINE, "):\n", sep = "")
print(as.data.frame(summ_chk), row.names = FALSE)
cat("\nRank correlations of 1 - WPE (m = ", M_HEADLINE, ") with forecast results and data diagnostics:\n", sep = "")
print(cors %>% mutate(across(c(spearman, p, spearman_unflagged), ~ round(.x, 3))), row.names = FALSE)
for (f in c(out_rds, out_csv, out_chk, out_chk_csv, out_cor, out_cor_csv, out_sum, out_sum2, out_png)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
end_log("finished")
cat("f03_permutation_entropy.R completed\n")
