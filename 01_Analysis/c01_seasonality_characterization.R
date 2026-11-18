# =============================================================================
# Script:  c01_seasonality_characterization.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Fits one harmonic seasonal curve (periods of 12, 6 and 4 months) to each
#   coastwide diversity series (2 sampling methods x q = 0, 1, 2) over the
#   whole record and reports how strong and how stable the seasonality is
#   (summarized in Online Resource 2, Table S4, by h14). These full-record
#   curves describe the series; the curves the forecasts use are refitted at
#   every forecast date in c02. Estuary-level curves are also fitted in c02.
#
# Run after: c00, b09, lib/seasonal.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   01_Analysis/Data/b9_forecast_settings/gamma_q012_filled.rds (b09, from b07)
#
# Outputs:
#   01_Analysis/Data/c1_seasonality_characterization/
#     seasonality_summary.csv  one row per sampling method and q (read by h14)
#     seasonality_fits.rds     named list BS__q0 ... BT__q2: model, metrics,
#                              seasonal_curve (read by c02, h14)
#   02_Results/01_Figures/c1_seasonality_characterization/
#     seasonal_cycle_curves__by_q_by_gear.pdf and .tiff
#   01_Analysis/Logs/c1_seasonality_characterization/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   - Model per series: q = intercept + sin1 + cos1 + sin2 + cos2 + sin3 +
#     cos3 on calendar month (lib/seasonal.R), natural scale, one lm() fit
#     over the whole record. The number of harmonics comes from lib/seasonal.R,
#     so c01, c02 and e03 always agree (c02 checks it). The summary keeps the
#     annual amplitude as seasonal_amplitude and adds the 6- and 4-month
#     amplitudes.
#   - No year or trend term: year-to-year variation and trend stay in the
#     departures from the seasonal curve, which are what the models forecast.
#     Trend is tested out of sample by the trend versions of the models.
#   - Natural scale, not log: the seasonal swing grows with the level by less
#     than one for one, so the log scale over-corrects and leaves more
#     seasonality in the departures.
#   - The curve is fitted on observed months only (filled == FALSE).
#   - Residual autocorrelation and the Ljung-Box test are computed on the
#     complete monthly grid with the filled months set to NA, so every lag is
#     a lag in calendar months.
#   - The split-half check compares the annual amplitude in the two halves of
#     the record.
#   - Fit settings are constants below (at least 60 months, Ljung-Box lag 12,
#     ACF to lag 24). A missing input stops the script, and so does any of
#     the 6 fits failing.
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
script_stem <- "c1_seasonality_characterization"   # output folder stem (short form c1_)

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_gamma    <- file.path(PATHS$data_dir, "b9_forecast_settings", "gamma_q012_filled.rds")

N_HARMONICS         <- SEASONAL_N_HARMONICS   # from lib/seasonal.R
MIN_N_MONTHS        <- 60     # fewer observed months: the fit fails and the script stops
LJUNGBOX_LAG        <- 12
ACF_MAX_LAG         <- 24
Q_LEVELS            <- 0:2

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
dirs <- script_dirs(script_stem, c("data", "logs", "figs"))
archive_previous(c(dirs$data, dirs$figs))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

# -----------------------------------------------------------------------------
# 4) Input
# -----------------------------------------------------------------------------
stop_if_missing(c(in_settings, in_gamma), "b09_forecast_settings.R")
settings <- load_settings(in_settings)
gamma <- readRDS(in_gamma)
expect_cols(gamma, c("GEAR", "YEAR", "MONTH", "q0", "q1", "q2", "filled"), "gamma_q012_filled")
if (!setequal(unique(gamma$GEAR), settings$gears)) {
  stop("Gears in gamma_q012_filled differ from settings: ", paste(sort(unique(gamma$GEAR)), collapse = ", "), call. = FALSE)
}
if (anyDuplicated(gamma[, c("GEAR", "YEAR", "MONTH")])) stop("Duplicate GEAR x YEAR x MONTH rows in gamma_q012_filled", call. = FALSE)

# Long: one row per gear, q and month, on the complete monthly series (c00
# checks it is complete). `filled` marks April and May 2020.
ts_dat <- gamma %>%
  mutate(GEAR = as.character(GEAR), YEAR = as.integer(YEAR), MONTH = as.integer(MONTH)) %>%
  pivot_longer(all_of(paste0("q", Q_LEVELS)), names_to = "q_raw", values_to = "value") %>%
  mutate(q = as.integer(sub("^q", "", q_raw)), value = as.numeric(value)) %>%
  select(GEAR, q, YEAR, MONTH, filled, value) %>%
  arrange(GEAR, q, YEAR, MONTH)
cat("\nCoastwide rows:", nrow(gamma), "| filled months excluded from the fits:", sum(gamma$filled), "\n")

# -----------------------------------------------------------------------------
# 5) Harmonic fit
# -----------------------------------------------------------------------------
q_label <- function(q) {
  c("0" = "q = 0 (richness)", "1" = "q = 1 (Shannon-weighted diversity)",
    "2" = "q = 2 (dominance-weighted diversity)")[[as.character(q)]]
}

# Adds the harmonic columns sin1, cos1, ... (lib/seasonal.R) to a data frame with MONTH.
make_harmonics <- function(df, K) {
  bind_cols(df, harmonic_X(df$MONTH, K)[, -1, drop = FALSE])
}

# df: one complete monthly series (GEAR, q, YEAR, MONTH, filled, value).
# The curve is fitted on the observed months; the residuals are then placed
# on the monthly grid with NA in the filled months, so the ACF and Ljung-Box
# lags are calendar months. Box.test passes NA through to acf().
fit_series <- function(df, K) {
  df  <- df %>% arrange(YEAR, MONTH)
  obs <- df %>% filter(!filled, !is.na(value)) %>% make_harmonics(K)
  n_obs <- nrow(obs)
  n_months <- n_distinct(paste(obs$YEAR, obs$MONTH, sep = "-"))
  if (n_months < MIN_N_MONTHS) {
    return(list(ok = FALSE, reason = paste0("Insufficient months: ", n_months, " < ", MIN_N_MONTHS),
                n_obs = n_obs, n_months = n_months))
  }

  fml <- stats::as.formula(paste("value ~", paste(names(harmonic_X(1, K))[-1], collapse = " + ")))
  mod <- stats::lm(fml, data = obs)
  co <- stats::coef(mod)
  if (any(!is.finite(co))) return(list(ok = FALSE, reason = "non-finite coefficients", n_obs = n_obs, n_months = n_months))
  amp_k <- function(k) if (k <= K) sqrt(co[[paste0("sin", k)]]^2 + co[[paste0("cos", k)]]^2) else NA_real_
  amp_annual <- amp_k(1)

  resid <- stats::residuals(mod)
  var_y <- stats::var(obs$value)
  var_expl <- if (is.na(var_y) || var_y == 0) NA_real_ else 1 - stats::var(resid) / var_y

  resid_grid <- rep(NA_real_, nrow(df))
  resid_grid[!df$filled & !is.na(df$value)] <- resid
  ac_vals <- as.numeric(stats::acf(resid_grid, plot = FALSE, lag.max = ACF_MAX_LAG, na.action = stats::na.pass)$acf)[-1]
  lb_p <- stats::Box.test(resid_grid, lag = LJUNGBOX_LAG, type = "Ljung-Box")$p.value

  months_df <- tibble::tibble(MONTH = 1:12) %>%
    make_harmonics(K) %>%
    mutate(value_hat = stats::predict(mod, newdata = .)) %>%
    mutate(seasonal_centered = value_hat - mean(value_hat))

  list(
    ok = TRUE,
    n_obs = n_obs,
    n_months = n_months,
    model = mod,
    metrics = tibble::tibble(
      seasonal_amplitude = amp_annual,
      seasonal_peak_to_trough = 2 * amp_annual,
      amp_semiannual = amp_k(2),
      amp_4month = amp_k(3),
      variance_explained_season = var_expl,
      resid_acf_lag1 = ac_vals[1],
      resid_acf_max_abs_lag1_24 = max(abs(ac_vals), na.rm = TRUE),
      ljungbox_p_lag12 = lb_p
    ),
    seasonal_curve = months_df
  )
}

# Amplitude in the first and second half of the years (split year included
# in the first half).
fit_split_half <- function(df, K) {
  yrs <- sort(unique(df$YEAR))
  split_year <- yrs[ceiling(length(yrs) / 2)]
  fe <- fit_series(filter(df, YEAR <= split_year), K)
  fl <- fit_series(filter(df, YEAR > split_year), K)
  ae <- if (isTRUE(fe$ok)) fe$metrics$seasonal_amplitude else NA_real_
  al <- if (isTRUE(fl$ok)) fl$metrics$seasonal_amplitude else NA_real_
  tibble::tibble(split_half_ok = isTRUE(fe$ok) & isTRUE(fl$ok), split_year = as.integer(split_year),
                 amp_early = ae, amp_late = al, amp_diff_late_minus_early = al - ae)
}

# -----------------------------------------------------------------------------
# 6) Fit every gear and q once
# -----------------------------------------------------------------------------
gq_list <- ts_dat %>% group_by(GEAR, q) %>% group_split()
names(gq_list) <- vapply(gq_list, function(x) paste0(x$GEAR[1], "__q", x$q[1]), character(1))
fits_full <- lapply(gq_list, fit_series, K = N_HARMONICS)

failed <- names(fits_full)[!vapply(fits_full, function(f) isTRUE(f$ok), logical(1))]
if (length(failed)) {
  stop("Seasonal fit failed: ", paste(paste0(failed, " (", vapply(fits_full[failed], `[[`, "", "reason"), ")"),
                                      collapse = "; "), call. = FALSE)
}
n_expected <- length(settings$gears) * length(Q_LEVELS)
if (length(fits_full) != n_expected) stop(length(fits_full), " fits; expected ", n_expected, call. = FALSE)

fits_tbl <- bind_rows(lapply(names(gq_list), function(nm) {
  df <- gq_list[[nm]]
  fit <- fits_full[[nm]]
  tibble::tibble(GEAR = df$GEAR[1], q = df$q[1], q_label = q_label(df$q[1]),
                 ok = TRUE, reason = NA_character_, n_obs = fit$n_obs, n_months = fit$n_months) %>%
    bind_cols(fit$metrics, fit_split_half(df, N_HARMONICS))
})) %>%
  arrange(GEAR, q)

cat("\n=== Seasonality summary (harmonic regression, natural scale, observed months) ===\n")
print(as.data.frame(fits_tbl %>% select(-q_label, -ok, -reason) %>%
                      mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7) Outputs
# -----------------------------------------------------------------------------
out_csv <- file.path(dirs$data, "seasonality_summary.csv")
out_rds <- file.path(dirs$data, "seasonality_fits.rds")
readr::write_csv(fits_tbl, out_csv, na = "NA")
saveRDS(fits_full, out_rds)

seasonal_curves <- bind_rows(lapply(names(fits_full), function(nm) {
  df <- gq_list[[nm]]
  fits_full[[nm]]$seasonal_curve %>%
    transmute(GEAR = df$GEAR[1], q = df$q[1], q_label = q_label(df$q[1]), MONTH, seasonal_centered)
}))
p <- ggplot2::ggplot(seasonal_curves, ggplot2::aes(x = MONTH, y = seasonal_centered, group = 1)) +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::geom_line() +
  ggplot2::scale_x_continuous(breaks = 1:12) +
  ggplot2::labs(x = "Month", y = "Seasonal component (centered)",
                title = "Seasonal cycle (harmonic regression; centered)",
                subtitle = paste(N_HARMONICS, "harmonics (12-, 6- and 4-month waves), full record")) +
  ggplot2::facet_grid(q_label ~ GEAR, scales = "free_y") +
  ggplot2::theme_bw()
fig_pdf  <- file.path(dirs$figs, "seasonal_cycle_curves__by_q_by_gear.pdf")
fig_tiff <- file.path(dirs$figs, "seasonal_cycle_curves__by_q_by_gear.tiff")
ggplot2::ggsave(fig_pdf,  plot = p, width = 12, height = 8, units = "in", dpi = 300)
ggplot2::ggsave(fig_tiff, plot = p, width = 12, height = 8, units = "in", dpi = 300, compression = "lzw")

for (f in c(out_csv, out_rds, fig_pdf, fig_tiff)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("c01_seasonality_characterization.R completed\n")
