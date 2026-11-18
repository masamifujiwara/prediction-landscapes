# =============================================================================
# Script:  f04_environment_trends.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Trends in water temperature and salinity for each estuary and sampling
#   method, and for the coast as a whole, for the Methods (study system),
#   Online Resource 1 and Online Resource 2 (Table S9, h17). Trends are
#   estimated estuary by estuary because a coastwide mean of the raw values
#   mixes estuaries with very different salinity, which can hide trends
#   within estuaries.
#
# Run after: b08, b09, lib/helpers.R, lib/forecast_rules.R, lib/seasonal.R
#   (independent of the forecasts)
#
# Inputs:
#   - 01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   - 01_Analysis/Data/b8_build_environment_timeseries/env_alpha_index.rds (b08)
#   - 01_Analysis/Data/b8_build_environment_timeseries/env_gamma_index.rds (b08)
#
# Outputs:
#   01_Analysis/Data/f4_environment_trends/:
#   - f4_environment_trends.rds / .csv   one row per variable and series
#   - f4_yearly_departures.rds           yearly mean departures, for a figure
#   02_Results/03_Diagnostics/f4_environment_trends/f4_summary.csv  counts and
#     ranges by variable, scale and sampling method
#   01_Analysis/Logs/f4_environment_trends/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   1. Estuary by estuary (seven estuaries) and coastwide, each sampling
#      method separately: 16 series per variable (TEMP__mean, SAL__mean of
#      b08).
#   2. Trend measured on the departures from the same three-harmonic seasonal
#      curve as the forecasts (lib/seasonal.R), fitted here to the whole
#      record of sampled months (description, as in f02). April and May 2020
#      (filled, not sampled) are left out.
#   3. Main test: Mann-Kendall test and Sen slope on the yearly means of the
#      departures (years with at least MIN_MONTHS_YEAR sampled months), as f02
#      does for the biological series. Yearly means remove most of the
#      month-to-month autocorrelation; the lag-1 correlation of the yearly
#      means is reported, and a trend-free prewhitened Mann-Kendall test (Yue
#      et al. 2002) is given as a check.
#   4. Within each variable and sampling method, the p-values of the seven
#      estuaries are also adjusted for multiple tests (Benjamini-Hochberg).
#   5. Check: the Sen slope of the monthly departures (all sampled months).
#   6. Check: a simpler calculation (Mann-Kendall test and Sen slope on
#      coastwide yearly means of the raw values, complete years only) is
#      repeated and printed to the log, to show how it differs from the
#      estuary-by-estuary result.
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
MIN_MONTHS_YEAR <- 6L        # sampled months for a yearly mean (as f02)
ALPHA           <- 0.05
VARIABLES       <- c(temperature = "TEMP__mean", salinity = "SAL__mean")
AREA_NAMES <- c(
  "2" = "Galveston Bay",
  "3" = "Matagorda Bay",     "4" = "San Antonio Bay",
  "5" = "Aransas Bay",       "6" = "Corpus Christi Bay",
  "7" = "Upper Laguna Madre", "8" = "Lower Laguna Madre")

script_stem <- "f4_environment_trends"          # a test copy uses its own stem, e.g. "f4_environment_trends_TEST"
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_env_a    <- file.path(PATHS$data_dir, "b8_build_environment_timeseries", "env_alpha_index.rds")
in_env_g    <- file.path(PATHS$data_dir, "b8_build_environment_timeseries", "env_gamma_index.rds")

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr")
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
stop_if_missing(c(in_settings, in_env_a, in_env_g),
                c("b09_forecast_settings.R", rep("b08_build_environment_timeseries.R", 2)))
settings <- load_settings(in_settings)

env_a <- readRDS(in_env_a); env_g <- readRDS(in_env_g)
expect_cols(env_a, c("BAY_ID", "GEAR", "YEAR", "MONTH", "DATE", "filled", unname(VARIABLES)), "env_alpha_index")
expect_cols(env_g, c("GEAR", "YEAR", "MONTH", "DATE", "filled", unname(VARIABLES)), "env_gamma_index")
env <- bind_rows(
  env_a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, YEAR = as.integer(YEAR), MONTH = as.integer(MONTH),
                      DATE, filled, TEMP__mean, SAL__mean),
  env_g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, YEAR = as.integer(YEAR), MONTH = as.integer(MONTH),
                      DATE, filled, TEMP__mean, SAL__mean)) %>%
  mutate(sampled = !filled & !(format(DATE, "%Y-%m") %in% format(settings$gap_months, "%Y-%m")),
         area = ifelse(scale == "bay", unname(AREA_NAMES[as.character(BAY_ID)]), "Coastwide"))
if (any(is.na(env$area))) stop("Unknown BAY_ID: ", paste(unique(env$BAY_ID[is.na(env$area)]), collapse = ", "), call. = FALSE)
if (any(env$filled & !(format(env$DATE, "%Y-%m") %in% format(settings$gap_months, "%Y-%m")))) {
  stop("Filled months other than the gap months", call. = FALSE)
}
KEYS <- c("scale", "BAY_ID", "area", "GEAR")
env$sid <- paste(env$scale, env$BAY_ID, env$GEAR, sep = "|")
n_series <- length(unique(env$sid))
if (n_series != 16L) stop(n_series, " series; expected 16 (7 estuaries and the coast, 2 sampling methods)", call. = FALSE)
cat("Series per variable:", n_series, "| rows:", nrow(env), "\n")

# -----------------------------------------------------------------------------
# 3) Helpers
# -----------------------------------------------------------------------------
# Mann-Kendall test (via Kendall's tau against time) and Sen slope per year.
# Same calculation as in f02.
mann_kendall <- function(y, t) {
  ok <- is.finite(y); y <- y[ok]; t <- t[ok]
  ct <- suppressWarnings(stats::cor.test(t, y, method = "kendall", exact = FALSE))
  ij <- utils::combn(length(y), 2)
  sen <- stats::median((y[ij[2, ]] - y[ij[1, ]]) / (t[ij[2, ]] - t[ij[1, ]]))
  c(tau = unname(ct$estimate), p = ct$p.value, sen_per_year = sen)
}
# Trend-free prewhitening (Yue et al. 2002): remove the Sen trend, remove the
# lag-1 correlation of the residuals if it is significant, add the trend back,
# and test again. Returns the lag-1 correlation and the new p-value.
mk_tfpw <- function(y, t) {
  ok <- is.finite(y); y <- y[ok]; t <- t[ok]
  b <- mann_kendall(y, t)[["sen_per_year"]]
  r <- y - b * t
  r1 <- stats::acf(r, lag.max = 1, plot = FALSE)$acf[2]
  if (abs(r1) <= stats::qnorm(1 - ALPHA / 2) / sqrt(length(r))) {
    return(c(r1 = r1, prewhitened = 0, p_tfpw = mann_kendall(y, t)[["p"]]))
  }
  rw <- r[-1] - r1 * r[-length(r)]
  yw <- rw + b * t[-1]
  c(r1 = r1, prewhitened = 1, p_tfpw = mann_kendall(yw, t[-1])[["p"]])
}

# -----------------------------------------------------------------------------
# 4) Per variable and series
# -----------------------------------------------------------------------------
res <- list(); yr_rows <- list()
for (v in names(VARIABLES)) {
  col <- VARIABLES[[v]]
  for (s in sort(unique(env$sid))) {
    d <- env[env$sid == s, ]; d <- d[order(d$DATE), ]
    if (any(diff(d$YEAR * 12L + d$MONTH) != 1L)) stop(s, ": not a complete monthly series", call. = FALSE)
    x <- d[[col]]
    use <- d$sampled & is.finite(x)
    if (sum(use) < 120L) stop(s, " ", v, ": fewer than 120 sampled months", call. = FALSE)
    fit <- fit_harmonic(d$MONTH[use], x[use])
    dep <- rep(NA_real_, nrow(d)); dep[use] <- x[use] - curve_from_fit(fit)[d$MONTH[use]]
    yr <- data.frame(YEAR = d$YEAR[use], dep = dep[use]) %>% group_by(YEAR) %>%
      summarise(n = n(), mean_dep = mean(dep), .groups = "drop") %>% filter(n >= MIN_MONTHS_YEAR)
    mk <- mann_kendall(yr$mean_dep, yr$YEAR)
    pw <- mk_tfpw(yr$mean_dep, yr$YEAR)
    t_year <- (d$YEAR + (d$MONTH - 0.5) / 12)[use]
    sen_monthly <- mann_kendall(dep[use], t_year)[["sen_per_year"]]
    res[[length(res) + 1L]] <- data.frame(variable = v, d[1, KEYS],
      first_month = min(d$DATE[use]), last_month = max(d$DATE[use]), n_sampled = sum(use), n_years = nrow(yr),
      median = stats::median(x[use]), q10 = unname(stats::quantile(x[use], 0.1)), q90 = unname(stats::quantile(x[use], 0.9)),
      seasonal_r2 = summary(fit)$r.squared, sd_departure = stats::sd(dep[use]),
      mk_tau = mk[["tau"]], mk_p = mk[["p"]], sen_per_decade = 10 * mk[["sen_per_year"]],
      r1_yearly = pw[["r1"]], prewhitened = pw[["prewhitened"]] == 1, p_tfpw = pw[["p_tfpw"]],
      sen_monthly_per_decade = 10 * sen_monthly, stringsAsFactors = FALSE, row.names = NULL)
    yr_rows[[length(yr_rows) + 1L]] <- data.frame(variable = v, d[1, KEYS], yr, row.names = NULL)
  }
}
tab <- do.call(rbind, res); yr_tab <- do.call(rbind, yr_rows)
tab <- tab %>% group_by(variable, scale, GEAR) %>%
  mutate(p_bh = ifelse(scale == "bay", stats::p.adjust(mk_p, method = "BH"), mk_p)) %>% ungroup() %>%
  mutate(direction = case_when(mk_p < ALPHA & mk_tau > 0 ~ "increase", mk_p < ALPHA & mk_tau < 0 ~ "decrease", TRUE ~ "no clear trend")) %>%
  as.data.frame()

# Checks: 32 rows, finite core values.
if (nrow(tab) != 32L || anyDuplicated(tab[, c("variable", "scale", "BAY_ID", "GEAR")])) stop("Expected 32 unique rows", call. = FALSE)
core <- c("mk_tau", "mk_p", "sen_per_decade", "p_tfpw", "sen_monthly_per_decade", "median")
bad <- sapply(core, function(k) sum(!is.finite(tab[[k]])))
if (any(bad > 0)) { print(bad[bad > 0]); stop("Non-finite values (listed above)", call. = FALSE) }

# -----------------------------------------------------------------------------
# 5) Check: the simpler calculation on coastwide yearly means of the raw
#    values (complete years)
# -----------------------------------------------------------------------------
check_old <- env %>% filter(scale == "coastwide") %>% group_by(GEAR, YEAR) %>%
  summarise(n = sum(sampled), TEMP = mean(TEMP__mean[sampled]), SAL = mean(SAL__mean[sampled]), .groups = "drop") %>%
  filter(n == 12L) %>% group_by(GEAR) %>%
  summarise(temp_sen_per_decade = 10 * mann_kendall(TEMP, YEAR)[["sen_per_year"]], temp_tau = mann_kendall(TEMP, YEAR)[["tau"]],
            sal_sen_per_decade = 10 * mann_kendall(SAL, YEAR)[["sen_per_year"]], sal_tau = mann_kendall(SAL, YEAR)[["tau"]],
            sal_p = mann_kendall(SAL, YEAR)[["p"]], .groups = "drop")

# -----------------------------------------------------------------------------
# 6) Summary, write
# -----------------------------------------------------------------------------
rng <- function(x) sprintf("%.2f to %.2f", min(x), max(x))
summ <- tab %>% group_by(variable, scale, GEAR) %>% summarise(
  series = n(), increase = sum(direction == "increase"), decrease = sum(direction == "decrease"),
  increase_bh = sum(p_bh < ALPHA & mk_tau > 0), decrease_bh = sum(p_bh < ALPHA & mk_tau < 0),
  increase_tfpw = sum(p_tfpw < ALPHA & mk_tau > 0), decrease_tfpw = sum(p_tfpw < ALPHA & mk_tau < 0),
  sen_per_decade_range = rng(sen_per_decade), sen_per_decade_median = stats::median(sen_per_decade), .groups = "drop") %>%
  as.data.frame()

out_rds <- file.path(dirs$data, "f4_environment_trends.rds"); out_csv <- file.path(dirs$data, "f4_environment_trends.csv")
out_yr  <- file.path(dirs$data, "f4_yearly_departures.rds"); out_sum <- file.path(dirs$diag, "f4_summary.csv")
old <- list.files(dirs$data, pattern = "[.](rds|csv)$", full.names = TRUE); if (length(old)) archive_files(old)
old_d <- list.files(dirs$diag, pattern = "[.](csv)$", full.names = TRUE); if (length(old_d)) archive_files(old_d)
saveRDS(tab, out_rds); readr::write_csv(tab, out_csv, na = "")
saveRDS(yr_tab, out_yr); readr::write_csv(summ, out_sum, na = "")

show <- tab %>% transmute(variable, area, GEAR, n_years, median = round(median, 1), sen_per_decade = round(sen_per_decade, 2),
                          tau = round(mk_tau, 2), p = signif(mk_p, 2), p_bh = signif(p_bh, 2), r1 = round(r1_yearly, 2),
                          p_tfpw = signif(p_tfpw, 2), sen_monthly = round(sen_monthly_per_decade, 2), direction)
cat("\nTrend of the yearly mean departures from the seasonal curve (Sen slope per decade; Mann-Kendall):\n")
print(show, row.names = FALSE)
cat("\nSummary by variable, scale and sampling method (counts at p < ", ALPHA, "):\n", sep = ""); print(summ, row.names = FALSE)
cat("\nCheck: hand calculation (coastwide yearly means of the raw values, complete years):\n")
print(as.data.frame(check_old), row.names = FALSE)
for (f in c(out_rds, out_csv, out_yr, out_sum)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("f04_environment_trends.R completed\n")
