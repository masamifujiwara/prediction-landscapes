# =============================================================================
# Script:  b09_forecast_settings.R
# Author:  Masami Fujiwara
#
# Purpose:
#   1. Writes settings.rds, the forecast settings that every later script
#      reads (with load_settings()) instead of setting them itself: lead
#      times, covariate lags, forecast dates, training start, tolerance and
#      the COVID rule.
#   2. Fills April and May 2020 (no sampling anywhere, COVID) in every
#      estuary and coastwide q0, q1, q2 diversity series and marks them
#      filled = TRUE. The b06 and b07 outputs stay observed only; every later
#      diversity script reads the filled series written here.
#
# Run after: b05, b06, b07, lib/forecast_rules.R
#
# Inputs:
#   01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds (b05)
#   01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds (b06)
#   01_Analysis/Data/01_TimeSeries/gamma_q012.rds (b07)
#
# Outputs:
#   01_Analysis/Data/b9_forecast_settings/
#     settings.rds            forecast settings
#     alpha_q012_filled.rds   b06 output without abund_vec, plus the filled
#                             months and a `filled` column
#     gamma_q012_filled.rds   b07 output, likewise
#   02_Results/03_Diagnostics/b9_forecast_settings/
#     b9_filled_target_values.csv  the filled values, with the observed
#                                  values in March and June 2020
#     b9_settings.txt              the settings in plain text
#   01_Analysis/Logs/b9_forecast_settings/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   1. April and May 2020 are the only months missing inside any diversity
#      series; the script stops if it finds any other gap.
#   2. The diversity values are filled by seasonal-anomaly interpolation on
#      the log scale: the departure from each calendar month's mean is
#      interpolated linearly across the gap and the mean is added back
#      (fill_gap_seasonal_anomaly() in lib/forecast_rules.R).
#   3. The calendar-month means use only months up to June 2020 (the first
#      observed month after the gap), so the fill uses no later data.
#   4. Forecasts whose target month is April or May 2020 are excluded from
#      every score (all models, all lead times, the references included)
#      through drop_unscored_targets() in lib/forecast_rules.R. Forecasts that
#      use those months as lags or inside the training window are kept.
#   5. April and May 2020 are not used as forecast dates
#      (is_allowed_origin()). The values there are fills that use June 2020
#      (for the covariates too, b08), and a forecast of last month's value
#      from them would start from a filled value.
#   6. After this fill, and with the series ends trimmed in b06 to b08, no
#      model series has a missing month, so every model sees the same
#      complete series.
#   7. The filled values must keep q0 >= q1 >= q2 (the script stops
#      otherwise).
#   8. Lead times: 1 to 24 months for the models without covariates and the
#      references (horizons); 1 to 12 months for the covariate models
#      (horizons_covariate). One fit per forecast date gives every lead time
#      for the models without covariates; leads beyond 12 months would need a
#      lag-24 covariate set for the covariate models.
#   9. Covariate sets: lags 1 and 12 for lead 1, lag 12 for leads 2 to 12 (two
#      fits per forecast date). No lag-24 covariate.
#  10. Forecast dates: every month (origin_step 1) from the first month with
#      min_train_months (120) sampled months, counted from the series start;
#      expanding window. The main analysis uses the same forecast dates at
#      every lead time: the last one is max(horizons) months before the
#      series end (origin_set = "common"). Leads 1 and 12 are also scored from
#      all possible forecast dates as a check (check_all_origins_horizons),
#      and the comparison over periods uses all possible forecast dates
#      (period_origins = "all"). No series is scored from a forecast date
#      before first_origin_date (April 1992), the first possible forecast
#      date of the series that start latest, so every diversity and abundance
#      series is scored over the same calendar months. Only the forecast
#      dates are aligned: the earlier data stay in every series (seasonal
#      curves, training, lags).
#  11. Training start: every model, with or without covariates, starts its
#      training at month training_start_index (13) of each series, the first
#      month with a lag-12 covariate, so paired models see the same years. The
#      seasonal curves are not affected (lib/seasonal.R fits them on all
#      sampled months up to the forecast date).
#  12. Tolerance for absolute accuracy: 30% of the series mean level for
#      diversity and 0.5 log(CPUE + 1) units for abundance, about twice the
#      haul-to-haul sampling uncertainty of a typical estuary-month value.
#      The range in tolerance_range is shown in the supplement.
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
# 1) Settings written to settings.rds
# -----------------------------------------------------------------------------
# Output folder stem (short form b9_, as b03 to b08).
script_stem <- "b9_forecast_settings"

in_keys  <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.rds")
in_alpha <- file.path(PATHS$data_dir, "b6_compute_hill_alpha_q012", "alpha_q012.rds")
in_gamma <- file.path(PATHS$data_dir, "01_TimeSeries", "gamma_q012.rds")

GAP_MONTHS <- as.Date(c("2020-04-01", "2020-05-01"))   # COVID: no sampling (Method note 1)
CLIM_END   <- as.Date("2020-06-01")                    # Method note 3

settings <- list(
  # Lead times (Method note 8).
  horizons               = 1:24,     # models without covariates and reference forecasts
  horizons_covariate     = 1:12,     # ARIMAX and state-space models with covariates
  # Covariate lags per lead time (Method note 9): list name = leads it applies to.
  covariate_sets         = list(lead_1 = c(1L, 12L), leads_2_12 = 12L),
  # Forecast dates and training (Method notes 10 and 11).
  min_train_months       = 120L,     # first origin: 120 sampled months from the series start
  origin_step            = 1L,
  origin_set             = "common", # same origins at every lead: last origin = series end - max(horizons)
  first_origin_date      = as.Date("1992-04-01"),  # common first forecast date for every series (Method note 10)
  check_all_origins_horizons = c(1L, 12L),  # also scored from all possible origins (check)
  period_origins         = "all",    # period comparison (1992-2007, 2008-2019, 2020-2024) uses all origins
  training_start_index   = 13L,      # every model trains from month 13 of its series
  min_train_non_na       = 24L,      # fit only with at least this many training values
  # Absolute accuracy (Method note 12).
  tolerance              = list(diversity_rel = 0.30, abundance_log = 0.5),
  tolerance_range        = list(diversity_rel = c(0.20, 0.60), abundance_log = c(0.3, 0.9)),
  pi_levels              = c(80, 95),
  freq                   = 12L,
  seed                   = 20261006L, # as c09; scripts using future.seed = TRUE are unchanged
  gears                  = c("BS", "BT"),
  bays                   = NULL,      # from index_keys.rds below
  # The COVID rule (Method notes 1 to 5).
  gap_months             = GAP_MONTHS,
  target_fill            = list(method = "seasonal anomaly, log scale", climatology_end = CLIM_END),
  unscored_target_months = GAP_MONTHS,
  skip_origin_months     = GAP_MONTHS,
  # No open settings. The seasonal curve (refitted at each forecast date)
  # and its harmonics are set in lib/seasonal.R.
  open_decisions         = character(0)
)

options(scipen = 999, warn = 1)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

# -----------------------------------------------------------------------------
# 4) Inputs
# -----------------------------------------------------------------------------
stop_if_missing(c(in_keys, in_alpha, in_gamma),
                c("b05_define_spatiotemporal_units.R", "b06_compute_hill_alpha_q012.R",
                  "b07_compute_hill_gamma_q012.R"))
index_keys <- readRDS(in_keys)
if (is.null(index_keys$keys$valid_bays)) stop("index_keys.rds has no keys$valid_bays", call. = FALSE)
settings$bays <- as.integer(index_keys$keys$valid_bays)

alpha <- readRDS(in_alpha)
gamma <- readRDS(in_gamma)
expect_cols(alpha, c("BAY_ID", "GEAR", "YEAR", "MONTH", "q0", "q1", "q2"), "alpha_q012")
expect_cols(gamma, c("GEAR", "YEAR", "MONTH", "q0", "q1", "q2"), "gamma_q012")
alpha$abund_vec <- NULL
gamma$abund_vec <- NULL

if (!setequal(unique(alpha$BAY_ID), settings$bays)) stop("Bays in alpha_q012 differ from keys$valid_bays", call. = FALSE)
if (!setequal(unique(alpha$GEAR), settings$gears) || !setequal(unique(gamma$GEAR), settings$gears)) {
  stop("Gears differ from ", paste(settings$gears, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# 5) Fill April and May 2020 in every series
# -----------------------------------------------------------------------------
# One series = one bay x gear (alpha) or one gear (gamma). The grid runs from
# the series' first to its last month; the missing months must be exactly
# GAP_MONTHS (Method note 1).
fill_series <- function(df) {
  t <- month_index(df$YEAR, df$MONTH)
  grid <- tibble::tibble(t = seq(min(t), max(t))) %>%
    mutate(YEAR = as.integer(t %/% 12L), MONTH = as.integer(t %% 12L + 1L))
  out <- grid %>% left_join(df, by = c("YEAR", "MONTH"))
  for (k in intersect(c("BAY_ID", "GEAR"), names(df))) out[[k]] <- df[[k]][1]
  gap <- is.na(out$q0)
  gap_dates <- month_date(out$YEAR[gap], out$MONTH[gap])
  if (!setequal(gap_dates, GAP_MONTHS) || anyDuplicated(gap_dates)) {
    stop("Series ", paste(unlist(df[1, intersect(c("BAY_ID", "GEAR"), names(df))]), collapse = " "),
         ": missing months are ", paste(format(gap_dates, "%Y-%m"), collapse = ", "),
         "; expected exactly ", paste(format(GAP_MONTHS, "%Y-%m"), collapse = ", "), call. = FALSE)
  }
  for (v in c("q0", "q1", "q2")) {
    out[[v]] <- fill_gap_seasonal_anomaly(out[[v]], out$YEAR, out$MONTH, fill = gap, clim_end = CLIM_END)
  }
  if ("n_samples" %in% names(out)) out$n_samples[gap] <- 0L
  if ("YEAR_MONTH" %in% names(out)) out$YEAR_MONTH[gap] <- sprintf("%04d-%02d", out$YEAR[gap], out$MONTH[gap])
  out$filled <- gap
  out %>% select(-t)
}

alpha_filled <- alpha %>%
  group_by(BAY_ID, GEAR) %>% group_split() %>% lapply(fill_series) %>% bind_rows() %>%
  select(all_of(names(alpha)), filled) %>% arrange(BAY_ID, GEAR, YEAR, MONTH)
gamma_filled <- gamma %>%
  group_by(GEAR) %>% group_split() %>% lapply(fill_series) %>% bind_rows() %>%
  select(all_of(names(gamma)), filled) %>% arrange(GEAR, YEAR, MONTH)

# -----------------------------------------------------------------------------
# 6) Checks
# -----------------------------------------------------------------------------
n_series <- nrow(distinct(alpha, BAY_ID, GEAR)) + nrow(distinct(gamma, GEAR))
n_filled <- sum(alpha_filled$filled) + sum(gamma_filled$filled)
if (n_filled != n_series * length(GAP_MONTHS)) stop("Filled ", n_filled, " months; expected ", n_series * length(GAP_MONTHS), call. = FALSE)
if (nrow(alpha_filled) != nrow(alpha) + sum(alpha_filled$filled)) stop("alpha row count changed beyond the filled months", call. = FALSE)
if (nrow(gamma_filled) != nrow(gamma) + sum(gamma_filled$filled)) stop("gamma row count changed beyond the filled months", call. = FALSE)
for (d in list(alpha_filled, gamma_filled)) {
  if (anyNA(d[, c("q0", "q1", "q2")])) stop("NA in q0, q1 or q2 after filling", call. = FALSE)
  if (any(d$q1 > d$q0 + 1e-9 | d$q2 > d$q1 + 1e-9)) stop("q0 >= q1 >= q2 fails after filling (method note 7)", call. = FALSE)
}
# Observed rows are unchanged.
chk <- function(orig, filled, keys) {
  m <- filled %>% filter(!filled) %>% select(-filled)
  isTRUE(all.equal(as.data.frame(arrange(orig, across(all_of(keys)))), as.data.frame(arrange(m, across(all_of(keys)))),
                   check.attributes = FALSE))
}
if (!chk(alpha, alpha_filled, c("BAY_ID", "GEAR", "YEAR", "MONTH"))) stop("Observed alpha rows changed", call. = FALSE)
if (!chk(gamma, gamma_filled, c("GEAR", "YEAR", "MONTH"))) stop("Observed gamma rows changed", call. = FALSE)
cat("Series:", n_series, "| months filled:", n_filled, "| observed rows unchanged\n")

# -----------------------------------------------------------------------------
# 7) Outputs
# -----------------------------------------------------------------------------
out_settings <- file.path(dirs$data, "settings.rds")
out_alpha    <- file.path(dirs$data, "alpha_q012_filled.rds")
out_gamma    <- file.path(dirs$data, "gamma_q012_filled.rds")
out_values   <- file.path(dirs$diag, "b9_filled_target_values.csv")
out_text     <- file.path(dirs$diag, "b9_settings.txt")

saveRDS(settings, out_settings)
saveRDS(alpha_filled, out_alpha)
saveRDS(gamma_filled, out_gamma)

# Filled values with the observed neighbors (March and June 2020).
near <- c(month_date(2020, 3), GAP_MONTHS, CLIM_END)
values <- bind_rows(
  alpha_filled %>% mutate(scale = "bay"),
  gamma_filled %>% mutate(scale = "coastwide", BAY_ID = NA_integer_)
) %>%
  filter(month_date(YEAR, MONTH) %in% near) %>%
  transmute(scale, BAY_ID, GEAR, YEAR, MONTH, filled,
            q0 = round(q0, 3), q1 = round(q1, 3), q2 = round(q2, 3)) %>%
  arrange(scale, BAY_ID, GEAR, YEAR, MONTH)
readr::write_csv(values, out_values, na = "")

writeLines(utils::capture.output(print_settings(settings)), out_text)

cat("\n=== Settings ===\n")
print_settings(settings)
cat("\n=== Filled months (q0, q1, q2) ===\n")
print(as.data.frame(values %>% filter(filled) %>% select(-filled)), row.names = FALSE)

for (p in c(out_settings, out_alpha, out_gamma, out_values, out_text)) cat("Saved:", p, "\n")

end_log("finished")
cat("b09_forecast_settings.R completed\n")
