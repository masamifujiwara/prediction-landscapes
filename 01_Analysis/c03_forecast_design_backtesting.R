# =============================================================================
# Script:  c03_forecast_design_backtesting.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Writes the forecast dates of every diversity series (42 estuary series
#   and 6 coastwide series: estuary or coastwide x sampling method x q). Every
#   diversity model and reference script (c07 to c10, c12, d01) reads this
#   table instead of building its own, so all models use exactly the same
#   forecast dates and training start.
#
# Run after: c02, b09, lib/forecast_rules.R
#
# Inputs:
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09)
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/seasonal_coefs_by_origin_scaled.rds (c02)
#
# Outputs:
#   01_Analysis/Data/c3_forecast_design_backtesting/forecast_origins_diversity.rds
#     scale, BAY_ID (NA coastwide), GEAR, q, origin_i, origin_date,
#     train_start_i, n_train_sampled, max_lead, in_common
#   02_Results/03_Diagnostics/c3_forecast_design_backtesting/c3_origins_by_series.csv
#   01_Analysis/Logs/c3_forecast_design_backtesting/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   - The rules are in forecast_origins() (lib/forecast_rules.R) and the
#     settings come from b09.
#   - Forecast dates: every month from first_origin_date (April 1992, the same
#     for every series) to the month before the series end, except April and
#     May 2020; expanding window; at least 120 sampled months from the series
#     start.
#   - Training starts at month 13 of each series for every model.
#   - in_common marks the common set of forecast dates used at every lead
#     time in the main analysis (at least 24 months before the series end);
#     the other forecast dates serve the checks at leads 1 and 12 and the
#     comparison over periods. max_lead is the longest lead time with a
#     target inside the record (covariate models use at most 12).
#   - Every forecast date must have a seasonal curve in the c02 per-forecast-
#     date coefficient file (the script stops otherwise); the model scripts
#     rebuild the departures from the seasonal curve for a forecast date with
#     anomalies_at_origin() (lib/seasonal.R).
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
script_stem <- "c3_forecast_design_backtesting"   # output folder stem (short form c3_)

in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
in_alpha    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds")
in_gamma    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds")
in_coefs    <- file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "seasonal_coefs_by_origin_scaled.rds")
producers   <- c("b09_forecast_settings.R", rep("c02_deseasonalize_alpha_gamma_scaled.R", 3))
out_name    <- "forecast_origins_diversity.rds"
diag_name   <- "c3_origins_by_series.csv"
KEY         <- "q"

options(scipen = 999, warn = 1, width = 200)
required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 2) Folders, archive, log
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

# -----------------------------------------------------------------------------
# 3) Inputs
# -----------------------------------------------------------------------------
stop_if_missing(c(in_settings, in_alpha, in_gamma, in_coefs), producers)
settings <- load_settings(in_settings)
a <- readRDS(in_alpha); g <- readRDS(in_gamma)
expect_cols(a, c("BAY_ID", "GEAR", "q", "date", "filled"), "alpha_q012_deseas_scaled")
expect_cols(g, c("GEAR", "q", "date", "filled"), "gamma_q012_deseas_scaled")
ts_all <- bind_rows(
  a %>% transmute(scale = "bay", BAY_ID = as.integer(BAY_ID), GEAR, q = as.integer(q), date, filled),
  g %>% transmute(scale = "coastwide", BAY_ID = NA_integer_, GEAR, q = as.integer(q), date, filled))
n_expected <- (length(settings$bays) + 1L) * length(settings$gears) * 3L
coefs <- readRDS(in_coefs)
keys <- c("scale", "BAY_ID", "GEAR", KEY)
series <- ts_all %>% distinct(across(all_of(keys))) %>% arrange(across(all_of(keys)))
if (nrow(series) != n_expected) stop(nrow(series), " series; expected ", n_expected, call. = FALSE)

# -----------------------------------------------------------------------------
# 4) Origins per series (forecast_origins(), lib/forecast_rules.R)
# -----------------------------------------------------------------------------
parts <- vector("list", nrow(series))
for (i in seq_len(nrow(series))) {
  k <- series[i, ]
  d <- ts_all %>% filter(scale == k$scale, GEAR == k$GEAR, .data[[KEY]] == k[[KEY]],
                         (is.na(BAY_ID) & is.na(k$BAY_ID)) | (BAY_ID %in% k$BAY_ID)) %>% arrange(date)
  parts[[i]] <- bind_cols(k, forecast_origins(d$date, d$filled, settings))
}
origins <- bind_rows(parts)

# -----------------------------------------------------------------------------
# 5) Checks
# -----------------------------------------------------------------------------
# Every origin has a seasonal curve at the same row and month, flagged allowed.
chk <- origins %>% left_join(coefs %>% select(all_of(keys), origin_i, coef_date = origin_date, allowed_origin),
                             by = c(keys, "origin_i"))
if (anyNA(chk$coef_date)) stop(sum(is.na(chk$coef_date)), " origins have no seasonal curve in ", basename(in_coefs), call. = FALSE)
if (any(chk$coef_date != chk$origin_date) || !all(chk$allowed_origin)) stop("Origin months or flags differ from the coefficient file", call. = FALSE)
if (any(origins$origin_date %in% settings$skip_origin_months)) stop("An origin falls in skip_origin_months", call. = FALSE)
if (any(origins$origin_date < settings$first_origin_date)) stop("An origin precedes first_origin_date", call. = FALSE)
by_series <- origins %>% group_by(across(all_of(keys))) %>%
  summarise(n_origins = n(), n_common = sum(in_common), first_origin = min(origin_date), last_origin = max(origin_date),
            last_common = max(origin_date[in_common]), min_n_train = min(n_train_sampled), .groups = "drop")
if (n_distinct(by_series$n_origins) != 1L || n_distinct(by_series$n_common) != 1L) {
  cat("Note: origin counts differ between series:\n"); print(as.data.frame(count(by_series, n_origins, n_common)))
}

# -----------------------------------------------------------------------------
# 6) Write
# -----------------------------------------------------------------------------
out_rds  <- file.path(dirs$data, out_name)
out_diag <- file.path(dirs$diag, diag_name)
saveRDS(origins, out_rds)
readr::write_csv(by_series, out_diag, na = "")

cat("\nSeries:", nrow(series), "| origins:", nrow(origins), "| per series", min(by_series$n_origins), "to", max(by_series$n_origins),
    "| common set", min(by_series$n_common), "to", max(by_series$n_common), "\n")
cat("Origins", format(min(origins$origin_date), "%Y-%m"), "to", format(max(origins$origin_date), "%Y-%m"),
    "| last common origin", format(max(by_series$last_common), "%Y-%m"), "| training from month", settings$training_start_index,
    "| fewest training months at the first origin", min(by_series$min_n_train), "\n")
for (f in c(out_rds, out_diag)) cat("Saved:", f, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("c03_forecast_design_backtesting.R completed\n")
