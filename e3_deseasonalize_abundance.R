# =============================================================================
# Script: e3_deseasonalize_abundance.R
#
# Purpose:
# - Load the aggregated, log-transformed abundance data from e2.
# - Quantify and remove strong seasonal migratory pulses for each of the 8 
#   target species using harmonic regression (annual + semi-annual).
# - Produce deseasonalized indices (`value_deseas`) that mirror the exact 
#   data structure used for the Alpha and Gamma diversity forecasting.
#
# Outputs:
# - 01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds
# - 01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds
# - 02_Results/03_Diagnostics/e3_deseasonalize_abundance/seasonality_summary.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS (Windows & Linux Safe)
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  # 2) Rscript / R CMD BATCH: --file=...
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  # 3) source("script.R")
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  # 4) fallback: assume you launched from 01_Analysis
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

options(scipen = 999, warn = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Packages
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) Setup Directories & Inputs
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e3_deseasonalize_abundance"
data_out_dir <- file.path(PATHS$analysis_dir, "Data", script_stem)
diag_out_dir <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)

dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_out_dir, recursive = TRUE, showWarnings = FALSE)

in_alpha <- file.path(PATHS$analysis_dir, "Data", "e2_aggregate_abundance", "alpha_abundance_targets.rds")
in_gamma <- file.path(PATHS$analysis_dir, "Data", "e2_aggregate_abundance", "gamma_abundance_targets.rds")

stopifnot(file.exists(in_alpha), file.exists(in_gamma))
alpha_abund <- readRDS(in_alpha)
gamma_abund <- readRDS(in_gamma)

# Rename log_cpue to 'value' to match the forecasting script expectations
alpha_abund <- alpha_abund %>% mutate(value = log_cpue)
gamma_abund <- gamma_abund %>% mutate(value = log_cpue)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Harmonic Regression & Deseasonalization Functions
# ─────────────────────────────────────────────────────────────────────────────
make_harmonics <- function(df) {
  df %>% mutate(
    m = as.integer(MONTH),
    sin1 = sin(2 * pi * m / 12),
    cos1 = cos(2 * pi * m / 12),
    sin2 = sin(4 * pi * m / 12),
    cos2 = cos(4 * pi * m / 12)
  )
}

# Function to fit harmonic model, center the seasonality to retain overall mean 
# abundance, and calculate residuals (value_deseas)
deseasonalize_group <- function(df) {
  df2 <- df %>% filter(!is.na(value)) %>% make_harmonics()
  
  # Minimum data check
  if (nrow(df2) < 24) {
    return(df %>% mutate(seasonal_component = NA_real_, residual_component = NA_real_, value_deseas = NA_real_, seasonal_amplitude = NA_real_))
  }
  
  mod <- stats::lm(value ~ sin1 + cos1 + sin2 + cos2, data = df2)
  
  # Extract amplitude for diagnostics
  co <- stats::coef(mod)
  amp_annual <- sqrt(unname(co["sin1"])^2 + unname(co["cos1"])^2)
  
  # Generate predictions for all 12 months (the seasonal map)
  month_grid <- data.frame(MONTH = 1:12) %>% make_harmonics()
  month_grid$value_hat <- stats::predict(mod, newdata = month_grid)
  month_grid$seasonal_centered <- month_grid$value_hat - mean(month_grid$value_hat)
  
  # Merge seasonal pattern back to the time series
  df <- df %>%
    left_join(month_grid %>% select(MONTH, seasonal_centered), by = "MONTH") %>%
    mutate(
      seasonal_component = seasonal_centered,
      residual_component = value - seasonal_component,
      value_deseas       = residual_component, # Primary target for prediction
      seasonal_amplitude = amp_annual          # Tracked for diagnostics
    )
  
  return(df)
}

# ─────────────────────────────────────────────────────────────────────────────
# 5) Apply Deseasonalization to Gamma (Coastwide)
# ─────────────────────────────────────────────────────────────────────────────
cat("Applying harmonic deseasonalization to Gamma (Coastwide) Abundance...\n")

gamma_deseas <- gamma_abund %>%
  group_by(GEAR, SPECIES_CODE, Common_Name) %>%
  group_modify(~ deseasonalize_group(.x)) %>%
  ungroup() %>%
  arrange(GEAR, SPECIES_CODE, date)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Apply Deseasonalization to Alpha (Bay-Level)
# ─────────────────────────────────────────────────────────────────────────────
cat("Applying harmonic deseasonalization to Alpha (Bay-level) Abundance...\n")

alpha_deseas <- alpha_abund %>%
  group_by(BAY_ID, GEAR, SPECIES_CODE, Common_Name) %>%
  group_modify(~ deseasonalize_group(.x)) %>%
  ungroup() %>%
  arrange(BAY_ID, GEAR, SPECIES_CODE, date)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Extract Seasonal Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
# Create a quick summary table of the seasonal amplitude per species and gear
gamma_diag <- gamma_deseas %>%
  distinct(GEAR, SPECIES_CODE, Common_Name, seasonal_amplitude) %>%
  arrange(desc(seasonal_amplitude))

alpha_diag <- alpha_deseas %>%
  distinct(BAY_ID, GEAR, SPECIES_CODE, Common_Name, seasonal_amplitude) %>%
  group_by(GEAR, SPECIES_CODE, Common_Name) %>%
  summarise(mean_alpha_amplitude = mean(seasonal_amplitude, na.rm = TRUE), .groups = "drop")

diag_summary <- gamma_diag %>%
  rename(gamma_amplitude = seasonal_amplitude) %>%
  left_join(alpha_diag, by = c("GEAR", "SPECIES_CODE", "Common_Name")) %>%
  arrange(GEAR, desc(gamma_amplitude))

# ─────────────────────────────────────────────────────────────────────────────
# 8) Save Derived Outputs & Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
alpha_out <- file.path(data_out_dir, "alpha_abundance_deseas.rds")
gamma_out <- file.path(data_out_dir, "gamma_abundance_deseas.rds")
diag_out  <- file.path(diag_out_dir, "seasonality_summary.csv")

saveRDS(alpha_deseas, alpha_out)
saveRDS(gamma_deseas, gamma_out)
write_csv(diag_summary, diag_out)

cat("\n✅ Deseasonalization complete.\n")
cat("💾 Saved alpha deseasonalized abundance to:", alpha_out, "\n")
cat("💾 Saved gamma deseasonalized abundance to:", gamma_out, "\n")
cat("🧾 Saved seasonality diagnostics to:", diag_out, "\n")