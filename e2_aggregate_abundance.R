# =============================================================================
# Script: e2_aggregate_abundance.R
#
# Purpose:
# - Extract catch data for 8 target species across Bag Seine and Bay Trawl.
# - Ensure zero-catch stations are correctly accounted for in effort.
# - Compute monthly mean CPUE at the bay level (Alpha).
# - Compute monthly mean CPUE at the coastwide level (Gamma).
# - Apply log(CPUE + 1) transformation to address skewness/zero-inflation.
#
# Outputs:
# - 01_Analysis/Data/e2_aggregate_abundance/alpha_abundance.rds
# - 01_Analysis/Data/e2_aggregate_abundance/gamma_abundance.rds
# - 02_Results/03_Diagnostics/e2_aggregate_abundance/abundance_summary.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
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
# 3) Setup Directories & Targets
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e2_aggregate_abundance"
data_out_dir <- file.path(PATHS$analysis_dir, "Data", script_stem)
diag_out_dir <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)

dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_out_dir, recursive = TRUE, showWarnings = FALSE)

# Define the 8 target species
target_species <- tibble::tribble(
  ~Common_Name,       ~SPECIES_CODE,
  "Atlantic croaker", "602",
  "Spot",             "608",
  "Pinfish",          "633",
  "Blue crab",        "9605",
  "Brown shrimp",     "9618",
  "Bay anchovy",      "601",
  "White shrimp",     "9600",
  "Gulf menhaden",    "604"
)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Load Input Data (Station Core + Organisms)
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(PATHS$clean_data_rdata))
load(PATHS$clean_data_rdata) # Loads bagseine_organisms, baytrawl_organisms

station_core_path <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
stopifnot(file.exists(station_core_path))
station_core <- readRDS(station_core_path)

# ─────────────────────────────────────────────────────────────────────────────
# 5) CPUE Calculation Function
# ─────────────────────────────────────────────────────────────────────────────
aggregate_cpue <- function(station_df, org_df, gear_label) {
  
  # Ensure we only use valid core bays (1-8) mapped to index units
  station_base <- station_df %>%
    filter(MAJOR_AREA %in% 1:8) %>%
    select(SAMPLE_ID, BAY_ID = MAJOR_AREA, YEAR = YEAR_i, MONTH = MONTH_i)
  
  # Cross-join stations with target species to create a complete zero-inflated grid
  grid_full <- tidyr::expand_grid(
    SAMPLE_ID = station_base$SAMPLE_ID,
    SPECIES_CODE = target_species$SPECIES_CODE
  ) %>%
    left_join(station_base, by = "SAMPLE_ID")
  
  # Filter organism data for targets and aggregate by station
  # (in case of multiple measurement records per station)
  org_targets <- org_df %>%
    filter(as.character(SPECIES_CODE) %in% target_species$SPECIES_CODE) %>%
    group_by(SAMPLE_ID, SPECIES_CODE = as.character(SPECIES_CODE)) %>%
    summarise(catch = sum(as.numeric(catch), na.rm = TRUE), .groups = "drop")
  
  # Merge catch into the base grid and set NA catch to 0
  catch_grid <- grid_full %>%
    left_join(org_targets, by = c("SAMPLE_ID", "SPECIES_CODE")) %>%
    mutate(catch = coalesce(catch, 0))
  
  # Calculate Alpha (Bay-level) CPUE
  alpha_df <- catch_grid %>%
    group_by(BAY_ID, YEAR, MONTH, SPECIES_CODE) %>%
    summarise(
      n_samples = n(),
      mean_cpue = mean(catch),
      log_cpue  = log1p(mean_cpue), # log(CPUE + 1)
      .groups = "drop"
    ) %>%
    mutate(GEAR = gear_label)
  
  # Calculate Gamma (Coastwide) CPUE
  gamma_df <- catch_grid %>%
    group_by(YEAR, MONTH, SPECIES_CODE) %>%
    summarise(
      n_samples = n(),
      mean_cpue = mean(catch),
      log_cpue  = log1p(mean_cpue), # log(CPUE + 1)
      .groups = "drop"
    ) %>%
    mutate(GEAR = gear_label)
  
  list(alpha = alpha_df, gamma = gamma_df)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6) Process BS and BT Data
# ─────────────────────────────────────────────────────────────────────────────
cat("Processing Bag Seine (BS) abundance...\n")
bs_results <- aggregate_cpue(station_core$station_bs_core, bagseine_organisms, "BS")

cat("Processing Bay Trawl (BT) abundance...\n")
bt_results <- aggregate_cpue(station_core$station_bt_core, baytrawl_organisms, "BT")

# Combine and attach common names
alpha_abundance <- bind_rows(bs_results$alpha, bt_results$alpha) %>%
  left_join(target_species, by = "SPECIES_CODE") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))) %>%
  select(BAY_ID, GEAR, YEAR, MONTH, date, SPECIES_CODE, Common_Name, n_samples, mean_cpue, log_cpue) %>%
  arrange(BAY_ID, GEAR, SPECIES_CODE, date)

gamma_abundance <- bind_rows(bs_results$gamma, bt_results$gamma) %>%
  left_join(target_species, by = "SPECIES_CODE") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))) %>%
  select(GEAR, YEAR, MONTH, date, SPECIES_CODE, Common_Name, n_samples, mean_cpue, log_cpue) %>%
  arrange(GEAR, SPECIES_CODE, date)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Save Derived Data & Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
alpha_out <- file.path(data_out_dir, "alpha_abundance_targets.rds")
gamma_out <- file.path(data_out_dir, "gamma_abundance_targets.rds")

saveRDS(alpha_abundance, alpha_out)
saveRDS(gamma_abundance, gamma_out)

# Create a diagnostic summary of overall means and sparsity per gear/species
diag_summary <- gamma_abundance %>%
  group_by(GEAR, Common_Name, SPECIES_CODE) %>%
  summarise(
    total_months = n(),
    months_with_zero_catch = sum(mean_cpue == 0),
    overall_mean_cpue = mean(mean_cpue),
    .groups = "drop"
  ) %>%
  arrange(GEAR, desc(overall_mean_cpue))

diag_csv <- file.path(diag_out_dir, "abundance_summary.csv")
write_csv(diag_summary, diag_csv)

cat("✅ Abundance calculation complete.\n")
cat("💾 Saved alpha abundance to:", alpha_out, "\n")
cat("💾 Saved gamma abundance to:", gamma_out, "\n")
cat("🧾 Diagnostics summary saved to:", diag_csv, "\n")