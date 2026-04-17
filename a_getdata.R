# =============================================================================
# Script: a_getdata.R
#
# Purpose:
#   Reference script documenting how raw TPWD organism and station data
#   were originally loaded, combined, and minimally harmonized for analysis.
#   This script is NOT intended to be rerun by external users.
#
#   The primary analysis workflow begins from the saved data object
#   "Data.RData", which is distributed separately (via GitHub or other means).
#   Users are expected to work directly with those prepared data objects
#   rather than executing this script.
#
# Notes:
#   • This script is retained for transparency and reproducibility only.
#   • File paths, directory structure, and local dependencies are specific
#     to the original development environment.
#   • Subsequent data filtering (e.g., gear-specific taxonomic exclusions)
#     is intentionally handled in later scripts, not here.
#
# Author: Masami Fujiwara
# Created: 2025/12/19
# Revised: 2026/01/16 (BS + BT only; exclude GN)
# Revised: 2026/01/16 (Filter to bays 1–8 in BOTH station and organisms)
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 0. Set working directory to this script's folder
# ─────────────────────────────────────────────────────────────────────────────
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else if (!is.null(sys.frames()[[1]]$ofile)) {
  script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
} else {
  stop("Cannot determine script directory. Are you running this interactively?")
}
setwd(script_dir)
cat("📁 Working directory:", getwd(), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load required packages
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 2. Define file paths relative to script_dir
# ─────────────────────────────────────────────────────────────────────────────
base_data_dir <- file.path(dirname(dirname(script_dir)), "0001 TPWD Data Updated")

# --- Bag Seine (BS)
file_bs_1 <- file.path(base_data_dir, "02 Bag Seine", "TPWD BS 05-24.xlsx")
file_bs_2 <- file.path(base_data_dir, "02 Bag Seine", "TPWD BS 82-04.xlsx")

# --- Bay Trawl (BT)
file_bt_1 <- file.path(base_data_dir, "03 Bay Trawl", "TPWD BT 03-24.xlsx")
file_bt_2 <- file.path(base_data_dir, "03 Bay Trawl", "TPWD BT 82-03.xlsx")

# --- Station data (BS + BT only)
file_station_bs <- file.path(
  base_data_dir, "04 Station Data",
  "TPWD BS 82-24 Corrected 24 hour format.xlsx"
)
file_station_bt <- file.path(
  base_data_dir, "04 Station Data",
  "TPWD BT 82-24 Hydro with time.xlsx"
)

# --- Species
file_species <- file.path(
  dirname(dirname(script_dir)),
  "0001 TPWD Data Updated", "06 Species Data",
  "TPWD_Spp_codes_MF.xlsx"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3. Load organism data and combine (stack by rows)
# ─────────────────────────────────────────────────────────────────────────────
cat("📥 Loading organism data (BS + BT only)...\n")

bagseine_organisms <- bind_rows(
  read_excel(file_bs_1),
  read_excel(file_bs_2)
)

baytrawl_organisms <- bind_rows(
  read_excel(file_bt_1),
  read_excel(file_bt_2)
)

# ────────────────────────────────────────────────────────────────
# Helper function: combine_species_615_212_333()
# Robust to SPECIES_CODE being character or numeric
# ────────────────────────────────────────────────────────────────
combine_species_615_212_333 <- function(df) {
  required_cols <- c("SAMPLE_ID", "SPECIES_CODE", "catch", "catch_per_hour", "length_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df |>
    dplyr::mutate(SPECIES_CODE = as.character(SPECIES_CODE))
  
  target_codes <- c("615", "212", "333")
  
  target_samples <- df |>
    dplyr::filter(SPECIES_CODE %in% target_codes) |>
    dplyr::pull(SAMPLE_ID) |>
    unique()
  
  combined <- df |>
    dplyr::filter(SAMPLE_ID %in% target_samples,
                  SPECIES_CODE %in% target_codes) |>
    dplyr::group_by(SAMPLE_ID) |>
    dplyr::summarize(
      MAJOR_AREA     = dplyr::first(MAJOR_AREA),
      YEAR           = dplyr::first(YEAR),
      MONTH          = dplyr::first(MONTH),
      LATITUDE       = dplyr::first(LATITUDE),
      LONGITUDE      = dplyr::first(LONGITUDE),
      DEPTH          = mean(DEPTH, na.rm = TRUE),
      TEMPERATURE    = mean(TEMPERATURE, na.rm = TRUE),
      SALINITY       = mean(SALINITY, na.rm = TRUE),
      DISS_OXYGEN    = mean(DISS_OXYGEN, na.rm = TRUE),
      TURBIDITY      = mean(TURBIDITY, na.rm = TRUE),
      catch          = sum(catch, na.rm = TRUE),
      catch_per_hour = mean(catch_per_hour, na.rm = TRUE),
      length_num     = mean(length_num, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(SPECIES_CODE = "333")
  
  df_out <- df |>
    dplyr::filter(!(SAMPLE_ID %in% target_samples & SPECIES_CODE %in% target_codes)) |>
    dplyr::bind_rows(combined) |>
    dplyr::arrange(SAMPLE_ID)
  
  return(df_out)
}

bagseine_organisms <- combine_species_615_212_333(bagseine_organisms)
baytrawl_organisms <- combine_species_615_212_333(baytrawl_organisms)

cat("✅ Organism data successfully combined (BS + BT).\n")

# ─────────────────────────────────────────────────────────────────────────────
# Keep only required columns in organism data (early)
# ─────────────────────────────────────────────────────────────────────────────
bagseine_organisms <- bagseine_organisms |>
  dplyr::select(SAMPLE_ID, SPECIES_CODE, catch)

baytrawl_organisms <- baytrawl_organisms |>
  dplyr::select(SAMPLE_ID, SPECIES_CODE, catch)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Load station data (BS + BT only)
# ─────────────────────────────────────────────────────────────────────────────
cat("📥 Loading station data (BS + BT only)...\n")

STATION_BS <- read_excel(file_station_bs)
STATION_BT <- read_excel(file_station_bt)

cat("✅ Station data loaded (BS + BT).\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4b. Filter BOTH station tables to bays 1–8, then filter organisms by SAMPLE_ID
# ─────────────────────────────────────────────────────────────────────────────
valid_bays <- 1:8

# Identify bay field in each station table
#   BS: MAJOR_AREA_CODE
#   BT: MAJOR_AREA
stopifnot("MAJOR_AREA_CODE" %in% names(STATION_BS))
stopifnot("MAJOR_AREA" %in% names(STATION_BT))

# Filter station tables
STATION_BS <- STATION_BS |>
  dplyr::mutate(MAJOR_AREA_CODE = as.integer(MAJOR_AREA_CODE)) |>
  dplyr::filter(MAJOR_AREA_CODE %in% valid_bays)%>%
  dplyr::rename(SAMPLE_ID = HYDRO_ID)

STATION_BT <- STATION_BT |>
  dplyr::mutate(MAJOR_AREA = as.integer(MAJOR_AREA)) |>
  dplyr::filter(MAJOR_AREA %in% valid_bays)%>%
  dplyr::rename(SAMPLE_ID = HYDRO_ID)

# Build valid SAMPLE_ID sets from filtered stations
bs_ids <- STATION_BS |>
  dplyr::pull(SAMPLE_ID) |>
  unique()

bt_ids <- STATION_BT |>
  dplyr::pull(SAMPLE_ID) |>
  unique()

# Filter organism tables strictly by those SAMPLE_IDs
bagseine_organisms <- bagseine_organisms |>
  dplyr::filter(SAMPLE_ID %in% bs_ids)

baytrawl_organisms <- baytrawl_organisms |>
  dplyr::filter(SAMPLE_ID %in% bt_ids)

cat("✅ Filtered stations and organisms to bays 1–8.\n")

# Optional: report how many samples remain
cat("────────────────────────────────────────────\n")
cat("Station BS (bays 1–8):", nrow(STATION_BS), "rows\n")
cat("Station BT (bays 1–8):", nrow(STATION_BT), "rows\n")
cat("Bag Seine organisms (bays 1–8 by station IDs):", nrow(bagseine_organisms), "rows\n")
cat("Bay Trawl organisms (bays 1–8 by station IDs):", nrow(baytrawl_organisms), "rows\n")
cat("────────────────────────────────────────────\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Load species data
# ─────────────────────────────────────────────────────────────────────────────
cat("📥 Reading species data from:", file_species, "\n")
species <- readxl::read_excel(file_species)
cat("✅ Species data loaded successfully —", nrow(species), "records\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Save output
# ─────────────────────────────────────────────────────────────────────────────
output_file <- file.path(script_dir, "Data.RData")

save(
  bagseine_organisms,
  baytrawl_organisms,
  STATION_BS,
  STATION_BT,
  species,
  file = output_file
)

cat("💾 All data frames saved successfully to:\n", output_file, "\n")
# =============================================================================
