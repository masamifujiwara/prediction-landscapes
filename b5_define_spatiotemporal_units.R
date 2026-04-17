# =============================================================================
# Script: b5_define_spatiotemporal_units.R
#
# Purpose (project directory convention):
#   Define the aggregation units ("index keys") for alpha (α) and gamma (γ) series.
#
# Recommended default:
#   - α units: BAY_ID × GEAR × YEAR × MONTH
#   - γ units: GEAR × YEAR × MONTH  (coastwide pooling across bays)
#
# This script:
#   - Loads station_core (derived data from b3)
#   - Standardizes BAY_ID, YEAR, MONTH (integer)
#   - Defines a complete YEAR×MONTH time grid (based on observed range)
#   - Summarizes sample counts for α cells and γ cells
#   - Flags sparse α cells (< N_min samples)
#   - Saves index_keys.rds (derived data) for downstream scripts (e.g., b6)
#
# Inputs (derived data):
#   - 01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds
#
# Outputs (derived data):
#   - 01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds
#
# Diagnostics:
#   - 02_Results/03_Diagnostics/b5_define_spatiotemporal_units/
#       * b5_counts_alpha_cells__bay_gear_year_month.csv
#       * b5_sparse_alpha_cells__lt_<N_min>_samples.csv
#       * b5_counts_gamma_cells__gear_year_month.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 2) Rscript (works with nohup)
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 3) Fallback: current working directory
  # (Safe for nohup if you cd into 01_Analysis first)
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

get_project_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  normalizePath(file.path(here, ".."), winslash = "/", mustWork = TRUE)
}
project_dir <- get_project_root()

paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)  # loads PATHS
cat("🧭 Project directory:", PATHS$project_dir, "\n")

options(scipen = 999, warn = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Packages (minimal)
# ─────────────────────────────────────────────────────────────────────────────
required_pkgs <- c("dplyr", "tidyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) Parameters
# ─────────────────────────────────────────────────────────────────────────────
valid_bays <- 1:8
N_min <- 5  # sparse-cell threshold (alpha)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Load station_core (derived data from b3)
# ─────────────────────────────────────────────────────────────────────────────
station_core_path <- file.path(
  PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds"
)
stopifnot(file.exists(station_core_path))

station_core <- readRDS(station_core_path)
stopifnot(all(c("station_bs_core", "station_bt_core") %in% names(station_core)))

station_all <- bind_rows(
  station_core$station_bs_core,
  station_core$station_bt_core
)

# Standardize keys
station_all <- station_all %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  ) %>%
  filter(BAY_ID %in% valid_bays)

# Sanity checks
stopifnot(!any(is.na(station_all$BAY_ID)))
stopifnot(!any(is.na(station_all$YEAR)))
stopifnot(!any(is.na(station_all$MONTH)))
stopifnot(all(station_all$MONTH >= 1 & station_all$MONTH <= 12))

# ─────────────────────────────────────────────────────────────────────────────
# 5) Define time grid (YEAR × MONTH) based on observed range
# ─────────────────────────────────────────────────────────────────────────────
year_min <- min(station_all$YEAR, na.rm = TRUE)
year_max <- max(station_all$YEAR, na.rm = TRUE)

time_grid <- tidyr::expand_grid(
  YEAR = year_min:year_max,
  MONTH = 1:12
) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

# ─────────────────────────────────────────────────────────────────────────────
# 6) Lookup tables: bay coverage over time
# ─────────────────────────────────────────────────────────────────────────────
bay_lookup <- tibble(BAY_ID = valid_bays)

bay_time_coverage <- station_all %>%
  distinct(BAY_ID, YEAR, MONTH) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH)) %>%
  group_by(BAY_ID) %>%
  summarise(
    first_year_month  = min(YEAR_MONTH),
    last_year_month   = max(YEAR_MONTH),
    n_months_observed = n(),
    .groups = "drop"
  ) %>%
  arrange(BAY_ID)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Counts for α cells: BAY_ID × GEAR × YEAR × MONTH
# ─────────────────────────────────────────────────────────────────────────────
alpha_counts <- station_all %>%
  count(BAY_ID, GEAR, YEAR, MONTH, name = "n_samples") %>%
  arrange(BAY_ID, GEAR, YEAR, MONTH)

alpha_dist <- alpha_counts %>%
  summarise(
    n_cells   = n(),
    min_n     = min(n_samples),
    p10       = as.numeric(stats::quantile(n_samples, 0.10)),
    median_n  = stats::median(n_samples),
    p90       = as.numeric(stats::quantile(n_samples, 0.90)),
    max_n     = max(n_samples)
  )

alpha_sparse <- alpha_counts %>%
  filter(n_samples < N_min) %>%
  arrange(BAY_ID, GEAR, YEAR, MONTH)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Counts for γ cells: GEAR × YEAR × MONTH (coastwide pooling)
# ─────────────────────────────────────────────────────────────────────────────
gamma_counts <- station_all %>%
  count(GEAR, YEAR, MONTH, name = "n_samples") %>%
  arrange(GEAR, YEAR, MONTH)

# ─────────────────────────────────────────────────────────────────────────────
# 9) Save diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
diag_dir <- PATHS$diag_dir

alpha_counts_path <- file.path(diag_dir, "b5_counts_alpha_cells__bay_gear_year_month.csv")
alpha_sparse_path <- file.path(diag_dir, paste0("b5_sparse_alpha_cells__lt_", N_min, "_samples.csv"))
gamma_counts_path <- file.path(diag_dir, "b5_counts_gamma_cells__gear_year_month.csv")

write_csv(alpha_counts, alpha_counts_path)
write_csv(alpha_sparse, alpha_sparse_path)
write_csv(gamma_counts, gamma_counts_path)

cat("🧾 Saved: ", alpha_counts_path, "\n", sep = "")
cat("🧾 Saved: ", alpha_sparse_path, "\n", sep = "")
cat("🧾 Saved: ", gamma_counts_path, "\n", sep = "")

cat("\n=== α cell sample count distribution (BAY_ID×GEAR×YEAR×MONTH) ===\n")
print(alpha_dist)

# ─────────────────────────────────────────────────────────────────────────────
# 10) Save index_keys (derived data) -> 01_Analysis/Data/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
out_dir <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

index_keys <- list(
  keys = list(
    spatial_key      = "BAY_ID",
    bay_field_source = "MAJOR_AREA",
    valid_bays       = valid_bays,
    alpha_unit       = c("BAY_ID", "GEAR", "YEAR", "MONTH"),
    gamma_unit       = c("GEAR", "YEAR", "MONTH")
  ),
  time_grid = time_grid,
  lookups = list(
    bay_lookup        = bay_lookup,
    bay_time_coverage = bay_time_coverage
  ),
  diagnostics = list(
    alpha_counts = alpha_counts,
    gamma_counts = gamma_counts,
    alpha_dist   = alpha_dist,
    N_min_sparse_flag = N_min
  )
)

out_rds <- file.path(out_dir, "index_keys.rds")
saveRDS(index_keys, out_rds)
cat("💾 Saved: ", out_rds, "\n", sep = "")

cat("✅ b5_define_spatiotemporal_units.R completed successfully\n")