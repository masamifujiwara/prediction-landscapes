# =============================================================================
# Script: b3_harmonize_station_metadata.R
#
# Purpose:
#   - Create harmonized "core" station metadata tables for BS and BT:
#       * numeric YEAR_i, MONTH_i (from YEAR/MONTH)
#       * DATE from COMPLETION_DTTM
#       * consistent MAJOR_AREA (integer) across gears
#       * standardized covariate names:
#           TEMP, SAL, DEPTH_M, DO, TURB, LAT, LON
#   - Save derived output:
#       01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds
#   - Save diagnostics:
#       02_Results/03_Diagnostics/b3_harmonize_station_metadata/
#         - b3_missingness_by_year__station_core.csv
#         - b3_station_coverage_counts__year_month.png
#
# Inputs:
#   - 01_Analysis/Data.RData (STATION_BS, STATION_BT)
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

# ─────────────────────────────────────────────────────────────────────────────
# 2) Options
# ─────────────────────────────────────────────────────────────────────────────
options(scipen = 999, warn = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Packages (minimal)
# ─────────────────────────────────────────────────────────────────────────────
required_pkgs <- c("dplyr", "ggplot2", "readr", "lubridate")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(lubridate)
})

# ─────────────────────────────────────────────────────────────────────────────
# 4) Load Data.RData
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(PATHS$clean_data_rdata))
load(PATHS$clean_data_rdata)

stopifnot(exists("STATION_BS"), exists("STATION_BT"))

# ─────────────────────────────────────────────────────────────────────────────
# 5) Harmonization utilities
# ─────────────────────────────────────────────────────────────────────────────
to_int_safely <- function(x) suppressWarnings(as.integer(x))

harmonize_station <- function(df, gear_label) {
  
  major_area_field <- if (gear_label == "BS") "MAJOR_AREA_CODE" else "MAJOR_AREA"
  
  required_fields <- c(
    "SAMPLE_ID",
    major_area_field,
    "COMPLETION_DTTM", "YEAR", "MONTH",
    "GRID", "LATITUDE", "LONGITUDE",
    "DEPTH", "TEMPERATURE", "SALINITY", "DISS_OXYGEN", "TURBIDITY"
  )
  
  missing_fields <- setdiff(required_fields, names(df))
  if (length(missing_fields) > 0) {
    stop(
      paste0("Missing required fields for ", gear_label, ":\n  - ",
             paste(missing_fields, collapse = "\n  - ")),
      call. = FALSE
    )
  }
  
  # --- timezone handling (UTC -> America/Chicago conversion) -----------------
  
  DATE_TIME_UTC <- df$COMPLETION_DTTM
  
  # Guardrail: ensure it is treated as UTC before conversion
  if (is.null(attr(DATE_TIME_UTC, "tzone")) ||
      !nzchar(attr(DATE_TIME_UTC, "tzone"))) {
    attr(DATE_TIME_UTC, "tzone") <- "UTC"
  }
  
  # Convert to Central time (DST-aware)
  DATE_TIME_CENTRAL <- lubridate::with_tz(
    DATE_TIME_UTC,
    tzone = "America/Chicago"
  )
  
  out <- df %>%
    mutate(
      GEAR = gear_label,
      MAJOR_AREA = to_int_safely(.data[[major_area_field]]),
      
      DATE_TIME = DATE_TIME_CENTRAL,                       # use precomputed vector
      DATE      = as.Date(DATE_TIME),
      YEAR_i    = as.integer(format(DATE_TIME, "%Y")),
      MONTH_i   = as.integer(format(DATE_TIME, "%m")),
      
      TEMP    = as.numeric(TEMPERATURE),
      SAL     = as.numeric(SALINITY),
      DEPTH_M = as.numeric(DEPTH),
      DO      = as.numeric(DISS_OXYGEN),
      TURB    = as.numeric(TURBIDITY),
      
      LAT = as.numeric(LATITUDE),
      LON = as.numeric(LONGITUDE)
    )
  
  # Hard stops: ID uniqueness + required derived fields present
  if (anyDuplicated(out$SAMPLE_ID) > 0) {
    dup_ids <- unique(out$SAMPLE_ID[duplicated(out$SAMPLE_ID)])
    stop(
      paste0("Duplicated SAMPLE_ID detected in ", gear_label, ". Example IDs:\n  - ",
             paste(head(dup_ids, 20), collapse = "\n  - ")),
      call. = FALSE
    )
  }
  
  if (any(is.na(out$DATE))) stop("Missing DATE after parsing COMPLETION_DTTM in ", gear_label, call. = FALSE)
  if (any(is.na(out$YEAR_i)) || any(is.na(out$MONTH_i))) stop("Missing YEAR_i or MONTH_i in ", gear_label, call. = FALSE)
  if (any(is.na(out$MAJOR_AREA))) stop("Missing MAJOR_AREA in ", gear_label, call. = FALSE)
  
  out
}

station_bs <- harmonize_station(STATION_BS, "BS")
station_bt <- harmonize_station(STATION_BT, "BT")

core_covars <- c("TEMP", "SAL", "DEPTH_M", "DO", "TURB")

station_bs_core <- station_bs %>%
  select(SAMPLE_ID, GEAR, MAJOR_AREA, DATE, YEAR_i, MONTH_i, LAT, LON, all_of(core_covars)) %>%
  filter(MAJOR_AREA %in% 1:8)

station_bt_core <- station_bt %>%
  select(SAMPLE_ID, GEAR, MAJOR_AREA, DATE, YEAR_i, MONTH_i, LAT, LON, all_of(core_covars)) %>%
  filter(MAJOR_AREA %in% 1:8)

stopifnot(!any(is.na(station_bs_core$SAMPLE_ID)))
stopifnot(!any(is.na(station_bt_core$SAMPLE_ID)))

cat("✅ Built station_bs_core and station_bt_core\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
diag_dir <- PATHS$diag_dir

missingness_by_year <- function(df_core) {
  df_core %>%
    group_by(GEAR, YEAR = YEAR_i) %>%
    summarise(
      n_samples = n(),
      across(all_of(core_covars),
             list(n_miss = ~sum(is.na(.x)),
                  prop_miss = ~mean(is.na(.x))),
             .names = "{.col}__{.fn}"),
      .groups = "drop"
    )
}

miss_all <- bind_rows(
  missingness_by_year(station_bs_core),
  missingness_by_year(station_bt_core)
) %>% arrange(GEAR, YEAR)

miss_csv <- file.path(diag_dir, "b3_missingness_by_year__station_core.csv")
readr::write_csv(miss_all, miss_csv)
cat("🧾 Missingness report saved:", miss_csv, "\n")

coverage_counts <- function(df_core) {
  df_core %>%
    mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR_i, MONTH_i)) %>%
    count(GEAR, YEAR_i, MONTH_i, YEAR_MONTH, name = "n_samples") %>%
    arrange(GEAR, YEAR_i, MONTH_i)
}

cov_all <- bind_rows(
  coverage_counts(station_bs_core),
  coverage_counts(station_bt_core)
)

p1 <- ggplot(cov_all, aes(x = as.Date(paste0(YEAR_MONTH, "-01")), y = n_samples)) +
  geom_line() +
  facet_wrap(~GEAR, scales = "free_y", ncol = 1) +
  labs(
    x = "Year-month",
    y = "Number of station samples",
    title = "Station metadata coverage over time (BS vs BT)"
  )

plot_path <- file.path(diag_dir, "b3_station_coverage_counts__year_month.png")
ggsave(plot_path, p1, width = 10, height = 7, dpi = 300)
cat("📈 Coverage plot saved:", plot_path, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Save derived data -> 01_Analysis/Data/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
out_dir <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir, "station_core.rds")
saveRDS(
  list(
    station_bs_core = station_bs_core,
    station_bt_core = station_bt_core,
    core_covars     = core_covars
  ),
  out_rds
)

cat("💾 Saved derived station core:", out_rds, "\n")
cat("✅ b3_harmonize_station_metadata.R completed successfully\n")