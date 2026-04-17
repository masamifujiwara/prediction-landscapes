# =============================================================================
# Script: b4_join_organisms_to_station.R
#
# Purpose (project directory convention):
#   - Join organisms ↔ station_core by SAMPLE_ID for BS and BT separately.
#   - Create analysis-ready “event-level” tables that RETAIN zero-catch stations:
#       * One row per observed species per sample (observed-only organisms table)
#       * PLUS one sentinel row per station sample with no organism rows:
#           SPECIES_CODE = NA, catch = 0
#
# Notes:
#   - Organism tables are observed-only: missing species in a sample imply catch = 0.
#   - We retain true zero-catch samples (no organism records) via sentinel rows.
#   - Sentinel rows must be excluded from diversity vectors (SPECIES_CODE is NA).
#
# Inputs:
#   - 01_Analysis/Data.RData
#       * bagseine_organisms, baytrawl_organisms
#   - 01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds
#
# Outputs (derived data):
#   - 01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds
#   - 01_Analysis/Data/b4_join_organisms_to_station/bt_events.rds
#
# Diagnostics:
#   - 02_Results/03_Diagnostics/b4_join_organisms_to_station/
#       * b4_join_diagnostics__BS_BT.csv
#       * b4_unmatched_org_SAMPLE_ID__BS.csv
#       * b4_unmatched_org_SAMPLE_ID__BT.csv
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
required_pkgs <- c("dplyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) Load inputs
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(PATHS$clean_data_rdata))
load(PATHS$clean_data_rdata)

stopifnot(exists("bagseine_organisms"), exists("baytrawl_organisms"))

station_core_path <- file.path(
  PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds"
)
stopifnot(file.exists(station_core_path))
station_core <- readRDS(station_core_path)

stopifnot(all(c("station_bs_core", "station_bt_core") %in% names(station_core)))
station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

# Optional guarantee (already done in b3): keep bays 1–8
station_bs_core <- station_bs_core %>% filter(MAJOR_AREA %in% 1:8)
station_bt_core <- station_bt_core %>% filter(MAJOR_AREA %in% 1:8)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Validate organism tables + standardize catch / species_code types
# ─────────────────────────────────────────────────────────────────────────────
validate_org_table <- function(org_df, label) {
  expected <- c("SAMPLE_ID", "SPECIES_CODE", "catch")
  if (!identical(sort(names(org_df)), sort(expected))) {
    stop(label, " has unexpected columns: ", paste(names(org_df), collapse = ", "), call. = FALSE)
  }
  stopifnot(is.numeric(org_df$SAMPLE_ID))
  stopifnot(is.character(org_df$SPECIES_CODE) || is.factor(org_df$SPECIES_CODE))
  stopifnot(is.numeric(org_df$catch) || is.integer(org_df$catch))
}

validate_org_table(bagseine_organisms, "bagseine_organisms")
validate_org_table(baytrawl_organisms, "baytrawl_organisms")

# Treat NA catch as 0 and ensure SPECIES_CODE is character
bagseine_organisms <- bagseine_organisms %>%
  mutate(
    SPECIES_CODE = as.character(SPECIES_CODE),
    catch = dplyr::if_else(is.na(catch), 0, as.numeric(catch))
  )

baytrawl_organisms <- baytrawl_organisms %>%
  mutate(
    SPECIES_CODE = as.character(SPECIES_CODE),
    catch = dplyr::if_else(is.na(catch), 0, as.numeric(catch))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Join stations -> organisms (retain sentinel rows for zero-catch stations)
# ─────────────────────────────────────────────────────────────────────────────
join_keep_zeros <- function(station_df, org_df, gear_label) {
  
  stopifnot("SAMPLE_ID" %in% names(station_df))
  stopifnot(!anyDuplicated(station_df$SAMPLE_ID))
  
  # Station-left join retains all station samples.
  # Samples with no organism records become a single sentinel row:
  # SPECIES_CODE = NA, catch = 0
  events <- station_df %>%
    left_join(org_df, by = "SAMPLE_ID") %>%
    mutate(
      GEAR = gear_label,
      catch = dplyr::if_else(is.na(catch), 0, as.numeric(catch)),
      SPECIES_CODE = as.character(SPECIES_CODE)  # NA stays NA
    )
  
  # Diagnostics ---------------------------------------------------------------
  n_station <- nrow(station_df)
  n_events  <- nrow(events)
  
  station_with_any <- events %>%
    group_by(SAMPLE_ID) %>%
    summarise(has_any_record = any(!is.na(SPECIES_CODE)), .groups = "drop") %>%
    summarise(
      n_with_any = sum(has_any_record),
      n_without_any = sum(!has_any_record),
      pct_without_any = 100 * mean(!has_any_record),
      .groups = "drop"
    )
  
  # Organism sample IDs that do not exist in station_core
  unmatched_org_ids <- org_df %>%
    anti_join(station_df %>% select(SAMPLE_ID), by = "SAMPLE_ID") %>%
    distinct(SAMPLE_ID) %>%
    arrange(SAMPLE_ID)
  
  # Sentinel row integrity: exactly 1 sentinel per zero-catch sample
  sentinel_counts <- events %>%
    group_by(SAMPLE_ID) %>%
    summarise(
      n_rows = n(),
      n_non_na_species = sum(!is.na(SPECIES_CODE)),
      n_sentinel = sum(is.na(SPECIES_CODE) & catch == 0),
      .groups = "drop"
    ) %>%
    mutate(is_zero_catch_sample = (n_non_na_species == 0))
  
  bad_sentinel <- sentinel_counts %>%
    filter(is_zero_catch_sample & n_sentinel != 1)
  
  # Catch sanity
  n_neg <- sum(events$catch < 0, na.rm = TRUE)
  q999  <- suppressWarnings(stats::quantile(events$catch, probs = 0.999, na.rm = TRUE))
  n_hi  <- sum(events$catch > q999, na.rm = TRUE)
  
  if (nrow(bad_sentinel) > 0) {
    stop(
      "Sentinel row integrity failed for ", gear_label, ": ",
      nrow(bad_sentinel), " samples do not have exactly 1 sentinel row.",
      call. = FALSE
    )
  }
  
  list(
    events = events,
    diag = list(
      n_station_samples = n_station,
      n_event_rows = n_events,
      station_with_any = station_with_any,
      n_unmatched_org_sample_ids = nrow(unmatched_org_ids),
      unmatched_org_sample_ids = unmatched_org_ids,
      n_negative_catch = n_neg,
      catch_q999 = as.numeric(q999),
      n_above_q999 = n_hi
    )
  )
}

bs_res <- join_keep_zeros(station_bs_core, bagseine_organisms, "BS")
bt_res <- join_keep_zeros(station_bt_core, baytrawl_organisms, "BT")

bs_events <- bs_res$events
bt_events <- bt_res$events

# ─────────────────────────────────────────────────────────────────────────────
# 6) Save derived outputs -> 01_Analysis/Data/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
out_dir <- file.path(PATHS$data_dir, "b4_join_organisms_to_station")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bs_out <- file.path(out_dir, "bs_events.rds")
bt_out <- file.path(out_dir, "bt_events.rds")

saveRDS(bs_events, bs_out)
saveRDS(bt_events, bt_out)

cat("💾 Saved: ", bs_out, "\n", sep = "")
cat("💾 Saved: ", bt_out, "\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Save diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
diag_dir <- PATHS$diag_dir

join_diag <- bind_rows(
  tibble(
    gear = "BS",
    n_station_samples = bs_res$diag$n_station_samples,
    n_event_rows = bs_res$diag$n_event_rows,
    n_samples_with_any_record = bs_res$diag$station_with_any$n_with_any,
    n_samples_without_any_record = bs_res$diag$station_with_any$n_without_any,
    pct_samples_without_any_record = bs_res$diag$station_with_any$pct_without_any,
    n_unmatched_org_sample_ids = bs_res$diag$n_unmatched_org_sample_ids,
    n_negative_catch = bs_res$diag$n_negative_catch,
    catch_q999 = bs_res$diag$catch_q999,
    n_above_q999 = bs_res$diag$n_above_q999
  ),
  tibble(
    gear = "BT",
    n_station_samples = bt_res$diag$n_station_samples,
    n_event_rows = bt_res$diag$n_event_rows,
    n_samples_with_any_record = bt_res$diag$station_with_any$n_with_any,
    n_samples_without_any_record = bt_res$diag$station_with_any$n_without_any,
    pct_samples_without_any_record = bt_res$diag$station_with_any$pct_without_any,
    n_unmatched_org_sample_ids = bt_res$diag$n_unmatched_org_sample_ids,
    n_negative_catch = bt_res$diag$n_negative_catch,
    catch_q999 = bt_res$diag$catch_q999,
    n_above_q999 = bt_res$diag$n_above_q999
  )
)

readr::write_csv(join_diag, file.path(diag_dir, "b4_join_diagnostics__BS_BT.csv"))
readr::write_csv(bs_res$diag$unmatched_org_sample_ids, file.path(diag_dir, "b4_unmatched_org_SAMPLE_ID__BS.csv"))
readr::write_csv(bt_res$diag$unmatched_org_sample_ids, file.path(diag_dir, "b4_unmatched_org_SAMPLE_ID__BT.csv"))

print(join_diag)
cat("🧾 Diagnostics saved under: ", diag_dir, "\n", sep = "")

cat("\nNOTE: Zero-catch stations are retained as sentinel rows (SPECIES_CODE = NA, catch = 0).\n")
cat("      Exclude sentinel rows when building abundance vectors for diversity.\n\n")

cat("✅ b4_join_organisms_to_station.R completed successfully\n")