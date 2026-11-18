# =============================================================================
# Script:  b04_join_organisms_to_station.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Joins the organism records to the station table by SAMPLE_ID, separately
#   for bag seine (BS) and bay trawl (BT), and keeps samples with no catch:
#     - one row per species recorded in a sample;
#     - plus one placeholder row (SPECIES_CODE = NA, catch = 0) for each
#       sample with no organism record.
#   Placeholder rows must be left out of the abundance vectors for diversity
#   (their SPECIES_CODE is NA).
#
# Run after: a00, b03.
#
# Inputs:
#   01_Analysis/Data.RData: bagseine_organisms, baytrawl_organisms, species (a00)
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#
# Outputs:
#   01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds, bt_events.rds
#     station columns from b03, plus SPECIES_CODE (text) and catch; read by
#     b06, b07, f07, h11
#   02_Results/03_Diagnostics/b4_join_organisms_to_station/
#     b4_join_diagnostics__BS_BT.csv
#     b4_zero_catch_rows_removed.csv
#     b4_unmatched_org_SAMPLE_ID__BS.csv, b4_unmatched_org_SAMPLE_ID__BT.csv
#   Run log: 01_Analysis/Logs/b4_join_organisms_to_station/run_<timestamp>.txt
#   Earlier outputs in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> first. Output folders use the stem
#   b4_join_organisms_to_station.
#
# Method notes:
#   - Organism tables list recorded species only, so a species absent from a
#     sample means a catch of 0.
#   - Organism rows recorded with a catch of 0 are treated as entry errors and
#     removed before the join. A recorded row then always means a catch above
#     0, and a sample whose rows were all 0 gets a placeholder row like any
#     other empty sample. Hill numbers and composition do not change, because
#     the later scripts use only positive catches. The removed rows are listed
#     in b4_zero_catch_rows_removed.csv.
#   - A missing catch, or an organism SAMPLE_ID with no station row, stops the
#     script: a00 guarantees that neither occurs, so either means the inputs
#     have changed.
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

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
# Stem used for the output folders.
script_stem <- "b4_join_organisms_to_station"
in_station  <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")

options(scipen = 999, warn = 1)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "readr")
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
# 4) Load inputs
# -----------------------------------------------------------------------------
stop_if_missing(c(PATHS$clean_data_rdata, in_station),
                c("a00_getdata.R", "b03_harmonize_station_metadata.R"))
load(PATHS$clean_data_rdata)
for (obj in c("bagseine_organisms", "baytrawl_organisms", "species")) {
  if (!exists(obj)) stop("Data.RData has no object ", obj, " (rerun a00_getdata.R)", call. = FALSE)
}

station_core <- readRDS(in_station)
for (el in c("station_bs_core", "station_bt_core")) {
  if (!el %in% names(station_core)) stop("station_core.rds has no element ", el, call. = FALSE)
}
station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

# -----------------------------------------------------------------------------
# 5) Check the organism tables and set column types
# -----------------------------------------------------------------------------
validate_org_table <- function(org_df, label) {
  expected <- c("SAMPLE_ID", "SPECIES_CODE", "catch")
  if (!identical(sort(names(org_df)), sort(expected))) {
    stop(label, " has unexpected columns: ", paste(names(org_df), collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(org_df$SAMPLE_ID)) stop(label, ": SAMPLE_ID is not numeric", call. = FALSE)
  if (!(is.character(org_df$SPECIES_CODE) || is.factor(org_df$SPECIES_CODE))) {
    stop(label, ": SPECIES_CODE is not text", call. = FALSE)
  }
  if (!is.numeric(org_df$catch)) stop(label, ": catch is not numeric", call. = FALSE)
  # a00 removes the presence-only codes, which were the only rows without a
  # catch, so a missing catch here means the input has changed.
  if (anyNA(org_df$catch)) {
    stop(label, ": ", sum(is.na(org_df$catch)), " rows with no catch (expected none after a00)", call. = FALSE)
  }
}
validate_org_table(bagseine_organisms, "bagseine_organisms")
validate_org_table(baytrawl_organisms, "baytrawl_organisms")

bagseine_organisms <- bagseine_organisms %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE), catch = as.numeric(catch))
baytrawl_organisms <- baytrawl_organisms %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE), catch = as.numeric(catch))

# A row with a catch of 0 is treated as an entry error (see header). Keep a
# list of the removed rows, then drop them so that every organism row means a
# catch above 0.
zero_rows <- bind_rows(
  bagseine_organisms %>% filter(catch == 0) %>% mutate(GEAR = "BS"),
  baytrawl_organisms %>% filter(catch == 0) %>% mutate(GEAR = "BT")
)
n_zero <- c(BS = sum(zero_rows$GEAR == "BS"), BT = sum(zero_rows$GEAR == "BT"))
bagseine_organisms <- bagseine_organisms %>% filter(catch != 0)
baytrawl_organisms <- baytrawl_organisms %>% filter(catch != 0)
cat("Organism rows with a catch of 0 removed: BS", n_zero[["BS"]], " BT", n_zero[["BT"]], "\n")

# -----------------------------------------------------------------------------
# 6) Join stations to organisms (one placeholder row per zero-catch sample)
# -----------------------------------------------------------------------------
join_keep_zeros <- function(station_df, org_df, gear_label) {

  if (!"SAMPLE_ID" %in% names(station_df)) stop(gear_label, ": station table has no SAMPLE_ID", call. = FALSE)
  if (anyDuplicated(station_df$SAMPLE_ID)) stop(gear_label, ": duplicated SAMPLE_ID in station table", call. = FALSE)

  # Station-left join keeps every station sample. A sample with no organism
  # record gets one placeholder row: SPECIES_CODE = NA, catch = 0.
  events <- station_df %>%
    left_join(org_df, by = "SAMPLE_ID") %>%
    mutate(
      GEAR = gear_label,
      catch = dplyr::if_else(is.na(catch), 0, as.numeric(catch)),
      SPECIES_CODE = as.character(SPECIES_CODE)  # NA stays NA
    )

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

  # Organism sample IDs with no station row. a00 keeps only organisms whose
  # sample is in the station table, so any here means the inputs disagree.
  unmatched_org_ids <- org_df %>%
    anti_join(station_df %>% select(SAMPLE_ID), by = "SAMPLE_ID") %>%
    distinct(SAMPLE_ID) %>%
    arrange(SAMPLE_ID)

  # Placeholder integrity: exactly one placeholder per zero-catch sample.
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

  n_neg <- sum(events$catch < 0, na.rm = TRUE)
  q999  <- suppressWarnings(stats::quantile(events$catch, probs = 0.999, na.rm = TRUE))
  n_hi  <- sum(events$catch > q999, na.rm = TRUE)

  if (nrow(bad_sentinel) > 0) {
    stop("Placeholder row check failed for ", gear_label, ": ",
         nrow(bad_sentinel), " samples do not have exactly 1 placeholder row.", call. = FALSE)
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

# -----------------------------------------------------------------------------
# 7) Diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# -----------------------------------------------------------------------------
diag_row <- function(res, gear_label) {
  tibble(
    GEAR = gear_label,
    n_station_samples = res$diag$n_station_samples,
    n_zero_catch_rows_removed = unname(n_zero[[gear_label]]),
    n_event_rows = res$diag$n_event_rows,
    n_samples_with_any_record = res$diag$station_with_any$n_with_any,
    n_samples_without_any_record = res$diag$station_with_any$n_without_any,
    pct_samples_without_any_record = res$diag$station_with_any$pct_without_any,
    n_unmatched_org_sample_ids = res$diag$n_unmatched_org_sample_ids,
    n_negative_catch = res$diag$n_negative_catch,
    catch_q999 = res$diag$catch_q999,
    n_above_q999 = res$diag$n_above_q999
  )
}
join_diag <- bind_rows(diag_row(bs_res, "BS"), diag_row(bt_res, "BT"))

readr::write_csv(join_diag, file.path(dirs$diag, "b4_join_diagnostics__BS_BT.csv"))

# Removed zero-catch rows, with the sample's bay, year and month and the
# species names from the species list (NA for codes not in the list).
sp_names <- species %>%
  transmute(SPECIES_CODE = as.character(.data[[names(species)[1]]]),
            latin_name = .data[[names(species)[2]]],
            common_name = .data[[names(species)[3]]]) %>%
  distinct(SPECIES_CODE, .keep_all = TRUE)
zero_list <- zero_rows %>%
  left_join(bind_rows(station_bs_core, station_bt_core) %>%
              select(SAMPLE_ID, GEAR, MAJOR_AREA, YEAR_i, MONTH_i),
            by = c("SAMPLE_ID", "GEAR")) %>%
  left_join(sp_names, by = "SPECIES_CODE") %>%
  select(GEAR, SAMPLE_ID, MAJOR_AREA, YEAR_i, MONTH_i, SPECIES_CODE, common_name, latin_name, catch) %>%
  arrange(GEAR, YEAR_i, MONTH_i, SAMPLE_ID, SPECIES_CODE)
readr::write_csv(zero_list, file.path(dirs$diag, "b4_zero_catch_rows_removed.csv"))
readr::write_csv(bs_res$diag$unmatched_org_sample_ids, file.path(dirs$diag, "b4_unmatched_org_SAMPLE_ID__BS.csv"))
readr::write_csv(bt_res$diag$unmatched_org_sample_ids, file.path(dirs$diag, "b4_unmatched_org_SAMPLE_ID__BT.csv"))
print(as.data.frame(join_diag), row.names = FALSE)
cat("Diagnostics saved under:", dirs$diag, "\n")

# -----------------------------------------------------------------------------
# 8) Checks before writing
# -----------------------------------------------------------------------------
if (any(join_diag$n_unmatched_org_sample_ids > 0)) {
  stop("Organism records with no station sample (see b4_unmatched_org_SAMPLE_ID__*.csv). ",
       "a00 should have removed them.", call. = FALSE)
}
if (any(join_diag$n_negative_catch > 0)) stop("Negative catch values found", call. = FALSE)
for (ev in list(bs_events, bt_events)) {
  if (nrow(ev) == 0) stop("An events table is empty", call. = FALSE)
  if (anyNA(ev$SAMPLE_ID) || anyNA(ev$MAJOR_AREA) || anyNA(ev$YEAR_i) || anyNA(ev$MONTH_i)) {
    stop("Missing key values (SAMPLE_ID, MAJOR_AREA, YEAR_i or MONTH_i) in the events", call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# 9) Save derived outputs -> 01_Analysis/Data/<script_stem>/
# -----------------------------------------------------------------------------
bs_out <- file.path(dirs$data, "bs_events.rds")
bt_out <- file.path(dirs$data, "bt_events.rds")
saveRDS(bs_events, bs_out)
saveRDS(bt_events, bt_out)
cat("Saved:", bs_out, "\n")
cat("Saved:", bt_out, "\n")

cat("\nNOTE: Zero-catch samples are kept as placeholder rows (SPECIES_CODE = NA, catch = 0).\n")
cat("      Leave placeholder rows out when building abundance vectors for diversity.\n")

end_log("finished")
cat("b04_join_organisms_to_station.R completed\n")
