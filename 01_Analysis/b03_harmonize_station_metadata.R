# =============================================================================
# Script:  b03_harmonize_station_metadata.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds the harmonized station ("core") tables for bag seine (BS) and bay
#   trawl (BT): one row per sample with the integer bay ID (MAJOR_AREA), year
#   (YEAR_i) and month (MONTH_i), and the two covariates that later scripts use,
#   water temperature (TEMP) and salinity (SAL). Other station columns (date,
#   position, depth, dissolved oxygen, turbidity) are not used downstream and
#   are not kept. Unusual values are kept as recorded and listed in a
#   diagnostics file.
#
# Run after: a00 (Data.RData), b01.
#
# Inputs:
#   01_Analysis/Data.RData: STATION_BS, STATION_BT (a00; estuaries 2 to 8 only)
#
# Outputs:
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds
#     list(station_bs_core, station_bt_core, core_covars); read by b04 to b08,
#     e02, f07, h12
#   02_Results/03_Diagnostics/b3_harmonize_station_metadata/
#     b3_missingness_by_year__station_core.csv
#     b3_unusual_values__station_core.csv
#     b3_station_coverage_counts__year_month.png
#   Run log: 01_Analysis/Logs/b3_harmonize_station_metadata/run_<timestamp>.txt
#   Earlier outputs in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> first. Output folders use the stem
#   b3_harmonize_station_metadata.
#
# Method notes:
#   - YEAR_i and MONTH_i come from the raw YEAR and MONTH columns. The stored
#     completion times look like local clock times, so they are not converted
#     from UTC; converting them would not change any sample's year or month.
#   - Salinity above 80 and water temperature of exactly 0 are listed as
#     unusual but kept: they were reviewed and judged to be real readings
#     (hypersaline lagoon conditions and freezes).
#   - The script does not filter estuaries; it only checks that each sampling
#     method has the estuaries in EXPECTED_BAYS, which a00 decides.
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
script_stem <- "b3_harmonize_station_metadata"

# Bays each gear must have. a00 decides which bays are kept; this is only a
# check, so a change there is noticed here instead of passing on silently.
EXPECTED_BAYS <- 2:8

# Ranges for the unusual-values list (report only; no value is changed).
# A value is listed if it is below the lower bound, above the upper bound, or
# equal to a bound marked as listed (water temperature of exactly 0).
UNUSUAL <- list(
  TEMP = list(lower = 0, upper = 40, list_at_lower = TRUE),
  SAL  = list(lower = 0, upper = 80, list_at_lower = FALSE)
)

options(scipen = 999, warn = 1)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

# -----------------------------------------------------------------------------
# 4) Load Data.RData
# -----------------------------------------------------------------------------
stop_if_missing(PATHS$clean_data_rdata, "a00_getdata.R")
load(PATHS$clean_data_rdata)
for (obj in c("STATION_BS", "STATION_BT")) {
  if (!exists(obj)) stop("Data.RData has no object ", obj, " (rerun a00_getdata.R)", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 5) Harmonization
# -----------------------------------------------------------------------------
to_int_safely <- function(x) suppressWarnings(as.integer(x))

harmonize_station <- function(df, gear_label) {

  major_area_field <- if (gear_label == "BS") "MAJOR_AREA_CODE" else "MAJOR_AREA"

  required_fields <- c("SAMPLE_ID", major_area_field, "YEAR", "MONTH",
                       "TEMPERATURE", "SALINITY")
  missing_fields <- setdiff(required_fields, names(df))
  if (length(missing_fields) > 0) {
    stop("Missing required fields for ", gear_label, ":\n  - ",
         paste(missing_fields, collapse = "\n  - "), call. = FALSE)
  }

  # YEAR and MONTH are stored as text ("1982", "01") in the raw tables.
  out <- df %>%
    mutate(
      GEAR       = gear_label,
      MAJOR_AREA = to_int_safely(.data[[major_area_field]]),
      YEAR_i     = to_int_safely(YEAR),
      MONTH_i    = to_int_safely(MONTH),
      TEMP       = as.numeric(TEMPERATURE),
      SAL        = as.numeric(SALINITY)
    )

  # Hard stops: unique sample IDs and no missing keys.
  if (anyNA(out$SAMPLE_ID)) stop("Missing SAMPLE_ID in ", gear_label, call. = FALSE)
  if (anyDuplicated(out$SAMPLE_ID) > 0) {
    dup_ids <- unique(out$SAMPLE_ID[duplicated(out$SAMPLE_ID)])
    stop("Duplicated SAMPLE_ID in ", gear_label, ". Example IDs:\n  - ",
         paste(head(dup_ids, 20), collapse = "\n  - "), call. = FALSE)
  }
  if (anyNA(out$YEAR_i) || anyNA(out$MONTH_i)) stop("Missing YEAR_i or MONTH_i in ", gear_label, call. = FALSE)
  if (anyNA(out$MAJOR_AREA)) stop("Missing or non-integer MAJOR_AREA in ", gear_label, call. = FALSE)

  out
}

station_bs <- harmonize_station(STATION_BS, "BS")
station_bt <- harmonize_station(STATION_BT, "BT")

core_covars <- c("TEMP", "SAL")
core_cols   <- c("SAMPLE_ID", "GEAR", "MAJOR_AREA", "YEAR_i", "MONTH_i", core_covars)

station_bs_core <- station_bs %>% select(all_of(core_cols))
station_bt_core <- station_bt %>% select(all_of(core_cols))

# -----------------------------------------------------------------------------
# 6) Checks before writing
# -----------------------------------------------------------------------------
for (g in list(list("BS", station_bs_core), list("BT", station_bt_core))) {
  d <- g[[2]]
  if (nrow(d) == 0) stop("No ", g[[1]], " samples", call. = FALSE)
  bays <- sort(unique(d$MAJOR_AREA))
  if (!identical(bays, as.integer(EXPECTED_BAYS))) {
    stop(g[[1]], " has bays ", paste(bays, collapse = " "), "; expected ",
         paste(EXPECTED_BAYS, collapse = " "), call. = FALSE)
  }
  if (any(d$MONTH_i < 1L | d$MONTH_i > 12L)) stop("MONTH_i outside 1 to 12 in ", g[[1]], call. = FALSE)
  cat(g[[1]], ": ", nrow(d), " samples, bays ", paste(bays, collapse = " "),
      ", years ", min(d$YEAR_i), " to ", max(d$YEAR_i), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# 7) Diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# -----------------------------------------------------------------------------
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

miss_csv <- file.path(dirs$diag, "b3_missingness_by_year__station_core.csv")
readr::write_csv(miss_all, miss_csv)
cat("Missingness report saved:", miss_csv, "\n")

# Unusual values, kept as recorded (see header).
unusual <- bind_rows(lapply(list(station_bs_core, station_bt_core), function(d) {
  bind_rows(lapply(names(UNUSUAL), function(v) {
    r <- UNUSUAL[[v]]
    x <- d[[v]]
    hit <- !is.na(x) & (x < r$lower | x > r$upper | (r$list_at_lower & x == r$lower))
    d[hit, c("GEAR", "SAMPLE_ID", "MAJOR_AREA", "YEAR_i", "MONTH_i")] %>%
      mutate(variable = v, value = x[hit])
  }))
})) %>% arrange(variable, GEAR, YEAR_i, MONTH_i, SAMPLE_ID)

unusual_csv <- file.path(dirs$diag, "b3_unusual_values__station_core.csv")
readr::write_csv(unusual, unusual_csv)
cat("Unusual values (kept as recorded):", nrow(unusual), "rows, saved:", unusual_csv, "\n")
print(as.data.frame(count(unusual, variable, GEAR)), row.names = FALSE)

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

plot_path <- file.path(dirs$diag, "b3_station_coverage_counts__year_month.png")
ggsave(plot_path, p1, width = 10, height = 7, dpi = 300)
cat("Coverage plot saved:", plot_path, "\n")

# -----------------------------------------------------------------------------
# 8) Save derived data -> 01_Analysis/Data/<script_stem>/
# -----------------------------------------------------------------------------
out_rds <- file.path(dirs$data, "station_core.rds")
saveRDS(
  list(
    station_bs_core = station_bs_core,
    station_bt_core = station_bt_core,
    core_covars     = core_covars
  ),
  out_rds
)
cat("Saved station core:", out_rds, "\n")

end_log("finished")
cat("b03_harmonize_station_metadata.R completed\n")
