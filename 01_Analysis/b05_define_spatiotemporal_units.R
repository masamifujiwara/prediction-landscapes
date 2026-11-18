# =============================================================================
# Script:  b05_define_spatiotemporal_units.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Defines the units the diversity and covariate series are built on:
#     - alpha (estuary scale): BAY_ID x GEAR x YEAR x MONTH
#     - gamma (coastwide):     GEAR x YEAR x MONTH (all estuaries pooled)
#   and the complete monthly time grid from January of the first year to
#   December of the last year sampled. Counts the samples (seine hauls or
#   trawl tows) in each cell and flags cells with fewer than N_min samples
#   (report only; nothing is removed). Records, per estuary and sampling method
#   and coastwide per method, which months of the grid have no samples.
#
# Run after: b03.
#
# Inputs:
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#
# Outputs:
#   01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds
#     keys (key names, valid_bays), time_grid, lookups (coverage per estuary
#     and sampling method), diagnostics; read by b06, b07, b08, b09
#   02_Results/03_Diagnostics/b5_define_spatiotemporal_units/
#     b5_counts_alpha_cells__bay_gear_year_month.csv
#     b5_sparse_alpha_cells__lt_<N_min>_samples.csv
#     b5_counts_gamma_cells__gear_year_month.csv
#     b5_coverage_by_bay_gear.csv
#   Run log: 01_Analysis/Logs/b5_define_spatiotemporal_units/run_<timestamp>.txt
#   Earlier outputs in the data and diagnostics folders are moved to
#   ARCHIVE_<timestamp> first. Output folders use the stem
#   b5_define_spatiotemporal_units.
#
# Method notes:
#   - valid_bays is taken from the estuaries present in station_core, not typed
#     here; a00 decides which estuaries are kept and b03 checks them.
#   - Months with no samples are split into those before the first sample
#     (leading), inside the record (interior) and after the last sample
#     (trailing), per estuary and sampling method and coastwide per method
#     (lookups$coverage and b5_coverage_by_bay_gear.csv).
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
script_stem <- "b5_define_spatiotemporal_units"
in_station  <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")

N_min <- 5  # a bay x gear x month cell with fewer samples than this is flagged (report only)

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
# 4) Load station_core (from b03)
# -----------------------------------------------------------------------------
stop_if_missing(in_station, "b03_harmonize_station_metadata.R")
station_core <- readRDS(in_station)
for (el in c("station_bs_core", "station_bt_core")) {
  if (!el %in% names(station_core)) stop("station_core.rds has no element ", el, call. = FALSE)
}

station_all <- bind_rows(
  station_core$station_bs_core,
  station_core$station_bt_core
)

# Standard key names and integer types
station_all <- station_all %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  )

for (k in c("BAY_ID", "YEAR", "MONTH")) {
  if (anyNA(station_all[[k]])) stop("Missing ", k, " in station_core", call. = FALSE)
}
if (any(station_all$MONTH < 1L | station_all$MONTH > 12L)) stop("MONTH outside 1 to 12", call. = FALSE)
if (!setequal(unique(station_all$GEAR), c("BS", "BT"))) {
  stop("Expected gears BS and BT; found ", paste(sort(unique(station_all$GEAR)), collapse = ", "), call. = FALSE)
}

# The bays are the ones present in the data. a00 decides which bays are kept
# and b03 checks them, so no bay list is typed here.
valid_bays <- sort(unique(station_all$BAY_ID))
cat("Bays:", paste(valid_bays, collapse = " "), "\n")

# -----------------------------------------------------------------------------
# 5) Time grid (YEAR x MONTH) over the observed years
# -----------------------------------------------------------------------------
year_min <- min(station_all$YEAR, na.rm = TRUE)
year_max <- max(station_all$YEAR, na.rm = TRUE)

time_grid <- tidyr::expand_grid(
  YEAR = year_min:year_max,
  MONTH = 1:12
) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

# -----------------------------------------------------------------------------
# 6) Coverage: which grid months have samples, per bay and gear and coastwide
# -----------------------------------------------------------------------------
# For each series, the months of time_grid with no sample are split into those
# before the first sample (leading), inside the record (interior) and after
# the last sample (trailing). Interior months are listed by name.
grid_ym <- time_grid$YEAR_MONTH

coverage_one <- function(ym_sampled) {
  have <- grid_ym %in% ym_sampled
  first <- min(which(have)); last <- max(which(have))
  interior_missing <- grid_ym[first:last][!have[first:last]]
  tibble(
    first_year_month          = grid_ym[first],
    last_year_month           = grid_ym[last],
    n_months_grid             = length(grid_ym),
    n_months_sampled          = sum(have),
    n_months_missing_leading  = first - 1L,
    n_months_missing_interior = length(interior_missing),
    n_months_missing_trailing = length(grid_ym) - last,
    missing_interior_months   = paste(interior_missing, collapse = "; ")
  )
}

sampled_months <- station_all %>%
  distinct(GEAR, BAY_ID, YEAR, MONTH) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

coverage <- bind_rows(
  sampled_months %>%
    group_by(scale = "alpha", GEAR, BAY_ID) %>%
    group_modify(~ coverage_one(.x$YEAR_MONTH)) %>%
    ungroup(),
  sampled_months %>%
    group_by(scale = "gamma", GEAR) %>%
    group_modify(~ coverage_one(.x$YEAR_MONTH)) %>%
    ungroup() %>%
    mutate(BAY_ID = NA_integer_)
) %>%
  select(scale, GEAR, BAY_ID, everything()) %>%
  arrange(scale, GEAR, BAY_ID)

# -----------------------------------------------------------------------------
# 7) Samples per alpha cell: BAY_ID x GEAR x YEAR x MONTH
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 8) Samples per gamma cell: GEAR x YEAR x MONTH (coastwide)
# -----------------------------------------------------------------------------
gamma_counts <- station_all %>%
  count(GEAR, YEAR, MONTH, name = "n_samples") %>%
  arrange(GEAR, YEAR, MONTH)

# -----------------------------------------------------------------------------
# 9) Diagnostics -> 02_Results/03_Diagnostics/<script_stem>/
# -----------------------------------------------------------------------------
alpha_counts_path <- file.path(dirs$diag, "b5_counts_alpha_cells__bay_gear_year_month.csv")
alpha_sparse_path <- file.path(dirs$diag, paste0("b5_sparse_alpha_cells__lt_", N_min, "_samples.csv"))
gamma_counts_path <- file.path(dirs$diag, "b5_counts_gamma_cells__gear_year_month.csv")
coverage_path     <- file.path(dirs$diag, "b5_coverage_by_bay_gear.csv")

readr::write_csv(alpha_counts, alpha_counts_path)
readr::write_csv(alpha_sparse, alpha_sparse_path)
readr::write_csv(gamma_counts, gamma_counts_path)
readr::write_csv(coverage, coverage_path, na = "")

cat("Saved:", alpha_counts_path, "\n")
cat("Saved:", alpha_sparse_path, "\n")
cat("Saved:", gamma_counts_path, "\n")
cat("Saved:", coverage_path, "\n")

cat("\n=== Months with no samples (grid of", length(grid_ym), "months) ===\n")
print(as.data.frame(coverage %>% select(scale, GEAR, BAY_ID, first_year_month, n_months_missing_leading,
                                        n_months_missing_interior, missing_interior_months,
                                        n_months_missing_trailing)), row.names = FALSE)

cat("\n=== Samples per alpha cell (BAY_ID x GEAR x YEAR x MONTH) ===\n")
print(as.data.frame(alpha_dist), row.names = FALSE)
cat("Cells with fewer than", N_min, "samples:", nrow(alpha_sparse), "\n")

# -----------------------------------------------------------------------------
# 10) Save index_keys -> 01_Analysis/Data/<script_stem>/
# -----------------------------------------------------------------------------
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
    coverage = coverage
  ),
  diagnostics = list(
    alpha_counts = alpha_counts,
    gamma_counts = gamma_counts,
    alpha_dist   = alpha_dist,
    N_min_sparse_flag = N_min
  )
)

out_rds <- file.path(dirs$data, "index_keys.rds")
saveRDS(index_keys, out_rds)
cat("Saved:", out_rds, "\n")

end_log("finished")
cat("b05_define_spatiotemporal_units.R completed\n")
