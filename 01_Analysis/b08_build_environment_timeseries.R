# =============================================================================
# Script:  b08_build_environment_timeseries.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds the monthly environmental covariate series used by the covariate
#   models: monthly mean water temperature (TEMP) and salinity (SAL) per
#   estuary and sampling method (bay level, alpha) and per sampling method
#   over the whole coast (coastwide, gamma), on the monthly grid from b05,
#   with lagged copies (1, 3, 6 and 12 months). Also writes a table of the
#   months dropped and filled in each series and the correlation between the
#   covariates.
#
# Run after: b03, b05
#
# Inputs:
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds (b05):
#     keys$valid_bays and time_grid
#
# Outputs:
#   01_Analysis/Data/b8_build_environment_timeseries/
#     env_alpha_index.rds, env_alpha_index.csv  read by c09, c10, e07, e08,
#                                               f04, g04, h13
#     env_gamma_index.rds, env_gamma_index.csv  read by c09, c10, e07, e08,
#                                               f04, h13
#   02_Results/03_Diagnostics/b8_build_environment_timeseries/
#     env_missingness.csv  per series: first and last month, months dropped
#                          at the start and end, months filled inside
#     cor_matrix__alpha_index__means.csv, cor_matrix__gamma_index__means.csv
#   01_Analysis/Logs/b8_build_environment_timeseries/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   - Only TEMP and SAL are built, the covariates the models use.
#   - Coastwide means are taken over all samples of the sampling method in
#     the month (effort weighted, as the coastwide diversity in b07). A
#     coastwide month counts as sampled only when every estuary was sampled,
#     as in b07; the script stops if a month inside the record misses an
#     estuary.
#   - Each series starts at its first sampled month and ends at its last.
#     Months before or after are not filled, so no value is carried backward
#     into the start of a series and from there into the 12-month lags.
#   - Months inside a series with no samples (April and May 2020, no sampling
#     anywhere because of COVID) are filled by linear interpolation between
#     the neighboring months and marked filled = TRUE. The fill uses June
#     2020, a small look-ahead that matters only for forecasts made at those
#     two forecast dates, which b09 excludes as forecast dates. Lags are
#     taken after filling.
#   - Missing months are counted before filling (env_missingness.csv).
#   - Lags 1, 3, 6 and 12 are written; the models use lags 1 and 12 (e07 and
#     e08 also the current month).
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

# -------------------------------------------------------------------------
# Folders, archive, log
# -------------------------------------------------------------------------
# Output folder stem (short form b8_).
script_stem <- "b8_build_environment_timeseries"
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
data_out_dir <- dirs$data
diag_out_dir <- dirs$diag
archive_previous(c(data_out_dir, diag_out_dir))
start_log(dirs$logs, script_stem)
cat("Project directory:", PATHS$project_dir, "\n\n")

# -------------------------------------------------------------------------
# Packages (minimal)
# -------------------------------------------------------------------------
required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages(library(dplyr))

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------
covars_core <- c("TEMP", "SAL")
lag_months  <- c(1, 3, 6, 12)
cor_method  <- "pearson"

# -------------------------------------------------------------------------
# Inputs (fixed paths; b03 and b05 keep older copies in ARCHIVE_ subfolders)
# -------------------------------------------------------------------------
station_core_path <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
index_keys_path   <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.rds")
stop_if_missing(c(station_core_path, index_keys_path),
                c("b03_harmonize_station_metadata.R", "b05_define_spatiotemporal_units.R"))

cat("Using station_core:", station_core_path, "\n")
cat("Using index_keys:", index_keys_path, "\n\n")

station_core <- readRDS(station_core_path)
index_keys   <- readRDS(index_keys_path)

# Expectation from upstream objects (b3/b5)
station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays
if (is.null(valid_bays) || length(valid_bays) == 0) {
  stop("index_keys.rds has no keys$valid_bays (rerun b05_define_spatiotemporal_units.R)", call. = FALSE)
}
time_grid  <- index_keys$time_grid
expect_cols(time_grid, c("YEAR", "MONTH", "YEAR_MONTH"), "index_keys$time_grid")

# -------------------------------------------------------------------------
# Unified station table
# -------------------------------------------------------------------------
station_all <- bind_rows(station_bs_core, station_bt_core) %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i),
    YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH)
  ) %>%
  filter(BAY_ID %in% valid_bays) %>%
  select(SAMPLE_ID, GEAR, BAY_ID, YEAR, MONTH, YEAR_MONTH, any_of(covars_core))

covars_present <- covars_core
expect_cols(station_all, covars_core, "station_core")
# b03 stops on duplicated SAMPLE_ID, so each row here is one sample.

# -------------------------------------------------------------------------
# Monthly means on the full grid, then trim, fill and lag
# -------------------------------------------------------------------------
mean_cols <- paste0(covars_present, "__mean")
grid_ym <- time_grid$YEAR_MONTH

monthly_means <- function(df, keys) {
  df %>%
    group_by(across(all_of(c(keys, "YEAR", "MONTH", "YEAR_MONTH")))) %>%
    summarise(
      n_samples = n(),
      n_bays    = n_distinct(BAY_ID),
      across(all_of(covars_present), ~ { m <- mean(.x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m },
             .names = "{.col}__mean"),
      .groups = "drop"
    )
}

# One series on the full grid -> trimmed to its sampled span, interior gaps
# filled linearly (and marked), lags added. `need_bays` > 0 applies the
# coastwide rule: a month counts as sampled only when all bays were sampled.
build_series <- function(sub, need_bays = 0L) {
  full <- tibble::tibble(YEAR_MONTH = grid_ym) %>%
    left_join(time_grid, by = "YEAR_MONTH") %>%
    left_join(sub, by = c("YEAR", "MONTH", "YEAR_MONTH")) %>%
    arrange(YEAR, MONTH)
  sampled <- !is.na(full$n_samples)
  partial <- if (need_bays > 0) sampled & full$n_bays < need_bays else rep(FALSE, nrow(full))
  complete <- sampled & !partial
  if (!any(complete)) stop("A series has no complete month", call. = FALSE)
  first <- min(which(complete)); last <- max(which(complete))
  inside <- seq(first, last)
  if (any(partial[inside])) {
    stop("Coastwide months inside the record with fewer than all bays sampled: ",
         paste(full$YEAR_MONTH[inside][partial[inside]], collapse = ", "), call. = FALSE)
  }
  out <- full[inside, ]
  # Months with no sample get no mean; months with a partial set of bays are
  # only possible at the start or end, and those are cut above.
  out$filled <- !stats::complete.cases(out[, mean_cols])
  for (mc in mean_cols) {
    if (any(is.na(out[[mc]]))) {
      x <- out[[mc]]; i <- which(!is.na(x))
      out[[mc]] <- stats::approx(i, x[i], xout = seq_along(x), rule = 1)$y
    }
  }
  for (mc in mean_cols) for (k in lag_months) out[[paste0(mc, "__lag", k)]] <- dplyr::lag(out[[mc]], k)
  attr(out, "coverage") <- tibble::tibble(
    first_year_month = full$YEAR_MONTH[first], last_year_month = full$YEAR_MONTH[last],
    n_months = length(inside),
    n_months_dropped_leading = first - 1L, n_months_dropped_trailing = nrow(full) - last,
    n_months_filled = sum(out$filled),
    filled_months = paste(out$YEAR_MONTH[out$filled], collapse = "; ")
  )
  out
}

# Alpha: one series per bay and gear
alpha_means <- monthly_means(station_all, c("BAY_ID", "GEAR"))
alpha_parts <- alpha_means %>% group_by(BAY_ID, GEAR) %>% group_split()
alpha_built <- lapply(alpha_parts, function(p) {
  o <- build_series(p %>% select(-BAY_ID, -GEAR))
  cov <- attr(o, "coverage")
  o <- o %>% mutate(BAY_ID = p$BAY_ID[1], GEAR = p$GEAR[1], .before = 1)
  list(data = o, coverage = bind_cols(tibble::tibble(scale = "alpha", GEAR = p$GEAR[1], BAY_ID = p$BAY_ID[1]), cov))
})
env_alpha_index <- bind_rows(lapply(alpha_built, `[[`, "data")) %>%
  mutate(DATE = as.Date(paste0(YEAR_MONTH, "-01"))) %>%
  select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH, DATE, n_samples, filled, all_of(mean_cols), contains("__lag")) %>%
  arrange(BAY_ID, GEAR, DATE)

# Gamma: one series per gear, all bays pooled; a month needs every bay
gamma_means <- monthly_means(station_all, "GEAR")
gamma_built <- lapply(split(gamma_means, gamma_means$GEAR), function(p) {
  o <- build_series(p %>% select(-GEAR), need_bays = length(valid_bays))
  cov <- attr(o, "coverage")
  o <- o %>% mutate(GEAR = p$GEAR[1], .before = 1)
  list(data = o, coverage = bind_cols(tibble::tibble(scale = "gamma", GEAR = p$GEAR[1], BAY_ID = NA_integer_), cov))
})
env_gamma_index <- bind_rows(unname(lapply(gamma_built, `[[`, "data"))) %>%
  mutate(DATE = as.Date(paste0(YEAR_MONTH, "-01"))) %>%
  select(GEAR, YEAR, MONTH, YEAR_MONTH, DATE, n_samples, n_bays, filled, all_of(mean_cols), contains("__lag")) %>%
  arrange(GEAR, DATE)

# -------------------------------------------------------------------------
# Checks
# -------------------------------------------------------------------------
if (nrow(env_gamma_index) != nrow(distinct(env_gamma_index, GEAR, YEAR_MONTH)) ||
    nrow(env_alpha_index) != nrow(distinct(env_alpha_index, BAY_ID, GEAR, YEAR_MONTH))) {
  stop("Duplicated months in the covariate series", call. = FALSE)
}
if (anyNA(env_alpha_index[, mean_cols]) || anyNA(env_gamma_index[, mean_cols])) {
  stop("Missing covariate means after filling", call. = FALSE)
}
for (g in c("BS", "BT")) {
  b <- sort(unique(env_alpha_index$BAY_ID[env_alpha_index$GEAR == g]))
  if (!identical(as.integer(b), as.integer(sort(valid_bays)))) stop("Gear ", g, " is missing bays", call. = FALSE)
}

# -------------------------------------------------------------------------
# Missing months, counted before filling
# -------------------------------------------------------------------------
env_missingness <- bind_rows(unname(c(lapply(alpha_built, `[[`, "coverage"),
                                      lapply(gamma_built, `[[`, "coverage")))) %>%
  arrange(scale, GEAR, BAY_ID)
miss_path <- file.path(diag_out_dir, "env_missingness.csv")
readr::write_csv(env_missingness, miss_path, na = "")
cat("\n=== Covariate series: months dropped and filled (grid of", length(grid_ym), "months) ===\n")
print(as.data.frame(env_missingness %>% select(scale, GEAR, BAY_ID, first_year_month, last_year_month,
                                               n_months, n_months_dropped_leading, n_months_filled, filled_months)),
      row.names = FALSE)
cat("Saved:", miss_path, "\n\n")

mean_cols_alpha <- mean_cols
mean_cols_gamma <- mean_cols

# -------------------------------------------------------------------------
# Correlation matrices
# -------------------------------------------------------------------------
write_cor_matrix <- function(df, cols, out_csv, method = "pearson") {
  sub <- df %>%
    select(all_of(cols)) %>%
    select(where(is.numeric)) %>%
    select(where(~ !all(is.na(.x))))
  
  if (ncol(sub) < 2) {
    warning("Not enough non-NA numeric columns for correlations: ", out_csv)
    return(invisible(NULL))
  }
  
  cor_mat <- suppressWarnings(stats::cor(sub, use = "pairwise.complete.obs", method = method))
  cor_df <- data.frame(covar = rownames(cor_mat), cor_mat, row.names = NULL, check.names = FALSE)
  
  readr::write_csv(cor_df, out_csv)
  cat("Saved:", out_csv, "\n")
}

write_cor_matrix(
  env_alpha_index, mean_cols_alpha,
  file.path(diag_out_dir, "cor_matrix__alpha_index__means.csv"),
  method = cor_method
)

write_cor_matrix(
  env_gamma_index, mean_cols_gamma,
  file.path(diag_out_dir, "cor_matrix__gamma_index__means.csv"),
  method = cor_method
)

# -------------------------------------------------------------------------
# Save outputs (derived data)
# -------------------------------------------------------------------------
out_alpha_rds <- file.path(data_out_dir, "env_alpha_index.rds")
out_gamma_rds <- file.path(data_out_dir, "env_gamma_index.rds")

saveRDS(env_alpha_index, out_alpha_rds)
saveRDS(env_gamma_index, out_gamma_rds)

readr::write_csv(env_alpha_index, file.path(data_out_dir, "env_alpha_index.csv"))
readr::write_csv(env_gamma_index, file.path(data_out_dir, "env_gamma_index.csv"))

end_log("finished")
cat("b08_build_environment_timeseries.R completed\n")
