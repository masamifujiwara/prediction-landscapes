# =============================================================================
# Script: b8_build_environment_timeseries.R
# Purpose: Build monthly environmental covariate time series (alpha + gamma grids)
# Outputs:
#   - Derived data: 01_Analysis/Data/b8_build_environment_timeseries/...
#   - Diagnostics:  02_Results/03_Diagnostics/b8_build_environment_timeseries/...
# =============================================================================

rm(list = ls())

# -------------------------------------------------------------------------
# Required bootstrap snippet (must include exactly)
# -------------------------------------------------------------------------
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

# -------------------------------------------------------------------------
# PATHS.RData usage (required)
# -------------------------------------------------------------------------
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file) # loads PATHS

stopifnot(exists("PATHS"))

# Validate PATHS project_dir matches computed project_dir
paths_project_dir <- tryCatch(PATHS$project_dir, error = function(e) NULL)
if (is.null(paths_project_dir) || !nzchar(paths_project_dir)) {
  stop("PATHS$project_dir is missing/empty. Rebuild PATHS in b1.")
}
paths_project_dir <- normalizePath(paths_project_dir, winslash = "/", mustWork = TRUE)
stopifnot(identical(paths_project_dir, project_dir))

# -------------------------------------------------------------------------
# Logging (to 01_Analysis/Logs/<script_stem>/...)
# -------------------------------------------------------------------------
script_stem <- "b8_build_environment_timeseries"

logs_dir <- file.path(PATHS$project_dir, "01_Analysis", "Logs", script_stem)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs_dir, paste0(script_stem, "__", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

# Write console output to both console and log
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message", append = TRUE)

cat("Script:", script_stem, "\n")
cat("Working directory:", getwd(), "\n")
cat("Project directory:", project_dir, "\n")
cat("Log file:", log_file, "\n\n")

on.exit({
  cat("\n[END]", script_stem, "completed.\n")
  sink(type = "message")
  sink()
  close(zz)
}, add = TRUE)

# -------------------------------------------------------------------------
# Packages (minimal)
# -------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forecast)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------------------------------------------------------------
# Parameters (scientific logic unchanged)
# -------------------------------------------------------------------------
covars_core <- c("TEMP", "SAL", "DEPTH_M", "DO", "TURB")
lag_months  <- c(1, 3, 6, 12)
cor_method  <- "pearson"

# If TRUE, enforce unique SAMPLE_ID within gear before monthly summaries.
# Turn OFF only if you *expect* multiple rows per sample_id and want them counted.
DEDUP_BY_SAMPLE <- TRUE

# -------------------------------------------------------------------------
# Output directories (must follow conventions)
# -------------------------------------------------------------------------
data_out_dir <- file.path(PATHS$project_dir, "01_Analysis", "Data", script_stem)
diag_out_dir <- file.path(PATHS$project_dir, "02_Results", "03_Diagnostics", script_stem)

dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Inputs: locate derived files created upstream (station_core.rds, index_keys.rds)
#   We search under 01_Analysis/Data/ to avoid hardcoding prior script stems.
# -------------------------------------------------------------------------
analysis_data_root <- file.path(PATHS$project_dir, "01_Analysis", "Data")
stopifnot(dir.exists(analysis_data_root))

find_single_file <- function(root_dir, filename) {
  hits <- list.files(root_dir, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("Input not found: ", filename, " under ", root_dir)
  if (length(hits) > 1) {
    stop(
      "Ambiguous input (multiple matches) for ", filename, " under ", root_dir, "\n",
      "Matches:\n", paste(" -", normalizePath(hits, winslash = "/"), collapse = "\n")
    )
  }
  normalizePath(hits, winslash = "/", mustWork = TRUE)
}

station_core_path <- find_single_file(analysis_data_root, "station_core.rds")
index_keys_path   <- find_single_file(analysis_data_root, "index_keys.rds")

cat("Using station_core:", station_core_path, "\n")
cat("Using index_keys:", index_keys_path, "\n\n")

station_core <- readRDS(station_core_path)
index_keys   <- readRDS(index_keys_path)

# Expectation from upstream objects (b3/b5)
station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays %||% 1:8
time_grid  <- index_keys$time_grid
stopifnot(all(c("YEAR", "MONTH", "YEAR_MONTH") %in% names(time_grid)))

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

covars_present <- covars_core[covars_core %in% names(station_all)]
covars_missing <- setdiff(covars_core, covars_present)
if (length(covars_missing) > 0) {
  warning("Missing covariates in station_core (skipping): ", paste(covars_missing, collapse = ", "))
}
stopifnot(length(covars_present) > 0)

if (isTRUE(DEDUP_BY_SAMPLE)) {
  dup_n <- station_all %>% count(GEAR, SAMPLE_ID) %>% filter(n > 1) %>% nrow()
  if (dup_n > 0) warning("Detected duplicated SAMPLE_ID rows within GEAR (n=", dup_n, "). Deduplicating.")
  station_all <- station_all %>% distinct(GEAR, SAMPLE_ID, .keep_all = TRUE)
}

# -------------------------------------------------------------------------
# α and γ index grids
# -------------------------------------------------------------------------
alpha_index <- tidyr::expand_grid(
  BAY_ID = valid_bays,
  GEAR   = sort(unique(station_all$GEAR)),
  time_grid
) %>%
  select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH) %>%
  arrange(BAY_ID, GEAR, YEAR, MONTH)

gamma_index <- tidyr::expand_grid(
  GEAR = sort(unique(station_all$GEAR)),
  time_grid
) %>%
  select(GEAR, YEAR, MONTH, YEAR_MONTH) %>%
  arrange(GEAR, YEAR, MONTH)

# -------------------------------------------------------------------------
# α summaries
# -------------------------------------------------------------------------
alpha_env_summary <- station_all %>%
  group_by(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH) %>%
  summarise(
    n_samples = n(),
    across(all_of(covars_present),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

env_alpha_index <- alpha_index %>%
  left_join(alpha_env_summary, by = c("BAY_ID", "GEAR", "YEAR", "MONTH", "YEAR_MONTH")) %>%
  mutate(DATE = as.Date(paste0(YEAR_MONTH, "-01"))) %>%
  arrange(BAY_ID, GEAR, DATE)

mean_cols_alpha <- paste0(covars_present, "__mean")

# -------------------------------------------------------------------------
# Impute missing environmental covariates before calculating lags
# -------------------------------------------------------------------------
env_alpha_index <- env_alpha_index %>%
  group_by(BAY_ID, GEAR) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(across(all_of(mean_cols_alpha), ~ as.numeric(forecast::na.interp(.x)))) %>%
  ungroup()
# -------------------------------------------------------------------------

if (length(lag_months) > 0) {
  env_alpha_index <- env_alpha_index %>%
    group_by(BAY_ID, GEAR) %>%
    arrange(DATE, .by_group = TRUE) %>%
    mutate(
      across(
        all_of(mean_cols_alpha),
        .fns = setNames(lapply(lag_months, function(k) function(x) dplyr::lag(x, k)),
                        paste0("lag", lag_months)),
        .names = "{.col}__{.fn}"
      )
    ) %>%
    ungroup()
}

# -------------------------------------------------------------------------
# γ summaries
# -------------------------------------------------------------------------
gamma_env_summary <- station_all %>%
  group_by(GEAR, YEAR, MONTH, YEAR_MONTH) %>%
  summarise(
    n_samples = n(),
    across(all_of(covars_present),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

env_gamma_index <- gamma_index %>%
  left_join(gamma_env_summary, by = c("GEAR", "YEAR", "MONTH", "YEAR_MONTH")) %>%
  mutate(DATE = as.Date(paste0(YEAR_MONTH, "-01"))) %>%
  arrange(GEAR, DATE)

mean_cols_gamma <- paste0(covars_present, "__mean")

# -------------------------------------------------------------------------
# NEW: Impute missing environmental covariates before calculating lags
# -------------------------------------------------------------------------
env_gamma_index <- env_gamma_index %>%
  group_by(GEAR) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(across(all_of(mean_cols_gamma), ~ as.numeric(forecast::na.interp(.x)))) %>%
  ungroup()
# -------------------------------------------------------------------------

if (length(lag_months) > 0) {
  env_gamma_index <- env_gamma_index %>%
    group_by(GEAR) %>%
    arrange(DATE, .by_group = TRUE) %>%
    mutate(
      across(
        all_of(mean_cols_gamma),
        .fns = setNames(lapply(lag_months, function(k) function(x) dplyr::lag(x, k)),
                        paste0("lag", lag_months)),
        .names = "{.col}__{.fn}"
      )
    ) %>%
    ungroup()
}

# -------------------------------------------------------------------------
# HARD uniqueness checks (unchanged)
# -------------------------------------------------------------------------
stopifnot(nrow(env_gamma_index) == nrow(distinct(env_gamma_index, GEAR, YEAR_MONTH)))
stopifnot(nrow(env_alpha_index) == nrow(distinct(env_alpha_index, BAY_ID, GEAR, YEAR_MONTH)))

# -------------------------------------------------------------------------
# Diagnostics (missingness heatmaps)
#   NOTE: Removed rlang dependency by creating an explicit facet label column.
#   Scientific logic unchanged; plotting-only change.
# -------------------------------------------------------------------------
make_missingness_heatmap <- function(df, facet_label_col, value_cols, title, out_png) {
  miss_rate <- df %>%
    mutate(
      DATE = as.Date(paste0(YEAR_MONTH, "-01")),
      facet_label = .data[[facet_label_col]]
    ) %>%
    select(facet_label, YEAR_MONTH, DATE, all_of(value_cols)) %>%
    pivot_longer(cols = all_of(value_cols), names_to = "covar", values_to = "value") %>%
    mutate(missing = is.na(value)) %>%
    group_by(facet_label, covar, YEAR_MONTH, DATE) %>%
    summarise(miss_rate = mean(missing), .groups = "drop")
  
  p <- ggplot(miss_rate, aes(x = DATE, y = covar, fill = miss_rate)) +
    geom_tile() +
    facet_wrap(~ facet_label, ncol = 1) +
    labs(x = "Time", y = NULL, title = title, fill = "Missing\nrate") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
  ggsave(out_png, p, width = 12, height = 8, dpi = 300)
  cat("Saved:", out_png, "\n")
}

env_alpha_index <- env_alpha_index %>%
  mutate(facet_alpha = paste0("BAY_ID=", BAY_ID, " | GEAR=", GEAR))

env_gamma_index <- env_gamma_index %>%
  mutate(facet_gamma = paste0("GEAR=", GEAR))

alpha_heat_path <- file.path(diag_out_dir, "missingness_heatmap__alpha_index__means.png")
gamma_heat_path <- file.path(diag_out_dir, "missingness_heatmap__gamma_index__means.png")

make_missingness_heatmap(
  env_alpha_index,
  facet_label_col = "facet_alpha",
  value_cols = mean_cols_alpha,
  title = "Missingness heatmap (α index): covariate means",
  out_png = alpha_heat_path
)

make_missingness_heatmap(
  env_gamma_index,
  facet_label_col = "facet_gamma",
  value_cols = mean_cols_gamma,
  title = "Missingness heatmap (γ index): covariate means",
  out_png = gamma_heat_path
)

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

cat("\nOK:", script_stem, "completed successfully\n")