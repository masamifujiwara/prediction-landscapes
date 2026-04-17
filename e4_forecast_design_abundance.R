# =============================================================================
# Script: e4_forecast_design_abundance.R
#
# Purpose:
# Define rolling-origin backtesting folds for forecasting population abundance 
# (per series × gear × species), using deseasonalized series from e3.
# 
# This strictly mirrors the `c3` script used for diversity to ensure a perfect 
# apples-to-apples predictability comparison between diversity and abundance.
#
# Outputs:
# - Derived: 01_Analysis/Data/e4_forecast_design_abundance/backtest_folds_abundance.rds
# - Diags:   02_Results/03_Diagnostics/e4_forecast_design_abundance/backtest_groups_check.csv
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS (Windows & Linux Safe)
# ─────────────────────────────────────────────────────────────────────────────
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
  ofile <- tryCatch(sys.frames()[[6]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  # 4) fallback: assume you launched from 01_Analysis
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
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
# 3) Script stem + directories
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e4_forecast_design_abundance"

dir_data_script <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_diag_script <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Inputs & Outputs
# ─────────────────────────────────────────────────────────────────────────────
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds")
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")

stopifnot(file.exists(IN_ALPHA), file.exists(IN_GAMMA))

OUT_FOLDS <- file.path(dir_data_script, "backtest_folds_abundance.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Forecast design parameters (Identical to Diversity `c3`)
# ─────────────────────────────────────────────────────────────────────────────
HORIZONS <- c(1, 12)
TRAINING_RULE <- "expanding"    
MIN_TRAIN_MONTHS <- 120         
FIXED_TRAIN_MONTHS <- 120       
ORIGIN_STEP <- 1
HOLDOUT_LAST_N <- 0
TARGET_COL <- "value_deseas"

# ─────────────────────────────────────────────────────────────────────────────
# 6) Helpers
# ─────────────────────────────────────────────────────────────────────────────
make_group_id <- function(series_name, bay_id, gear, sp) {
  paste(series_name, bay_id, gear, paste0("sp", sp), sep = "__")
}

complete_monthly <- function(df) {
  df <- df %>% arrange(date)
  d0 <- min(df$date, na.rm = TRUE)
  d1 <- max(df$date, na.rm = TRUE)
  grid <- tibble(date = seq.Date(d0, d1, by = "month"))
  grid %>% left_join(df, by = "date") %>% arrange(date)
}

build_folds_one_series <- function(df, series_name, horizons, training_rule, min_train, fixed_train, step, holdout) {
  df <- df %>% arrange(date)
  n <- nrow(df)
  
  n_effective <- if (holdout > 0) max(0, n - holdout) else n
  if (n_effective < (min_train + max(horizons) + 1)) {
    return(list(ok = FALSE, reason = "Too short for requested min_train + max_horizon", folds = list(), n = n, n_effective = n_effective))
  }
  
  folds_out <- list()
  for (h in horizons) {
    origins <- seq(from = min_train, to = (n_effective - h), by = step)
    folds_h <- lapply(origins, function(t0) {
      train_idx <- if (training_rule == "expanding") 1:t0 else max(1, t0 - fixed_train + 1):t0
      list(horizon = h, origin_i = t0, test_i = t0 + h, train_idx = train_idx)
    })
    folds_out[[paste0("h", h)]] <- folds_h
  }
  
  list(ok = TRUE, reason = NA_character_, folds = folds_out, n = n, n_effective = n_effective)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7) Load and Prepare Time Series
# ─────────────────────────────────────────────────────────────────────────────
alpha <- readRDS(IN_ALPHA)
gamma <- readRDS(IN_GAMMA)

### Alpha: Retain distinct BAY_IDs and exclude Sabine Lake (Bay 1)
alpha2 <- alpha %>% transmute(
  series = "alpha",
  BAY_ID = as.character(BAY_ID),
  GEAR = as.character(GEAR),
  SPECIES_CODE = as.character(SPECIES_CODE),
  date = as.Date(date),
  y = as.numeric(.data[[TARGET_COL]])
) %>% filter(BAY_ID != "1") 

### Gamma: Assign placeholder BAY_ID for coast-wide scale
gamma2 <- gamma %>% transmute(
  series = "gamma",
  BAY_ID = "ALL", 
  GEAR = as.character(GEAR),
  SPECIES_CODE = as.character(SPECIES_CODE),
  date = as.Date(date),
  y = as.numeric(.data[[TARGET_COL]])
)

ts_all <- bind_rows(alpha2, gamma2) %>% 
  filter(!is.na(date)) %>% 
  arrange(series, BAY_ID, GEAR, SPECIES_CODE, date)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Build Folds per (series × gear × species) on complete monthly grids
# ─────────────────────────────────────────────────────────────────────────────
groups_split <- ts_all %>% group_by(series, BAY_ID, GEAR, SPECIES_CODE) %>% group_split()
group_ids <- vapply(
  groups_split,
  function(x) make_group_id(unique(x$series), unique(x$BAY_ID), unique(x$GEAR), unique(x$SPECIES_CODE)),
  character(1)
)
groups_named <- stats::setNames(groups_split, group_ids)

cat("Building rolling-origin folds for", length(groups_named), "groups...\n")

groups_folds <- lapply(names(groups_named), function(gid) {
  df_g <- groups_named[[gid]]
  series_name <- unique(df_g$series)
  bay <- unique(df_g$BAY_ID)
  gear <- unique(df_g$GEAR)
  sp <- unique(df_g$SPECIES_CODE)
  
  df_g2 <- df_g %>% select(date, y) %>% complete_monthly()
  
  res <- build_folds_one_series(
    df = df_g2, 
    series_name = series_name, 
    horizons = HORIZONS, 
    training_rule = TRAINING_RULE, 
    min_train = MIN_TRAIN_MONTHS, 
    fixed_train = FIXED_TRAIN_MONTHS, 
    step = ORIGIN_STEP, 
    holdout = HOLDOUT_LAST_N
  )
  
  list(
    ok = res$ok,
    reason = res$reason,
    series = series_name,
    BAY_ID = bay,      
    GEAR = gear,
    SPECIES_CODE = sp,
    group_id = gid,
    dates = df_g2$date,
    y = df_g2$y,
    n = res$n,
    n_effective = res$n_effective,
    folds = res$folds
  )
})
names(groups_folds) <- names(groups_named)

# ─────────────────────────────────────────────────────────────────────────────
# 9) Save Fold Output & Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
backtest_folds <- list(
  created = Sys.time(),
  target_col = TARGET_COL,
  horizons = HORIZONS,
  training_rule = TRAINING_RULE,
  min_train_months = MIN_TRAIN_MONTHS,
  fixed_train_months = if (TRAINING_RULE == "fixed") FIXED_TRAIN_MONTHS else NA_integer_,
  origin_step = ORIGIN_STEP,
  holdout_last_n = HOLDOUT_LAST_N,
  groups = groups_folds
)

saveRDS(backtest_folds, OUT_FOLDS)

# Generate diagnostic check table
check_tbl <- dplyr::bind_rows(lapply(names(backtest_folds$groups), function(gid) {
  g <- backtest_folds$groups[[gid]]
  if (!isTRUE(g$ok)) {
    return(tibble(
      group_id = gid, series = g$series, BAY_ID = g$BAY_ID, GEAR = g$GEAR, SPECIES_CODE = g$SPECIES_CODE,
      ok = FALSE, reason = g$reason, n = g$n, n_effective = g$n_effective, n_folds = 0L
    ))
  }
  n_folds <- sum(vapply(g$folds, length, integer(1)))
  tibble(
    group_id = gid, series = g$series, BAY_ID = g$BAY_ID, GEAR = g$GEAR, SPECIES_CODE = g$SPECIES_CODE,
    ok = TRUE, reason = NA_character_, n = g$n, n_effective = g$n_effective, n_folds = n_folds
  )
}))

diag_out <- file.path(dir_diag_script, "backtest_groups_check.csv")
readr::write_csv(check_tbl, diag_out)

cat("\n✅ Forecast design complete.\n")
cat("💾 Saved abundance folds to:", OUT_FOLDS, "\n")
cat("🧾 Diagnostics saved to:", diag_out, "\n")