### =============================================================================
### Script: d1_compare_model_metrics_dictionary.R  (REVISED)
### Purpose:
### Build a SINGLE, enforced “metrics dictionary” and compute comparable forecast
### metrics across models, for:
### - targets: alpha AND gamma
### - q: 0, 1, 2
### - gear: bagseine (BS) and bay trawl (BT) separately
### - horizons: h = 1 and 12
### - folds: rolling-origin on COMPLETE monthly grid (assumed by upstream c* scripts)
### Key revision:
### - Added c10 and c11 script directories to capture SSM-XREG models.
### - Syntax fix: wrapped custom infix operator %||% in backticks.
### =============================================================================
rm(list = ls())

### -----------------------------
### Required bootstrap snippet
### -----------------------------
get_script_dir <- function() {
  ### 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) { 
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL) 
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/"))) 
  }
  ### 2) Rscript / R CMD BATCH: --file=...
  cmd <- commandArgs(trailingOnly = FALSE) 
  file_arg <- grep("^--file=", cmd, value = TRUE) 
  if (length(file_arg) == 1) { 
    p <- sub("^--file=", "", file_arg) 
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/"))) 
  }
  ### 3) source("script.R")
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL) 
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  ### 4) fallback: assume you launched from 01_Analysis
  normalizePath(getwd(), winslash = "/") 
}
script_dir <- get_script_dir() 
setwd(script_dir)

find_paths_file <- function(start_dir, max_up = 8) { 
  d <- normalizePath(start_dir, winslash = "/", mustWork = TRUE) 
  for (i in 0:max_up) { 
    candidate <- file.path(d, "01_Analysis", "PATHS.RData") 
    if (file.exists(candidate)) return(candidate) 
    d2 <- dirname(d); if (identical(d2, d)) break 
    d <- d2 
  } 
  NA_character_ 
}

paths_file <- find_paths_file(script_dir, max_up = 8) 
if (is.na(paths_file)) { stop("Could not find 01_Analysis/PATHS.RData by walking up from: ", script_dir) } 
load(paths_file) # loads PATHS 
stopifnot(exists("PATHS"))

### Optional integrity check if PATHS$project_dir exists
if (!is.null(PATHS$project_dir)) {
  paths_proj <- normalizePath(PATHS$project_dir, winslash = "/", mustWork = TRUE) 
}

### -----------------------------
### Dependencies
### -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "purrr") 
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)] 
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({ 
  library(dplyr) 
  library(tidyr) 
  library(readr) 
  library(stringr) 
  library(purrr) 
})

### FIXED: added backticks around %||%
`%||%` <- function(a, b) if (!is.null(a)) a else b

### -----------------------------
### Script stem + directories
### -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "d1_compare_model_metrics_dictionary.R")) 
if (!nzchar(script_stem) || script_stem == "NULL") script_stem <- "d1_compare_model_metrics_dictionary"

dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem) 
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE) 
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE) 
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)

OUT_DICT_CSV_D <- file.path(dir_data_script,  "metrics_dictionary_alpha_gamma.csv") 
OUT_DICT_RDS_D <- file.path(dir_data_script,  "metrics_dictionary_alpha_gamma.rds") 
OUT_DICT_CSV_R <- file.path(dir_model_script, "metrics_dictionary_alpha_gamma.csv") 
OUT_DICT_RDS_R <- file.path(dir_model_script, "metrics_dictionary_alpha_gamma.rds") 
OUT_DISC_LOG   <- file.path(dir_diag_script,  "discovery_log.csv")

### -----------------------------
### Config
### -----------------------------
HORIZONS <- c(1, 12) 
Q_LEVELS <- c("0", "1", "2") 
GEARS_KEEP <- c("BS", "BT") 
TARGETS_KEEP <- c("alpha", "gamma")

BASELINE_MODEL_PATTERNS <- c("seasonal_naive", "seasonal-naive", "seasonal naive", "snaive", "naive_seasonal")

### Canonical forecast-producing scripts (REVISED to include c10 and c11)
FORECAST_SCRIPT_DIRS <- c( 
  "c4p_baseline_forecasts_gamma_parallel", 
  "c5_state_space_models_forecasts", 
  "c6_gamma_arimax_from_env_index", 
  "c7_alpha_forecasting_multilevel_or_bay_replication", 
  "c8_alpha_state_space_models_forecasts_parallel", 
  "c9_alpha_arimax_from_env_index_parallel",
  "c10_alpha_state_space_with_env_covariates_parallel",
  "c11_gamma_state_space_with_env_covariates_parallel"
)

### -----------------------------
### Metric utilities (dictionary)
### -----------------------------
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE)) 
mae  <- function(e) mean(abs(e), na.rm = TRUE)

coverage_rate <- function(y, lo, hi) { 
  ok <- !is.na(y) & !is.na(lo) & !is.na(hi) 
  if (!any(ok)) return(NA_real_) 
  mean(y[ok] >= lo[ok] & y[ok] <= hi[ok]) 
}

### -----------------------------
### Forecast-table standardization
### -----------------------------
standardize_forecast_tbl <- function(x, source_path) { 
  if (!is.data.frame(x)) return(NULL) 
  nm <- names(x)
  
  ### truth columns 
  if (!("y_true" %in% nm)) { 
    if ("truth" %in% nm)  x <- dplyr::rename(x, y_true = truth) 
    if ("actual" %in% nm) x <- dplyr::rename(x, y_true = actual) 
    if ("y" %in% nm)      x <- dplyr::rename(x, y_true = y) 
    if ("y_obs" %in% nm)  x <- dplyr::rename(x, y_true = y_obs) 
    if ("obs" %in% nm)    x <- dplyr::rename(x, y_true = obs) 
  }
  ### prediction columns 
  if (!("y_pred" %in% nm)) { 
    if ("pred" %in% nm)     x <- dplyr::rename(x, y_pred = pred) 
    if ("forecast" %in% nm) x <- dplyr::rename(x, y_pred = forecast) 
    if ("mean" %in% nm)     x <- dplyr::rename(x, y_pred = mean) 
    if ("y_hat" %in% nm)    x <- dplyr::rename(x, y_pred = y_hat) 
  }
  ### horizon
  if (!("horizon" %in% names(x))) { 
    if ("h" %in% names(x)) x <- dplyr::rename(x, horizon = h) 
    if ("H" %in% names(x)) x <- dplyr::rename(x, horizon = H) 
  }
  ### date
  if (!("test_date" %in% names(x))) { 
    if ("date" %in% names(x)) x <- dplyr::rename(x, test_date = date) 
    if ("DATE" %in% names(x)) x <- dplyr::rename(x, test_date = DATE) 
  }
  ### gear
  if (!("GEAR" %in% names(x))) { 
    gcol <- intersect(names(x), c("gear", "Gear")) 
    if (length(gcol) >= 1) x <- dplyr::rename(x, GEAR = all_of(gcol)) 
  }
  ### q
  if (!("q" %in% names(x))) { 
    qcol <- intersect(names(x), c("Q", "hill_q", "order_q")) 
    if (length(qcol) >= 1) x <- dplyr::rename(x, q = all_of(qcol)) 
  }
  ### bay optional
  if (!("BAY" %in% names(x))) { 
    bcol <- intersect(names(x), c("BAY_ID", "bay", "Bay", "BAYID")) 
    if (length(bcol) >= 1) x <- dplyr::rename(x, BAY = all_of(bcol)) 
  }
  ### ok flag optional
  if (!("ok" %in% names(x))) { 
    if ("success" %in% names(x)) x <- dplyr::rename(x, ok = success) 
  }
  ### model
  if (!("model" %in% names(x))) { 
    inferred <- basename(source_path) 
    inferred <- str_replace_all(inferred, "\\.rds$", "")
    x$model <- inferred 
  }
  ### PI cols optional
  if (!("lower80" %in% names(x))) x$lower80 <- NA_real_
  if (!("upper80" %in% names(x))) x$upper80 <- NA_real_ 
  if (!("lower95" %in% names(x))) x$lower95 <- NA_real_
  if (!("upper95" %in% names(x))) x$upper95 <- NA_real_
  
  need <- c("y_true", "y_pred", "horizon", "GEAR", "q", "model") 
  if (any(!need %in% names(x))) return(NULL)
  
  x <- x %>% mutate( 
    GEAR = as.character(GEAR), 
    q = as.character(q), 
    horizon = as.integer(horizon), 
    test_date = if ("test_date" %in% names(.)) as.Date(test_date) else as.Date(NA) 
  ) %>% filter(horizon %in% HORIZONS, q %in% Q_LEVELS, GEAR %in% GEARS_KEEP)
  
  if (nrow(x) == 0) return(NULL)
  
  ### infer target from path/model
  path_l <- tolower(source_path) 
  model_l <- tolower(paste(unique(x$model), collapse = " "))
  target <- dplyr::case_when(
    str_detect(path_l, "alpha") | str_detect(model_l, "alpha") ~ "alpha",
    str_detect(path_l, "gamma") | str_detect(model_l, "gamma") ~ "gamma",
    TRUE ~ NA_character_
  )
  x$target <- target
  
  ### gamma sometimes has no bay; fill ALL
  if (!("BAY" %in% names(x))) x$BAY <- "ALL"
  x$BAY <- as.character(x$BAY)
  x$source_path <- source_path 
  x 
}

is_baseline_model <- function(model_name) { 
  ml <- tolower(model_name) 
  any(vapply(BASELINE_MODEL_PATTERNS, function(p) str_detect(ml, fixed(p)), logical(1))) || 
    (str_detect(ml, "seasonal") && str_detect(ml, "naive")) 
}

### -----------------------------
### Discover forecast RDS files (canonical locations only)
### -----------------------------
discover_rds_in_dir <- function(root) { 
  if (!dir.exists(root)) return(character(0)) 
  list.files(root, recursive = TRUE, pattern = "\\.rds$", full.names = TRUE) 
}

roots <- c(
  ### Results-side model outputs
  file.path(PATHS$results_dir, "04_Models", FORECAST_SCRIPT_DIRS),
  ### Analysis-side derived outputs
  file.path(PATHS$analysis_dir, "Data", FORECAST_SCRIPT_DIRS) 
)

rds_files <- unique(unlist(lapply(roots, discover_rds_in_dir))) 
rds_files <- rds_files[str_detect(tolower(rds_files), "forecast")]  

if (length(rds_files) == 0) { 
  stop( "No forecast .rds files found in canonical locations.\n", 
        "Checked roots:\n  - ", paste(roots, collapse = "\n  - ") ) 
}

disc_rows <- list() 
fc_list <- list()

for (fp in rds_files) { 
  obj <- tryCatch(readRDS(fp), error = function(e) e) 
  if (inherits(obj, "error")) { 
    disc_rows[[length(disc_rows) + 1]] <- tibble(file = fp, status = "read_error", note = conditionMessage(obj)) 
    next 
  }
  
  df <- obj
  if (is.list(obj) && !is.data.frame(obj)) { 
    if ("predictions" %in% names(obj) && is.data.frame(obj$predictions)) {
      df <- obj$predictions 
    } else { 
      idx <- which(vapply(obj, is.data.frame, logical(1))) 
      if (length(idx) >= 1) df <- obj[[idx[1]]] 
    } 
  }
  
  std <- standardize_forecast_tbl(df, fp) 
  if (is.null(std)) { 
    disc_rows[[length(disc_rows) + 1]] <- tibble(file = fp, status = "skipped", note = "Not a recognized forecast table.") 
    next 
  }
  
  disc_rows[[length(disc_rows) + 1]] <- tibble( 
    file = fp, 
    status = "loaded", 
    note = paste0("n=", nrow(std), "; target=", std$target %||% "NA") 
  ) 
  fc_list[[length(fc_list) + 1]] <- std 
}

discovery_log <- bind_rows(disc_rows) 
readr::write_csv(discovery_log, OUT_DISC_LOG)

fc_all <- bind_rows(fc_list) %>% filter(target %in% TARGETS_KEEP)

if (nrow(fc_all) == 0) { 
  stop("No usable alpha/gamma forecast tables were discovered. See: ", OUT_DISC_LOG) 
}

### -----------------------------
### Baseline metrics per stratum
### -----------------------------
fc_all <- fc_all %>% 
  mutate(
    model_clean = as.character(model), 
    model_is_baseline = vapply(model_clean, is_baseline_model, logical(1))
  )

baseline_tbl <- fc_all %>% 
  filter(model_is_baseline) %>% 
  group_by(target, GEAR, q, BAY, horizon) %>% 
  summarise( 
    n_base = sum(!is.na(y_true) & !is.na(y_pred)), 
    rmse_base = rmse(y_pred - y_true), 
    mae_base  = mae(y_pred - y_true), 
    .groups = "drop" 
  )

##### -----------------------------
##### Metrics dictionary for every model
##### -----------------------------
metrics_dict <- fc_all %>% 
  group_by(target, GEAR, q, BAY, horizon, model_clean) %>% 
  summarise( 
    n = sum(!is.na(y_true) & !is.na(y_pred)), 
    rmse = rmse(y_pred - y_true), 
    mae  = mae(y_pred - y_true), 
    coverage80 = coverage_rate(y_true, lower80, upper80), 
    coverage95 = coverage_rate(y_true, lower95, upper95), 
    frac_failed = if ("ok" %in% names(pick(everything()))) mean(!ok, na.rm = TRUE) else NA_real_, 
    .groups = "drop" 
  ) %>% 
  left_join(baseline_tbl, by = c("target", "GEAR", "q", "BAY", "horizon")) %>% 
  mutate( 
    skill_rmse_ratio = ifelse(is.na(rmse_base) | rmse_base == 0, NA_real_, rmse / rmse_base), 
    skill_mae_ratio  = ifelse(is.na(mae_base)  | mae_base  == 0, NA_real_, mae  / mae_base), 
    skill_rmse_improve_pct = ifelse(is.na(skill_rmse_ratio), NA_real_, 100 * (1 - skill_rmse_ratio)), 
    skill_mae_improve_pct  = ifelse(is.na(skill_mae_ratio),  NA_real_, 100 * (1 - skill_mae_ratio)), 
    metrics_version = "d1_v1", 
    horizons = paste(HORIZONS, collapse = ","), 
    target_definition = "deseasonalized (value_deseas/residual_component upstream)", 
    folds_definition  = "rolling-origin on complete monthly grid (upstream)", 
    baseline_definition = "seasonal-naive (identified by model name patterns)" 
  ) %>% 
  arrange(target, GEAR, q, BAY, horizon, model_clean)

readr::write_csv(metrics_dict, OUT_DICT_CSV_D)
readr::write_csv(metrics_dict, OUT_DICT_CSV_R)
saveRDS(metrics_dict, OUT_DICT_RDS_D)
saveRDS(metrics_dict, OUT_DICT_RDS_R)

# --- REVISION: Save fold-level paired errors for bootstrap CIs and effect sizes ---
# Recompute fold-level errors for every baseline row, then deduplicate to exactly
# one baseline observation per (target, GEAR, q, BAY, horizon, test_date).
# Multiple rows can exist if the same forecast date appears in more than one
# source file — we keep the first occurrence after deduplication.
baseline_folds <- fc_all %>%
  filter(model_is_baseline, !is.na(y_true), !is.na(y_pred)) %>%
  mutate(
    error_base     = y_pred - y_true,
    abs_error_base = abs(y_pred - y_true)
  ) %>%
  select(target, GEAR, q, BAY, horizon, test_date, error_base, abs_error_base) %>%
  # One row per join key — resolves many-to-many caused by duplicate source files
  distinct(target, GEAR, q, BAY, horizon, test_date, .keep_all = TRUE)

# Join onto advanced models at the fold level
fc_joined <- fc_all %>%
  filter(!model_is_baseline, !is.na(y_true), !is.na(y_pred)) %>%
  mutate(
    error     = y_pred - y_true,
    abs_error = abs(y_pred - y_true)
  ) %>%
  left_join(baseline_folds,
            by = c("target", "GEAR", "q", "BAY", "horizon", "test_date"))

OUT_JOINED_RDS_D <- file.path(dir_data_script,  "fc_joined_diversity.rds")
OUT_JOINED_RDS_R <- file.path(dir_model_script, "fc_joined_diversity.rds")
saveRDS(fc_joined, OUT_JOINED_RDS_D)
saveRDS(fc_joined, OUT_JOINED_RDS_R)

# Diagnostic: confirm join quality
n_adv        <- fc_all %>% filter(!model_is_baseline, !is.na(y_true), !is.na(y_pred)) %>% nrow()
n_joined     <- nrow(fc_joined)
n_has_base   <- sum(!is.na(fc_joined$error_base))
pct_matched  <- round(100 * n_has_base / n_joined, 1)
cat(sprintf(
  "Fold-level join: %d advanced rows, %d joined rows, %d (%.1f%%) have a matched baseline error.\n",
  n_adv, n_joined, n_has_base, pct_matched
))
if (pct_matched < 80) warning("Less than 80% of advanced rows matched a baseline — check test_date alignment.")
cat("Saved fold-level paired diversity errors to:", OUT_JOINED_RDS_D, "\n")
# ---------------------------------------------------------------------------------

cat("\nCompleted d1 consolidation successfully. Check:", OUT_DICT_CSV_D, "\n")