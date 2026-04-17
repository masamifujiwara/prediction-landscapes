### =============================================================================
### Script: f1_export_rds_to_csv.R
### Purpose: Load necessary intermediate and final .rds files and export them 
###          as .csv files. Handles standard data frames, forecast list objects, 
###          and extracts the R2 table from the seasonality_fits list.
### =============================================================================
rm(list = ls())

### ─────────────────────────────────────────────────────────────────────────────
### 0) Bootstrap: setwd() to script folder + load PATHS
### ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  ofile <- tryCatch(sys.frames()[]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

# Find project root and load PATHS
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)  # loads PATHS
stopifnot(exists("PATHS"))

### ─────────────────────────────────────────────────────────────────────────────
### 1) Packages & Output Directory
### ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

script_stem <- "f1_export_rds_to_csv"
data_out_dir <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

### ─────────────────────────────────────────────────────────────────────────────
### 2) Export Standard Data Frames
### ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Exporting Standard Data Frames ---\n")

simple_files <- list(
  alpha_div_deseas   = file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds"),
  gamma_div_deseas   = file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "gamma_q012_deseas.rds"),
  seasonal_map       = file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "seasonal_map.rds"),
  alpha_abund_deseas = file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds"),
  gamma_abund_deseas = file.path(PATHS$analysis_dir, "Data", "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds")
)

for (name in names(simple_files)) {
  if (file.exists(simple_files[[name]])) {
    df <- readRDS(simple_files[[name]])
    out_csv <- file.path(data_out_dir, paste0(name, ".csv"))
    write_csv(df, out_csv)
    cat("✅ Saved:", out_csv, "\n")
  } else {
    cat("⚠️ Missing:", simple_files[[name]], "\n")
  }
}

### ─────────────────────────────────────────────────────────────────────────────
### 3) Extract R² from Seasonality Fits (figure_S1_v7 logic)
### ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Extracting Seasonality R-Squared Table ---\n")

seas_fits_file <- file.path(PATHS$analysis_dir, "Data", "c1_seasonality_characterization", "seasonality_fits.RDS")

if (file.exists(seas_fits_file)) {
  seas_fits <- readRDS(seas_fits_file)
  
  extract_r2 <- function(fits) {
    if (is.data.frame(fits)) {
      r2_col <- intersect(c("var_expl", "r2", "R2", "r_squared"), names(fits))
      q_col  <- intersect(c("q", "Q"), names(fits))
      if (length(r2_col) > 0 && length(q_col) > 0) {
        return(fits %>% select(q = all_of(q_col), r2 = all_of(r2_col)))
      }
    }
    if (is.list(fits) && !is.data.frame(fits)) {
      rows <- lapply(names(fits), function(nm) {
        el <- fits[[nm]]
        r2 <- NULL
        if (is.list(el)) {
          if ("metrics" %in% names(el) && is.data.frame(el$metrics)) {
            r2 <- el$metrics$var_expl
          } else if ("var_expl" %in% names(el)) {
            r2 <- el$var_expl
          }
        }
        if (is.null(r2)) return(NULL)
        data.frame(key = nm, r2 = r2, stringsAsFactors = FALSE)
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows) > 0) return(bind_rows(rows))
    }
    return(NULL)
  }
  
  r2_table <- extract_r2(seas_fits)
  if (!is.null(r2_table)) {
    out_csv <- file.path(data_out_dir, "seasonality_r2_values.csv")
    write_csv(r2_table, out_csv)
    cat("✅ Saved:", out_csv, "\n")
  } else {
    cat("⚠️ Could not extract R2 table from structure.\n")
  }
} else {
  cat("⚠️ Missing:", seas_fits_file, "\n")
}

### ─────────────────────────────────────────────────────────────────────────────
### 4) Extract Prediction Data Frames from Model Object Lists
### ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Extracting Forecast Predictions ---\n")

model_files <- list(
  abund_preds_baseline = file.path(PATHS$analysis_dir, "Data", "e5_baseline_forecasts_abundance", "baseline_forecasts_abundance.rds"),
  abund_preds_ssm      = file.path(PATHS$analysis_dir, "Data", "e6_state_space_forecasts_abundance", "ssm_forecasts_abundance.rds"),
  abund_preds_arimax   = file.path(PATHS$analysis_dir, "Data", "e7_arimax_env_forecasts_abundance", "arimax_forecasts_abundance.rds"),
  abund_preds_ssmxreg  = file.path(PATHS$analysis_dir, "Data", "e8_ssm_xreg_forecasts_abundance", "ssm_xreg_forecasts_abundance.rds")
)

for (name in names(model_files)) {
  if (file.exists(model_files[[name]])) {
    obj <- readRDS(model_files[[name]])
    # The models save a list, the data frame is named "predictions"
    if ("predictions" %in% names(obj) && is.data.frame(obj$predictions)) {
      out_csv <- file.path(data_out_dir, paste0(name, ".csv"))
      write_csv(obj$predictions, out_csv)
      cat("✅ Saved:", out_csv, "\n")
    } else {
      cat("⚠️ File loaded but no '$predictions' data.frame found:", model_files[[name]], "\n")
    }
  } else {
    cat("⚠️ Missing:", model_files[[name]], "\n")
  }
}

cat("\n🎉 f1_export_rds_to_csv.R completed successfully. All available CSVs are in:", data_out_dir, "\n")