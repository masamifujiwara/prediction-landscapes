# =============================================================================
# Script: tabS4_mae_results.R
#
# Appendix S20 — Complete MAE Results: All Model Classes, All Targets
#
# Excel output with two sheets:
#   Sheet 1 "About": Table title, column descriptions, abbreviations, notes
#   Sheet 2 "Data":  Full MAE table with frozen header, alternating shading
#
# Requires: openxlsx package
# Output: 01_Analysis/tables/Table_S4_MAE.xlsx
# =============================================================================

rm(list = ls())

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
  normalizePath(getwd(), winslash = "/")
}
script_dir <- get_script_dir()
setwd(script_dir)
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
load(file.path(project_dir, "01_Analysis", "PATHS.RData"))

need <- c("dplyr", "tidyr", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(openxlsx)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Lookups
# ─────────────────────────────────────────────────────────────────────────────
bay_labels <- c(
  "2"="Galveston Bay", "3"="Matagorda Bay", "4"="San Antonio Bay",
  "5"="Aransas Bay",   "6"="Corpus Christi Bay",
  "7"="Upper Laguna Madre", "8"="Lower Laguna Madre",
  "ALL"="Coastwide (Gamma)")

sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  stringsAsFactors = FALSE)

model_labels_div <- c(
  "sarima"                    = "SARIMA",
  "ssm"                       = "SSM",
  "ssm_ll_seas12"             = "SSM",
  "sarimax"                   = "ARIMAX",
  "arimax_seasonal_autoarima" = "ARIMAX",
  "ssm_ll_seas12_xreg"        = "SSM-XREG",
  "seasonal_naive"            = "Baseline")

model_labels_abd <- c(
  "sarima"             = "SARIMA",
  "ssm_ll_seas12"      = "SSM",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12_xreg" = "SSM-XREG",
  "seasonal_naive"     = "Baseline")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Load data
# ─────────────────────────────────────────────────────────────────────────────
div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Diversity RMSE — pivot to wide
# ─────────────────────────────────────────────────────────────────────────────
div_models <- c("seasonal_naive","sarima","ssm","ssm_ll_seas12",
                "sarimax","arimax_seasonal_autoarima","ssm_ll_seas12_xreg")

div_wide <- div_raw %>%
  filter(model_clean %in% div_models) %>%
  mutate(
    model_display = recode(model_clean, !!!model_labels_div),
    scale   = ifelse(BAY == "ALL", "Gamma", "Alpha"),
    estuary = recode(as.character(BAY), !!!bay_labels),
    q_label = paste0("q = ", q)
  ) %>%
  group_by(GEAR, q_label, scale, estuary, horizon, model_display) %>%
  summarise(mae_val = mean(mae, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = model_display, values_from = mae_val) %>%
  # Rename Baseline column if it exists
  rename_with(~ "Baseline_MAE", .cols = any_of("Baseline")) %>%
  mutate(
    `Target type` = "Community diversity",
    Gear    = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    Target  = q_label,
    Scale   = scale,
    Estuary = estuary,
    Horizon = paste0("h = ", horizon)
  ) %>%
  select(`Target type`, Gear, Scale, Target, Estuary, Horizon,
         `Baseline MAE` = Baseline_MAE,
         SARIMA, SSM, ARIMAX, `SSM-XREG`) %>%
  arrange(Gear, Target, Scale, Estuary, Horizon)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Abundance RMSE — pivot to wide
# ─────────────────────────────────────────────────────────────────────────────
abd_wide <- abd_raw %>%
  filter(model %in% names(model_labels_abd)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    model_display = recode(model, !!!model_labels_abd),
    scale   = ifelse(series == "gamma", "Gamma", "Alpha"),
    estuary = recode(as.character(BAY), !!!bay_labels)
  ) %>%
  group_by(GEAR, common_name, scale, estuary, horizon, model_display) %>%
  summarise(mae_val = mean(mae, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = model_display, values_from = mae_val) %>%
  rename_with(~ "Baseline_MAE", .cols = any_of("Baseline")) %>%
  mutate(
    `Target type` = "Population abundance",
    Gear    = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    Target  = common_name,
    Scale   = scale,
    Estuary = estuary,
    Horizon = paste0("h = ", horizon)
  ) %>%
  select(`Target type`, Gear, Scale, Target, Estuary, Horizon,
         `Baseline MAE` = Baseline_MAE,
         SARIMA, SSM, ARIMAX, `SSM-XREG`) %>%
  arrange(Gear, Target, Scale, Estuary, Horizon)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Combine and round
# ─────────────────────────────────────────────────────────────────────────────
tab_data <- bind_rows(div_wide, abd_wide) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("Data dimensions:", nrow(tab_data), "rows x", ncol(tab_data), "cols\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Sheet 1 — About this table
# ─────────────────────────────────────────────────────────────────────────────
about_data <- data.frame(
  Item = c(
    "Table title",
    "",
    "COLUMN DESCRIPTIONS",
    "Target type",
    "Gear",
    "Scale",
    "Target",
    "Estuary",
    "Horizon",
    "Baseline MAE",
    "SARIMA",
    "SSM",
    "ARIMAX",
    "SSM-XREG",
    "",
    "ABBREVIATIONS",
    "MAE",
    "SARIMA",
    "SSM",
    "ARIMAX",
    "SSM-XREG",
    "BS",
    "BT",
    "q",
    "",
    "NOTES",
    "Units (diversity)",
    "Units (abundance)",
    "NA values",
    "Rounding"
  ),
  Description = c(
    paste0("Appendix S20. Complete MAE results for all model classes and all ",
           "forecasting targets at h = 1 month and h = 12 months. MAE values ",
           "are in the original units of each target (effective species for ",
           "diversity; log CPUE for abundance)."),
    "",
    "",
    "Biological level of organization: Community diversity or Population abundance",
    "Sampling gear: Bag Seine (littoral, active) or Bay Trawl (demersal, active)",
    "Spatial scale: Alpha = individual estuary; Gamma = coastwide aggregate",
    "Specific forecasting target: diversity order (q = 0, 1, 2) or species common name",
    "Estuary name for alpha-scale targets; 'Coastwide (Gamma)' for gamma-scale targets",
    "Forecast horizon: h = 1 (one month ahead) or h = 12 (twelve months ahead)",
    "Root mean square error of the seasonal naive baseline model",
    "Root mean square error of the SARIMA model",
    "Root mean square error of the univariate state-space model (SSM)",
    "Root mean square error of the environment-driven ARIMAX model",
    "Root mean square error of the SSM with exogenous regressors (SSM-XREG)",
    "",
    "",
    "Mean Absolute Error",
    "Seasonal Autoregressive Integrated Moving Average",
    "Univariate State-Space Model (local level + seasonal structure)",
    "Autoregressive Integrated Moving Average with eXogenous variables (temperature, salinity)",
    "State-Space Model with eXogenous Regressors (temperature, salinity)",
    "Bag Seine gear type",
    "Bay Trawl gear type",
    "Hill diversity order (q = 0: species richness; q = 1: Shannon; q = 2: dominance-weighted)",
    "",
    "",
    "Effective species (Hill numbers)",
    "Log-transformed catch per unit effort (log CPUE)",
    "NA indicates model did not converge or insufficient data for that stratum",
    "All MAE values rounded to 3 decimal places"
  ),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build workbook
# ─────────────────────────────────────────────────────────────────────────────
wb <- createWorkbook()

# ── Sheet 1: About ──────────────────────────────────────────────────────────
addWorksheet(wb, "About")
writeData(wb, "About", about_data, colNames = FALSE, startRow = 1)

# Styles for About sheet
title_style <- createStyle(
  fontSize = 11, fontName = "Times New Roman", fontColour = "#000000",
  textDecoration = "bold", wrapText = TRUE, valign = "top")
section_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman", fontColour = "#FFFFFF",
  textDecoration = "bold", fgFill = "#2F5496", halign = "left")
normal_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman", wrapText = TRUE, valign = "top")
label_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman",
  textDecoration = "bold", valign = "top")

# Apply styles
addStyle(wb, "About", title_style,  rows = 1, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", section_style, rows = c(3, 16), cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", section_style, rows = 26, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", label_style,
         rows = c(4:14, 17:24, 27:30), cols = 1, gridExpand = TRUE)
addStyle(wb, "About", normal_style,
         rows = c(4:14, 17:24, 27:30), cols = 2, gridExpand = TRUE)

# Row heights and column widths
setRowHeights(wb, "About", rows = 1, heights = 80)
setColWidths(wb, "About", cols = 1, widths = 22)
setColWidths(wb, "About", cols = 2, widths = 80)

# ── Sheet 2: Data ───────────────────────────────────────────────────────────
addWorksheet(wb, "Data")
writeData(wb, "Data", tab_data, startRow = 1, colNames = TRUE)

# Header style
hdr_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman",
  textDecoration = "bold", fgFill = "#2F5496",
  fontColour = "#FFFFFF", halign = "center",
  border = "Bottom", borderColour = "#000000")
addStyle(wb, "Data", hdr_style, rows = 1, cols = 1:ncol(tab_data),
         gridExpand = TRUE)

# Alternating row shading by Target type group
# Find row indices for each group
n_rows <- nrow(tab_data)
div_rows <- which(tab_data$`Target type` == "Community diversity") + 1
abd_rows <- which(tab_data$`Target type` == "Population abundance") + 1

shade_style <- createStyle(fgFill = "#EEF2F7",
                            fontSize = 9, fontName = "Times New Roman")
plain_style <- createStyle(fgFill = "#FFFFFF",
                            fontSize = 9, fontName = "Times New Roman")
center_style_shade <- createStyle(fgFill = "#EEF2F7", halign = "center",
                                   fontSize = 9, fontName = "Times New Roman")
center_style_plain <- createStyle(fgFill = "#FFFFFF", halign = "center",
                                   fontSize = 9, fontName = "Times New Roman")

# Shade diversity rows (light blue)
if (length(div_rows) > 0) {
  addStyle(wb, "Data", shade_style, rows = div_rows, cols = 1:4,
           gridExpand = TRUE)
  addStyle(wb, "Data", center_style_shade, rows = div_rows, cols = 5:ncol(tab_data),
           gridExpand = TRUE)
}
# Plain for abundance rows
if (length(abd_rows) > 0) {
  addStyle(wb, "Data", plain_style, rows = abd_rows, cols = 1:4,
           gridExpand = TRUE)
  addStyle(wb, "Data", center_style_plain, rows = abd_rows, cols = 5:ncol(tab_data),
           gridExpand = TRUE)
}

# Column widths
setColWidths(wb, "Data", cols = 1,  widths = 20)  # Target type
setColWidths(wb, "Data", cols = 2,  widths = 11)  # Gear
setColWidths(wb, "Data", cols = 3,  widths = 9)   # Scale
setColWidths(wb, "Data", cols = 4,  widths = 20)  # Target
setColWidths(wb, "Data", cols = 5,  widths = 20)  # Estuary
setColWidths(wb, "Data", cols = 6,  widths = 9)   # Horizon
setColWidths(wb, "Data", cols = 7,  widths = 12)  # Baseline MAE
setColWidths(wb, "Data", cols = 8,  widths = 10)  # SARIMA
setColWidths(wb, "Data", cols = 9,  widths = 10)  # SSM
setColWidths(wb, "Data", cols = 10, widths = 10)  # ARIMAX
setColWidths(wb, "Data", cols = 11, widths = 10)  # SSM-XREG

# Freeze header row
freezePane(wb, "Data", firstRow = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Save
# ─────────────────────────────────────────────────────────────────────────────
out_xlsx <- file.path(tab_dir, "Table_S4_MAE.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("\u2705 Table S3 saved:", out_xlsx, "\n")
cat("Rows:", nrow(tab_data), "| Cols:", ncol(tab_data), "\n")
