# =============================================================================
# Script: tabS5_dm_results.R
#
# Appendix S21 — Diebold-Mariano Test Results: Population Abundance Targets
#
# DM tests compare each advanced model vs seasonal naive baseline.
# Available for: SARIMA, SARIMAX, SSM, SSM-XREG
# Targets: all abundance strata (8 species x BS/BT x Alpha 7 bays + Gamma)
# Horizons: h=1 and h=12
#
# Excel output with two sheets:
#   Sheet 1 "About": Table title, column descriptions, abbreviations, notes
#   Sheet 2 "Data":  DM test results table
#
# Output: 01_Analysis/tables/Table_S5_DM.xlsx
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

need <- c("dplyr", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(openxlsx)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Lookups
# ─────────────────────────────────────────────────────────────────────────────
bay_labels <- c(
  "2"="Galveston Bay",     "3"="Matagorda Bay",
  "4"="San Antonio Bay",   "5"="Aransas Bay",
  "6"="Corpus Christi Bay","7"="Upper Laguna Madre",
  "8"="Lower Laguna Madre","ALL"="Coastwide (Gamma)")

sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  stringsAsFactors = FALSE)

model_labels <- c(
  "sarima"             = "SARIMA",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12"      = "SSM",
  "ssm_ll_seas12_xreg" = "SSM-XREG")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Load and prepare data
# ─────────────────────────────────────────────────────────────────────────────
abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

tab_data <- abd_raw %>%
  filter(model %in% names(model_labels),
         !is.na(dm_statistic)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    Model   = recode(model, !!!model_labels),
    Gear    = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    Scale   = ifelse(series == "gamma", "Gamma", "Alpha"),
    Estuary = recode(as.character(BAY), !!!bay_labels),
    Horizon = paste0("h = ", horizon),
    `DM statistic`  = round(dm_statistic, 3),
    `p-value`       = ifelse(dm_p_value < 0.001, "< 0.001",
                             sprintf("%.4f", dm_p_value)),
    `Significant (p < 0.05)` = ifelse(is_significantly_better, "Yes", "No"),
    `Cohen's d`     = round(cohens_d, 3),
    `Skill (%ΔRMSE)` = round(skill_rmse_improve_pct, 2)
  ) %>%
  select(
    Model, Gear, Scale, Species = common_name, Estuary, Horizon,
    `Skill (%ΔRMSE)`,
    `DM statistic`, `p-value`,
    `Significant (p < 0.05)`,
    `Cohen's d`
  ) %>%
  arrange(Model, Gear, Scale, Species, Estuary, Horizon)

cat("Data dimensions:", nrow(tab_data), "rows x", ncol(tab_data), "cols\n")
cat("Models:", paste(unique(tab_data$Model), collapse=", "), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Sheet 1 — About
# ─────────────────────────────────────────────────────────────────────────────
about_data <- data.frame(
  Item = c(
    "Table title",
    "",
    "COLUMN DESCRIPTIONS",
    "Model",
    "Gear",
    "Scale",
    "Species",
    "Estuary",
    "Horizon",
    "Skill (%ΔRMSE)",
    "DM statistic",
    "p-value",
    "Significant (p < 0.05)",
    "Cohen's d",
    "",
    "ABBREVIATIONS",
    "DM test",
    "SARIMA",
    "SSM",
    "ARIMAX",
    "SSM-XREG",
    "RMSE",
    "",
    "NOTES",
    "Alternative hypothesis",
    "DM test implementation",
    "Cohen's d interpretation",
    "NA values",
    "Diversity targets"
  ),
  Description = c(
    paste0("Appendix S21. Diebold-Mariano (DM) paired forecast accuracy test ",
           "results for all population abundance forecasting strata. Each row ",
           "represents one model compared against the seasonal naive baseline ",
           "for one species-gear-scale-estuary-horizon combination. The DM test ",
           "evaluates whether the advanced model produces statistically ",
           "significantly smaller squared forecast errors than the baseline."),
    "",
    "",
    "Advanced forecasting model compared against the seasonal naive baseline",
    "Sampling gear type: Bag Seine or Bay Trawl",
    "Spatial scale: Alpha = individual estuary; Gamma = coastwide aggregate",
    "Common name of the target species",
    "Estuary name for alpha-scale; 'Coastwide (Gamma)' for gamma-scale",
    "Forecast horizon: h = 1 (one month ahead) or h = 12 (twelve months ahead)",
    "Percentage improvement in RMSE over the seasonal naive baseline",
    "Diebold-Mariano test statistic (negative = advanced model better)",
    "One-sided p-value testing whether advanced model error < baseline error",
    "Whether DM test p-value < 0.05 (Yes = statistically significant improvement)",
    "Cohen's d effect size computed on paired fold-level error differences",
    "",
    "",
    "Diebold-Mariano paired forecast accuracy test (Diebold & Mariano 1995)",
    "Seasonal Autoregressive Integrated Moving Average",
    "Univariate State-Space Model",
    "Autoregressive Integrated Moving Average with eXogenous variables",
    "State-Space Model with eXogenous Regressors",
    "Root Mean Square Error",
    "",
    "",
    "alternative = 'less' (tests whether advanced model error is less than baseline)",
    "Implemented via forecast package in R; power = 2 (squared error loss)",
    "Small: |d| < 0.2; Medium: 0.2-0.5; Large: |d| > 0.5",
    "NA indicates insufficient data or model did not converge for that stratum",
    "DM tests were not computed for community diversity targets (Wilcoxon test used instead)"
  ),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Build workbook
# ─────────────────────────────────────────────────────────────────────────────
wb <- createWorkbook()

# ── Sheet 1: About ──────────────────────────────────────────────────────────
addWorksheet(wb, "About")
writeData(wb, "About", about_data, colNames = FALSE, startRow = 1)

title_style <- createStyle(
  fontSize = 11, fontName = "Times New Roman",
  textDecoration = "bold", wrapText = TRUE, valign = "top")
section_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman",
  fontColour = "#FFFFFF", textDecoration = "bold",
  fgFill = "#2F5496", halign = "left")
label_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman",
  textDecoration = "bold", valign = "top")
normal_style <- createStyle(
  fontSize = 10, fontName = "Times New Roman",
  wrapText = TRUE, valign = "top")

addStyle(wb, "About", title_style,   rows = 1,  cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", section_style, rows = 3,  cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", section_style, rows = 16, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", section_style, rows = 24, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "About", label_style,
         rows = c(4:14, 17:22, 25:29), cols = 1, gridExpand = TRUE)
addStyle(wb, "About", normal_style,
         rows = c(4:14, 17:22, 25:29), cols = 2, gridExpand = TRUE)

setRowHeights(wb, "About", rows = 1, heights = 90)
setColWidths(wb, "About", cols = 1, widths = 25)
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

# Row shading — alternate by model
models <- unique(tab_data$Model)
shade_cols <- c("#EEF2F7", "#FFFFFF", "#E8F4E8", "#FFF8E8")
names(shade_cols) <- models

for (i in seq_along(models)) {
  row_idx <- which(tab_data$Model == models[i]) + 1
  row_style <- createStyle(
    fgFill = shade_cols[i],
    fontSize = 9, fontName = "Times New Roman")
  row_style_ctr <- createStyle(
    fgFill = shade_cols[i], halign = "center",
    fontSize = 9, fontName = "Times New Roman")
  addStyle(wb, "Data", row_style,     rows = row_idx, cols = 1:5,
           gridExpand = TRUE)
  addStyle(wb, "Data", row_style_ctr, rows = row_idx, cols = 6:ncol(tab_data),
           gridExpand = TRUE)
}

# Column widths
setColWidths(wb, "Data", cols = 1,  widths = 11)  # Model
setColWidths(wb, "Data", cols = 2,  widths = 11)  # Gear
setColWidths(wb, "Data", cols = 3,  widths = 9)   # Scale
setColWidths(wb, "Data", cols = 4,  widths = 18)  # Species
setColWidths(wb, "Data", cols = 5,  widths = 20)  # Estuary
setColWidths(wb, "Data", cols = 6,  widths = 9)   # Horizon
setColWidths(wb, "Data", cols = 7,  widths = 12)  # Skill
setColWidths(wb, "Data", cols = 8,  widths = 12)  # DM statistic
setColWidths(wb, "Data", cols = 9,  widths = 10)  # p-value
setColWidths(wb, "Data", cols = 10, widths = 18)  # Significant
setColWidths(wb, "Data", cols = 11, widths = 10)  # Cohen's d

freezePane(wb, "Data", firstRow = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Save
# ─────────────────────────────────────────────────────────────────────────────
out_xlsx <- file.path(tab_dir, "Table_S5_DM.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("\u2705 Table S5 saved:", out_xlsx, "\n")
cat("Rows:", nrow(tab_data), "| Cols:", ncol(tab_data), "\n")
