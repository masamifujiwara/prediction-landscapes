# =============================================================================
# Script: tabS6_pi_coverage.R
#
# Appendix S22 (Table S6) -- Empirical Prediction Interval Coverage
#
# Output: Word document (.docx) formatted for PDF conversion
#   - Appendix S22 as section title (bold 12pt)
#   - Table S6 as table caption (12pt double-spaced)
#   - Footnote below table (10pt)
#   - Landscape orientation via body_end_section_landscape()
#   - Table body 8pt, header 9pt (consistent with Tables S3-S5)
#
# Requires: dplyr, tidyr, officer, flextable
# Output: 01_Analysis/tables/Table_S6_Coverage.docx
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

need <- c("dplyr", "tidyr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(officer); library(flextable)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Lookups
# =============================================================================
sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  stringsAsFactors = FALSE)

model_labels_div <- c(
  "sarima"             = "SARIMA",
  "ssm"                = "SSM",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12_xreg" = "SSM-XREG",
  "seasonal_naive"     = "Baseline")

model_labels_abd <- c(
  "sarima"             = "SARIMA",
  "ssm_ll_seas12"      = "SSM",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12_xreg" = "SSM-XREG",
  "seasonal_naive"     = "Baseline")

model_order <- c("Baseline","SARIMA","SSM","ARIMAX","SSM-XREG")

# =============================================================================
# 2) Load data
# =============================================================================
div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

# =============================================================================
# 3) Diversity coverage (Gamma)
# =============================================================================
div_cov <- div_raw %>%
  filter(BAY == "ALL",
         model_clean %in% names(model_labels_div),
         !is.na(coverage80)) %>%
  mutate(
    Model       = recode(model_clean, !!!model_labels_div),
    Target_type = "Community diversity",
    Target      = paste0("q = ", q),
    Gear        = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    Horizon     = paste0("h = ", horizon),
    cov80       = round(coverage80 * 100, 1),
    cov95       = round(coverage95 * 100, 1)
  ) %>%
  select(Target_type, Target, Gear, Horizon, Model, cov80, cov95)

# =============================================================================
# 4) Abundance coverage (Gamma)
# =============================================================================
abd_cov <- abd_raw %>%
  filter(series == "gamma", BAY == "ALL",
         model %in% names(model_labels_abd),
         !is.na(coverage80)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    Model       = recode(model, !!!model_labels_abd),
    Target_type = "Population abundance",
    Target      = common_name,
    Gear        = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    Horizon     = paste0("h = ", horizon),
    cov80       = round(coverage80 * 100, 1),
    cov95       = round(coverage95 * 100, 1)
  ) %>%
  select(Target_type, Target, Gear, Horizon, Model, cov80, cov95)

# =============================================================================
# 5) Combine and pivot wide
# =============================================================================
all_cov <- bind_rows(div_cov, abd_cov) %>%
  mutate(Model = factor(Model, levels = model_order)) %>%
  arrange(Target_type, Gear, Target, Horizon, Model)

cov_wide <- all_cov %>%
  pivot_wider(
    names_from  = Model,
    values_from = c(cov80, cov95),
    names_glue  = "{Model}_{.value}"
  ) %>%
  select(Target_type, Target, Gear, Horizon,
         ends_with("_cov80"), ends_with("_cov95")) %>%
  rename_with(~ gsub("_cov80$", " (80%)", .x)) %>%
  rename_with(~ gsub("_cov95$", " (95%)", .x)) %>%
  arrange(Target_type, Gear, Target, Horizon)

cov_wide <- cov_wide %>%
  mutate(across(where(is.numeric), ~ paste0(.x, "%")))

names(cov_wide)[1:4] <- c("Type", "Target", "Gear", "Horizon")

cat("Dims:", nrow(cov_wide), "x", ncol(cov_wide), "\n")

# =============================================================================
# 6) Build flextable
# =============================================================================
n_tot <- ncol(cov_wide)

ft <- flextable(as.data.frame(cov_wide)) %>%
  # ── Font and size ─────────────────────────────────────────────────────────
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 8, part = "body") %>%
  fontsize(size = 9, part = "header") %>%
  bold(part = "header") %>%
  # ── Alignment ─────────────────────────────────────────────────────────────
  align(j = 1:4,     align = "left",   part = "all") %>%
  align(j = 5:n_tot, align = "center", part = "all") %>%
  # ── Borders ───────────────────────────────────────────────────────────────
  border_remove() %>%
  hline_top(border    = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.0), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  # ── Shading: community diversity rows ─────────────────────────────────────
  bg(i = which(cov_wide$Type == "Community diversity"),
     bg = "#F0F0F0", part = "body") %>%
  bg(i = which(cov_wide$Type == "Population abundance"),
     bg = "#FFFFFF", part = "body") %>%
  # ── Merge Type column ─────────────────────────────────────────────────────
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top", part = "body") %>%
  # ── Column widths ─────────────────────────────────────────────────────────
  width(j = 1,       width = 1.20) %>%
  width(j = 2,       width = 1.30) %>%
  width(j = 3,       width = 0.85) %>%
  width(j = 4,       width = 0.65) %>%
  width(j = 5:n_tot, width = 0.62) %>%
  # ── Row height and padding ────────────────────────────────────────────────
  height_all(height = 0.20) %>%
  hrule(rule = "atleast") %>%
  padding(padding.top = 1, padding.bottom = 1,
          padding.left = 3, padding.right = 3, part = "all") %>%
  line_spacing(space = 1, part = "all")

# =============================================================================
# 7) Caption and footnote text
# =============================================================================
caption_txt <- paste0(
  "Table S6. Empirical prediction interval (PI) coverage for all forecasting ",
  "models at the gamma (coastwide) scale. Values are the proportion of observed ",
  "values falling within the stated PI, expressed as a percentage. Nominal ",
  "coverage levels are 80% and 95%. ",
  "Type: biological level of organization (Community diversity or Population abundance). ",
  "Target: diversity order (q = 0, 1, 2) or species common name. ",
  "Gear: sampling gear (Bag Seine = littoral; Bay Trawl = demersal). ",
  "Horizon: forecast horizon (h = 1, one month ahead; h = 12, twelve months ahead). ",
  "Each model is shown with its 80% and 95% PI coverage in parentheses. ",
  "Community diversity rows are shaded grey; population abundance rows are white."
)

footnote_txt <- paste0(
  "Abbreviations: PI = Prediction Interval; ",
  "SARIMA = Seasonal Autoregressive Integrated Moving Average; ",
  "SSM = Univariate State-Space Model; ",
  "ARIMAX = Autoregressive Integrated Moving Average with eXogenous variables; ",
  "SSM-XREG = State-Space Model with eXogenous Regressors; ",
  "Baseline = seasonal naive model. ",
  "Notes: a well-calibrated model should achieve approximately 80% coverage at ",
  "the 80% nominal level and 95% coverage at the 95% nominal level."
)

# =============================================================================
# 8) Assemble Word document (landscape)
# =============================================================================
landscape_section <- block_section(
  prop_section(
    page_size = page_size(orient = "landscape"),
    type = "continuous"
  )
)

doc <- read_docx() %>%
  # Appendix title
  # body_add_fpar(
  #   fpar(
  #     ftext(
  #       "Appendix S22",
  #       fp_text(
  #         font.size = 12,
  #         bold = TRUE,
  #         font.family = "Times New Roman"
  #       )
  #     )
  #   ),
  #   style = "Normal"
  # ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Caption
  body_add_fpar(
    fpar(
      ftext(
        caption_txt,
        fp_text(
          font.size = 12,
          font.family = "Times New Roman"
        )
      ),
      fp_p = fp_par(line_spacing = 2)
    ),
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Table
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>%
  
  # Footnote
  body_add_fpar(
    fpar(
      ftext(
        footnote_txt,
        fp_text(
          font.size = 10,
          font.family = "Times New Roman"
        )
      )
    ),
    style = "Normal"
  ) %>%
  
  # End block as one continuous landscape section
  body_end_block_section(value = landscape_section)

# =============================================================================
# 9) Save
# =============================================================================
out_docx <- file.path(tab_dir, "Table_S6_Coverage.docx")
print(doc, target = out_docx)
cat("\u2705 Table S6 (Appendix S22) saved:", out_docx, "\n")
