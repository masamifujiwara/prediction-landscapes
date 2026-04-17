# =============================================================================
# Script: tabS2_harmonic_regression_fits.R
#
# Appendix S18 (Table S2) -- Harmonic Regression Model Fits
#
# Output: Word document (.docx) formatted for PDF conversion
#   - Appendix S18 as section title (bold 12pt)
#   - Table S2 as table caption (12pt double-spaced)
#   - Footnote below table (10pt)
#   - Landscape orientation via body_end_section_landscape()
#   - Table body 8pt, header 9pt (consistent with Tables S3-S5)
#
# Requires: dplyr, officer, flextable
# Output: 01_Analysis/tables/Table_S2.docx
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

need <- c("dplyr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(officer); library(flextable)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Load seasonality summary
# =============================================================================
seas <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "c1_seasonality_characterization", "seasonality_summary.csv"),
  stringsAsFactors = FALSE)

# =============================================================================
# 2) Format table
# =============================================================================
tab <- seas %>%
  arrange(GEAR, q) %>%
  mutate(
    Gear                    = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    `Diversity order`       = q_label,
    `N observations`        = formatC(n_obs, format = "d", big.mark = ","),
    `Seasonal amplitude`    = round(seasonal_amplitude, 3),
    `Peak-to-trough range`  = round(seasonal_peak_to_trough, 3),
    `Variance explained (%)` = round(variance_explained_season * 100, 1),
    `Residual ACF (lag 1)`  = round(resid_acf_lag1, 3),
    `Ljung-Box p (lag 12)`  = ifelse(
      ljungbox_p_lag12 < 0.001, "< 0.001",
      sprintf("%.3f", ljungbox_p_lag12))
  ) %>%
  select(Gear, `Diversity order`, `N observations`,
         `Seasonal amplitude`, `Peak-to-trough range`,
         `Variance explained (%)`,
         `Residual ACF (lag 1)`,
         `Ljung-Box p (lag 12)`)

# =============================================================================
# 3) Build flextable
# =============================================================================
ft <- flextable(tab) %>%
  # ── Column widths ─────────────────────────────────────────────────────────
  width(j = 1, width = 0.90) %>%
  width(j = 2, width = 1.00) %>%
  width(j = 3, width = 0.90) %>%
  width(j = 4, width = 0.95) %>%
  width(j = 5, width = 1.00) %>%
  width(j = 6, width = 1.00) %>%
  width(j = 7, width = 0.95) %>%
  width(j = 8, width = 0.95) %>%
  # ── Font and size ─────────────────────────────────────────────────────────
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 8, part = "body") %>%
  fontsize(size = 9, part = "header") %>%
  bold(part = "header") %>%
  # ── Alignment ─────────────────────────────────────────────────────────────
  align(j = 1:2, align = "left",   part = "all") %>%
  align(j = 3:8, align = "center", part = "all") %>%
  # ── Borders ───────────────────────────────────────────────────────────────
  border_remove() %>%
  hline_top(border    = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.0), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  # ── Merge gear column ─────────────────────────────────────────────────────
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top", part = "body") %>%
  # ── Shading for Bag Seine rows ────────────────────────────────────────────
  bg(i = 1:3, bg = "#F0F0F0", part = "body") %>%
  # ── Row height and padding ────────────────────────────────────────────────
  height_all(height = 0.20) %>%
  hrule(rule = "atleast") %>%
  padding(padding.top = 1, padding.bottom = 1,
          padding.left = 3, padding.right = 3, part = "all") %>%
  line_spacing(space = 1, part = "all")

# =============================================================================
# 4) Caption and footnote text
# =============================================================================
caption_txt <- paste0(
  "Table S2. Harmonic regression deseasonalization model fits for community ",
  "diversity targets. Models incorporate annual and semi-annual harmonic ",
  "components fitted to the full time series for each gear type and diversity ",
  "order. Seasonal amplitude is the fitted annual harmonic amplitude in units ",
  "of effective species. Peak-to-trough range is the difference between the ",
  "maximum and minimum fitted seasonal values. Variance explained is the ",
  "proportion of total time-series variance attributable to the fitted seasonal ",
  "component. Residual ACF at lag 1 is the first-order autocorrelation of the ",
  "deseasonalized residuals. Ljung-Box p-value is computed at lag 12 and tests ",
  "for residual serial correlation. All models were fitted coastwide (gamma ",
  "scale); abundance targets were deseasonalized using bay-specific harmonic ",
  "fits (not shown). Bag Seine rows are shaded grey; Bay Trawl rows are white."
)

footnote_txt <- paste0(
  "Abbreviations: BS = Bag Seine; BT = Bay Trawl; ",
  "ACF = autocorrelation function; ",
  "q = Hill diversity order (q = 0: species richness; q = 1: Shannon diversity; ",
  "q = 2: dominance-weighted diversity)."
)

# =============================================================================
# 5) Assemble Word document (landscape)
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
  #       "Appendix S18",
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
  
  # End this block as one landscape continuous section
  body_end_block_section(value = landscape_section)
# =============================================================================
# 6) Save
# =============================================================================
out_docx <- file.path(tab_dir, "Table_S2.docx")
print(doc, target = out_docx)
cat("\u2705 Table S2 (Appendix S18) saved:", out_docx, "\n")

cat("\n--- Table S2 preview ---\n")
print(as.data.frame(tab), row.names = FALSE)
