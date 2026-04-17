# =============================================================================
# Script: make_individual_figure_docx.R
#
# Purpose:
#   For Figure_S1.tiff through Figure_S14.tiff, create one Word file per figure:
#   Figure_S1.docx, Figure_S2.docx, ..., Figure_S14.docx
#   Each Word file contains the embedded TIFF and its caption below.
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 0) Bootstrap
# ─────────────────────────────────────────────────────────────────────────────
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
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

# Figure directory from your plotting scripts
fig_dir <- file.path(PATHS$analysis_dir, "figures")
stopifnot(dir.exists(fig_dir))

# ─────────────────────────────────────────────────────────────────────────────
# 1) Packages
# ─────────────────────────────────────────────────────────────────────────────
need <- c("officer")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

# magick is optional; used only to read image aspect ratio
has_magick <- requireNamespace("magick", quietly = TRUE)

suppressPackageStartupMessages({
  library(officer)
})

# ─────────────────────────────────────────────────────────────────────────────
# 2) Captions (from Figure_Captions_Revised.docx)
# ─────────────────────────────────────────────────────────────────────────────
captions <- c(
  S1  = "Figure S1. Representative harmonic regression deseasonalization fits for community diversity and population abundance time series. Left column (panels A-D): bag seine community diversity at gamma (coastwide) scale for q = 0 (Panel A), q = 1 (Panel B), and q = 2 (Panel C), and bag seine bay anchovy log-CPUE at alpha (Galveston Bay) scale (Panel D). Right column (panels E-H): bay trawl community diversity at gamma scale for q = 0 (Panel E), q = 1 (Panel F), and q = 2 (Panel G), and bag seine bay anchovy log-CPUE at gamma (coastwide) scale (Panel H; note that Panel H shows bag seine gear, matching Panel D, rather than bay trawl gear). Within each panel, the top sub-panel shows the observed time series (grey line) and the fitted seasonal curve from harmonic regression (colored line); the bottom sub-panel shows the deseasonalized residuals used as model inputs. Orange indicates bag seine gear; blue indicates bay trawl gear.",
  
  S2  = "Figure S2. Rolling-origin cross-validation design. Six representative folds are shown as horizontal bars. Dark blue bars represent the expanding training window for each fold; light grey bars represent withheld future data. Open diamonds mark the h = 1 month test observation and filled circles mark the h = 12 month test observation for each fold origin. The dashed vertical line indicates the minimum required training window of 120 months. The rolling origin advances one month per fold.",
  
  S3  = "Figure S3. SARIMA, SSM, ARIMAX, and SSM-XREG forecasting skill (%DRMSE over seasonal naive baseline) for all community diversity strata at h = 1 month (Panel A) and h = 12 months (Panel B). Each row represents one spatial stratum (gamma scale or individual alpha-scale estuary). Rows are ordered from bottom to top: Galveston, Matagorda, San Antonio, Aransas, Corpus Christi, Upper Laguna Madre, Lower Laguna Madre, and Gamma (coastwide). Panels are faceted by diversity order (q = 0, q = 1, q = 2; rows) and gear type (Bag Seine, Bay Trawl; columns). The dotted vertical line marks zero skill improvement. Values below the clip boundary (-25%) are shown as text labels at the left panel edge. Colors and shapes distinguish model classes: SARIMA (orange circle), SSM (blue triangle), ARIMAX (green square), SSM-XREG (pink diamond).",
  
  S4  = "Figure S4. SARIMA, SSM, ARIMAX, and SSM-XREG forecasting skill (%DRMSE over seasonal naive baseline) for all population abundance targets at gamma (coastwide) scale at h = 1 month (Panel A) and h = 12 months (Panel B). Each row represents one species, ordered from bottom to top by ascending mean SARIMA skill at h = 1. Panels are faceted by gear type (Bag Seine, Bay Trawl). The dotted vertical line marks zero skill improvement. Values below the clip boundary (-25%) are shown as text labels at the left panel edge. Colors and shapes distinguish model classes as in Figure S3.",
  
  S5  = "Figure S5. SARIMA, SSM, ARIMAX, and SSM-XREG forecasting skill (%DMAE over seasonal naive baseline) for all community diversity strata at h = 1 month (Panel A) and h = 12 months (Panel B). Layout, ordering, and encoding are identical to Figure S3.",
  
  S6  = "Figure S6. SARIMA, SSM, ARIMAX, and SSM-XREG forecasting skill (%DMAE over seasonal naive baseline) for all population abundance targets at gamma (coastwide) scale at h = 1 month (Panel A) and h = 12 months (Panel B). Layout, ordering, and encoding are identical to Figure S4.",
  
  S7  = "Figure S7. Empirical prediction interval coverage for all forecasting models across diversity and abundance targets at gamma (coastwide) scale. Panel A shows 80% prediction interval coverage; Panel B shows 95% prediction interval coverage. Within each panel, results are faceted by gear type (Bag Seine, Bay Trawl; columns) and forecast horizon (h = 1 month, h = 12 months; rows). Targets on the x-axis include community diversity (q = 0, q = 1, q = 2; grey background) and eight population abundance species, ordered from left to right by ascending mean SARIMA skill at h = 1. The horizontal dashed line marks the nominal coverage level (80% in Panel A, 95% in Panel B). Colors and shapes distinguish model classes as in Figure S3.",
  
  S8  = "Figure S8. SARIMA forecasting skill (%DRMSE) at h = 1 month and h = 12 months for all community diversity strata. Each panel shows one diversity order (q = 0, q = 1, q = 2; rows) and gear type (Bag Seine, Bay Trawl; columns). Within each panel, thin grey lines with open circles represent individual alpha-scale estuaries (7 bays); the thick dark line with filled circles represents the gamma-scale (coastwide) value. Text labels indicate the skill value at h = 1 and h = 12 for the gamma-scale line, and the proportional skill retention (h = 12 skill as a percentage of h = 1 skill) is shown at the midpoint of the gamma-scale line.",
  
  S9  = "Figure S9. SARIMA forecasting skill (%DRMSE) at h = 1 month and h = 12 months for all eight population abundance species at gamma (coastwide) scale, for both gear types. Each line connects the h = 1 and h = 12 skill values for one species-gear combination. Line color distinguishes species (eight colors) and linetype distinguishes gear (solid = bag seine, dashed = bay trawl). Point shape also encodes gear (circle = bag seine, triangle = bay trawl). Species names and proportional skill retention (h = 12 as a percentage of h = 1) are labeled at the right side of each line.",
  
  S10 = "Figure S10. SARIMA forecasting skill (%DRMSE) at h = 1 month and h = 12 months for eight dominant species at gamma (coastwide) scale, faceted by gear type (Bag Seine, Bay Trawl). Each row represents one species, ordered from bottom to top by ascending mean h = 1 SARIMA skill. Blue filled circles indicate h = 1 skill; red filled triangles indicate h = 12 skill. Thin grey lines connect the two horizon values for each species. The dotted vertical line marks zero skill.",
  
  S11 = "Figure S11. Comparison of SARIMA forecasting skill (%DRMSE at h = 1 month) between bag seine (x-axis) and bay trawl (y-axis) gear types. Panel A shows community diversity targets: open circles represent individual alpha-scale estuaries and filled circles represent gamma-scale (coastwide) values, with color indicating diversity order (q = 0, q = 1, q = 2). Panel B shows population abundance targets at gamma scale: filled circles represent individual species, with color and label indicating species identity. The dashed diagonal line is the 1:1 reference; points above the diagonal indicate higher bay trawl skill and points below indicate higher bag seine skill.",
  
  S12 = "Figure S12. SARIMA forecasting skill (%DRMSE at h = 1 month) for all eight population abundance species across seven individual estuaries (alpha scale) and both gear types. Panels are faceted by gear type (Bag Seine, Bay Trawl; rows) and estuary (columns: Galveston, Matagorda, San Antonio, Aransas, Corpus Christi, Upper Laguna Madre, Lower Laguna Madre). Each row within a panel represents one species, ordered from bottom to top by ascending mean alpha-scale skill. Horizontal error bars show bootstrap 95% confidence intervals (B = 2,000 resamples). Blue filled circles indicate statistically significant improvement over the seasonal naive baseline (Diebold-Mariano test, p < 0.05); grey open circles indicate improvement that is not statistically significant (p >= 0.05). The red vertical bar (|) marks the corresponding gamma-scale (coastwide) SARIMA skill for each species and gear. The dotted vertical line marks zero skill improvement.",
  
  S13 = "Figure S13. SARIMA forecasting skill (%DRMSE at h = 12 months) for all eight population abundance species across seven individual estuaries (alpha scale) and both gear types. Layout, ordering, and encoding are identical to Figure S12, except that grey open circles indicate skill values that are not statistically significant (p >= 0.05), which at this horizon includes both near-zero and negative skill values.",
  
  S14 = "Figure S14. Comparison of SARIMA and SSM-XREG forecasting skill (%DRMSE at h = 1 month) for eight population abundance species at gamma (coastwide) scale, for bag seine (Panel A) and bay trawl (Panel B) gear types. Each row represents one species, ordered from bottom to top by ascending mean SARIMA skill. Orange filled circles indicate SARIMA skill; pink filled diamonds indicate SSM-XREG skill. The dotted vertical line marks zero skill improvement. Values below the clip boundary (-60%, wider than in Figures S3-S6 to accommodate extreme SSM-XREG failures) are shown as text labels. Pink-shaded rows in Panel A highlight the three bag seine transient littoral taxa (white shrimp, blue crab, pinfish) for which SSM-XREG skill collapses below -25%."
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Utility function: get image size in inches
# ─────────────────────────────────────────────────────────────────────────────
get_img_dims_in <- function(img_file,
                            max_width_in = 6.5,
                            max_height_in = 8.0) {
  # Default if magick is unavailable
  if (!has_magick) {
    return(list(width = max_width_in, height = max_height_in))
  }
  
  info <- magick::image_info(magick::image_read(img_file))
  img_w <- info$width[1]
  img_h <- info$height[1]
  
  aspect <- img_h / img_w
  
  width_in  <- max_width_in
  height_in <- width_in * aspect
  
  if (height_in > max_height_in) {
    height_in <- max_height_in
    width_in  <- height_in / aspect
  }
  
  list(width = width_in, height = height_in)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4) Text formatting
# ─────────────────────────────────────────────────────────────────────────────
caption_fp <- fp_text(
  font.size = 12,
  font.family = "Times New Roman"
)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Build one Word file per figure
# ─────────────────────────────────────────────────────────────────────────────
for (i in 1:14) {
  
  fig_tag   <- paste0("S", i)
  fig_name  <- paste0("Figure_", fig_tag, ".tiff")
  docx_name <- paste0("Figure_", fig_tag, ".docx")
  
  fig_file  <- file.path(fig_dir, fig_name)
  out_docx  <- file.path(fig_dir, docx_name)
  
  if (!file.exists(fig_file)) {
    warning("Skipping missing figure: ", fig_file)
    next
  }
  
  if (!fig_tag %in% names(captions)) {
    warning("Skipping because caption not found for: ", fig_tag)
    next
  }
  
  dims <- get_img_dims_in(
    img_file = fig_file,
    max_width_in = 6.5,
    max_height_in = 8.0
  )
  
  doc <- read_docx()
  
  doc <- doc |>
    body_add_par("", style = "Normal") |>
    body_add_img(
      src = fig_file,
      width = dims$width,
      height = dims$height,
      style = "centered"
    ) |>
    body_add_par("", style = "Normal") |>
    body_add_fpar(
      fpar(
        ftext(captions[[fig_tag]], prop = caption_fp)
      )
    )
  
  print(doc, target = out_docx)
  message("Saved: ", out_docx)
}

cat("\nDone. Individual Word files were saved in:\n", fig_dir, "\n")