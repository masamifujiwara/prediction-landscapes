# =============================================================================
# File:    lib/figures.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Shared look and output for the figure scripts: one copy of the theme,
#   colors, labels and the save step, so every figure matches.
#
# Used by: g02 to g08, g11 to g17, h01, h02, h15
#
# Functions:
#   PAL, GEAR_COLS, SCALE_COLS        color palette; colors by sampling method
#                                     and by scale
#   GEAR_LAB, SCALE_LAB, CLASS_LAB    labels for sampling method, scale and
#                                     target class
#   FONT_SIZE                         base font size (points)
#   theme_paper()                     ggplot2 theme of every figure
#   save_figure()                     writes one figure as TIFF and PDF
#
# Outputs:
#   save_figure() writes 01_Analysis/figures/<name>.tiff and <name>.pdf.
#   Earlier copies are moved to 01_Analysis/figures/ARCHIVE_<timestamp>.
#
# Method notes:
#   Colors: three categorical colors checked for color-blind separation
#   (blue, orange, aqua), plus greys. Sampling method: bag seine blue, bay
#   trawl orange. Estuaries against the coast: estuaries blue, coastwide
#   orange. Text is never colored.
#   Labels follow the manuscript: "Estuaries" for the bay scale and
#   "Biodiversity" for the diversity class.
#   Canvas sizes are in millimeters, with height/width at or below
#   8.0 / 6.5 = 1.231 (checked when saving), so g20 can place every figure at
#   full text width.
#   TIFF at 600 dpi with LZW compression; PDF with cairo when available.
#   All text, including axis tick labels and legend text, is set at the base
#   size of 8 pt (journal guideline: lettering of about 8 to 12 pt).
#   Axis bounds come from the data (no fixed clips).
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

PAL <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a",
         grey_dark = "#52514e", grey_mid = "#9a9993", grey_light = "#d6d5d0")
GEAR_COLS  <- c("Bag seine" = PAL[["blue"]], "Bay trawl" = PAL[["orange"]])
SCALE_COLS <- c("Estuaries" = PAL[["blue"]], "Coastwide" = PAL[["orange"]])
GEAR_LAB   <- c(BS = "Bag seine", BT = "Bay trawl")
SCALE_LAB  <- c(bay = "Estuaries", coastwide = "Coastwide")
CLASS_LAB  <- c(diversity = "Biodiversity", abundance = "Abundance")
FONT_SIZE  <- 8

theme_paper <- function(base = FONT_SIZE) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.25),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70", linewidth = 0.3),
                   strip.text = ggplot2::element_text(size = base),
                   axis.text = ggplot2::element_text(size = base),
                   legend.text = ggplot2::element_text(size = base),
                   legend.position = "bottom", legend.key.size = ggplot2::unit(3.5, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0), plot.margin = ggplot2::margin(3, 4, 2, 2))
}

# Writes p as TIFF and PDF under 01_Analysis/figures/<name>.*, archiving any
# earlier copies first. Stops if the canvas is taller than 1.231 x its width.
save_figure <- function(p, name, w_mm, h_mm, fig_dir = file.path(PATHS$analysis_dir, "figures")) {
  if (h_mm / w_mm > 8.0 / 6.5) stop(name, ": height/width ", round(h_mm / w_mm, 3), " is above 1.231", call. = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(fig_dir, paste0(name, c(".tiff", ".pdf")))
  old <- out[file.exists(out)]
  if (length(old)) {
    arch <- file.path(fig_dir, paste0("ARCHIVE_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(arch, showWarnings = FALSE)
    if (!all(file.rename(old, file.path(arch, basename(old))))) stop("Could not archive ", paste(old, collapse = ", "), call. = FALSE)
  }
  grDevices::tiff(out[1], width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 600, compression = "lzw")
  print(p); grDevices::dev.off()
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  pdf_dev(out[2], width = w_mm / 25.4, height = h_mm / 25.4)
  print(p); grDevices::dev.off()
  for (f in out) cat("Saved:", f, "\n")
  invisible(out)
}
