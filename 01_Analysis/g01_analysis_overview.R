# =============================================================================
# Script:  g01_analysis_overview.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Renders Figure 1, the overview of the analysis (data, seasonal curve
#   refitted at each forecast date, forecast models and simple references,
#   out-of-sample evaluation, data diagnostics). The figure is a drawn diagram
#   kept as an SVG file; this script only converts it to the TIFF (600 dpi)
#   and PDF that the Word builder (g20) expects, at a width of 170 mm.
#
# Run after: none (the SVG source is drawn by hand)
#
# Inputs:
#   01_Analysis/figures/source/Figure_1_analysis_overview.svg (drawn by hand,
#     kept in the repository)
#
# Outputs:
#   01_Analysis/figures/Figure_1_analysis_overview.tiff and .pdf
#   01_Analysis/Logs/g1_analysis_overview/run_<timestamp>.txt
#   Earlier copies of the TIFF and PDF are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The figure sits in the Methods next to Table 1 (the models); the model box
#   names the model groups and refers to Table 1 for the full list. To change
#   the figure, edit the diagram in the SVG file, then rerun this script. The
#   canvas height follows the SVG viewBox, and the script stops if the height
#   is more than 1.231 times the width. Needs the package rsvg
#   (install.packages("rsvg")).
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

# ---- Bootstrap: find the project folder and build PATHS (same in every script)
# Looks for the marker file .project_root above, in order: the file being
# source()d (RStudio Source, source()), the file given to Rscript / R -f /
# R CMD BATCH (also under nohup), the file open in RStudio, the working folder.
local({
  a <- commandArgs(FALSE)
  here <- c(unlist(lapply(rev(sys.frames()), function(f) get0("ofile", envir = f, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", a, value = TRUE)),
            a[which(a == "-f") + 1L],
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
              tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
            getwd())
  here <- gsub("~+~", " ", here[!is.na(here) & nzchar(here)], fixed = TRUE)  # Rscript on Linux writes spaces as ~+~
  for (h in here) {
    d <- normalizePath(h, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) d <- dirname(d)
    repeat {
      if (file.exists(file.path(d, ".project_root"))) {
        assign(".PROJECT_ROOT", d, envir = globalenv())
        return(invisible())
      }
      if (dirname(d) == d) break
      d <- dirname(d)
    }
  }
  stop("Project folder not found: no .project_root above ", paste(unique(here), collapse = " or "))
})
source(file.path(.PROJECT_ROOT, "01_Analysis", "lib", "paths.R"))
PATHS <- get_paths(.PROJECT_ROOT)
setwd(PATHS$analysis_dir)
# -----------------------------------------------------------------------------
source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
if (!requireNamespace("rsvg", quietly = TRUE)) stop("Package rsvg is needed: install.packages(\"rsvg\")", call. = FALSE)
script_stem <- "g1_analysis_overview"
in_svg <- file.path(PATHS$analysis_dir, "figures", "source", "Figure_1_analysis_overview.svg")
W_MM <- 170

dirs <- script_dirs(script_stem, c("logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_svg, "the diagram source (kept by hand, not made by a script)")

# Canvas: the SVG's own width and height, scaled to W_MM.
svg_txt <- paste(readLines(in_svg, warn = FALSE), collapse = " ")
dims <- as.numeric(regmatches(svg_txt, regexec('viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_txt))[[1]][2:3])
if (length(dims) != 2L || anyNA(dims)) stop("Could not read the viewBox of ", in_svg, call. = FALSE)
H_MM <- W_MM * dims[2] / dims[1]
cat("Canvas:", W_MM, "x", round(H_MM, 1), "mm\n")

fig_dir <- file.path(PATHS$analysis_dir, "figures")
out <- file.path(fig_dir, paste0("Figure_1_analysis_overview", c(".tiff", ".pdf")))
if (H_MM / W_MM > 8.0 / 6.5) stop("Diagram height/width above 1.231", call. = FALSE)
old <- out[file.exists(out)]
if (length(old)) archive_files(old)
bm <- rsvg::rsvg(in_svg, width = round(W_MM / 25.4 * 600))
grDevices::tiff(out[1], width = W_MM / 25.4, height = H_MM / 25.4, units = "in", res = 600, compression = "lzw")
grid::grid.newpage(); grid::grid.raster(bm, interpolate = TRUE); grDevices::dev.off()
rsvg::rsvg_pdf(in_svg, out[2], width = W_MM / 25.4 * 72, height = H_MM / 25.4 * 72)
for (f in out) cat("Saved:", f, "\n")
end_log("finished")
cat("g01_analysis_overview.R completed\n")
