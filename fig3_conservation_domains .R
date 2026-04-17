# =============================================================================
# Script: fig3_conservation_domains.R  (v3)
#
# Key change: all cell text trimmed to max 3 short lines per cell.
# Cell heights increased slightly. Figure height 160 mm.
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

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
library(ggplot2)

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Layout
# ─────────────────────────────────────────────────────────────────────────────
X_LAB <- 0.60;  W_LAB <- 1.10
X_D1  <- 2.20;  X_D2  <- 3.80;  X_D3  <- 5.40;  W_D <- 1.45

H_HDR <- 0.72   # header row
H_R1  <- 1.30   # characteristics row
H_R2  <- 1.10   # conservation strategy
H_R3  <- 0.90   # management philosophy (short — 2 lines)
H_R4  <- 1.10   # decision tools

# y-centres (header at top, rows below)
Y_HDR <- 5.20
Y_R1  <- Y_HDR - H_HDR/2 - H_R1/2 - 0.06
Y_R2  <- Y_R1  - H_R1/2  - H_R2/2 - 0.06
Y_R3  <- Y_R2  - H_R2/2  - H_R3/2 - 0.06
Y_R4  <- Y_R3  - H_R3/2  - H_R4/2 - 0.06

# ─────────────────────────────────────────────────────────────────────────────
# Colors
# ─────────────────────────────────────────────────────────────────────────────
col_D1_hdr <- "#C0392B"; col_D1_fill <- "#FADBD8"
col_D2_hdr <- "#D68910"; col_D2_fill <- "#FCF3CF"
col_D3_hdr <- "#1E8449"; col_D3_fill <- "#D5F5E3"
col_lab_hdr <- "#808B96"; col_lab_fill <- "#EBF5FB"
txt_hdr  <- "white"
txt_body <- "grey15"
txt_lab  <- "grey20"

# ─────────────────────────────────────────────────────────────────────────────
# Cell content — MAX 3 short lines per cell
# ─────────────────────────────────────────────────────────────────────────────

# Row 1: Characteristics
r1_d1 <- "Low predictability\nRare species (q = 0)\nLocal scale"
r1_d2 <- "Moderate predictability\nMixed community (q = 1)\nLocal to regional scale"
r1_d3 <- "High predictability\nDominant species (q = 2)\nRegional scale"

# Row 2: Conservation strategy
r2_d1 <- "Endangered species mgmt.\nHotspot conservation\nLocal richness monitoring"
r2_d2 <- "Adaptive monitoring\nScenario-based planning\nIterative model updating"
r2_d3 <- "Ecosystem-based mgmt.\nHarvest management\nClimate adaptation"

# Row 3: Management philosophy (2 lines only)
r3_d1 <- "Precautionary & diagnostic\n(protect, buffer, characterize)"
r3_d2 <- "Adaptive & uncertainty-aware\n(monitor, update, scenario-plan)"
r3_d3 <- "Forecast-based & quantitative\n(predict, target, optimize)"

# Row 4: Decision tools
r4_d1 <- "Habitat protection\nPopulation buffering\nLTRE / demographic analysis"
r4_d2 <- "Systematic monitoring\nForecast uncertainty bounds\nScenario projections"
r4_d3 <- "Predictive time-series models\nBiologically defensible targets\nSpatial planning tools"

# Row labels
lab_r1 <- "Characteristics"
lab_r2 <- "Conservation\nstrategy"
lab_r3 <- "Management\nphilosophy"
lab_r4 <- "Decision\ntools"

# ─────────────────────────────────────────────────────────────────────────────
# Build tile and text data frames
# ─────────────────────────────────────────────────────────────────────────────
mk_tile <- function(x, y, w, h, fill)
  data.frame(x=x, y=y, w=w, h=h, fill=fill, stringsAsFactors=FALSE)

tiles <- rbind(
  # Header
  mk_tile(X_LAB, Y_HDR, W_LAB, H_HDR, col_lab_hdr),
  mk_tile(X_D1,  Y_HDR, W_D,   H_HDR, col_D1_hdr),
  mk_tile(X_D2,  Y_HDR, W_D,   H_HDR, col_D2_hdr),
  mk_tile(X_D3,  Y_HDR, W_D,   H_HDR, col_D3_hdr),
  # Row 1
  mk_tile(X_LAB, Y_R1, W_LAB, H_R1, col_lab_fill),
  mk_tile(X_D1,  Y_R1, W_D,   H_R1, col_D1_fill),
  mk_tile(X_D2,  Y_R1, W_D,   H_R1, col_D2_fill),
  mk_tile(X_D3,  Y_R1, W_D,   H_R1, col_D3_fill),
  # Row 2
  mk_tile(X_LAB, Y_R2, W_LAB, H_R2, col_lab_fill),
  mk_tile(X_D1,  Y_R2, W_D,   H_R2, col_D1_fill),
  mk_tile(X_D2,  Y_R2, W_D,   H_R2, col_D2_fill),
  mk_tile(X_D3,  Y_R2, W_D,   H_R2, col_D3_fill),
  # Row 3
  mk_tile(X_LAB, Y_R3, W_LAB, H_R3, col_lab_fill),
  mk_tile(X_D1,  Y_R3, W_D,   H_R3, col_D1_fill),
  mk_tile(X_D2,  Y_R3, W_D,   H_R3, col_D2_fill),
  mk_tile(X_D3,  Y_R3, W_D,   H_R3, col_D3_fill),
  # Row 4
  mk_tile(X_LAB, Y_R4, W_LAB, H_R4, col_lab_fill),
  mk_tile(X_D1,  Y_R4, W_D,   H_R4, col_D1_fill),
  mk_tile(X_D2,  Y_R4, W_D,   H_R4, col_D2_fill),
  mk_tile(X_D3,  Y_R4, W_D,   H_R4, col_D3_fill)
)

mk_txt <- function(x, y, label, size, fontface, colour)
  data.frame(x=x, y=y, label=label, size=size,
             fontface=fontface, colour=colour, stringsAsFactors=FALSE)

SZ_HDR  <- 2.8
SZ_LAB  <- 2.3
SZ_BODY <- 2.15

txt <- rbind(
  # Headers
  mk_txt(X_LAB, Y_HDR, "",                     SZ_HDR, "bold",  txt_hdr),
  mk_txt(X_D1,  Y_HDR, "Domain I\nStochastic",      SZ_HDR, "bold",  txt_hdr),
  mk_txt(X_D2,  Y_HDR, "Domain II\nTransitional",   SZ_HDR, "bold",  txt_hdr),
  mk_txt(X_D3,  Y_HDR, "Domain III\nStructured", SZ_HDR, "bold",  txt_hdr),
  # Row labels
  mk_txt(X_LAB, Y_R1, lab_r1, SZ_LAB, "bold", txt_lab),
  mk_txt(X_LAB, Y_R2, lab_r2, SZ_LAB, "bold", txt_lab),
  mk_txt(X_LAB, Y_R3, lab_r3, SZ_LAB, "bold", txt_lab),
  mk_txt(X_LAB, Y_R4, lab_r4, SZ_LAB, "bold", txt_lab),
  # Content
  mk_txt(X_D1, Y_R1, r1_d1, SZ_BODY, "plain", txt_body),
  mk_txt(X_D2, Y_R1, r1_d2, SZ_BODY, "plain", txt_body),
  mk_txt(X_D3, Y_R1, r1_d3, SZ_BODY, "plain", txt_body),
  mk_txt(X_D1, Y_R2, r2_d1, SZ_BODY, "plain", txt_body),
  mk_txt(X_D2, Y_R2, r2_d2, SZ_BODY, "plain", txt_body),
  mk_txt(X_D3, Y_R2, r2_d3, SZ_BODY, "plain", txt_body),
  mk_txt(X_D1, Y_R3, r3_d1, SZ_BODY, "plain", txt_body),
  mk_txt(X_D2, Y_R3, r3_d2, SZ_BODY, "plain", txt_body),
  mk_txt(X_D3, Y_R3, r3_d3, SZ_BODY, "plain", txt_body),
  mk_txt(X_D1, Y_R4, r4_d1, SZ_BODY, "plain", txt_body),
  mk_txt(X_D2, Y_R4, r4_d2, SZ_BODY, "plain", txt_body),
  mk_txt(X_D3, Y_R4, r4_d3, SZ_BODY, "plain", txt_body)
)

# ─────────────────────────────────────────────────────────────────────────────
# Arrow position
# ─────────────────────────────────────────────────────────────────────────────
arrow_y <- Y_R4 - H_R4/2 - 0.18

# ─────────────────────────────────────────────────────────────────────────────
# Build plot
# ─────────────────────────────────────────────────────────────────────────────
fig3 <- ggplot() +
  geom_tile(data = tiles,
            aes(x = x, y = y, width = w, height = h, fill = fill),
            colour = "white", linewidth = 0.6) +
  scale_fill_identity() +
  geom_text(data = txt,
            aes(x = x, y = y, label = label,
                size = size, fontface = fontface, colour = colour),
            lineheight = 1.05, vjust = 0.5, hjust = 0.5) +
  scale_colour_identity() +
  scale_size_identity() +
  # Arrow
  annotate("segment",
           x = X_D1 - W_D/2 + 0.05, xend = X_D3 + W_D/2 - 0.05,
           y = arrow_y, yend = arrow_y,
           colour = "grey35", linewidth = 0.7,
           arrow = arrow(length = unit(0.20, "cm"),
                         type = "closed", ends = "last")) +
  annotate("text",
           x = (X_D1 + X_D3) / 2, y = arrow_y - 0.19,
           label = "Increasing predictability",
           size = 2.4, colour = "grey30", fontface = "italic", hjust = 0.5) +
  coord_cartesian(
    xlim = c(0.02, 6.18),
    ylim = c(arrow_y - 0.36, Y_HDR + H_HDR/2 + 0.08),
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(3, 3, 6, 3, "mm"))

# ─────────────────────────────────────────────────────────────────────────────
# Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 150

out_tiff <- file.path(fig_dir, "Figure_3.tiff")
out_pdf  <- file.path(fig_dir, "Figure_3.pdf")

tiff(out_tiff, width = mm_to_in(W), height = mm_to_in(H),
     units = "in", res = 600, compression = "lzw")
print(fig3); dev.off()

pdf(out_pdf, width = mm_to_in(W), height = mm_to_in(H))
print(fig3); dev.off()

cat("\u2705 Figure 3 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
