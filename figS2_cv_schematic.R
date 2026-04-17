# =============================================================================
# Script: figS2_cv_schematic.R
#
# Figure S2 — Rolling-Origin Cross-Validation Schematic
#
# Conceptual diagram showing:
#   - Multiple rolling-origin folds stacked vertically
#   - Each fold: expanding training window (dark) + test points at h=1 and h=12
#   - Minimum training window = 120 months annotated
#   - Arrow showing direction of rolling (left to right)
#
# Built entirely with ggplot2 geometry — no data files needed.
# Width: 140 mm (between single and double column)
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
# Layout constants
# ─────────────────────────────────────────────────────────────────────────────
# x-axis: time (arbitrary units, 0-160)
# y-axis: fold number (1 = bottom, N = top)
# Each fold is a horizontal bar

N_FOLDS    <- 6      # number of folds to show (representative subset)
TOTAL_T    <- 155    # total time units shown
MIN_TRAIN  <- 80     # minimum training window (represents 120 months)
FOLD_STEP  <- 8      # how much origin advances per fold
BAR_H      <- 0.55   # height of each bar
GAP        <- 0.45   # gap between bars (bar_h + gap = 1.0 = fold spacing)

# Fold origins (training end points)
origins <- MIN_TRAIN + (0:(N_FOLDS-1)) * FOLD_STEP   # 80, 88, 96, 104, 112, 120

# h=1 and h=12 test points
h1_pts  <- origins + 1
h12_pts <- origins + 12

# ─────────────────────────────────────────────────────────────────────────────
# Colors
# ─────────────────────────────────────────────────────────────────────────────
col_train  <- "#4A90D9"    # blue — training window
col_future <- "#D9D9D9"    # light grey — future (not yet used)
col_h1     <- "#E69F00"    # orange — h=1 test point
col_h12    <- "#CC3311"    # red    — h=12 test point
col_minwin <- "#2C6E9E"    # darker blue for minimum window annotation

# ─────────────────────────────────────────────────────────────────────────────
# Build data frames for rectangles and points
# ─────────────────────────────────────────────────────────────────────────────
# Training bars (expand with each fold)
train_bars <- data.frame(
  xmin = 0,
  xmax = origins,
  ymin = (1:N_FOLDS) - BAR_H/2,
  ymax = (1:N_FOLDS) + BAR_H/2,
  fold = 1:N_FOLDS
)

# Future / unused bars (from origin to end of series)
future_bars <- data.frame(
  xmin = origins + 1,
  xmax = TOTAL_T,
  ymin = (1:N_FOLDS) - BAR_H/2,
  ymax = (1:N_FOLDS) + BAR_H/2,
  fold = 1:N_FOLDS
)

# Minimum training window bracket (shown on fold 1)
min_bracket <- data.frame(
  x = c(0, MIN_TRAIN),
  y = rep(0.28, 2)
)

# Test points
test_h1 <- data.frame(
  x    = h1_pts,
  y    = 1:N_FOLDS,
  fold = 1:N_FOLDS
)

test_h12 <- data.frame(
  x    = h12_pts,
  y    = 1:N_FOLDS,
  fold = 1:N_FOLDS
)

# ─────────────────────────────────────────────────────────────────────────────
# Build plot
# ─────────────────────────────────────────────────────────────────────────────
legend_x1 <- 150
legend_x2 <- 160
legend_tx <- 162
legend_pt <- 155

fig_S2 <- ggplot() +

  # Future bars (background — drawn first)
  geom_rect(data = future_bars,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = col_future, colour = "white", linewidth = 0.3) +

  # Training bars
  geom_rect(data = train_bars,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = col_train, colour = "white", linewidth = 0.3) +

  # Minimum training window bracket on lowest fold
  annotate("segment",
           x = 0, xend = MIN_TRAIN,
           y = 0.28, yend = 0.28,
           colour = col_minwin, linewidth = 0.6) +
  annotate("segment", x = 0,         xend = 0,
           y = 0.28, yend = 0.38, colour = col_minwin, linewidth = 0.6) +
  annotate("segment", x = MIN_TRAIN, xend = MIN_TRAIN,
           y = 0.28, yend = 0.38, colour = col_minwin, linewidth = 0.6) +
  annotate("text", x = MIN_TRAIN/2, y = -0.02,
           label = "Minimum training\nwindow (120 months)",
           size = 2.0, colour = col_minwin, hjust = 0.5,
           fontface = "italic", lineheight = 0.95) +

  # Vertical dashed line at minimum window boundary
  annotate("segment",
           x = MIN_TRAIN, xend = MIN_TRAIN,
           y = 0.45, yend = N_FOLDS + 0.5,
           colour = col_minwin, linewidth = 0.4, linetype = "dashed") +

  # h=1 test points
  geom_point(data = test_h1,
             aes(x=x, y=y),
             shape = 23, fill = col_h1, colour = "white",
             size = 3.0, stroke = 0.3) +

  # h=12 test points
  geom_point(data = test_h12,
             aes(x=x, y=y),
             shape = 21, fill = col_h12, colour = "white",
             size = 3.0, stroke = 0.3) +

  # Rolling arrow below bars
  annotate("segment",
           x = 5, xend = TOTAL_T - 5,
           y = -0.42, yend = -0.42,
           colour = "grey40", linewidth = 0.6,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("text", x = TOTAL_T/2, y = -0.68,
           label = "Time -> (rolling origin advances one month per fold)",
           size = 2.1, colour = "grey30", hjust = 0.5, fontface = "italic") +

  # Fold labels on left
  annotate("text",
           x = rep(-3, N_FOLDS), y = 1:N_FOLDS,
           label = paste0("Fold ", 1:N_FOLDS),
           size = 2.1, hjust = 1, colour = "grey25") +

  # "Origin" tick marks and labels above each fold
  annotate("segment",
           x = origins, xend = origins,
           y = (1:N_FOLDS) + BAR_H/2,
           yend = (1:N_FOLDS) + BAR_H/2 + 0.15,
           colour = "grey50", linewidth = 0.4) +

  # Legend in far-right margin area
  # Legend background box (lower placement)
  annotate("rect",
           xmin = 148, xmax = 195,
           ymin = 1.6, ymax = 3.6,
           fill = "white", colour = "grey70", linewidth = 0.3) +
  
  # Training window
  annotate("rect", xmin = legend_x1, xmax = legend_x2,
           ymin = 3.15, ymax = 3.45,
           fill = col_train, colour = "white") +
  annotate("text", x = legend_tx, y = 3.30,
           label = "Training window",
           size = 2.1, hjust = 0, colour = "grey15") +
  
  # Withheld data
  annotate("rect", xmin = legend_x1, xmax = legend_x2,
           ymin = 2.70, ymax = 3.00,
           fill = col_future, colour = "grey70", linewidth = 0.3) +
  annotate("text", x = legend_tx, y = 2.85,
           label = "Withheld data",
           size = 2.1, hjust = 0, colour = "grey15") +
  
  # Test: h = 1
  geom_point(data = data.frame(x = legend_pt, y = 2.40),
             aes(x = x, y = y),
             shape = 23, fill = col_h1, colour = "white",
             size = 2.8, stroke = 0.3, inherit.aes = FALSE) +
  annotate("text", x = legend_tx, y = 2.40,
           label = "Test: h = 1",
           size = 2.1, hjust = 0, colour = "grey15") +
  
  # Test: h = 12
  geom_point(data = data.frame(x = legend_pt, y = 1.95),
             aes(x = x, y = y),
             shape = 21, fill = col_h12, colour = "white",
             size = 2.8, stroke = 0.3, inherit.aes = FALSE) +
  annotate("text", x = legend_tx, y = 1.95,
           label = "Test: h = 12",
           size = 2.1, hjust = 0, colour = "grey15") +

  # Axis labels
  labs(x = NULL, y = NULL,
       title = NULL) +

  coord_cartesian(
    xlim = c(-18, TOTAL_T + 45),
    ylim = c(-0.95, N_FOLDS + 0.9),
    clip = "off"
  ) +

  theme_void(base_size = 8) +
  theme(
    plot.margin = margin(4, 18, 8, 10, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 140; H <- 80

out_tiff <- file.path(fig_dir, "Figure_S2.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S2.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S2); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S2); dev.off()

cat("\u2705 Figure S2 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
