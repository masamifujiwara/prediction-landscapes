# =============================================================================
# Script: figS11_gear_comparison.R
#
# Figure S11 — Gear Comparison: Bag Seine vs Bay Trawl (h=1)
#
# SARIMA skill (% RMSE improvement) at h=1 comparing BS and BT for:
#   - Diversity targets: q=0, q=1, q=2 at Alpha (median across bays) and Gamma
#   - Abundance targets: 8 species at Gamma scale
#
# Layout: two panels side by side
#   Panel A: Diversity targets
#   Panel B: Abundance targets (Gamma)
#
# Encoding:
#   X-axis: Bag Seine skill
#   Y-axis: Bay Trawl skill
#   Diagonal line: 1:1 (equal skill)
#   Points above diagonal: BT > BS; below: BS > BT
#   Color: spatial scale (Alpha = open, Gamma = filled) for diversity;
#          species for abundance
#
# Width: 170 mm; Height: 90 mm
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

need <- c("ggplot2", "dplyr", "tidyr", "patchwork", "ggrepel")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(ggrepel)
})

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load data
# ─────────────────────────────────────────────────────────────────────────────
div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay Anchovy","Atlantic Croaker","Gulf Menhaden","Spot",
                   "Pinfish","White Shrimp","Blue Crab","Brown Shrimp"),
  stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Diversity: pivot BS vs BT into wide format
# ─────────────────────────────────────────────────────────────────────────────
div_wide <- div_raw %>%
  filter(model_clean == "sarima", horizon == 1) %>%
  mutate(scale = ifelse(BAY == "ALL", "Gamma", "Alpha"),
         q_label = paste0("q = ", q)) %>%
  select(scale, q_label, BAY, GEAR, skill_rmse_improve_pct) %>%
  pivot_wider(names_from = GEAR, values_from = skill_rmse_improve_pct) %>%
  filter(!is.na(BS), !is.na(BT))

# ─────────────────────────────────────────────────────────────────────────────
# 3) Abundance: pivot BS vs BT into wide format (Gamma only)
# ─────────────────────────────────────────────────────────────────────────────
abd_wide <- abd_raw %>%
  filter(model == "sarima", horizon == 1,
         series == "gamma", BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  select(common_name, GEAR, skill_rmse_improve_pct) %>%
  pivot_wider(names_from = GEAR, values_from = skill_rmse_improve_pct) %>%
  filter(!is.na(BS), !is.na(BT))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
q_cols <- c("q = 0" = "#D73027",
            "q = 1" = "#FC8D59",
            "q = 2" = "#4575B4")

sp_cols <- c(
  "Spot"             = "#E69F00",
  "Pinfish"          = "#56B4E9",
  "Atlantic Croaker" = "#009E73",
  "Bay Anchovy"      = "#CC79A7",
  "Gulf Menhaden"    = "#0072B2",
  "Blue Crab"        = "#D55E00",
  "Brown Shrimp"     = "#F0E442",
  "White Shrimp"     = "#999999"
)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey91", linewidth = 0.3),
    axis.text         = element_text(size = 6.5),
    axis.title        = element_text(size = 7),
    legend.text       = element_text(size = 6.5),
    legend.title      = element_text(size = 7, face = "bold"),
    legend.key.size   = unit(0.38, "cm"),
    plot.tag          = element_text(size = 8, face = "bold"),
    plot.margin       = margin(2, 3, 2, 3, "mm")
  )

# Axis range for both panels (shared)
all_skills <- c(div_wide$BS, div_wide$BT, abd_wide$BS, abd_wide$BT)
ax_lim <- c(floor(min(all_skills, na.rm=TRUE) / 5) * 5,
            ceiling(max(all_skills, na.rm=TRUE) / 5) * 5)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Panel A — Diversity gear comparison
# ─────────────────────────────────────────────────────────────────────────────
pA <- ggplot(div_wide, aes(x = BS, y = BT)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.5) +
  geom_point(data = filter(div_wide, scale == "Alpha"),
             aes(colour = q_label),
             shape = 1, size = 2.0, stroke = 0.7, alpha = 0.6) +
  geom_point(data = filter(div_wide, scale == "Gamma"),
             aes(colour = q_label),
             shape = 16, size = 3.0) +
  geom_text_repel(data = filter(div_wide, scale == "Gamma"),
                  aes(label = q_label, colour = q_label),
                  size = 2.0, show.legend = FALSE,
                  nudge_y = 1.5, segment.size = 0.3,
                  seed = 42) +
  scale_colour_manual(values = q_cols, name = "Diversity order") +
  scale_x_continuous(limits = ax_lim) +
  scale_y_continuous(limits = ax_lim) +
  annotate("text", x = ax_lim[1] + 1, y = ax_lim[2] - 1,
           label = "BT > BS", hjust = 0, vjust = 1,
           size = 2.0, colour = "grey40", fontface = "italic") +
  annotate("text", x = ax_lim[2] - 1, y = ax_lim[1] + 1,
           label = "BS > BT", hjust = 1, vjust = 0,
           size = 2.0, colour = "grey40", fontface = "italic") +
  labs(x = bquote("Bag Seine skill (%" * Delta * "RMSE)"),
       y = bquote("Bay Trawl skill (%" * Delta * "RMSE)"),
       tag = "A", title = "Diversity") +
  coord_equal() +
  base_th +
  theme(plot.title    = element_text(size = 7, face = "bold", hjust = 0.5),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1))

# ─────────────────────────────────────────────────────────────────────────────
# 6) Panel B — Abundance gear comparison (Gamma)
# ─────────────────────────────────────────────────────────────────────────────
pB <- ggplot(abd_wide, aes(x = BS, y = BT, colour = common_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.5) +
  geom_point(shape = 16, size = 3.0) +
  geom_text_repel(aes(label = common_name),
                  size = 1.9, show.legend = FALSE,
                  box.padding = 0.3, segment.size = 0.3,
                  seed = 42) +
  scale_colour_manual(values = sp_cols, name = "Species") +
  scale_x_continuous(limits = ax_lim) +
  scale_y_continuous(limits = ax_lim) +
  annotate("text", x = ax_lim[1] + 1, y = ax_lim[2] - 1,
           label = "BT > BS", hjust = 0, vjust = 1,
           size = 2.0, colour = "grey40", fontface = "italic") +
  annotate("text", x = ax_lim[2] - 1, y = ax_lim[1] + 1,
           label = "BS > BT", hjust = 1, vjust = 0,
           size = 2.0, colour = "grey40", fontface = "italic") +
  labs(x = bquote("Bag Seine skill (%" * Delta * "RMSE)"),
       y = bquote("Bay Trawl skill (%" * Delta * "RMSE)"),
       tag = "B", title = "Abundance (Gamma scale)") +
  coord_equal() +
  base_th +
  theme(plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
        legend.position = "none")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Assemble
# ─────────────────────────────────────────────────────────────────────────────
fig_S11 <- (pA | pB) +
  plot_layout(guides = "keep")

# ─────────────────────────────────────────────────────────────────────────────
# 8) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 95

out_tiff <- file.path(fig_dir, "Figure_S11.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S11.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S11); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S11); dev.off()

cat("\u2705 Figure S11 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
