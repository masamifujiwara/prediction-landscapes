# =============================================================================
# Script: figS9_horizon_decay_abundance.R
#
# Figure S9 — Horizon Decay Curves: Abundance Targets (Gamma scale)
#
# SARIMA skill (% RMSE improvement) at h=1 and h=12 for all 8 species
# at Gamma scale, both gears (BS and BT).
#
# Layout: single panel, species as lines, gear as linetype
# Color: species (8 colors from colorblind-safe palette)
# Linetype: BS = solid, BT = dashed
# Labels: h=1 and h=12 skill values + retention % on each line
#
# Width: 170 mm; Height: 130 mm
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

need <- c("ggplot2", "dplyr", "tidyr", "ggrepel")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggrepel)
})

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load data
# ─────────────────────────────────────────────────────────────────────────────
abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay Anchovy","Atlantic Croaker","Gulf Menhaden","Spot",
                   "Pinfish","White Shrimp","Blue Crab","Brown Shrimp"),
  stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Prepare data
# ─────────────────────────────────────────────────────────────────────────────
dat <- abd_raw %>%
  filter(model == "sarima", series == "gamma", BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    gear_label = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    group_id   = paste0(common_name, "_", GEAR)
  )

# Species order by descending mean h=1 skill
sp_order <- dat %>%
  filter(horizon == 1) %>%
  group_by(common_name) %>%
  summarise(m = mean(skill_rmse_improve_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(common_name)

dat <- dat %>%
  mutate(common_name = factor(common_name, levels = sp_order))

# ─────────────────────────────────────────────────────────────────────────────
# 3) Label data frames
# ─────────────────────────────────────────────────────────────────────────────
retention_df <- dat %>%
  select(common_name, gear_label, group_id, horizon,
         skill_rmse_improve_pct) %>%
  pivot_wider(names_from = horizon, values_from = skill_rmse_improve_pct,
              names_prefix = "h") %>%
  mutate(
    retention = round(100 * h12 / h1, 1),
    ret_label = paste0(retention, "%"),
    x = 12.4,
    y = h12
  )

# Right-side labels at h=12 (species name + retention)
right_labels <- retention_df %>%
  mutate(label = paste0(common_name, " (", ret_label, ")"))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
# 8 colorblind-friendly colors (Okabe-Ito extended)
sp_colors <- c(
  "Spot"             = "#E69F00",
  "Pinfish"          = "#56B4E9",
  "Atlantic Croaker" = "#009E73",
  "Bay Anchovy"      = "#F0E442",
  "Gulf Menhaden"    = "#0072B2",
  "Blue Crab"        = "#D55E00",
  "Brown Shrimp"     = "#CC79A7",
  "White Shrimp"     = "#999999"
)

gear_lty <- c("Bag Seine" = "solid", "Bay Trawl" = "dashed")

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(colour = "grey91", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.text          = element_text(size = 6),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.margin        = margin(2, 45, 2, 3, "mm")  # wide right margin for labels
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Build figure
# ─────────────────────────────────────────────────────────────────────────────
fig_S9 <- ggplot(dat,
                 aes(x = horizon, y = skill_rmse_improve_pct,
                     colour = common_name, linetype = gear_label,
                     group = group_id)) +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = gear_label), size = 2.2) +
  # Right-side species + retention labels — repelled to avoid overlap
  geom_text_repel(data = right_labels,
                  aes(x = 12, y = y, label = label,
                      colour = common_name),
                  hjust        = 0,
                  direction    = "y",
                  nudge_x      = 0.6,
                  segment.size = 0.3,
                  segment.colour = "grey60",
                  size         = 1.9,
                  show.legend  = FALSE,
                  box.padding  = 0.15,
                  force        = 1.5,
                  xlim         = c(12.5, 20),
                  seed         = 42) +
  scale_colour_manual(values = sp_colors, name = "Species") +
  scale_linetype_manual(values = gear_lty, name = "Gear") +
  scale_shape_manual(values = c("Bag Seine" = 16, "Bay Trawl" = 17),
                     name = "Gear") +
  scale_x_continuous(breaks = c(1, 12),
                     labels = c("h = 1", "h = 12"),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Forecast horizon",
       y = bquote("SARIMA skill (%" * Delta * "RMSE)")) +
  base_th +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "vertical") +
  guides(
    colour   = guide_legend(title = "Species", nrow = 2,
                            override.aes = list(linetype = "solid",
                                                shape = NA,
                                                linewidth = 1.0)),
    linetype = guide_legend(title = "Gear", order = 1,
                            override.aes = list(colour = "grey30")),
    shape    = guide_legend(title = "Gear", order = 1,
                            override.aes = list(colour = "grey30"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 120

out_tiff <- file.path(fig_dir, "Figure_S9.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S9.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S9); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S9); dev.off()

cat("\u2705 Figure S9 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
