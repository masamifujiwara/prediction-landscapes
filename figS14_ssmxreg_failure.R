# =============================================================================
# Script: figS14_ssmxreg_failure.R
#
# Figure S14 — SSM-XREG Failure: Gear-Specific Collapse
#
# Shows SARIMA vs SSM-XREG skill (% RMSE improvement) at h=1 for all 8
# species at Gamma scale, both gears. The figure highlights the catastrophic
# SSM-XREG failure for Bag Seine transient taxa (Blue Crab, White Shrimp,
# Pinfish) while Bay Trawl shows positive skill for the same species.
#
# Layout: two panels side by side (Bag Seine | Bay Trawl)
# Y-axis: 8 species ordered by SARIMA skill (descending, top = highest)
# X-axis: skill (% RMSE improvement), clipped at -160%
# Models: SARIMA (orange circle) and SSM-XREG (pink diamond) only
# Reference line at 0
# Extreme values annotated with text labels
#
# Width: 170 mm; Height: 110 mm
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

need <- c("ggplot2", "dplyr", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
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
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Prepare data — SARIMA and SSM-XREG only, Gamma, h=1
# ─────────────────────────────────────────────────────────────────────────────
model_keep <- c("sarima", "ssm_ll_seas12_xreg")
model_labels <- c("sarima"             = "SARIMA",
                  "ssm_ll_seas12_xreg" = "SSM-XREG")

dat <- abd_raw %>%
  filter(model %in% model_keep,
         series == "gamma", BAY == "ALL", horizon == 1) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    model_label = recode(model, !!!model_labels),
    gear_label  = factor(ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
                         levels = c("Bag Seine", "Bay Trawl")),
    model_label = factor(model_label, levels = c("SARIMA", "SSM-XREG"))
  )

# Species order: descending SARIMA skill (mean across gears), highest at top
sp_order <- dat %>%
  filter(model_label == "SARIMA") %>%
  group_by(common_name) %>%
  summarise(m = mean(skill_rmse_improve_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>% pull(common_name)   # ascending = lowest at bottom

dat <- dat %>%
  mutate(common_name = factor(common_name, levels = sp_order))

# ─────────────────────────────────────────────────────────────────────────────
# 3) Clip and annotate extreme values
# ─────────────────────────────────────────────────────────────────────────────
x_min <- -60   # clip boundary

outliers <- dat %>%
  filter(skill_rmse_improve_pct < x_min) %>%
  mutate(label = sprintf("%.0f%%", skill_rmse_improve_pct))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Highlight species with BS SSM-XREG failure
# ─────────────────────────────────────────────────────────────────────────────
failure_species <- c("White shrimp", "Blue crab", "Pinfish")

# Shading data for failure species rows
shade_dat <- dat %>%
  filter(common_name %in% failure_species) %>%
  distinct(common_name, gear_label) %>%
  mutate(fill = "Transient taxa\n(BS failure)")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
model_cols   <- c("SARIMA"    = "#E69F00",
                  "SSM-XREG"  = "#CC79A7")
model_shapes <- c("SARIMA"    = 16,
                  "SSM-XREG"  = 18)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 7, face = "bold"),
    axis.text.y        = element_text(size = 6.5),
    axis.text.x        = element_text(size = 6),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.40, "cm"),
    plot.tag           = element_text(size = 8, face = "bold"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Build panel builder
# ─────────────────────────────────────────────────────────────────────────────
build_panel <- function(gear_val, tag_letter) {
  d        <- dat      %>% filter(gear_label == gear_val)
  out      <- outliers %>% filter(gear_label == gear_val)
  shade    <- shade_dat %>% filter(gear_label == gear_val)

  ggplot(d, aes(x = skill_rmse_improve_pct, y = common_name,
                colour = model_label, shape = model_label)) +
    # Highlight failure species rows — only in Bag Seine panel
    { if (gear_val == "Bag Seine")
        geom_tile(data = data.frame(
                    common_name = factor(failure_species, levels = sp_order),
                    x_centre    = (x_min + 55) / 2),
                  aes(x = x_centre, y = common_name,
                      width = 55 - x_min, height = 0.9),
                  fill = "#FFE8E8", colour = NA,
                  inherit.aes = FALSE)
    } +
    # Zero reference
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey40", linewidth = 0.4) +
    # Points
    geom_point(size = 2.5,
               position = position_dodge(width = 0.5)) +
    # Outlier labels inside panel
    geom_label(data = out,
               aes(x = x_min + 1, y = common_name,
                   label = label, colour = model_label),
               hjust = 0, size = 1.8, fill = "white",
               label.size = 0, label.padding = unit(0.08, "lines"),
               show.legend = FALSE) +
    scale_colour_manual(values = model_cols, name = "Model") +
    scale_shape_manual( values = model_shapes, name = "Model") +
    coord_cartesian(xlim = c(x_min, 55)) +
    labs(x = bquote("RMSE skill (%" * Delta * "RMSE), h=1"),
         y = NULL,
         tag = tag_letter,
         title = gear_val) +
    base_th +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}

pA <- build_panel("Bag Seine",  "A")
pB <- build_panel("Bay Trawl",  "B")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Assemble
# ─────────────────────────────────────────────────────────────────────────────
fig_S14 <- (pA | pB) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 8) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 110

out_tiff <- file.path(fig_dir, "Figure_S14.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S14.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S14); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S14); dev.off()

cat("\u2705 Figure S14 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
