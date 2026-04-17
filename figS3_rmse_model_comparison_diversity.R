# =============================================================================
# Script: figS3_rmse_model_comparison_diversity.R
#
# Figure S3 — Full Model Comparison: RMSE Improvement, Diversity Targets
#
# Shows SARIMA, SSM, ARIMAX, and SSM-XREG skill (% RMSE improvement over
# seasonal naive) for all diversity strata:
#   q = 0, 1, 2  x  BS/BT  x  Alpha (7 bays) + Gamma
# Panels: h=1 (left) and h=12 (right)
# Each stratum is a row; models are colored dots.
#
# Width: 170 mm (double column); Height: ~200 mm
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

need <- c("ggplot2", "dplyr", "tidyr", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load data
# ─────────────────────────────────────────────────────────────────────────────
div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Filter and label
# ─────────────────────────────────────────────────────────────────────────────
# Models to show (exclude baseline)
model_keep <- c("sarima", "ssm_ll_seas12", "arimax_seasonal_autoarima",
                "ssm_ll_seas12_xreg", "sarimax")

# Clean model labels
model_labels <- c(
  "sarima"                   = "SARIMA",
  "ssm_ll_seas12"            = "SSM",
  "arimax_seasonal_autoarima"= "ARIMAX",
  "sarimax"                  = "ARIMAX",        # treat sarimax same as arimax
  "ssm_ll_seas12_xreg"       = "SSM-XREG"
)

# Bay labels
bay_labels <- c("2"="Galveston", "3"="Matagorda", "4"="San Antonio",
                "5"="Aransas",   "6"="Corpus Christi",
                "7"="Upper Laguna Madre", "8"="Lower Laguna Madre",
                "ALL"="Gamma (coastwide)")

dat <- div_raw %>%
  filter(model_clean %in% model_keep) %>%
  mutate(
    model_label = recode(model_clean, !!!model_labels),
    scale       = ifelse(BAY == "ALL", "Gamma", "Alpha"),
    bay_label   = recode(as.character(BAY), !!!bay_labels),
    q_label     = paste0("q = ", q),
    gear_label  = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    # Y-axis label: "q=0 | Galveston" etc.
    strat_label = paste0(bay_label)
  ) %>%
  # Collapse sarimax and arimax into one "ARIMAX" model — keep higher skill if both
  group_by(target, GEAR, q, BAY, horizon, model_label, scale, bay_label,
           q_label, gear_label, strat_label) %>%
  summarise(skill = max(skill_rmse_improve_pct, na.rm = TRUE),
            .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Y-axis ordering
#    Within each q x gear group: Gamma at top, then Alpha bays ordered
#    south to north (Lower LM=8, Upper LM=7, CC=6, Aransas=5,
#    San Antonio=4, Matagorda=3, Galveston=2)
# ─────────────────────────────────────────────────────────────────────────────
bay_order <- c("Gamma (coastwide)",
               "Galveston", "Matagorda", "San Antonio",
               "Aransas", "Corpus Christi",
               "Upper Laguna Madre", "Lower Laguna Madre")

# Full y-axis: for each q x gear combination, list bays
# We'll facet by q_label x gear_label and let y = bay_label within each facet
dat <- dat %>%
  mutate(bay_label = factor(bay_label, levels = rev(bay_order)))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
model_cols <- c(
  "SARIMA"    = "#E69F00",
  "SSM"       = "#56B4E9",
  "ARIMAX"    = "#009E73",
  "SSM-XREG"  = "#CC79A7"
)
model_shapes <- c(
  "SARIMA"    = 16,
  "SSM"       = 17,
  "ARIMAX"    = 15,
  "SSM-XREG"  = 18
)

dat <- dat %>%
  mutate(model_label = factor(model_label,
                               levels = c("SARIMA","SSM","ARIMAX","SSM-XREG")))

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.text.y        = element_text(size = 6),
    axis.text.x        = element_text(size = 6),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.tag           = element_text(size = 8, face = "bold"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Build one panel per horizon
# ─────────────────────────────────────────────────────────────────────────────
build_panel <- function(h_val, tag_letter) {
  d <- dat %>% filter(horizon == h_val)

  # Clip boundary
  x_min <- -25

  # Identify extreme outliers (skill < x_min) for annotation
  outliers <- d %>%
    filter(skill < x_min) %>%
    mutate(label = sprintf("%.0f%%", skill))

  ggplot(d, aes(x = skill, y = bay_label,
                colour = model_label, shape = model_label)) +
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey45", linewidth = 0.4) +
    geom_point(size = 2.0, alpha = 0.85,
               position = position_dodge(width = 0.6)) +
    # Outlier labels inside panel at left edge — no dodge, stacked by model
    geom_text(data = outliers,
              aes(x = x_min + 0.3, y = bay_label,
                  label = label, colour = model_label),
              size = 1.7, hjust = 0,
              show.legend = FALSE) +
    scale_colour_manual(values = model_cols, name = "Model") +
    scale_shape_manual( values = model_shapes, name = "Model") +
    # Use coord_cartesian to clip without removing data from legend etc.
    coord_cartesian(xlim = c(x_min, NA)) +
    facet_grid(q_label ~ gear_label, scales = "free_y", space = "free_y") +
    labs(x = bquote("RMSE skill (%" * Delta * "RMSE)"),
         y = NULL,
         tag = tag_letter,
         title = paste0("h = ", h_val, " month", ifelse(h_val > 1, "s", ""))) +
    base_th +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}

pA <- build_panel(1,  "A")
pB <- build_panel(12, "B")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Assemble: side by side, shared legend
# ─────────────────────────────────────────────────────────────────────────────
fig_S3 <- (pA / pB) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 260

out_tiff <- file.path(fig_dir, "Figure_S3.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S3.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S3); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S3); dev.off()

cat("\u2705 Figure S3 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
