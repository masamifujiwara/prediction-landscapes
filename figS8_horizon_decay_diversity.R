# =============================================================================
# Script: figS8_horizon_decay_diversity.R  (v2)
#
# Figure S8 — Horizon Decay Curves: Diversity Targets
# All labels built as proper data frames with facet columns so they
# render inside the correct facet panels.
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

need <- c("ggplot2", "dplyr", "tidyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load and prepare data
# ─────────────────────────────────────────────────────────────────────────────
div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

bay_labels <- c("2"="Galveston", "3"="Matagorda", "4"="San Antonio",
                "5"="Aransas",   "6"="Corpus Christi",
                "7"="Upper Laguna Madre", "8"="Lower Laguna Madre",
                "ALL"="Gamma (coastwide)")

dat <- div_raw %>%
  filter(model_clean == "sarima") %>%
  mutate(
    scale      = ifelse(BAY == "ALL", "Gamma", "Alpha"),
    bay_label  = recode(as.character(BAY), !!!bay_labels),
    q_label    = factor(paste0("q = ", q), levels = c("q = 0","q = 1","q = 2")),
    gear_label = factor(ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
                        levels = c("Bag Seine","Bay Trawl")),
    group_id   = paste0(bay_label, "_", GEAR, "_q", q)
  )

dat_gamma <- dat %>% filter(scale == "Gamma")
dat_alpha <- dat %>% filter(scale == "Alpha")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Build label data frames — must contain q_label and gear_label
#    so geom_text places them in the correct facet
# ─────────────────────────────────────────────────────────────────────────────

# Skill retention label (midpoint between h=1 and h=12)
retention_labels <- dat_gamma %>%
  select(q_label, gear_label, horizon, skill_rmse_improve_pct) %>%
  pivot_wider(names_from = horizon, values_from = skill_rmse_improve_pct,
              names_prefix = "h") %>%
  mutate(
    x     = 6.5,
    y     = (h1 + h12) / 2 + 1.5,
    label = paste0(round(100 * h12 / h1, 1), "% retention")
  )

# h=1 point labels (above each Gamma point)
h1_labels <- dat_gamma %>%
  filter(horizon == 1) %>%
  mutate(
    x     = 1,
    y     = skill_rmse_improve_pct + 1.8,
    label = paste0(round(skill_rmse_improve_pct, 1), "%")
  )

# h=12 point labels (above each Gamma point)
h12_labels <- dat_gamma %>%
  filter(horizon == 12) %>%
  mutate(
    x     = 12,
    y     = skill_rmse_improve_pct + 1.8,
    label = paste0(round(skill_rmse_improve_pct, 1), "%")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 3) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
col_gamma <- "grey15"
col_alpha <- "grey65"

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
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4) Build figure
# ─────────────────────────────────────────────────────────────────────────────
fig_S8 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  # Alpha bay lines
  geom_line(data = dat_alpha,
            aes(x = horizon, y = skill_rmse_improve_pct,
                group = group_id),
            colour = col_alpha, linewidth = 0.45, alpha = 0.7) +
  geom_point(data = dat_alpha,
             aes(x = horizon, y = skill_rmse_improve_pct,
                 group = group_id),
             colour = col_alpha, shape = 1, size = 1.5, stroke = 0.5) +
  # Gamma line
  geom_line(data = dat_gamma,
            aes(x = horizon, y = skill_rmse_improve_pct,
                group = group_id),
            colour = col_gamma, linewidth = 1.1) +
  geom_point(data = dat_gamma,
             aes(x = horizon, y = skill_rmse_improve_pct),
             colour = col_gamma, shape = 16, size = 2.5) +
  # h=1 labels above Gamma points
  geom_label(data = h1_labels,
             aes(x = x, y = y, label = label,
                 q_label = q_label, gear_label = gear_label),
             hjust = 0.5, size = 1.9, colour = col_gamma,
             fill = "white", label.size = 0, label.padding = unit(0.10, "lines")) +
  # h=12 labels above Gamma points
  geom_label(data = h12_labels,
             aes(x = x, y = y, label = label,
                 q_label = q_label, gear_label = gear_label),
             hjust = 0.5, size = 1.9, colour = col_gamma,
             fill = "white", label.size = 0, label.padding = unit(0.10, "lines")) +
  # Retention label at midpoint
  geom_label(data = retention_labels,
             aes(x = x, y = y, label = label,
                 q_label = q_label, gear_label = gear_label),
             hjust = 0.5, size = 1.8, colour = "grey40",
             fill = "white", label.size = 0, label.padding = unit(0.10, "lines"),
             fontface = "italic") +
  # Legend via dummy layers
  geom_line(data = data.frame(x=1, y=1,
                              scale="Gamma (coastwide)",
                              q_label=factor("q = 0", levels=levels(dat$q_label)),
                              gear_label=factor("Bag Seine",
                                                levels=levels(dat$gear_label))),
            aes(x=x, y=y, colour=scale, linewidth=scale),
            inherit.aes=FALSE) +
  geom_line(data = data.frame(x=1, y=1,
                              scale="Alpha (individual bays)",
                              q_label=factor("q = 0", levels=levels(dat$q_label)),
                              gear_label=factor("Bag Seine",
                                                levels=levels(dat$gear_label))),
            aes(x=x, y=y, colour=scale, linewidth=scale),
            inherit.aes=FALSE) +
  scale_colour_manual(
    values = c("Gamma (coastwide)"     = col_gamma,
               "Alpha (individual bays)" = col_alpha),
    name = "Spatial scale") +
  scale_linewidth_manual(
    values = c("Gamma (coastwide)"     = 1.1,
               "Alpha (individual bays)" = 0.45),
    name = "Spatial scale") +
  scale_x_continuous(breaks = c(1, 12),
                     labels = c("h = 1", "h = 12"),
                     expand = expansion(mult = c(0.08, 0.08))) +
  facet_grid(q_label ~ gear_label) +
  labs(x = "Forecast horizon",
       y = bquote("SARIMA skill (%" * Delta * "RMSE)")) +
  base_th +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 180

out_tiff <- file.path(fig_dir, "Figure_S8.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S8.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S8); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S8); dev.off()

cat("\u2705 Figure S8 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
