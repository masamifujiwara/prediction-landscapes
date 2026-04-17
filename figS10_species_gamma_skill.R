# =============================================================================
# Script: figS10_species_gamma_skill.R
#
# Figure S10 — Species-Level SARIMA Skill: Gamma Scale
#
# SARIMA skill (% RMSE improvement) at h=1 and h=12 for all 8 species
# at Gamma scale, both gears. h=1 and h=12 shown as different colors/shapes,
# connected by a thin line per species to show horizon decay.
#
# Faceted by gear (Bag Seine | Bay Trawl)
# Species ordered by descending h=1 mean skill (highest at top)
#
# Width: 120 mm (single column); Height: 100 mm
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

need <- c("ggplot2", "dplyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
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
    gear_label    = factor(ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
                           levels = c("Bag Seine", "Bay Trawl")),
    horizon_label = factor(paste0("h = ", horizon),
                           levels = c("h = 1", "h = 12"))
  )

# Order species by ascending h=1 mean skill (lowest at bottom, highest at top)
sp_order <- dat %>%
  filter(horizon == 1) %>%
  group_by(common_name) %>%
  summarise(m = mean(skill_rmse_improve_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>% pull(common_name)

dat <- dat %>%
  mutate(common_name = factor(common_name, levels = sp_order))

# ─────────────────────────────────────────────────────────────────────────────
# 3) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
h_cols   <- c("h = 1"  = "#2166AC",
              "h = 12" = "#D73027")
h_shapes <- c("h = 1"  = 16,
              "h = 12" = 17)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.text          = element_text(size = 6.5),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4) Build figure
# ─────────────────────────────────────────────────────────────────────────────
fig_S10 <- ggplot(dat,
                  aes(x = skill_rmse_improve_pct, y = common_name)) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  # Thin connecting line between h=1 and h=12 per species
  geom_line(aes(group = common_name),
            colour = "grey70", linewidth = 0.5) +
  # Points colored and shaped by horizon
  geom_point(aes(colour = horizon_label, shape = horizon_label),
             size = 2.4) +
  scale_colour_manual(values = h_cols,   name = "Horizon") +
  scale_shape_manual( values = h_shapes, name = "Horizon") +
  facet_grid(. ~ gear_label) +
  labs(x = bquote("SARIMA skill (%" * Delta * "RMSE)"),
       y = NULL) +
  base_th +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 120; H <- 100

out_tiff <- file.path(fig_dir, "Figure_S10.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S10.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S10); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S10); dev.off()

cat("\u2705 Figure S10 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
