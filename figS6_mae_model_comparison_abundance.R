# =============================================================================
# Script: figS6_mae_model_comparison_abundance.R
#
# Figure S6 — Full Model Comparison: MAE Improvement, Abundance Targets
#             (Gamma scale, all 8 species, both gears, h=1 and h=12)
#
# Mirrors Figure S3 design but uses abundance metrics dictionary.
# Rows: 8 species ordered by descending mean SARIMA skill
# Columns: Bag Seine | Bay Trawl
# Panels: h=1 (A, top) and h=12 (B, bottom)
# Models: SARIMA, SSM, ARIMAX, SSM-XREG as colored/shaped dots
# X-axis clipped at -25%; extreme values annotated with text labels
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
  common_name  = c("Bay Anchovy","Atlantic Croaker","Gulf Menhaden","Spot",
                   "Pinfish","White Shrimp","Blue Crab","Brown Shrimp"),
  stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Filter and label
# ─────────────────────────────────────────────────────────────────────────────
model_keep <- c("sarima","ssm_ll_seas12","sarimax","ssm_ll_seas12_xreg")

model_labels <- c(
  "sarima"                = "SARIMA",
  "ssm_ll_seas12"         = "SSM",
  "sarimax"               = "ARIMAX",
  "ssm_ll_seas12_xreg"    = "SSM-XREG"
)

dat <- abd_raw %>%
  filter(model %in% model_keep,
         series == "gamma",
         BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    model_label = recode(model, !!!model_labels),
    gear_label  = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    skill       = skill_mae_improve_pct
  ) %>%
  # Collapse sarimax duplicates — keep max skill per stratum
  group_by(common_name, GEAR, gear_label, horizon, model_label) %>%
  summarise(skill = max(skill, na.rm = TRUE), .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Species ordering: descending mean SARIMA skill at h=1 (top = highest)
# ─────────────────────────────────────────────────────────────────────────────
sp_order <- dat %>%
  filter(model_label == "SARIMA", horizon == 1) %>%
  group_by(common_name) %>%
  summarise(m = mean(skill, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>%   # ascending = lowest at bottom of y-axis
  pull(common_name)

dat <- dat %>%
  mutate(
    common_name  = factor(common_name, levels = sp_order),
    model_label  = factor(model_label,
                          levels = c("SARIMA","SSM","ARIMAX","SSM-XREG"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4) Aesthetics (same as Figure S3)
# ─────────────────────────────────────────────────────────────────────────────
model_cols <- c(
  "SARIMA"   = "#E69F00",
  "SSM"      = "#56B4E9",
  "ARIMAX"   = "#009E73",
  "SSM-XREG" = "#CC79A7"
)
model_shapes <- c(
  "SARIMA"   = 16,
  "SSM"      = 17,
  "ARIMAX"   = 15,
  "SSM-XREG" = 18
)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.text.y        = element_text(size = 6.5),
    axis.text.x        = element_text(size = 6),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.tag           = element_text(size = 8, face = "bold"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Panel builder
# ─────────────────────────────────────────────────────────────────────────────
x_min <- -25

build_panel <- function(h_val, tag_letter) {
  d <- dat %>% filter(horizon == h_val)

  outliers <- d %>%
    filter(skill < x_min) %>%
    mutate(label = sprintf("%.0f%%", skill))

  ggplot(d, aes(x = skill, y = common_name,
                colour = model_label, shape = model_label)) +
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey45", linewidth = 0.4) +
    geom_point(size = 2.2, alpha = 0.85,
               position = position_dodge(width = 0.6)) +
    geom_text(data = outliers,
              aes(x = x_min + 0.3, y = common_name,
                  label = label, colour = model_label),
              size = 1.8, hjust = 0,
              show.legend = FALSE) +
    scale_colour_manual(values = model_cols, name = "Model") +
    scale_shape_manual( values = model_shapes, name = "Model") +
    coord_cartesian(xlim = c(x_min, NA)) +
    facet_grid(. ~ gear_label) +
    labs(x = bquote("MAE skill (%" * Delta * "MAE)"),
         y = NULL,
         tag = tag_letter,
         title = paste0("h = ", h_val, " month",
                        ifelse(h_val > 1, "s", ""))) +
    base_th +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}

pA <- build_panel(1,  "A")
pB <- build_panel(12, "B")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Assemble
# ─────────────────────────────────────────────────────────────────────────────
fig_S4 <- (pA / pB) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 160

out_tiff <- file.path(fig_dir, "Figure_S6.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S6.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S4); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S4); dev.off()

cat("\u2705 Figure S6 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
