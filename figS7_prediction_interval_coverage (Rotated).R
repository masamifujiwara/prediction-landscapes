# =============================================================================
# Script: figS7_prediction_interval_coverage.R  (v2)
#
# Figure S7 — Prediction Interval Coverage
#
# Model name mapping (differs between Gamma and Alpha diversity):
#   Display    | Gamma div      | Alpha div                  | Abundance
#   SARIMA     | sarima         | sarima                     | sarima
#   SSM        | ssm            | ssm_ll_seas12              | ssm_ll_seas12
#   ARIMAX     | sarimax        | arimax_seasonal_autoarima  | sarimax
#   SSM-XREG   | ssm_ll_seas12_xreg | ssm_ll_seas12_xreg    | ssm_ll_seas12_xreg
#
# Layout: 2 rows x 2 cols
#   Row 1: h=1  | Row 2: h=12
#   Col 1: 80%  | Col 2: 95%
# Each panel: x = target (diversity q=0,1,2 + 8 species), faceted by gear
# Reference lines at nominal coverage (0.80, 0.95)
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
# 2) Diversity coverage — Gamma scale only
#    Model name mapping specific to Gamma diversity
# ─────────────────────────────────────────────────────────────────────────────
gamma_div_models <- c(
  "sarima"             = "SARIMA",
  "ssm"                = "SSM",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12_xreg" = "SSM-XREG"
)

div_cov <- div_raw %>%
  filter(BAY == "ALL",
         model_clean %in% names(gamma_div_models),
         !is.na(coverage80)) %>%
  mutate(
    model_label  = recode(model_clean, !!!gamma_div_models),
    target_label = paste0("q = ", q),
    gear_label   = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    target_type  = "Diversity"
  ) %>%
  select(target_type, target_label, gear_label, horizon,
         model_label, coverage80, coverage95)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Abundance coverage — Gamma scale only
# ─────────────────────────────────────────────────────────────────────────────
abd_models <- c(
  "sarima"             = "SARIMA",
  "ssm_ll_seas12"      = "SSM",
  "sarimax"            = "ARIMAX",
  "ssm_ll_seas12_xreg" = "SSM-XREG"
)

abd_cov <- abd_raw %>%
  filter(series == "gamma", BAY == "ALL",
         model %in% names(abd_models),
         !is.na(coverage80)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    model_label  = recode(model, !!!abd_models),
    target_label = common_name,
    gear_label   = ifelse(GEAR == "BS", "Bag Seine", "Bay Trawl"),
    target_type  = "Abundance"
  ) %>%
  select(target_type, target_label, gear_label, horizon,
         model_label, coverage80, coverage95)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Combine and set factor levels
# ─────────────────────────────────────────────────────────────────────────────
# Species ordered by descending mean SARIMA h=1 skill for consistency with
# other figures
sp_order <- abd_raw %>%
  filter(model == "sarima", horizon == 1, series == "gamma", BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  group_by(common_name) %>%
  summarise(m = mean(skill_rmse_improve_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>% pull(common_name)

div_order <- c("q = 0", "q = 1", "q = 2")
x_levels  <- c(div_order, sp_order)

all_cov <- bind_rows(div_cov, abd_cov) %>%
  mutate(
    target_label = factor(target_label, levels = x_levels),
    model_label  = factor(model_label,
                          levels = c("SARIMA","SSM","ARIMAX","SSM-XREG")),
    gear_label   = factor(gear_label, levels = c("Bag Seine","Bay Trawl"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Aesthetics
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

n_div <- length(div_order)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.text.x        = element_text(size = 6, angle = 38, hjust = 1),
    axis.text.y        = element_text(size = 6),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.tag           = element_text(size = 8, face = "bold"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Panel builder
# ─────────────────────────────────────────────────────────────────────────────
build_panel <- function(h_val, cov_col, nominal, tag_letter) {
  d <- all_cov %>%
    filter(horizon == h_val) %>%
    mutate(cov_val = .data[[cov_col]])

  ytitle <- paste0(nominal * 100, "% PI coverage")

  ggplot(d, aes(x = target_label, y = cov_val,
                colour = model_label, shape = model_label)) +
    # Shading for diversity group
    annotate("rect",
             xmin = 0.5, xmax = n_div + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "grey95", colour = NA) +
    # Separator between diversity and abundance
    geom_vline(xintercept = n_div + 0.5,
               colour = "grey55", linewidth = 0.5) +
    # Nominal reference line
    geom_hline(yintercept = nominal, linetype = "dashed",
               colour = "grey25", linewidth = 0.5) +
    geom_point(size = 2.0, alpha = 0.85,
               position = position_dodge(width = 0.65)) +
    # Group labels
    annotate("text", x = n_div / 2 + 0.5, y = Inf,
             label = "Diversity", vjust = 1.4,
             size = 2.0, fontface = "bold", colour = "grey30") +
    annotate("text", x = n_div + length(sp_order) / 2 + 0.5, y = Inf,
             label = "Abundance", vjust = 1.4,
             size = 2.0, fontface = "bold", colour = "grey30") +
    scale_colour_manual(values = model_cols, name = "Model") +
    scale_shape_manual( values = model_shapes, name = "Model") +
    scale_y_continuous(
      labels = function(x) paste0(round(x * 100), "%"),
      limits = c(0.50, 1.02),
      breaks = seq(0.50, 1.00, 0.10)
    ) +
    facet_grid(. ~ gear_label) +
    labs(x = NULL, y = ytitle,
         tag = tag_letter,
         title = paste0("h = ", h_val, " month",
                        ifelse(h_val > 1, "s", ""))) +
    base_th +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build panels
# ─────────────────────────────────────────────────────────────────────────────
pA <- build_panel(1,  "coverage80", 0.80, "A")
pB <- build_panel(1,  "coverage95", 0.95, "B")
pC <- build_panel(12, "coverage80", 0.80, "C")
pD <- build_panel(12, "coverage95", 0.95, "D")

# ─────────────────────────────────────────────────────────────────────────────
# 8) Assemble
# ─────────────────────────────────────────────────────────────────────────────
fig_S7 <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# ─────────────────────────────────────────────────────────────────────────────
# 9) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 180

out_tiff <- file.path(fig_dir, "Figure_S7.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S7.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S7); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S7); dev.off()

cat("\u2705 Figure S7 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
