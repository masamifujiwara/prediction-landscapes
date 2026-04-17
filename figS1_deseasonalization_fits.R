# =============================================================================
# Script: figS1_deseasonalization_fits.R  (v4)
#
# Title fix: place title inside the TOP data panel using labs(title=...).
# Add generous plot.margin(top) and title margin to push panel content down.
# Tag letter placed via labs(tag=...) on the top panel only.
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
  library(ggplot2); library(dplyr); library(patchwork); library(grid); library(gridtext)
})

fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────────
alpha_div <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds"))
gamma_div <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "c2_deseasonalize_alpha_gamma", "gamma_q012_deseas.rds"))
alpha_abd <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds"))
gamma_abd <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds"))
seas_fits <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "c1_seasonality_characterization", "seasonality_fits.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
col_obs  <- "grey60"
col_BS   <- "#E69F00"
col_BT   <- "#56B4E9"
col_zero <- "grey70"

# Base theme for data panels — NO title styling here
base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    axis.title        = element_text(size = 6.5),
    axis.text         = element_text(size = 6),
    plot.margin       = margin(1, 3, 1, 3, "mm")
  )

# Theme for the TOP panel that carries the title — extra top margin
top_th_extra <- theme(
  plot.title        = element_text(size = 7, face = "bold", hjust = 0,
                                   margin = margin(b = 4)),
  plot.tag          = element_text(size = 8, face = "bold"),
  plot.margin       = margin(6, 3, 0, 3, "mm")   # 6mm top gives title room
)

# ─────────────────────────────────────────────────────────────────────────────
# Diversity compound panel
# ─────────────────────────────────────────────────────────────────────────────
build_div_panel <- function(gear, q_val, scale, bay_id = NULL,
                            tag, title_text, gear_col) {
  if (scale == "Gamma") {
    ts_dat <- gamma_div %>% filter(GEAR == gear, q == q_val) %>% arrange(date)
  } else {
    ts_dat <- alpha_div %>%
      filter(GEAR == gear, q == q_val, BAY_ID == bay_id) %>% arrange(date)
  }
  if (nrow(ts_dat) == 0) return(NULL)

  fits_key <- paste0(gear, "__q", q_val)
  if (fits_key %in% names(seas_fits)) {
    sc <- seas_fits[[fits_key]]$seasonal_curve %>% select(MONTH, seasonal_centered)
    ts_dat <- ts_dat %>%
      left_join(sc, by = "MONTH") %>%
      mutate(fitted_curve = mean(value, na.rm = TRUE) + seasonal_centered)
  } else {
    ts_dat <- ts_dat %>%
      mutate(fitted_curve = value - value_deseas + mean(value, na.rm = TRUE))
  }

  # Top panel carries title + tag
  p_top <- ggplot(ts_dat, aes(x = date)) +
    geom_line(aes(y = value),        colour = col_obs,  linewidth = 0.35, alpha = 0.8) +
    geom_line(aes(y = fitted_curve), colour = gear_col, linewidth = 0.8) +
    labs(x = NULL, y = "Diversity\n(eff. species)",
         title = title_text, tag = tag) +
    base_th + top_th_extra +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_bot <- ggplot(ts_dat, aes(x = date)) +
    geom_hline(yintercept = 0, colour = col_zero, linewidth = 0.4) +
    geom_line(aes(y = value_deseas), colour = gear_col,
              linewidth = 0.35, alpha = 0.7) +
    labs(x = "Year", y = "Deseasonalized") +
    base_th +
    theme(plot.margin = margin(0, 3, 2, 3, "mm"))

  p_top / p_bot + plot_layout(heights = c(1, 0.45))
}

# ─────────────────────────────────────────────────────────────────────────────
# Abundance compound panel
# ─────────────────────────────────────────────────────────────────────────────
build_abd_panel <- function(gear, species_code, scale, bay_id = NULL,
                            tag, title_text, gear_col) {
  if (scale == "Gamma") {
    ts_dat <- gamma_abd %>%
      filter(GEAR == gear, SPECIES_CODE == species_code) %>% arrange(date)
  } else {
    ts_dat <- alpha_abd %>%
      filter(GEAR == gear, SPECIES_CODE == species_code, BAY_ID == bay_id) %>%
      arrange(date)
  }
  if (nrow(ts_dat) == 0) return(NULL)

  p_top <- ggplot(ts_dat, aes(x = date)) +
    geom_line(aes(y = value),
              colour = col_obs, linewidth = 0.35, alpha = 0.8) +
    geom_line(aes(y = seasonal_component + mean(value, na.rm = TRUE)),
              colour = gear_col, linewidth = 0.8) +
    labs(x = NULL, y = "log(CPUE + 1)",
         title = title_text, tag = tag) +
    base_th + top_th_extra +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_bot <- ggplot(ts_dat, aes(x = date)) +
    geom_hline(yintercept = 0, colour = col_zero, linewidth = 0.4) +
    geom_line(aes(y = value_deseas), colour = gear_col,
              linewidth = 0.35, alpha = 0.7) +
    labs(x = "Year", y = "Deseasonalized") +
    base_th +
    theme(plot.margin = margin(0, 3, 2, 3, "mm"))

  p_top / p_bot + plot_layout(heights = c(1, 0.45))
}

# ─────────────────────────────────────────────────────────────────────────────
# Build panels
# ─────────────────────────────────────────────────────────────────────────────
pA <- build_div_panel("BS", 0, "Gamma", tag="A",
      title_text="q = 0 (Richness) — Bag Seine, Gamma",   gear_col=col_BS)
pB <- build_div_panel("BS", 1, "Gamma", tag="B",
      title_text="q = 1 (Shannon) — Bag Seine, Gamma",    gear_col=col_BS)
pC <- build_div_panel("BS", 2, "Gamma", tag="C",
      title_text="q = 2 (Dominance) — Bag Seine, Gamma",  gear_col=col_BS)
pD <- build_abd_panel("BS", "601", "Alpha", bay_id=2, tag="D",
      title_text="Bay Anchovy — Bag Seine, Alpha (Galveston Bay)", gear_col=col_BS)

pE <- build_div_panel("BT", 0, "Gamma", tag="E",
      title_text="q = 0 (Richness) — Bay Trawl, Gamma",  gear_col=col_BT)
pF <- build_div_panel("BT", 1, "Gamma", tag="F",
      title_text="q = 1 (Shannon) — Bay Trawl, Gamma",   gear_col=col_BT)
pG <- build_div_panel("BT", 2, "Gamma", tag="G",
      title_text="q = 2 (Dominance) — Bay Trawl, Gamma", gear_col=col_BT)
pH <- build_abd_panel("BS", "601", "Gamma", tag="H",
      title_text="Bay Anchovy — Bag Seine, Gamma (Coastwide)", gear_col=col_BS)

# ─────────────────────────────────────────────────────────────────────────────
# Assemble
# ─────────────────────────────────────────────────────────────────────────────
make_hdr <- function(txt) {
  ggplot() +
annotate("text", x = 0.5, y = 0.35, label = txt,
         fontface = "bold", size = 3.0, hjust = 0.5, vjust = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(4, 0, 4, 0, "mm")
    )
}
# 
fig_S2 <-
  (make_hdr("Community Diversity") | make_hdr("Population Abundance")) /
  (pA | pE) /
  (pB | pF) /
  (pC | pG) /
  (pD | pH) +
  plot_layout(heights = c(0.055, 1, 1, 1, 1))

# Caption text
caption_txt <- paste(
  "Figure S1. Representative harmonic regression deseasonalization fits for community diversity and population abundance time series.",
  "Left column (panels A–D): bag seine community diversity at gamma (coastwide) scale for q = 0 (Panel A), q = 1 (Panel B), and q = 2 (Panel C), and bag seine bay anchovy log-CPUE at alpha (Galveston Bay) scale (Panel D).",
  "Right column (panels E–H): bay trawl community diversity at gamma scale for q = 0 (Panel E), q = 1 (Panel F), and q = 2 (Panel G), and bag seine bay anchovy log-CPUE at gamma (coastwide) scale (Panel H; note that Panel H shows bag seine gear, matching Panel D, rather than bay trawl gear).",
  "Within each panel, the top sub-panel shows the observed time series (grey line) and the fitted seasonal curve from harmonic regression (colored line); the bottom sub-panel shows the deseasonalized residuals used as model inputs.",
  "Orange indicates bag seine gear; blue indicates bay trawl gear."
)

# Main multi-panel figure
fig_core <-
  (make_hdr("Community Diversity") | make_hdr("Population Abundance")) /
  (pA | pE) /
  (pB | pF) /
  (pC | pG) /
  (pD | pH) +
  plot_layout(heights = c(0.055, 1, 1, 1, 1))

# Freeze main figure as a grob so patchwork won't try to re-layout it
fig_core_grob <- grid::grid.grabExpr(print(fig_core))
fig_core_panel <- patchwork::wrap_elements(full = fig_core_grob)

# Caption grob
caption_grob <- gridtext::textbox_grob(
  caption_txt,
  x = unit(0.5, "npc"),
  y = unit(1, "npc"),
  width = unit(0.82, "npc"),   # try 0.78 to 0.86 if needed
  hjust = 0.5,                 # center the caption block
  vjust = 1,
  halign = 0,                  # left-justify text within the block
  gp = grid::gpar(fontsize = 12),
  padding = unit(c(0, 0, 0, 0), "pt"),
  margin  = unit(c(0, 0, 0, 0), "pt"),
  box_gp  = grid::gpar(col = NA),
  r = unit(0, "pt")
)

caption_panel <- patchwork::wrap_elements(full = caption_grob)

# Final figure
fig_S1 <-
  fig_core_panel / caption_panel +
  plot_layout(heights = c(1, 0.24))

# ─────────────────────────────────────────────────────────────────────────────
# Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H1 <- 360; H2 <-280

out_tiff <- file.path(fig_dir, "Figure_S1.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S1.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H2),
     units="in", res=600, compression="lzw")
print(fig_S2); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H1))
print(fig_S1); dev.off()

cat("\u2705 Figure S1 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
