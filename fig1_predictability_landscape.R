# =============================================================================
# Script: fig1_predictability_landscape.R  (v14)
#
# Legend fix: keep BOTH colour and shape guides merged into one legend block
# by giving them the same title. shape = "none" suppressed the open symbols
# so we restore it. The key insight: to show open shapes in the legend,
# we need shape in aes() AND must NOT suppress the shape guide.
# We merge colour+shape by giving both the same name.
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 0) Bootstrap
# ─────────────────────────────────────────────────────────────────────────────
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
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

# ─────────────────────────────────────────────────────────────────────────────
# 1) Packages
# ─────────────────────────────────────────────────────────────────────────────
need <- c("ggplot2", "dplyr", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

# ─────────────────────────────────────────────────────────────────────────────
# 2) Data
# ─────────────────────────────────────────────────────────────────────────────
fig_dir <- file.path(PATHS$analysis_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

div_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3",
  "metrics_dictionary_alpha_gamma.csv"), stringsAsFactors = FALSE)

abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "metrics_dictionary_abundance.csv"), stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Lookups
# ─────────────────────────────────────────────────────────────────────────────
sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay Anchovy","Atlantic Croaker","Gulf Menhaden","Spot",
                   "Pinfish","White Shrimp","Blue Crab","Brown Shrimp"),
  stringsAsFactors = FALSE)

sp_order <- abd_raw %>%
  filter(model == "seasonal_naive", horizon == 1, BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  group_by(common_name) %>%
  summarise(m = mean(rmse_base, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>% pull(common_name)

q_levels    <- c("q = 0", "q = 1", "q = 2")
bay_offsets <- setNames(seq(-0.30, 0.30, length.out = 7),
                        as.character(c(2,3,4,5,6,7,8)))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
col_BS <- "#E69F00"
col_BT <- "#56B4E9"

gs_levels <- c("BS \u2013 Gamma", "BT \u2013 Gamma",
               "BS \u2013 Alpha", "BT \u2013 Alpha")
gs_cols   <- setNames(c(col_BS, col_BT, col_BS, col_BT), gs_levels)
gs_shapes <- setNames(c(16, 17, 1, 2),                   gs_levels)
# stroke for each entry: filled shapes need stroke=0, open shapes need stroke>0
gs_stroke <- setNames(c(0, 0, 0.9, 0.9), gs_levels)
gs_fill   <- setNames(c(col_BS, col_BT, NA, NA), gs_levels)

gear_lab <- c("BS" = "Bag Seine", "BT" = "Bay Trawl")
q_lty    <- c("q = 0" = "solid",  "q = 1" = "dashed", "q = 2" = "dotdash")
q_shp    <- c("q = 0" = 16,       "q = 1" = 17,       "q = 2" = 15)
sp_shp   <- setNames(c(16,17,15,18,1,2,0,8), sp_order)

base_th <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey91", linewidth = 0.3),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6.5),
    legend.title      = element_text(size = 7, face = "bold"),
    legend.text       = element_text(size = 6.5),
    legend.key.size   = unit(0.40, "cm"),
    plot.tag          = element_text(size = 9, face = "bold"),
    plot.margin       = margin(2, 2, 2, 2, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Horizontal legend row for panels A-D
#    Both colour AND shape mapped to the same variable with the same name
#    so ggplot2 merges them into one legend block.
#    override.aes explicitly sets shape and stroke per entry so open shapes
#    render correctly.
# ─────────────────────────────────────────────────────────────────────────────
dummy4 <- data.frame(
  x  = 1:4, y = rep(1, 4),
  gs = factor(gs_levels, levels = gs_levels))

legend_row <- ggplot(dummy4, aes(x = x, y = y, colour = gs, shape = gs)) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = gs_cols,   name = NULL, drop = FALSE) +
  scale_shape_manual( values = gs_shapes, name = NULL, drop = FALSE) +
  # Merge colour + shape into one legend by using same name (NULL = same)
  # override.aes sets correct size, stroke, and fill for each of the 4 rows
  guides(
    colour = guide_legend(
      nrow = 1,
      override.aes = list(
        shape  = unname(gs_shapes),   # 16, 17, 1, 2
        size   = rep(2.4, 4),
        stroke = unname(gs_stroke),   # 0, 0, 0.9, 0.9
        alpha  = rep(1, 4)
      )
    ),
    shape = "none"   # suppressed — colour guide already shows correct shapes
  ) +
  coord_cartesian(xlim = c(100, 101), ylim = c(100, 101)) +
  theme_void(base_size = 8) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.42, "cm"),
    legend.key        = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(1, 0, 1, 0, "mm"),
    plot.margin       = margin(0, 0, 0, 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Diversity panels A & C
# ─────────────────────────────────────────────────────────────────────────────
build_div <- function(mdl, val_col, ylab, tag, hline = FALSE) {

  da <- div_raw %>%
    filter(model_clean == mdl, horizon == 1, BAY != "ALL") %>%
    mutate(q_num  = as.integer(as.character(q)) + 1L,
           x_plot = q_num + bay_offsets[as.character(BAY)],
           value  = .data[[val_col]])

  dg <- div_raw %>%
    filter(model_clean == mdl, horizon == 1, BAY == "ALL") %>%
    mutate(q_num = as.integer(as.character(q)) + 1L,
           value = .data[[val_col]])

  ggplot() +
    { if (hline) geom_hline(yintercept = 0, linetype = "dotted",
                             colour = "grey45", linewidth = 0.4) } +
    geom_point(data = filter(da, GEAR == "BS"),
               aes(x = x_plot, y = value),
               shape = 1, colour = col_BS,
               size = 2.0, stroke = 0.8, alpha = 0.65) +
    geom_point(data = filter(da, GEAR == "BT"),
               aes(x = x_plot, y = value),
               shape = 2, colour = col_BT,
               size = 2.0, stroke = 0.8, alpha = 0.65) +
    geom_line(data = filter(dg, GEAR == "BS"),
              aes(x = q_num, y = value, group = 1),
              colour = col_BS, linewidth = 1.0) +
    geom_line(data = filter(dg, GEAR == "BT"),
              aes(x = q_num, y = value, group = 1),
              colour = col_BT, linewidth = 1.0) +
    geom_point(data = filter(dg, GEAR == "BS"),
               aes(x = q_num, y = value),
               shape = 16, colour = col_BS, size = 3.0) +
    geom_point(data = filter(dg, GEAR == "BT"),
               aes(x = q_num, y = value),
               shape = 17, colour = col_BT, size = 3.0) +
    scale_x_continuous(breaks = 1:3, labels = q_levels,
                       expand = expansion(mult = 0.22)) +
    labs(x = "Diversity order", y = ylab, tag = tag) +
    base_th +
    theme(legend.position = "none",
          plot.margin     = margin(2, 2, 0, 2, "mm"))  # reduce bottom margin
}

panel_A <- build_div("seasonal_naive", "rmse_base",
                     "Baseline RMSE\n(effective species)", "A")
panel_C <- build_div("sarima", "skill_rmse_improve_pct",
                     bquote("SARIMA skill, h=1 (%" * Delta * "RMSE)"), "C",
                     hline = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Abundance panels B & D
# ─────────────────────────────────────────────────────────────────────────────
build_abd <- function(mdl, val_col, ylab, tag, hline = FALSE) {

  da <- abd_raw %>%
    filter(model == mdl, horizon == 1, BAY != "ALL") %>%
    mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
    left_join(sp_lut, by = "SPECIES_CODE") %>%
    mutate(common_name = factor(common_name, levels = sp_order),
           value = .data[[val_col]])

  dg <- abd_raw %>%
    filter(model == mdl, horizon == 1, BAY == "ALL") %>%
    mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
    left_join(sp_lut, by = "SPECIES_CODE") %>%
    mutate(common_name = factor(common_name, levels = sp_order),
           value = .data[[val_col]])

  ggplot() +
    { if (hline) geom_hline(yintercept = 0, linetype = "dotted",
                             colour = "grey45", linewidth = 0.4) } +
    geom_point(data = filter(da, GEAR == "BS"),
               aes(x = common_name, y = value),
               shape = 1, colour = col_BS,
               size = 2.0, stroke = 0.8, alpha = 0.65,
               position = position_jitter(width = 0.10, seed = 42)) +
    geom_point(data = filter(da, GEAR == "BT"),
               aes(x = common_name, y = value),
               shape = 2, colour = col_BT,
               size = 2.0, stroke = 0.8, alpha = 0.65,
               position = position_jitter(width = 0.10, seed = 43)) +
    geom_point(data = filter(dg, GEAR == "BS"),
               aes(x = common_name, y = value),
               shape = 16, colour = col_BS, size = 3.0,
               position = position_nudge(x = -0.18)) +
    geom_point(data = filter(dg, GEAR == "BT"),
               aes(x = common_name, y = value),
               shape = 17, colour = col_BT, size = 3.0,
               position = position_nudge(x = +0.18)) +
    labs(x = "Species", y = ylab, tag = tag) +
    base_th +
    theme(axis.text.x     = element_text(angle = 38, hjust = 1, size = 6),
          legend.position = "none",
          plot.margin     = margin(2, 2, 0, 2, "mm"))  # reduce bottom margin
}

panel_B <- build_abd("seasonal_naive", "rmse_base",
                     "Baseline RMSE\n(log CPUE)", "B")
panel_D <- build_abd("sarima", "skill_rmse_improve_pct",
                     bquote("SARIMA skill, h=1 (%" * Delta * "RMSE)"), "D",
                     hline = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Panel E
# ─────────────────────────────────────────────────────────────────────────────
dE <- div_raw %>%
  filter(model_clean == "sarima", BAY == "ALL") %>%
  mutate(q_lab = factor(paste0("q = ", q), levels = q_levels))

panel_E <- ggplot(dE,
                  aes(x = horizon, y = skill_rmse_improve_pct,
                      colour = GEAR, shape = q_lab, linetype = q_lab,
                      group = interaction(GEAR, q_lab))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1,12), labels = c("h = 1","h = 12")) +
  scale_colour_manual(  values = c("BS" = col_BS, "BT" = col_BT),
                        labels = gear_lab, name = "Gear") +
  scale_shape_manual(   values = q_shp,  name = "Diversity\norder") +
  scale_linetype_manual(values = q_lty,  name = "Diversity\norder") +
  labs(x = "Forecast horizon",
       y = bquote("SARIMA skill (%" * Delta * "RMSE)"), tag = "E") +
  base_th +
  theme(legend.position  = "right",
        legend.box       = "vertical",
        legend.spacing.y = unit(0.0, "cm"))

# ─────────────────────────────────────────────────────────────────────────────
# 9) Panel F
# ─────────────────────────────────────────────────────────────────────────────
dF <- abd_raw %>%
  filter(model == "sarima", BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(common_name = factor(common_name, levels = sp_order))

panel_F <- ggplot(dF,
                  aes(x = horizon, y = skill_rmse_improve_pct,
                      group = interaction(GEAR, common_name))) +
  geom_line(data = filter(dF, GEAR == "BS"),
            aes(group = common_name),
            colour = col_BS, linetype = "solid",  linewidth = 0.7) +
  geom_line(data = filter(dF, GEAR == "BT"),
            aes(group = common_name),
            colour = col_BT, linetype = "dashed", linewidth = 0.7) +
  geom_point(aes(colour = GEAR, shape = common_name), size = 2.2) +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  scale_x_continuous(breaks = c(1,12), labels = c("h = 1","h = 12")) +
  scale_colour_manual(values = c("BS" = col_BS, "BT" = col_BT),
                      labels = gear_lab, name = "Gear") +
  scale_shape_manual(values = sp_shp, name = "Species") +
  labs(x = "Forecast horizon",
       y = bquote("SARIMA skill (%" * Delta * "RMSE)"), tag = "F") +
  base_th +
  theme(legend.position  = "right",
        legend.box       = "vertical",
        legend.spacing.y = unit(0.05, "cm"),
        legend.text      = element_text(size = 6.0),
        legend.key.size  = unit(0.32, "cm")) +
  guides(
    colour = guide_legend(title = "Gear", order = 1,
                          override.aes = list(
                            shape     = NA,
                            linewidth = 1.1,
                            linetype  = c("solid","dashed"))),
    shape  = guide_legend(title = "Species", order = 2,
                          override.aes = list(colour = "grey25", size = 2.0))
  )

# ─────────────────────────────────────────────────────────────────────────────
# 10) Assemble
# ─────────────────────────────────────────────────────────────────────────────
make_hdr <- function(txt)
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             fontface = "bold", size = 3.2, hjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(2, 0, 0, 0, "mm"))

fig1 <-
  (make_hdr("Community Diversity") | make_hdr("Population Abundance")) /
  (panel_A | panel_B) /
  (panel_C | panel_D) /
  legend_row /
  (panel_E | panel_F) +
  plot_layout(heights = c(0.05, 1, 1, 0.11, 1))

# ─────────────────────────────────────────────────────────────────────────────
# 11) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 222

out_tiff <- file.path(fig_dir, "Figure_1.tiff")
out_pdf  <- file.path(fig_dir, "Figure_1.pdf")

tiff(out_tiff, width = mm_to_in(W), height = mm_to_in(H),
     units = "in", res = 600, compression = "lzw")
print(fig1); dev.off()

pdf(out_pdf, width = mm_to_in(W), height = mm_to_in(H))
print(fig1); dev.off()

cat("\u2705 Figure 1 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
