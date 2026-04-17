# =============================================================================
# Script: figS12_bay_alpha_abundance_h1.R  (v2)
#
# Figure S12 — Bay-Level Alpha Abundance Results: h=1
#
# Additions from v1:
#   1. Bootstrap 95% CIs on each dot (error bars)
#   2. Vertical dashed line per facet showing Gamma-scale skill value
#      (illustrates portfolio effect — coastwide vs local)
#
# Bootstrap: B=2000, seed=42, fold-level resampling from pred_joined_abundance
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

abd_folds <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance",
  "pred_joined_abundance.rds"))

sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  stringsAsFactors = FALSE)

bay_labels <- c("2"="Galveston", "3"="Matagorda", "4"="San Antonio",
                "5"="Aransas",   "6"="Corpus Christi",
                "7"="Upper\nLaguna Madre", "8"="Lower\nLaguna Madre")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Bootstrap CI function
# ─────────────────────────────────────────────────────────────────────────────
boot_skill_ci <- function(error_adv, error_base, B = 2000, alpha = 0.05) {
  ok <- !is.na(error_adv) & !is.na(error_base)
  ea <- error_adv[ok];  eb <- error_base[ok];  n <- length(ea)
  if (n < 5) return(c(skill=NA, lo=NA, hi=NA))
  pt <- 100 * (1 - sqrt(mean(ea^2)) / sqrt(mean(eb^2)))
  bs <- vapply(seq_len(B), function(i) {
    idx <- sample(n, n, replace=TRUE)
    100 * (1 - sqrt(mean(ea[idx]^2)) / sqrt(mean(eb[idx]^2)))
  }, numeric(1))
  c(skill = pt,
    lo    = unname(quantile(bs, alpha/2)),
    hi    = unname(quantile(bs, 1-alpha/2)))
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) Compute bootstrap CIs for all 112 alpha strata
# ─────────────────────────────────────────────────────────────────────────────
cat("Computing bootstrap CIs for 112 alpha strata...\n")
set.seed(42)

bays    <- as.character(c(2,3,4,5,6,7,8))
gears   <- c("BS","BT")
species <- sp_lut$SPECIES_CODE

ci_rows <- list()
for (g in gears) {
  for (sp in species) {
    for (b in bays) {
      d <- abd_folds %>%
        filter(model == "sarima", horizon == 1,
               series == "alpha", GEAR == g,
               SPECIES_CODE == sp, BAY == b,
               !is.na(error_base))
      ci <- boot_skill_ci(d$error, d$error_base)
      ci_rows[[length(ci_rows)+1]] <- data.frame(
        GEAR=g, SPECIES_CODE=sp, BAY=b,
        skill=ci["skill"], lo=ci["lo"], hi=ci["hi"],
        stringsAsFactors=FALSE)
    }
  }
}
ci_df <- bind_rows(ci_rows)
cat("Done.\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Gamma-scale skill reference values (one per species x gear)
# ─────────────────────────────────────────────────────────────────────────────
gamma_ref <- abd_raw %>%
  filter(model == "sarima", horizon == 1,
         series == "gamma", BAY == "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  select(GEAR, SPECIES_CODE, gamma_skill = skill_rmse_improve_pct)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Merge everything
# ─────────────────────────────────────────────────────────────────────────────
dat <- abd_raw %>%
  filter(model == "sarima", horizon == 1,
         series == "alpha", BAY != "ALL") %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  left_join(ci_df, by = c("GEAR","SPECIES_CODE","BAY")) %>%
  left_join(gamma_ref, by = c("GEAR","SPECIES_CODE")) %>%
  mutate(
    bay_label  = factor(recode(as.character(BAY), !!!bay_labels),
                        levels = unname(bay_labels)),
    gear_label = factor(ifelse(GEAR=="BS","Bag Seine","Bay Trawl"),
                        levels = c("Bag Seine","Bay Trawl")),
    sig_label  = ifelse(is_significantly_better, "p < 0.05", "n.s.")
  )

# Species order: ascending mean skill (highest at top)
sp_order <- dat %>%
  group_by(common_name) %>%
  summarise(m = mean(skill_rmse_improve_pct, na.rm=TRUE), .groups="drop") %>%
  arrange(m) %>% pull(common_name)

dat <- dat %>%
  mutate(common_name = factor(common_name, levels=sp_order))

# Gamma reference — needs same faceting columns
gamma_lines <- dat %>%
  distinct(gear_label, bay_label, common_name, gamma_skill)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Aesthetics
# ─────────────────────────────────────────────────────────────────────────────
sig_cols   <- c("p < 0.05" = "#2166AC", "n.s." = "grey70")
sig_shapes <- c("p < 0.05" = 16,        "n.s." = 1)

base_th <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(size = 5.5, face = "bold"),
    axis.text.y        = element_text(size = 6),
    axis.text.x        = element_text(size = 5.5),
    axis.title         = element_text(size = 7),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.38, "cm"),
    plot.margin        = margin(2, 3, 2, 3, "mm")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build figure
# ─────────────────────────────────────────────────────────────────────────────
fig_S12 <- ggplot(dat,
                  aes(x = skill_rmse_improve_pct, y = common_name)) +
  # Zero reference
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  # Gamma-scale reference line — use geom_point with a different shape
  # to mark the gamma value at each species row without conflicting with
  # the discrete y-axis scale
  geom_point(data = gamma_lines,
             aes(x = gamma_skill, y = common_name),
             shape = 124,   # vertical bar character
             colour = "#D73027", size = 4.5, alpha = 0.8) +
  # Error bars (bootstrap 95% CI)
  geom_errorbarh(aes(xmin = lo, xmax = hi,
                     colour = sig_label),
                 height = 0.25, linewidth = 0.4) +
  # Points
  geom_point(aes(colour = sig_label, shape = sig_label),
             size = 1.8) +
  scale_colour_manual(values = sig_cols,   name = "DM test") +
  scale_shape_manual( values = sig_shapes, name = "DM test") +
  facet_grid(gear_label ~ bay_label) +
  labs(x = bquote("SARIMA skill, h=1 (%" * Delta * "RMSE)"),
       y = NULL,
       caption = "Red | = Gamma-scale (coastwide) SARIMA skill") +
  base_th +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        plot.caption     = element_text(size = 6, colour = "#D73027",
                                        hjust = 0))

# ─────────────────────────────────────────────────────────────────────────────
# 8) Export
# ─────────────────────────────────────────────────────────────────────────────
mm_to_in <- function(mm) mm / 25.4
W <- 170; H <- 150

out_tiff <- file.path(fig_dir, "Figure_S12.tiff")
out_pdf  <- file.path(fig_dir, "Figure_S12.pdf")

tiff(out_tiff, width=mm_to_in(W), height=mm_to_in(H),
     units="in", res=600, compression="lzw")
print(fig_S12); dev.off()

pdf(out_pdf, width=mm_to_in(W), height=mm_to_in(H))
print(fig_S12); dev.off()

cat("\u2705 Figure S12 saved:\n  TIFF:", out_tiff, "\n  PDF:", out_pdf, "\n")
