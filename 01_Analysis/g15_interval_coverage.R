# =============================================================================
# Script:  g15_interval_coverage.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S6 (Online Resource 3): how often the observed value fell inside
#   the 80% (top) and 95% (bottom) prediction intervals of SARIMA (automatic),
#   by months ahead (1 to 24). Lines: median over the series of each group,
#   with the middle 50% of the estuary series shaded; coastwide diversity
#   drawn as its three orders, labeled. Dashed gray line: the nominal level.
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_S6_interval_coverage.tiff and .pdf
#   01_Analysis/Data/g15_interval_coverage/g15_coverage_by_group.csv
#   01_Analysis/Logs/g15_interval_coverage/run_<timestamp>.txt
#   Earlier copies of the figure and CSV file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Coverage from f01: over all scored forecast dates for each series and
#     lead time; all series, estuaries and coastwide, at every lead time.
#   - The figure shows SARIMA (automatic), the base model; the values file
#     also holds the local level, ARIMAX and state-space with covariates
#     models (MAIN).
#   - Coastwide diversity is drawn as its three orders, because a median of
#     three would hide one series.
#   - Gear by color and symbol, estuaries with solid lines and filled symbols,
#     coastwide with dashed lines and open symbols, as in Figs. S2 and S3.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

# ---- Bootstrap: find the project folder and build PATHS (same in every script)
# Looks for the marker file .project_root above, in order: the file being
# source()d (RStudio Source, source()), the file given to Rscript / R -f /
# R CMD BATCH (also under nohup), the file open in RStudio, the working folder.
local({
  a <- commandArgs(FALSE)
  here <- c(unlist(lapply(rev(sys.frames()), function(f) get0("ofile", envir = f, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", a, value = TRUE)),
            a[which(a == "-f") + 1L],
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
              tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
            getwd())
  here <- gsub("~+~", " ", here[!is.na(here) & nzchar(here)], fixed = TRUE)  # Rscript on Linux writes spaces as ~+~
  for (h in here) {
    d <- normalizePath(h, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) d <- dirname(d)
    repeat {
      if (file.exists(file.path(d, ".project_root"))) {
        assign(".PROJECT_ROOT", d, envir = globalenv())
        return(invisible())
      }
      if (dirname(d) == d) break
      d <- dirname(d)
    }
  }
  stop("Project folder not found: no .project_root above ", paste(unique(here), collapse = " or "))
})
source(file.path(.PROJECT_ROOT, "01_Analysis", "lib", "paths.R"))
PATHS <- get_paths(.PROJECT_ROOT)
setwd(PATHS$analysis_dir)
# -----------------------------------------------------------------------------
source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
source(file.path(PATHS$analysis_dir, "lib", "figures.R"))
need <- c("dplyr", "tidyr", "ggplot2", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})

FIGURE_LABEL <- "Figure_S6"   # Online Resource 3 figure number (Fig. S6)
script_stem <- "g15_interval_coverage"
MODEL <- "sarima_auto"; MAIN <- c("sarima_auto", "ssm_level", "arimax", "ssm_xreg"); W_MM <- 170; H_MM <- 130
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_lead, "f01_forecast_scoring.R")
out_csv <- file.path(dirs$data, "g15_coverage_by_group.csv"); if (file.exists(out_csv)) archive_files(out_csv)

bl <- readRDS(in_lead) %>% filter(model %in% MAIN)
if (!all(MAIN %in% bl$model)) stop("Main models missing from f01", call. = FALSE)
long <- bl %>% select(class, scale, GEAR, model, lead, coverage80, coverage95) %>%
  pivot_longer(c(coverage80, coverage95), names_to = "level", values_to = "coverage") %>%
  mutate(nominal = ifelse(level == "coverage80", 0.8, 0.95))
grp <- long %>% filter(!is.na(coverage)) %>% group_by(class, scale, GEAR, model, lead, level, nominal) %>%
  summarise(series = n(), med = median(coverage), lo = quantile(coverage, 0.25), hi = quantile(coverage, 0.75), .groups = "drop")
readr::write_csv(grp, out_csv)
chk <- grp %>% filter(model == MODEL)
if (n_distinct(chk$lead) != 24L || anyNA(chk$med)) stop("Expected SARIMA coverage at 24 leads", call. = FALSE)
print(as.data.frame(chk %>% filter(lead %in% c(1, 12, 24)) %>% select(class, scale, GEAR, lead, level, med)), row.names = FALSE)

GROUPS <- c("Bag seine, estuaries", "Bag seine, coastwide", "Bay trawl, estuaries", "Bay trawl, coastwide")
G_COL <- setNames(GEAR_COLS[c("Bag seine", "Bag seine", "Bay trawl", "Bay trawl")], GROUPS)
G_LT  <- setNames(c("solid", "22", "solid", "22"), GROUPS)
G_SH  <- setNames(c(16, 1, 17, 2), GROUPS)
MARK <- c(1, 6, 12, 18, 24)
lab_lvl <- function(l) factor(ifelse(l == "coverage80", "80% interval", "95% interval"), c("80% interval", "95% interval"))
grp_of <- function(G, sc) factor(paste0(GEAR_LAB[G], ", ", tolower(SCALE_LAB[sc])), GROUPS)
d <- chk %>% filter(!(class == "diversity" & scale == "coastwide")) %>%
  mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Group = grp_of(GEAR, scale), Gear = factor(GEAR_LAB[GEAR], GEAR_LAB), Level = lab_lvl(level))
# Coastwide diversity: each order separately (a median of three hides one series; as in Figs. 5 and S3).
cwd <- bl %>% filter(model == MODEL, class == "diversity", scale == "coastwide") %>%
  select(class, scale, GEAR, target, lead, coverage80, coverage95) %>%
  pivot_longer(c(coverage80, coverage95), names_to = "level", values_to = "coverage") %>%
  mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Group = grp_of(GEAR, scale), Level = lab_lvl(level),
         qlab = paste0("q = ", sub("^q", "", target)))
if (nrow(distinct(cwd, GEAR, target)) != 6L) stop("Expected 2 gears x 3 orders coastwide", call. = FALSE)
nom <- d %>% distinct(Level, nominal)
last <- cwd %>% filter(lead == 24)
p <- ggplot(d, aes(lead, med, colour = Group, linetype = Group, shape = Group)) +
  geom_hline(data = nom, aes(yintercept = nominal), colour = PAL[["grey_mid"]], linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
  geom_ribbon(data = d %>% filter(scale == "bay"), mapping = aes(x = lead, ymin = lo, ymax = hi, fill = Gear),
              inherit.aes = FALSE, colour = NA, alpha = 0.15, show.legend = FALSE) +
  geom_line(linewidth = 0.5) + geom_point(data = d %>% filter(lead %in% MARK), size = 1.5, stroke = 0.6) +
  geom_line(data = cwd, aes(lead, coverage, group = interaction(GEAR, target)), linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = cwd %>% filter(lead %in% MARK), aes(lead, coverage), size = 1.3, stroke = 0.5, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = last, aes(x = lead, y = coverage, label = qlab, colour = Group), inherit.aes = FALSE,
                           size = (FONT_SIZE - 2) / .pt, show.legend = FALSE, direction = "y", hjust = 0, nudge_x = 1,
                           min.segment.length = 0.1, segment.size = 0.2, box.padding = 0.1, seed = 1) +
  facet_grid(Level ~ Class, scales = "free_y") +
  scale_colour_manual(values = G_COL, name = NULL) + scale_linetype_manual(values = G_LT, name = NULL) +
  scale_shape_manual(values = G_SH, name = NULL) + scale_fill_manual(values = GEAR_COLS, guide = "none") +
  scale_x_continuous(breaks = MARK, expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_continuous(labels = function(z) paste0(round(100 * z), "%")) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Months ahead", y = "Observed values inside the interval (%)") + theme_paper() + theme(legend.position = "bottom")
save_figure(p, paste0(FIGURE_LABEL, "_interval_coverage"), W_MM, H_MM)
cat("Saved:", out_csv, "\n")
end_log("finished")
cat("g15_interval_coverage.R completed\n")
