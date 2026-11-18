# =============================================================================
# Script:  g14_horizon_tolerance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S2 (Online Resource 3): how the accuracy horizon of SARIMA
#   (automatic) depends on the tolerance. Percentage of the series of each
#   group whose forecasts stay within the tolerance (top) at least 1 month
#   ahead and (bottom) at all 24 months, by the range rule, against the
#   tolerance: 20 to 60% of the series mean for diversity, 0.3 to 0.9 log
#   units for abundance. Dashed line: the tolerance used in the rest of the
#   paper, read from the forecast settings.
#
# Run after: b09, f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_horizon_by_tolerance.rds (f01)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09, tolerance)
#
# Outputs:
#   01_Analysis/figures/Figure_S2_horizon_tolerance.tiff and .pdf
#   01_Analysis/Data/g14_horizon_tolerance/g14_horizon_tolerance_values.csv
#   01_Analysis/Logs/g14_horizon_tolerance/run_<timestamp>.txt
#   Earlier copies of the figure and CSV file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - The tolerance is the same at every lead time; this figure shows a range
#     of tolerances for readers with other needs.
#   - The range rule is shown; value-only percentages are in the values file.
#   - Gear by color and symbol, estuaries with solid lines and filled symbols,
#     coastwide with dashed lines and open symbols, as in Fig. S3.
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
source(file.path(PATHS$analysis_dir, "lib", "forecast_rules.R"))
need <- c("dplyr", "tidyr", "ggplot2", "patchwork", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})

FIGURE_LABEL <- "Figure_S2"   # Online Resource 3 figure number (Fig. S2)
script_stem <- "g14_horizon_tolerance"
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 120
in_t <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_horizon_by_tolerance.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_t, in_settings), c("f01_forecast_scoring.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
tol <- c(diversity = settings$tolerance$diversity_rel, abundance = settings$tolerance$abundance_log)
out_csv <- file.path(dirs$data, "g14_horizon_tolerance_values.csv"); if (file.exists(out_csv)) archive_files(out_csv)

x <- readRDS(in_t) %>% filter(model == MODEL) %>% mutate(tolerance = round(tolerance, 2))
v <- x %>% group_by(class, scale, GEAR, tolerance) %>%
  summarise(series = n(), at_least_1 = mean(abs_horizon >= 1), all_24 = mean(abs_horizon == 24),
            at_least_1_value = mean(abs_horizon_value >= 1), all_24_value = mean(abs_horizon_value == 24), .groups = "drop")
if (n_distinct(v$series[v$class == "diversity"]) > 2 || !all(round(tol, 2) %in% v$tolerance)) stop("Unexpected tolerance table", call. = FALSE)
readr::write_csv(v, out_csv)
print(as.data.frame(v %>% filter(tolerance %in% round(tol, 2))), row.names = FALSE)

GROUPS <- c("Bag seine, estuaries", "Bag seine, coastwide", "Bay trawl, estuaries", "Bay trawl, coastwide")
G_COL <- setNames(GEAR_COLS[c("Bag seine", "Bag seine", "Bay trawl", "Bay trawl")], GROUPS)
G_LT  <- setNames(c("solid", "22", "solid", "22"), GROUPS)
G_SH  <- setNames(c(16, 1, 17, 2), GROUPS)
d <- v %>% pivot_longer(c(at_least_1, all_24), names_to = "what", values_to = "share") %>%
  mutate(What = factor(ifelse(what == "at_least_1", "At least 1 month ahead", "All 24 months"), c("At least 1 month ahead", "All 24 months")),
         Group = factor(paste0(GEAR_LAB[GEAR], ", ", tolower(SCALE_LAB[scale])), GROUPS))
pan <- function(cls) {
  dd <- d %>% filter(class == cls)
  xl <- if (cls == "diversity") "Tolerance (% of the series mean)" else "Tolerance (log units)"
  xf <- if (cls == "diversity") function(z) paste0(100 * z, "%") else function(z) format(z)
  ggplot(dd, aes(tolerance, share, colour = Group, linetype = Group, shape = Group)) +
    geom_vline(xintercept = tol[[cls]], colour = PAL[["grey_mid"]], linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.5, stroke = 0.6) +
    facet_wrap(~ What, ncol = 1) +
    scale_colour_manual(values = G_COL, name = NULL) + scale_linetype_manual(values = G_LT, name = NULL) +
    scale_shape_manual(values = G_SH, name = NULL) +
    scale_x_continuous(breaks = sort(unique(dd$tolerance)), labels = xf) +
    scale_y_continuous(limits = c(0, 1), labels = function(z) paste0(100 * z, "%")) +
    labs(title = CLASS_LAB[[cls]], x = xl, y = "Percentage of series") + theme_paper()
}
fig <- (pan("diversity") | pan("abundance")) + plot_layout(guides = "collect") & guides(colour = guide_legend(nrow = 1)) &
  theme(legend.position = "bottom")
save_figure(fig, paste0(FIGURE_LABEL, "_horizon_tolerance"), W_MM, H_MM)
cat("Saved:", out_csv, "\n")
end_log("finished")
cat("g14_horizon_tolerance.R completed\n")
