# =============================================================================
# Script:  g11_seasonal_fits.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Online Resource 3, Figure S1: the seasonal curves and the departures
#   from the seasonal curve, for eight example series. Each panel: top, the
#   observed monthly values (grey) with the seasonal curve (color); bottom,
#   the departures from the curve, which are what the forecast models
#   forecast. Left column bag seine, right column bay trawl, all coastwide.
#   Rows 1 to 3: diversity, q = 0, 1, 2. Row 4: brown shrimp, log(CPUE + 1).
#
# Run after: c02, e03, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds (c02)
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds (e03)
#   01_Analysis/Data/e3_deseasonalize_abundance/gamma_abundance_deseas.rds (e03)
#
# Outputs:
#   01_Analysis/figures/Figure_S1_seasonal_fits.tiff and .pdf (name set by
#     FIGURE_LABEL)
#   01_Analysis/Data/g11_seasonal_fits/g11_seasonal_fits_values.csv
#   01_Analysis/Logs/g11_seasonal_fits/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The curve shown is the whole-record curve (three harmonics: periods of
#   12, 6 and 4 months), for description only. The forecasts use a curve
#   refitted at each forecast date on the months up to that date; the caption
#   says so. Curve and departures are read from the c02 and e03 series files
#   (value = seasonal_component + value_deseas, checked here), so the figure
#   shows exactly the series the other scripts use. April and May 2020 (not
#   sampled, filled for model fitting) are left as a gap in every line,
#   because they are not observations. Columns are sampling methods in every
#   row and every panel is coastwide; the seasonal fits for the estuaries are
#   in Online Resource 2, Table S4 (h14). Brown shrimp is the abundance
#   example because it has a strong, regular seasonal cycle in both sampling
#   methods. Look, colors and saving come from lib/figures.R (bag seine blue,
#   bay trawl orange).
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
need <- c("dplyr", "ggplot2", "patchwork", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})

FIGURE_LABEL <- "Figure_S1"   # figure number in Online Resource 3
script_stem <- "g11_seasonal_fits"
D <- PATHS$data_dir
inputs <- c(div_a = file.path(D, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds"),
            div_g = file.path(D, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds"),
            abd_a = file.path(D, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds"),
            abd_g = file.path(D, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds"))
W_MM <- 170
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, c(rep("c02_deseasonalize_alpha_gamma_scaled.R", 2), rep("e03_deseasonalize_abundance.R", 2)))
out_vals <- file.path(dirs$data, "g11_seasonal_fits_values.csv")
if (file.exists(out_vals)) archive_files(out_vals)

# -----------------------------------------------------------------------------
# 1) The eight series
# -----------------------------------------------------------------------------
div_a <- readRDS(inputs[["div_a"]]); div_g <- readRDS(inputs[["div_g"]])
abd_a <- readRDS(inputs[["abd_a"]]); abd_g <- readRDS(inputs[["abd_g"]])
Q_LAB <- c("0" = "q = 0", "1" = "q = 1", "2" = "q = 2")   # plain labels
# Abundance examples: coastwide, both sampling methods (columns are sampling methods in every row).
ABD_TAXA <- c("9618" = "Brown shrimp")   # strong, regular seasonal cycle in both sampling methods
panels <- bind_rows(
  lapply(c("BS", "BT"), function(g) lapply(0:2, function(qq)
    div_g %>% filter(GEAR == g, q == qq) %>%
      mutate(panel = paste0(Q_LAB[[as.character(qq)]], "\n", GEAR_LAB[[g]], ", coastwide"), y_lab = "Biodiversity", row = qq + 1L))) %>% unlist(recursive = FALSE),
  lapply(c("BS", "BT"), function(g) lapply(seq_along(ABD_TAXA), function(k)
    abd_g %>% filter(GEAR == g, SPECIES_CODE == names(ABD_TAXA)[k]) %>%
      mutate(panel = paste0(ABD_TAXA[[k]], ", log(CPUE + 1)\n", GEAR_LAB[[g]], ", coastwide"), y_lab = "log(CPUE + 1)", row = 3L + k))) %>%
    unlist(recursive = FALSE)) %>%
  transmute(panel, row, y_lab, GEAR, date, filled, value, seasonal_component, value_deseas)
N_ROWS <- 3L + length(ABD_TAXA)
n_by <- table(panels$panel)
if (length(n_by) != 2L * N_ROWS || any(n_by < 500)) stop("Expected ", 2L * N_ROWS, " series of at least 500 months", call. = FALSE)
if (max(abs(panels$value - panels$seasonal_component - panels$value_deseas)) > 1e-8)
  stop("value is not seasonal_component + value_deseas", call. = FALSE)
# April and May 2020 are filled values: drawn as a gap
panels <- panels %>% mutate(across(c(value, value_deseas), ~ ifelse(filled, NA_real_, .x)))
readr::write_csv(panels, out_vals)

# -----------------------------------------------------------------------------
# 2) One panel: observed with curve on top, departures below
# -----------------------------------------------------------------------------
one_panel <- function(nm, tag, x_title) {
  d <- panels %>% filter(panel == nm); col <- GEAR_COLS[[GEAR_LAB[[d$GEAR[1]]]]]
  top <- ggplot(d, aes(date)) +
    geom_line(aes(y = value), colour = PAL[["grey_mid"]], linewidth = 0.3) +
    geom_line(aes(y = seasonal_component), colour = col, linewidth = 0.35) +
    labs(title = nm, tag = tag, x = NULL, y = d$y_lab[1]) + theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = FONT_SIZE, lineheight = 0.95, margin = margin(b = 2)),
          plot.tag = element_text(size = FONT_SIZE + 1, face = "bold"), plot.margin = margin(4, 4, 1, 2))
  bot <- ggplot(d, aes(date)) +
    geom_hline(yintercept = 0, colour = PAL[["grey_light"]], linewidth = 0.4) +
    geom_line(aes(y = value_deseas), colour = col, linewidth = 0.3) +
    scale_y_continuous(breaks = function(l) pretty(l, 3)) +
    labs(x = if (x_title) "Year" else NULL, y = "Departure") + theme_paper() + theme(plot.margin = margin(1, 4, 2, 2))
  top / bot + plot_layout(heights = c(1, 0.6))
}
key <- panels %>% distinct(panel, row, GEAR) %>% arrange(GEAR, row)
order_left <- key$panel[key$GEAR == "BS"]; order_right <- key$panel[key$GEAR == "BT"]
last <- seq_len(N_ROWS) == N_ROWS
P <- c(mapply(one_panel, order_left, letters[seq_len(N_ROWS)], last, SIMPLIFY = FALSE),
       mapply(one_panel, order_right, letters[N_ROWS + seq_len(N_ROWS)], last, SIMPLIFY = FALSE))
fig <- wrap_plots(wrap_plots(P[seq_len(N_ROWS)], ncol = 1), wrap_plots(P[N_ROWS + seq_len(N_ROWS)], ncol = 1), ncol = 2)
H_MM <- min(50 * N_ROWS, 205)   # 200 mm for four rows; capped by the Online Resources page (lib/figures.R)
save_figure(fig, paste0(FIGURE_LABEL, "_seasonal_fits"), W_MM, H_MM)
cat("Saved:", out_vals, "\n")
end_log("finished")
cat("g11_seasonal_fits.R completed\n")
