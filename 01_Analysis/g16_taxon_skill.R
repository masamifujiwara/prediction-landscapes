# =============================================================================
# Script:  g16_taxon_skill.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S8 (Online Resource 3): skill of SARIMA (automatic) abundance
#   forecasts against the cumulative average for each taxon, estuary and
#   sampling method, 1 month (circle) and 12 months (triangle) ahead, joined
#   by a line; filled when the bootstrap range lies above zero (clearly better
#   than the cumulative average), open otherwise. Panels: sampling method
#   (rows) by estuary, north to south, plus coastwide (columns). Taxa are
#   ordered by their median one-month skill over all panels.
#
# Run after: b09, f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09, estuaries)
#
# Outputs:
#   01_Analysis/figures/Figure_S8_taxon_skill.tiff and .pdf
#   01_Analysis/Data/g16_taxon_skill/g16_taxon_skill_values.csv
#   01_Analysis/Logs/g16_taxon_skill/run_<timestamp>.txt
#   Earlier copies of the figure and CSV file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Skill against the cumulative average, the paper's relative measure,
#     rather than RMSE divided by the standard deviation around the long-run
#     mean, which uses the future level.
#   - Scores from f01 (common forecast dates); coastwide is the last column.
#   - Filled or open from the bootstrap range (the paper's range rule).
#   - Estuary names are set in BAY_NAMES; the script stops if an estuary in
#     the settings has no name there.
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
need <- c("dplyr", "ggplot2", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

FIGURE_LABEL <- "Figure_S8"   # Online Resource 3 figure number (Fig. S8)
script_stem <- "g16_taxon_skill"
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 120
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_lead, in_settings), c("f01_forecast_scoring.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
out_csv <- file.path(dirs$data, "g16_taxon_skill_values.csv"); if (file.exists(out_csv)) archive_files(out_csv)

BAY_NAMES <- c("2" = "Galveston Bay", "3" = "Matagorda Bay", "4" = "San Antonio Bay", "5" = "Aransas Bay",
               "6" = "Corpus Christi Bay", "7" = "Upper Laguna Madre", "8" = "Lower Laguna Madre")
BAYS <- as.character(settings$bays)
if (!all(BAYS %in% names(BAY_NAMES))) stop("A bay in the settings has no name in BAY_NAMES", call. = FALSE)
COLS <- c(unname(BAY_NAMES[BAYS]), "Coastwide")

d <- readRDS(in_lead) %>% filter(class == "abundance", model == MODEL, lead %in% c(1L, 12L)) %>%
  transmute(GEAR, Common_Name, lead, skill = skill_past, lo = skill_past_lo, hi = skill_past_hi,
            where = ifelse(scale == "coastwide", "Coastwide", unname(BAY_NAMES[as.character(BAY_ID)])))
if (nrow(d) != 128L * 2L || anyNA(d$skill) || anyNA(d$where)) stop("Expected 128 abundance series at leads 1 and 12", call. = FALSE)
ord <- d %>% filter(lead == 1L) %>% group_by(Common_Name) %>% summarise(m = median(skill), .groups = "drop") %>% arrange(m) %>% pull(Common_Name)
d <- d %>% mutate(clear = lo > 0, Gear = factor(GEAR_LAB[GEAR], GEAR_LAB), Where = factor(where, COLS),
                  Taxon = factor(Common_Name, ord), Lead = factor(ifelse(lead == 1L, "1 month", "12 months"), c("1 month", "12 months")),
                  fill_key = ifelse(clear, as.character(Gear), "open"))
readr::write_csv(d %>% select(GEAR, where, Common_Name, lead, skill, lo, hi, clear), out_csv)
print(as.data.frame(d %>% group_by(GEAR, lead, scale = Where == "Coastwide") %>% summarise(series = n(), clear = sum(clear), median_skill = median(skill), .groups = "drop")), row.names = FALSE)

segs <- d %>% select(Gear, Where, Taxon, lead, skill) %>% tidyr::pivot_wider(names_from = lead, values_from = skill, names_prefix = "h")
p <- ggplot(d, aes(skill, Taxon)) +
  geom_vline(xintercept = 0, colour = PAL[["grey_mid"]], linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = segs, aes(x = h1, xend = h12, y = Taxon, yend = Taxon, colour = Gear), linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(shape = Lead, colour = Gear, fill = fill_key), size = 1.5, stroke = 0.5) +
  facet_grid(Gear ~ Where, labeller = labeller(Where = label_wrap_gen(10))) +
  scale_shape_manual(values = c("1 month" = 21, "12 months" = 24), name = NULL) +
  scale_colour_manual(values = GEAR_COLS, guide = "none") +
  scale_fill_manual(values = c(GEAR_COLS, open = "white"), breaks = c("Bag seine", "open"),
                    labels = c("Filled: 95% range above zero", "Open: 95% range includes zero"), name = NULL) +
  scale_x_continuous(breaks = c(-0.4, 0, 0.4)) +
  labs(x = "Skill against the cumulative average (proportion of squared error removed)", y = NULL) + theme_paper() +
  theme(strip.text.x = element_text(size = FONT_SIZE - 1), panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey40", colour = "grey40")),
         fill = guide_legend(override.aes = list(shape = 21, colour = "grey40", fill = c("grey40", "white"))))
save_figure(p, paste0(FIGURE_LABEL, "_taxon_skill"), W_MM, H_MM)
cat("Saved:", out_csv, "\n")
end_log("finished")
cat("g16_taxon_skill.R completed\n")
