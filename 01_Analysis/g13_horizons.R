# =============================================================================
# Script:  g13_horizons.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S5 (Online Resource 3): forecast horizons of SARIMA (automatic),
#   the number of months ahead for which the forecasts stay (top) within the
#   tolerance, (middle) better than the cumulative average and (bottom) better
#   than the 12-month average. Bars: percentage of the series of each group in
#   five horizon classes (0, 1 to 3, 4 to 11, 12 to 23, all 24 months); in each
#   pair, the range rule (left) and the value-only rule (right).
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_horizons.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_S5_horizons.tiff and .pdf
#   01_Analysis/Data/g13_horizons/g13_horizon_classes.csv
#   01_Analysis/Logs/g13_horizons/run_<timestamp>.txt
#   Earlier copies of the figure and CSV file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Horizons as computed in f01: counted from 1 month ahead to the first
#     month that fails (0 if 1 month ahead fails). Skill that is lost and later
#     returns (return of skill) is not part of the horizon; it is counted
#     separately in the values file and the caption.
#   - Range rule: the whole bootstrap range must meet the condition; value
#     rule: the central estimate alone. Accuracy: RMSE within the tolerance
#     (30% of the series mean for diversity, 0.5 log units for abundance).
#   - Percentages rather than counts, because groups differ in size; the
#     numbers of series are in the panel labels.
#   - The return-of-skill counts are checked against the numbers quoted in the
#     Fig. S5 caption (g20); the script stops if they differ.
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

FIGURE_LABEL <- "Figure_S5"   # Online Resource 3 figure number (Fig. S5)
script_stem <- "g13_horizons"
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 175
in_h <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_horizons.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_h, "f01_forecast_scoring.R")
out_csv <- file.path(dirs$data, "g13_horizon_classes.csv"); if (file.exists(out_csv)) archive_files(out_csv)

h <- readRDS(in_h) %>% filter(model == MODEL)
if (nrow(h) != 176L || any(h$max_lead != 24L)) stop("Expected 176 SARIMA series with 24 leads", call. = FALSE)
TYPES <- c(abs = "Accurate enough (within tolerance)", past = "Better than the cumulative average", m12 = "Better than the 12-month average")
CLS <- c("0", "1 to 3", "4 to 11", "12 to 23", "24")
cut_h <- function(x) factor(as.character(cut(x, c(-Inf, 0, 3, 11, 23, Inf), labels = CLS)), CLS)
long <- h %>% transmute(class, scale, GEAR,
                        abs_range = abs_horizon, abs_value = abs_horizon_value,
                        past_range = rel_horizon_past, past_value = rel_horizon_past_value,
                        m12_range = rel_horizon_12m, m12_value = rel_horizon_12m_value) %>%
  pivot_longer(-c(class, scale, GEAR), names_to = c("type", "rule"), names_sep = "_", values_to = "months")
if (anyNA(long$months)) stop("Missing horizons", call. = FALSE)
grp_n <- h %>% count(class, scale, GEAR, name = "series")
tab <- long %>% mutate(hclass = cut_h(months)) %>% count(class, scale, GEAR, type, rule, hclass, name = "n", .drop = FALSE) %>%
  filter(!is.na(class)) %>% left_join(grp_n, by = c("class", "scale", "GEAR")) %>% mutate(share = n / series)
ret <- h %>% group_by(class, scale, GEAR) %>% summarise(return_of_skill = sum(!is.na(return_of_skill_lead)), .groups = "drop")
readr::write_csv(tab %>% left_join(ret, by = c("class", "scale", "GEAR")), out_csv)
print(as.data.frame(ret), row.names = FALSE)
# Numbers quoted in the Fig. S5 caption (g20): series whose skill against the cumulative average returns.
ret_q <- ret %>% group_by(class, scale) %>% summarise(n = sum(return_of_skill), .groups = "drop")
ret_q <- setNames(ret_q$n, paste(ret_q$class, ret_q$scale, sep = "_"))
CAPTION_RETURN <- c(diversity_bay = 9L, abundance_bay = 11L, diversity_coastwide = 1L, abundance_coastwide = 2L)
if (!identical(unname(ret_q[names(CAPTION_RETURN)]), unname(CAPTION_RETURN)))
  stop("Return-of-skill counts differ from the Figure S5 caption in g20: ", paste(names(ret_q), ret_q, collapse = ", "), call. = FALSE)

d <- tab %>% mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Type = factor(TYPES[type], TYPES),
                    Rule = factor(ifelse(rule == "range", "Range rule", "Value only"), c("Range rule", "Value only")),
                    Group = factor(paste0(SCALE_LAB[scale], ",\n", tolower(GEAR_LAB[GEAR]), " (", series, ")")),
                    x = interaction(Rule, Group, sep = " | ", lex.order = TRUE))
grp_lev <- d %>% distinct(class, scale, GEAR, Group) %>% arrange(scale, GEAR) %>% pull(Group) %>% as.character() %>% unique()
d <- d %>% mutate(Group = factor(as.character(Group), grp_lev))
HCOLS <- c("0" = "#e6e5e1", "1 to 3" = "#b9d3f2", "4 to 11" = "#7fb0e8", "12 to 23" = "#2a78d6", "24" = "#173f73")
p <- ggplot(d, aes(x = Rule, y = share, fill = hclass)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2, position = position_stack(reverse = TRUE)) +
  facet_grid(Type ~ Class + Group, labeller = labeller(Type = label_wrap_gen(18))) +
  scale_fill_manual(values = HCOLS, name = "Horizon (months)") +
  scale_x_discrete(labels = c("Range rule" = "Range", "Value only" = "Value")) +
  scale_y_continuous(labels = function(v) paste0(100 * v, "%"), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Percentage of series") + theme_paper() +
  theme(panel.spacing.x = unit(c(1, 1, 1, 4, 1, 1, 1), "mm"), strip.text.x = element_text(size = FONT_SIZE - 1),
        axis.text.x = element_text(size = FONT_SIZE - 2), panel.grid.major.x = element_blank())
save_figure(p, paste0(FIGURE_LABEL, "_horizons"), W_MM, H_MM)
cat("Saved:", out_csv, "\n")
end_log("finished")
cat("g13_horizons.R completed\n")
