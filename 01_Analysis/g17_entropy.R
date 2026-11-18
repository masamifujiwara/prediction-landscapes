# =============================================================================
# Script:  g17_entropy.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S4 (Online Resource 3): intrinsic predictability from weighted
#   permutation entropy (pattern length 4; 1 minus the entropy) against the
#   skill of SARIMA (automatic) against the cumulative average, one point per
#   series. Top: the one-month step against skill 1 month ahead; bottom: the
#   twelve-month step (the same calendar month across years) against skill
#   12 months ahead. Diversity left, abundance right. Dashed line: median of
#   the shuffled copies (no order at all). Panel text: rank correlations for
#   estuaries and coastwide, with the number of series.
#
# Run after: f01, f03
#
# Inputs:
#   01_Analysis/Data/f3_permutation_entropy/f3_entropy_by_series.rds,
#     f3_series_checks.rds, f3_correlations.rds (f03)
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_S4_entropy.tiff and .pdf
#   01_Analysis/Data/g17_entropy/g17_entropy_points.csv,
#     g17_entropy_correlations.csv
#   01_Analysis/Logs/g17_entropy/run_<timestamp>.txt
#   Earlier copies of the figure and CSV files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Pattern length 4; both steps are reported. Permutation entropy is not
#     presented as evidence of chaos or density dependence (caption).
#   - Rank correlations are taken from f03, not recomputed here.
#   - Estuary series flagged for runs of zero catches are kept and not marked
#     in the figure; the correlations without them are in the values file and
#     the caption.
#   - Gear by color and symbol (circles bag seine, triangles bay trawl),
#     estuaries filled, coastwide open.
#   - The numbers quoted in the Fig. S4 caption (g20) are checked against the
#     current outputs; the script stops if they differ.
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
need <- c("dplyr", "ggplot2", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

FIGURE_LABEL <- "Figure_S4"   # Online Resource 3 figure number (Fig. S4)
script_stem <- "g17_entropy"
MODEL <- "sarima_auto"; M <- 4L; W_MM <- 170; H_MM <- 150
F3 <- file.path(PATHS$data_dir, "f3_permutation_entropy")
inputs <- c(ent = file.path(F3, "f3_entropy_by_series.rds"), chk = file.path(F3, "f3_series_checks.rds"),
            cor = file.path(F3, "f3_correlations.rds"), lead = file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds"))
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, c(rep("f03_permutation_entropy.R", 3), "f01_forecast_scoring.R"))
out <- file.path(dirs$data, c("g17_entropy_points.csv", "g17_entropy_correlations.csv"))
old <- out[file.exists(out)]; if (length(old)) archive_files(old)

keys <- c("class", "scale", "BAY_ID", "GEAR", "target")
ent <- readRDS(inputs[["ent"]]) %>% filter(m == M)
chk <- readRDS(inputs[["chk"]]) %>% select(all_of(keys), zero_flag)
sk <- readRDS(inputs[["lead"]]) %>% filter(model == MODEL, lead %in% c(1L, 12L)) %>% select(all_of(keys), lead, skill_past)
pts <- ent %>% mutate(lead = as.integer(step)) %>% inner_join(sk, by = c(keys, "lead")) %>% left_join(chk, by = keys) %>%
  mutate(open = zero_flag & step == 1L & scale == "bay")
if (nrow(pts) != 176L * 2L || anyNA(pts$skill_past) || anyNA(pts$predictability)) stop("Expected 176 series at both steps with skill", call. = FALSE)
cr <- readRDS(inputs[["cor"]]) %>% filter(m == M, (step == 1L & result == "skill_past_1") | (step == 12L & result == "skill_past_12"))
if (nrow(cr) != 8L) stop("Expected 8 rank correlations from f03", call. = FALSE)
readr::write_csv(pts %>% select(all_of(keys), step, wpe, predictability, wpe_shuffle_median, skill_past, zero_flag), out[1])
readr::write_csv(cr, out[2])
print(as.data.frame(cr %>% select(class, scale, step, n, spearman, spearman_unflagged)), row.names = FALSE)

GROUPS <- c("Bag seine, estuaries", "Bag seine, coastwide", "Bay trawl, estuaries", "Bay trawl, coastwide")
G_COL <- setNames(GEAR_COLS[c("Bag seine", "Bag seine", "Bay trawl", "Bay trawl")], GROUPS)
G_SH  <- setNames(c(16, 1, 17, 2), GROUPS)
# Numbers quoted in the Fig. S4 caption (g20): zero-catch bay series and correlations without them.
n_zero <- sum(pts$zero_flag[pts$step == 1L & pts$scale == "bay" & pts$class == "abundance"])
unf <- cr %>% filter(class == "abundance", scale == "bay") %>% arrange(step) %>% pull(spearman_unflagged)
if (n_zero != 24L || !isTRUE(all.equal(round(unf, 2), c(0.54, 0.31))))
  stop("Zero-catch figures differ from the Figure S4 caption in g20: ", n_zero, ", ", paste(round(unf, 2), collapse = ", "), call. = FALSE)
q90 <- pts %>% group_by(class, step) %>% summarise(q = stats::quantile(predictability, 0.9), .groups = "drop")
if (any(q90$q >= 0.06)) stop("Not nine in ten series below 0.06 in every class and step; the Figure S4 caption needs changing", call. = FALSE)
STEP_LAB <- c("1" = "One-month step, skill 1 month ahead", "12" = "Twelve-month step, skill 12 months ahead")
d <- pts %>% mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Gear = factor(GEAR_LAB[GEAR], GEAR_LAB),
                    Scale = factor(SCALE_LAB[scale], SCALE_LAB), Step = factor(STEP_LAB[as.character(step)], STEP_LAB),
                    Group = factor(paste0(GEAR_LAB[GEAR], ", ", tolower(SCALE_LAB[scale])), GROUPS))
ref <- d %>% group_by(Class, Step) %>% summarise(x = median(1 - wpe_shuffle_median), .groups = "drop")
lab <- cr %>% mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Step = factor(STEP_LAB[as.character(step)], STEP_LAB)) %>%
  group_by(Class, Step) %>% arrange(scale, .by_group = TRUE) %>%
  summarise(txt = paste0("Spearman: ", paste0(ifelse(scale == "bay", "estuaries ", "coastwide "),
                                                       sub("^-", "\u2212", sprintf("%.2f", spearman)), " (n = ", n, ")", collapse = "; ")), .groups = "drop")
p <- ggplot(d, aes(predictability, skill_past)) +
  geom_hline(yintercept = 0, colour = PAL[["grey_light"]], linewidth = 0.4) +
  geom_vline(data = ref, aes(xintercept = x), colour = PAL[["grey_mid"]], linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(colour = Group, shape = Group), size = 1.5, stroke = 0.6) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = txt), hjust = 1.03, vjust = 1.5, size = (FONT_SIZE - 1.5) / .pt, colour = PAL[["grey_dark"]]) +
  facet_grid(Step ~ Class, scales = "free", labeller = labeller(Step = label_wrap_gen(22))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +   # room for the correlation label
  scale_colour_manual(values = G_COL, name = NULL) + scale_shape_manual(values = G_SH, name = NULL) +
  labs(x = "Intrinsic predictability (1 \u2212 weighted permutation entropy)",
       y = "Skill against the cumulative average\n(proportion of squared error removed)") + theme_paper() + theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(5, "mm"))   # keeps the end labels of the two columns apart
save_figure(p, paste0(FIGURE_LABEL, "_entropy"), W_MM, H_MM)
for (f in out) cat("Saved:", f, "\n")
end_log("finished")
cat("g17_entropy.R completed\n")
