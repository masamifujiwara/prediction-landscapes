# =============================================================================
# Script:  h01_models.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table 1 of the main text: every forecasting model and simple
#   reference (rows, with its role: base model, covariates, trend test,
#   reference) against the eight data groups (columns: biodiversity and
#   abundance, estuaries and coastwide, bag seine and bay trawl). One measure
#   is shown, the median skill against the cumulative average, in two blocks:
#   one month and 12 months ahead.
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#
# Outputs:
#   01_Analysis/tables/Table_1_models.csv (values, rounded to three decimals)
#   01_Analysis/tables/Table_1_models.docx (caption above the table)
#   01_Analysis/Logs/h1_models/run_<timestamp>.txt
#   An earlier Table_1_models.docx is moved to ARCHIVE_<timestamp> first.
#
# Method notes:
#   Skill is used because it is comparable between biodiversity and abundance
#   (RMSE units differ). The full results for all models and measures are in
#   Tables S7 and S8 (h15). The caption explains how to read the skill values
#   (0 = no better than the cumulative average, negative = worse) and points to
#   the Methods for the definition of skill and the forecast design.
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
source(file.path(PATHS$analysis_dir, "lib", "figures.R"))   # labels
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(officer); library(flextable)})
script_stem <- "h1_models"
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
LEADS <- c(1L, 12L)
MODELS <- c(sarima_auto = "SARIMA (automatic)", ssm_level = "State-space",
            arimax = "ARIMAX (lagged temperature, salinity)", ssm_xreg = "State-space with lagged temperature, salinity",
            sarima_notrend = "SARIMA, no trend", sarima_lintrend = "SARIMA, linear trend", damped_trend = "Damped trend",
            ssm_lintrend = "State-space with trend",
            ref_naive = "Last month's value", ref_12month = "12-month average")
ROLE <- c(sarima_auto = "Base model", ssm_level = "Base model", arimax = "Covariates", ssm_xreg = "Covariates",
          sarima_notrend = "Trend test", sarima_lintrend = "Trend test", damped_trend = "Trend test",
          ssm_lintrend = "Trend test", ref_naive = "Reference", ref_12month = "Reference")
TABLE_LABEL <- "Table 1"   # main text
CAPTION <- paste(
  "Skill of each forecasting model and reference against the cumulative average, one and 12 months ahead",
  "(median over the series of each group; 0 = no better than the cumulative average, negative = worse;",
  "Methods: Forecast accuracy and skill). Number of series per column: biodiversity 21 in estuaries and 3 coastwide;",
  "abundance 56 in estuaries and 8 coastwide. The models are described in Online Resource 1 (Sect. 1.9); full results",
  "with uncertainty ranges are in Tables S7 and S8 (Online Resource 2).")

dirs <- script_dirs(script_stem, c("logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_lead, "f01_forecast_scoring.R")
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(tab_dir, "Table_1_models.csv"); out_docx <- file.path(tab_dir, "Table_1_models.docx")

v <- readRDS(in_lead) %>% filter(model %in% names(MODELS), lead %in% LEADS) %>%
  group_by(class, scale, GEAR, model, lead) %>% summarise(series = n(), skill = stats::median(skill_past), .groups = "drop") %>%
  mutate(col = paste0(c(diversity = "Div", abundance = "Abu")[class], "_", scale, "_", GEAR))   # fixed keys, independent of the display labels in CLASS_LAB
cols <- as.vector(outer(c("bay_BS", "bay_BT", "coastwide_BS", "coastwide_BT"), c("Div", "Abu"), function(a, b) paste0(b, "_", a)))
wide <- v %>% select(model, lead, col, skill) %>% pivot_wider(names_from = col, values_from = skill) %>%
  mutate(model = factor(model, names(MODELS))) %>% arrange(lead, model) %>% select(lead, model, all_of(cols))
if (nrow(wide) != length(MODELS) * length(LEADS)) stop("Unexpected number of rows", call. = FALSE)
readr::write_csv(wide %>% mutate(across(where(is.double), ~ round(.x, 3))), out_csv)

f2 <- function(x) ifelse(is.na(x), "", sub("^-0[.]00$", "0.00", sprintf("%.2f", x)))
body <- bind_rows(lapply(LEADS, function(h) {
  d <- wide %>% filter(lead == h)
  bind_rows(data.frame(Model = ifelse(h == 1L, "One month ahead", "12 months ahead")),
            data.frame(Model = MODELS[as.character(d$model)], Role = ROLE[as.character(d$model)],
                       d[, cols] %>% mutate(across(everything(), f2)), check.names = FALSE))
})) %>% mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
lead_rows <- which(body$Model %in% c("One month ahead", "12 months ahead"))
ft <- flextable(body) %>%
  set_header_labels(values = setNames(as.list(c("Model", "Role", rep(c("Bag seine", "Bay trawl"), 4))), names(body))) %>%
  add_header_row(values = c("", "", "Estuaries", "Coastwide", "Estuaries", "Coastwide"), colwidths = c(1, 1, 2, 2, 2, 2)) %>%
  add_header_row(values = c("", "", "Biodiversity", "Abundance"), colwidths = c(1, 1, 4, 4)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9, part = "all") %>% bold(part = "header") %>%
  bold(i = lead_rows, part = "body") %>%
  align(j = 2:ncol(body), align = "center", part = "all") %>% border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "header") %>% hline_bottom(border = fp_border(width = 0.75), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  padding(padding.top = 0, padding.bottom = 0, padding.left = 2, padding.right = 2, part = "all")
for (i in lead_rows) ft <- merge_at(ft, i = i, j = seq_len(ncol(body)), part = "body")
ft <- width(ft, width = c(1.95, 0.75, rep(0.475, 8)))
if (file.exists(out_docx)) archive_files(c(out_docx))
doc <- read_docx() %>% body_set_default_section(prop_section(page_size = page_size(width = 8.5, height = 11, orient = "portrait"),
    page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1, gutter = 0))) %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)), ftext(CAPTION, fp_text(font.family = "Times New Roman", font.size = 10)))) %>%
  body_add_flextable(ft)
print(doc, target = out_docx)
cat("Saved:", out_csv, "\nSaved:", out_docx, "\n")
print(as.data.frame(wide %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)
end_log("finished")
cat("h01_models.R completed\n")
