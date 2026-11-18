# =============================================================================
# Script:  h16_published_measure.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S5 (Online Resource 2): accuracy of SARIMA forecasts on a
#   second measure, RMSE / SD, 1 and 12 months ahead, for every sampling
#   method and target. RMSE / SD is the RMSE of the out-of-sample prediction
#   errors divided by the standard deviation of the departures from the
#   seasonal curve over the same target months (below 1 = more accurate than
#   the mean of those months, which is not available at the forecast date).
#   Estuaries: median over the seven estuaries; coastwide: the single series.
#   It also writes a values file that compares the revised values with those
#   of the submitted manuscript, for the response letter only.
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01; the table)
#   01_Analysis/Data/f1_forecast_scoring/f1_published_check.rds (f01; values file only)
#   01_Analysis/Data/reference_submitted/submitted_table_1_accuracy.csv,
#     submitted_table_s4_bay_abundance_accuracy.csv,
#     submitted_table_s5_skill_naive.csv (values copied from the tables of the
#     submitted version, kept in the repository; values file only)
#
# Outputs:
#   01_Analysis/tables/Table_S5_second_measure.docx
#   01_Analysis/Data/h16_published_measure/published_measure.csv (revised
#     against submitted values)
#   01_Analysis/Logs/h16_published_measure/run_<timestamp>.txt
#   Earlier copies of these files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The table is calculated on the common forecast dates (f1_by_lead:
#   rmse_scaled / sd_obs_scaled per series, then the median), so it matches
#   every other table and figure.
#   The values file uses f01's published check instead (all forecast dates,
#   as in the submitted measure) and the SARIMA rows of the submitted tables,
#   kept as CSV files so the script does not depend on where old outputs are.
#   Skill against last month's value is converted to the submitted
#   definition, 1 minus the RMSE ratio (f01 stores 1 minus the MSE ratio),
#   per series before the median.
#   The menhaden target differs: submitted Gulf menhaden (code 604), revised
#   Menhaden, Brevoortia spp. (604, 658, 161); the rows are matched by name.
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
need <- c("dplyr", "tidyr", "readr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(officer); library(flextable)})

TABLE_LABEL <- "Table S5"   # Online Resource 2 number
script_stem <- "h16_published_measure"
REF <- file.path(PATHS$data_dir, "reference_submitted")
inputs <- c(lead = file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds"),
            pub = file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_published_check.rds"),
            t1 = file.path(REF, "submitted_table_1_accuracy.csv"), s4 = file.path(REF, "submitted_table_s4_bay_abundance_accuracy.csv"),
            s5 = file.path(REF, "submitted_table_s5_skill_naive.csv"))
dirs <- script_dirs(script_stem, c("data", "logs"))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, c("f01_forecast_scoring.R", "f01_forecast_scoring.R", rep("(copied from the submitted version's tables)", 3)))
out_docx <- file.path(tab_dir, "Table_S5_second_measure.docx"); out_csv <- file.path(dirs$data, "published_measure.csv")
old <- c(out_docx, out_csv)[file.exists(c(out_docx, out_csv))]; if (length(old)) archive_files(old)
rd <- function(f) readr::read_csv(f, show_col_types = FALSE, progress = FALSE, col_types = readr::cols(.default = "c"))
num <- function(x) suppressWarnings(as.numeric(x))

GEARS <- c(BS = "Bag seine", BT = "Bay trawl")
Q_LAB <- c("0" = "q = 0", "1" = "q = 1", "2" = "q = 2")   # plain labels
# target key: q for diversity, common name for abundance (submitted Gulf menhaden = revised Menhaden)
key_of <- function(class, q, name) ifelse(class == "diversity", paste0("q", q), ifelse(name %in% c("Gulf menhaden", "Menhaden"), "Menhaden", name))

# -----------------------------------------------------------------------------
# 1) Revised values (SARIMA automatic, all forecast dates)
# -----------------------------------------------------------------------------
pub <- readRDS(inputs[["pub"]]) %>% filter(model == "sarima_auto", lead %in% c(1L, 12L)) %>%
  mutate(key = key_of(class, q, Common_Name),
         skill_naive = 100 * (1 - sqrt(1 - skill_naive)))   # the submitted definition: 1 - RMSE / RMSE of last month's value (f01: 1 - MSE / MSE)
if (nrow(pub) != 176L * 2L) stop("Expected 176 series at leads 1 and 12", call. = FALSE)
new <- pub %>% group_by(class, scale, GEAR, key, lead) %>%
  summarise(series = n(), rmse_sd = median(rmse_over_sd), skill = median(skill_naive), .groups = "drop") %>% mutate(source = "revised")

# -----------------------------------------------------------------------------
# 2) Submitted values (SARIMA rows)
# -----------------------------------------------------------------------------
t1 <- rd(inputs[["t1"]]); s4 <- rd(inputs[["s4"]]); s5 <- rd(inputs[["s5"]])
sub_rmse <- bind_rows(t1 %>% mutate(scale = ifelse(scale == "Alpha", "bay", "coastwide")),
                      s4 %>% mutate(scale = "bay")) %>%
  transmute(class, scale, GEAR, key = key_of(class, target, target_label), `1` = num(rmse_sd_h1), `12` = num(rmse_sd_h12)) %>%
  pivot_longer(c(`1`, `12`), names_to = "lead", values_to = "rmse_sd") %>% mutate(lead = as.integer(lead))
sub_skill <- s5 %>% filter(model_class == "SARIMA") %>%
  transmute(class, scale = ifelse(scale == "Alpha", "bay", "coastwide"), GEAR, key = key_of(class, target, target_label),
            lead = as.integer(horizon), skill = num(value))
sub <- full_join(sub_rmse, sub_skill, by = c("class", "scale", "GEAR", "key", "lead")) %>% mutate(source = "submitted")
if (anyDuplicated(sub[c("class", "scale", "GEAR", "key", "lead")])) stop("Duplicate rows in the submitted tables", call. = FALSE)
both <- bind_rows(new %>% select(-series), sub) %>%
  pivot_wider(id_cols = c(class, scale, GEAR, key, lead), names_from = source, values_from = c(rmse_sd, skill))
miss_sub <- both %>% filter(is.na(rmse_sd_submitted) | is.na(skill_submitted))
if (nrow(miss_sub)) { print(as.data.frame(miss_sub)); stop("Rows without a submitted value (listed above)", call. = FALSE) }
if (nrow(both) != 22L * 2L * 2L) stop("Expected 22 gear-target rows x 2 scales x 2 leads, found ", nrow(both), call. = FALSE)
readr::write_csv(both %>% mutate(across(where(is.double), ~ round(.x, 4))), out_csv)
cmp <- both %>% group_by(scale, lead) %>% summarise(rmse_sd_change = median(rmse_sd_revised - rmse_sd_submitted),
                                                     skill_change = median(skill_revised - skill_submitted), .groups = "drop")
print(as.data.frame(cmp), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3) Word table: RMSE / SD on the common forecast dates (f1_by_lead)
# -----------------------------------------------------------------------------
cm <- readRDS(inputs[["lead"]]) %>% filter(model == "sarima_auto", lead %in% c(1L, 12L)) %>%
  mutate(key = key_of(class, q, Common_Name), rmse_sd = rmse_scaled / sd_obs_scaled)
if (nrow(cm) != 176L * 2L || anyNA(cm$rmse_sd)) stop("Expected 176 series at leads 1 and 12 with RMSE / SD", call. = FALSE)
cm <- cm %>% group_by(class, scale, GEAR, key, lead) %>% summarise(series = n(), rmse_sd = median(rmse_sd), .groups = "drop")
TAXA <- cm %>% filter(class == "abundance") %>% distinct(key) %>% pull(key) %>% sort()
row_order <- c(paste0("q", 0:2), TAXA)
label_of <- function(k) ifelse(grepl("^q[0-2]$", k), Q_LAB[sub("q", "", k)], k)
tab_a <- cm %>% transmute(class, GEAR, key, col = paste(scale, lead, sep = "_"), txt = formatC(rmse_sd, format = "f", digits = 2)) %>%
  pivot_wider(names_from = col, values_from = txt) %>%
  mutate(key = factor(key, row_order), class = factor(class, c("diversity", "abundance"))) %>% arrange(GEAR, class, key) %>%
  transmute(`Sampling method` = GEARS[GEAR], Target = label_of(as.character(key)),
            `Estuaries, 1 month` = bay_1, `Estuaries, 12 months` = bay_12, `Coastwide, 1 month` = coastwide_1, `Coastwide, 12 months` = coastwide_12)
if (nrow(tab_a) != 22L || anyNA(tab_a)) stop("Expected 22 complete rows (2 sampling methods x 11 targets)", call. = FALSE)
cat("RMSE / SD, common forecast dates:\n"); print(as.data.frame(tab_a), row.names = FALSE)
style_ft <- function(df) {
  flextable(as.data.frame(df)) %>% font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8, part = "all") %>%
    bold(part = "header") %>% align(j = 3:6, align = "center", part = "all") %>% border_remove() %>%
    hline_top(border = fp_border(width = 1.5), part = "header") %>% hline_bottom(border = fp_border(width = 0.75), part = "header") %>%
    hline_bottom(border = fp_border(width = 1.5), part = "body") %>% merge_v(j = 1) %>% valign(j = 1, valign = "top", part = "body") %>%
    padding(padding.top = 0, padding.bottom = 0, padding.left = 3, padding.right = 3, part = "all") %>%
    width(width = c(1.0, 1.5, 1.0, 1.0, 1.0, 1.0))
}
caption <- paste0(
  "Accuracy of SARIMA forecasts on a second measure, 1 and 12 months ahead: the root mean squared error (RMSE) of ",
  "the out-of-sample prediction errors divided by the standard deviation of the departures from the seasonal curve over the same ",
  "target months. Values below 1 mean that the forecasts were more accurate than the mean of the target months, which is not ",
  "available at the forecast date. Calculated on the 367 forecast dates up to December 2022, as in the main text. Estuaries: median ",
  "over the seven estuaries; coastwide: the single series. Biodiversity: Hill numbers; abundance: log(CPUE + 1), CPUE being catch per ",
  "unit effort.")
note_a <- "Menhaden is the genus "
note_i <- "Brevoortia"
note_b <- " (Gulf menhaden, finescale menhaden, and records identified only to the genus, combined), as in all abundance results."
ptxt <- function(s, size = 10, bold = FALSE, keep = FALSE) fpar(ftext(s, fp_text(font.size = size, bold = bold, font.family = "Times New Roman")),
                                                                fp_p = fp_par(keep_with_next = keep))
doc <- read_docx() %>%
  body_set_default_section(prop_section(page_size = page_size(orient = "portrait", width = 8.5, height = 11),
                                        page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75))) %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman")),
                     ftext(caption, fp_text(font.size = 10, font.family = "Times New Roman")),
                     fp_p = fp_par(keep_with_next = FALSE))) %>%
  body_add_flextable(style_ft(tab_a)) %>%
  body_add_fpar(fpar(ftext(note_a, fp_text(font.size = 9, font.family = "Times New Roman")),
                     ftext(note_i, fp_text(font.size = 9, italic = TRUE, font.family = "Times New Roman")),
                     ftext(note_b, fp_text(font.size = 9, font.family = "Times New Roman"))))
print(doc, target = out_docx)
for (f in c(out_docx, out_csv)) cat("Saved:", f, "\n")
end_log("finished")
cat("h16_published_measure.R completed\n")
