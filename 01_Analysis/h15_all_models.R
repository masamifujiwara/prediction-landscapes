# =============================================================================
# Script:  h15_all_models.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Tables S7 (biodiversity) and S8 (abundance) of Online Resource 2:
#   every forecasting model and simple reference, by scale and sampling
#   method. For each model: median over the series of the RMSE of the
#   out-of-sample prediction errors 1, 12 and 24 months ahead (biodiversity:
#   proportion of the series mean; abundance: log units); median skill
#   against the cumulative average at the same lead times; for the trend and
#   covariate versions, median skill against the base model and the numbers
#   of series clearly better and clearly worse than the base model (bootstrap
#   range of the paired skill above or below zero) 1 and 12 months ahead; and
#   the number of forecast cases lost to failed fits, with their percentage.
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds, f1_failed_cases.rds (f01)
#
# Outputs:
#   01_Analysis/tables/Table_S7_all_models_diversity.docx
#   01_Analysis/tables/Table_S8_all_models_abundance.docx
#   01_Analysis/Data/h15_all_models/all_models.csv (all values, both classes)
#   01_Analysis/Logs/h15_all_models/run_<timestamp>.txt
#   Earlier copies of these files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The base model of a trend or covariate version is the model without the
#   feature tested: SARIMA (automatic) for ARIMAX and the SARIMA and damped
#   trend versions; the state-space model for the state-space versions with
#   covariates or a trend (PARTNER in the code).
#   Covariate models reach 12 months only (n/a at 24). The two covariate sets
#   behind each combined covariate model (lags 1 and 12 at lead 1; lag 12 at
#   leads 2 to 12) are not shown separately.
#   The cumulative average is the reference for skill; its own skill is 0 by
#   definition and is not shown.
#   Biodiversity and abundance are two separate tables rather than two panels
#   of one table.
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
script_stem <- "h15_all_models"
TABLE_LABEL <- c(diversity = "Table S7", abundance = "Table S8")   # Online Resources numbers
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
in_fail <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_failed_cases.rds")
LEADS <- c(1L, 12L, 24L)
MODELS <- c(sarima_auto = "SARIMA (automatic)", ssm_level = "State-space",
            arimax = "ARIMAX (temperature, salinity)", ssm_xreg = "State-space with covariates",
            sarima_notrend = "SARIMA, no trend", sarima_lintrend = "SARIMA, linear trend", damped_trend = "Damped trend",
            ssm_lintrend = "State-space with trend",
            ref_naive = "Last month's value (naive)", ref_12month = "12-month average", ref_past_average = "Cumulative average")
ROLE <- c(sarima_auto = "Base model", ssm_level = "Base model", arimax = "Covariates", ssm_xreg = "Covariates",
          sarima_notrend = "Trend test", sarima_lintrend = "Trend test", damped_trend = "Trend test", ssm_lintrend = "Trend test",
          ref_naive = "Reference", ref_12month = "Reference", ref_past_average = "Reference")
PARTNER <- c(arimax = "sarima_auto", ssm_xreg = "ssm_level", sarima_notrend = "sarima_auto", sarima_lintrend = "sarima_auto",
             damped_trend = "sarima_auto", ssm_lintrend = "ssm_level")

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_lead, in_fail), rep("f01_forecast_scoring.R", 2))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(dirs$data, "all_models.csv")
out_docx <- c(diversity = file.path(tab_dir, "Table_S7_all_models_diversity.docx"),
              abundance = file.path(tab_dir, "Table_S8_all_models_abundance.docx"))
if (file.exists(out_csv)) archive_files(out_csv)

bl <- readRDS(in_lead) %>% filter(model %in% names(MODELS))
miss_m <- setdiff(names(MODELS), unique(bl$model)); if (length(miss_m)) stop("Models missing from f01: ", paste(miss_m, collapse = ", "), call. = FALSE)
long <- bl %>% filter(lead %in% LEADS) %>% group_by(class, scale, GEAR, model, lead) %>%
  summarise(series = n(), rmse = stats::median(rmse_scaled), skill_past = stats::median(skill_past),
            skill_partner = if (model[1] %in% names(PARTNER)) stats::median(skill_partner) else NA_real_,
            better = if (model[1] %in% names(PARTNER)) sum(skill_partner_lo > 0, na.rm = TRUE) else NA_integer_,
            worse = if (model[1] %in% names(PARTNER)) sum(skill_partner_hi < 0, na.rm = TRUE) else NA_integer_, .groups = "drop")
# check the partner column really compares with the intended partner: SARIMA's own partner skill must be absent
fail <- readRDS(in_fail) %>% filter(model %in% names(MODELS)) %>% group_by(class, scale, GEAR, model) %>%
  summarise(failed = sum(cases_failed), .groups = "drop")
# Scored cases a model did forecast (all leads, common forecast dates): n in f01 excludes failed cases,
# so failed / (scored + failed) is the share lost.
scored <- bl %>% group_by(class, scale, GEAR, model) %>% summarise(scored = sum(n), .groups = "drop")
fail <- fail %>% left_join(scored, by = c("class", "scale", "GEAR", "model")) %>% mutate(failed_pct = 100 * failed / (scored + failed))
if (anyNA(fail$failed_pct)) stop("Scored-case count missing for a model", call. = FALSE)
wide <- long %>% pivot_wider(id_cols = c(class, scale, GEAR, model, series), names_from = lead,
                             values_from = c(rmse, skill_past, skill_partner, better, worse)) %>%
  left_join(fail, by = c("class", "scale", "GEAR", "model")) %>%
  mutate(failed = coalesce(failed, 0), failed_pct = coalesce(failed_pct, 0), class = factor(class, c("diversity", "abundance")), scale = factor(scale, c("bay", "coastwide")),
         model = factor(model, names(MODELS))) %>% arrange(class, scale, GEAR, model)
if (nrow(wide) != 8L * length(MODELS)) stop("Expected ", 8L * length(MODELS), " rows", call. = FALSE)
readr::write_csv(wide %>% mutate(across(where(is.double), ~ round(.x, 3))), out_csv)

f2 <- function(x) ifelse(is.na(x), "n/a", sub("^-0[.]00$", "0.00", sprintf("%.2f", x)))
trio <- function(a, b, c) paste(f2(a), f2(b), f2(c), sep = " / ")
disp <- wide %>% transmute(class = as.character(class), grp = paste0(SCALE_LAB[as.character(scale)], ", ", tolower(GEAR_LAB[GEAR])),
  Model = MODELS[as.character(model)], Role = ROLE[as.character(model)],
  `RMSE 1 / 12 / 24` = trio(rmse_1, rmse_12, rmse_24),
  `Skill vs cumulative average 1 / 12 / 24` = ifelse(model == "ref_past_average", "(reference)", trio(skill_past_1, skill_past_12, skill_past_24)),
  `Skill vs base model 1 / 12 / 24` = ifelse(is.na(skill_partner_1), "", trio(skill_partner_1, skill_partner_12, skill_partner_24)),
  `Better / worse than base model at 1; at 12` = ifelse(is.na(better_1), "", paste0(better_1, " / ", worse_1, "; ", better_12, " / ", worse_12)),
  `Failed forecasts (% of evaluated)` = ifelse(failed == 0, "0", paste0(formatC(failed, format = "d", big.mark = ","), " (", ifelse(failed_pct < 0.005, "<0.01", sprintf("%.2f", failed_pct)), "%)")))
build <- function(k) {
  d <- disp %>% filter(class == k)
  rows <- bind_rows(lapply(unique(d$grp), function(g) {
    hdr <- as.data.frame(as.list(setNames(c(g, rep("", 6)), names(d)[-(1:2)])), check.names = FALSE)
    bind_rows(hdr, d %>% filter(grp == g) %>% select(-class, -grp))
  }))
  hi <- which(rows$Model %in% unique(d$grp))
  ft <- flextable(rows) %>% font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8, part = "all") %>%
    bold(part = "header") %>% bold(i = hi, part = "body") %>% align(j = 2:ncol(rows), align = "center", part = "all") %>%
    border_remove() %>% hline_top(border = fp_border(width = 1.5), part = "header") %>%
    hline_bottom(border = fp_border(width = 0.75), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
    padding(padding.top = 0, padding.bottom = 0, padding.left = 2, padding.right = 2, part = "all") %>%
    line_spacing(space = 1, part = "all") %>% height_all(height = 0.17, part = "body") %>% hrule(rule = "exact", part = "body")
  for (i in hi) ft <- merge_at(ft, i = i, j = seq_len(ncol(rows)), part = "body")
  w <- c(1.9, 0.8, 1.25, 1.35, 1.35, 1.25, 0.95); width(ft, width = w * 9 / sum(w))
}
caption <- function(k) paste0(
  "Results of every forecasting model and simple reference for ",
  tolower(CLASS_LAB[[k]]), ifelse(k == "diversity", " (Hill numbers q = 0, 1, 2)",
                                  " (log(CPUE + 1) of eight taxa, CPUE being catch per unit effort)"),
  ", by scale and sampling method: median over the series of each group (",
  ifelse(k == "diversity", "estuaries: 21 series per sampling method, 7 estuaries \u00d7 q = 0, 1, 2; coastwide: 3",
                           "estuaries: 56 series per sampling method, 7 estuaries \u00d7 8 taxa; coastwide: 8"),
  "), calculated on the 367 forecast dates up to December 2022. All models forecast the monthly departures from the seasonal ",
  "curve, which was refitted at each forecast date using only the data up to that date. RMSE: root mean squared out-of-sample ",
  "prediction error ", ifelse(k == "diversity", "divided by the series mean (a relative error)", "in log units"),
  ", 1, 12 and 24 months ahead. Cumulative average: the seasonal curve plus the mean departure of all earlier months; 12-month ",
  "average: the same with the last 12 months only. Skill: 1 minus the ratio of the mean squared prediction error to that ",
  "of the cumulative average (0 = no better than the cumulative average). Base model: the model without the feature being tested ",
  "(covariates or trend): SARIMA (automatic) for ARIMAX and the SARIMA and damped trend versions; the state-space model for the ",
  "state-space models with covariates or a trend. Skill against the base model measures what that feature adds (0 = no gain). ",
  "Better / worse: numbers of series whose bootstrap range of skill against the base model lies above / below zero, 1 and 12 months ",
  "ahead. Failed forecasts: cases (series \u00d7 forecast date \u00d7 months ahead) lost because a model could not be fitted at a ",
  "forecast date, with their percentage of the model's evaluated cases; these cases were left out of the comparison with the base model ",
  "for both models. Models with covariates were run to 12 months (n/a at 24).")
for (k in c("diversity", "abundance")) {
  if (file.exists(out_docx[[k]])) archive_files(out_docx[[k]])
  doc <- read_docx() %>%
    body_set_default_section(prop_section(page_size = page_size(orient = "landscape", width = 8.5, height = 11),
                                          page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75))) %>%
    body_add_fpar(fpar(ftext(paste0(TABLE_LABEL[[k]], " "), fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)), ftext(caption(k), fp_text(font.family = "Times New Roman", font.size = 10)))) %>%
    body_add_flextable(build(k))
  print(doc, target = out_docx[[k]])
}
cat("Saved:", out_csv, "\n"); for (f in out_docx) cat("Saved:", f, "\n")
print(as.data.frame(disp %>% select(-class)), row.names = FALSE)
end_log("finished")
cat("h15_all_models.R completed\n")
