# =============================================================================
# Script: tabS1_sampling_design.R  (revised)
#
# Appendix S17 (Table S1) -- Sampling Design + Community Diversity Summary
#
# Output: Word document (.docx) formatted for PDF conversion
#   - Appendix S17 as section title (bold 12pt)
#   - Table S1 as table caption (12pt double-spaced)
#   - Footnote below table (10pt)
#   - Landscape orientation via body_end_section_landscape()
#   - Table body 8pt, header 9pt (consistent with Tables S3-S5)
#
# Requires: dplyr, tidyr, officer, flextable
# Output: 01_Analysis/tables/Table_S1.docx
# =============================================================================

rm(list = ls())
set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# 0) Helpers
# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
# 1) Packages
# ─────────────────────────────────────────────────────────────────────────────
need <- c("dplyr", "tidyr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(officer); library(flextable)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Bay labels
# ─────────────────────────────────────────────────────────────────────────────
bay_labels <- c(
  "1" = "Sabine Lake",
  "2" = "Galveston Bay",
  "3" = "Matagorda Bay",
  "4" = "San Antonio Bay",
  "5" = "Aransas Bay",
  "6" = "Corpus Christi Bay",
  "7" = "Upper Laguna Madre",
  "8" = "Lower Laguna Madre"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Load data
# ─────────────────────────────────────────────────────────────────────────────
station  <- readRDS(file.path(PATHS$analysis_dir, "Data",
                              "b3_harmonize_station_metadata", "station_core.rds"))
alpha_raw <- readRDS(file.path(PATHS$analysis_dir, "Data",
                               "b6_compute_hill_alpha_q012", "alpha_q012.rds"))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Sampling design summary
# ─────────────────────────────────────────────────────────────────────────────
summarise_gear <- function(df, gear_name) {
  df %>%
    mutate(YM = paste0(YEAR_i, "-", sprintf("%02d", MONTH_i))) %>%
    group_by(MAJOR_AREA) %>%
    summarise(
      Gear           = gear_name,
      Year_start     = min(YEAR_i,  na.rm = TRUE),
      Year_end       = max(YEAR_i,  na.rm = TRUE),
      N_months       = n_distinct(YM),
      N_samples      = n(),
      Mean_stn_month = round(n() / n_distinct(YM), 1),
      .groups = "drop"
    ) %>%
    mutate(
      Years      = paste0(Year_start, "\u2013", Year_end),
      MAJOR_AREA = as.character(MAJOR_AREA)
    )
}

bs_sum <- summarise_gear(station$station_bs_core, "Bag Seine")
bt_sum <- summarise_gear(station$station_bt_core, "Bay Trawl")

design_sum <- bind_rows(bs_sum, bt_sum) %>%
  filter(MAJOR_AREA %in% names(bay_labels)) %>%
  filter(MAJOR_AREA != "1") %>%
  mutate(
    Estuary    = recode(MAJOR_AREA, !!!bay_labels),
    Estuary    = factor(Estuary, levels = unname(bay_labels)[-1]),
    GEAR_clean = ifelse(Gear == "Bag Seine", "BS", "BT")
  ) %>%
  arrange(Estuary, Gear) %>%
  mutate(
    `Total station visits`    = formatC(N_samples,      format = "d", big.mark = ","),
    `Mean stations per month` = formatC(Mean_stn_month, format = "f", digits  = 1)
  ) %>%
  select(MAJOR_AREA, GEAR_clean, Estuary, Gear, Years,
         `Sampling months`        = N_months,
         `Total station visits`,
         `Mean stations per month`)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Prepare alpha diversity data
# ─────────────────────────────────────────────────────────────────────────────
alpha <- alpha_raw %>%
  filter(
    !is.na(q0), !is.na(q1), !is.na(q2),
    N_total > 0,
    coverage_fail == FALSE
  ) %>%
  mutate(
    BAY_ID = as.character(BAY_ID),
    GEAR   = as.character(GEAR)
  ) %>%
  filter(BAY_ID %in% names(bay_labels), BAY_ID != "1")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Median + IQR per bay x gear x q
# ─────────────────────────────────────────────────────────────────────────────
div_summary <- alpha %>%
  group_by(BAY_ID, GEAR) %>%
  summarise(
    med_q0 = median(q0, na.rm = TRUE),
    iqr_q0 = IQR(q0,    na.rm = TRUE),
    med_q1 = median(q1, na.rm = TRUE),
    iqr_q1 = IQR(q1,    na.rm = TRUE),
    med_q2 = median(q2, na.rm = TRUE),
    iqr_q2 = IQR(q2,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Median q=0 (IQR)` = sprintf("%.2f (%.2f)", med_q0, iqr_q0),
    `Median q=1 (IQR)` = sprintf("%.2f (%.2f)", med_q1, iqr_q1),
    `Median q=2 (IQR)` = sprintf("%.2f (%.2f)", med_q2, iqr_q2)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7) Cohen's d (bay vs. coastwide pooled) + bootstrap 95% CI
# ─────────────────────────────────────────────────────────────────────────────
B <- 2000

cohens_d_boot <- function(x_focal, x_ref, B = 2000) {
  x_focal <- x_focal[!is.na(x_focal)]
  x_ref   <- x_ref[!is.na(x_ref)]
  n1 <- length(x_focal)
  n2 <- length(x_ref)
  if (n1 < 3 || n2 < 3) return(list(d = NA_real_, lo = NA_real_, hi = NA_real_))

  pooled_sd <- function(a, b) {
    s1 <- var(a); s2 <- var(b)
    na <- length(a); nb <- length(b)
    denom <- na + nb - 2
    if (denom <= 0) return(NA_real_)
    sqrt(((na - 1) * s1 + (nb - 1) * s2) / denom)
  }

  pd <- pooled_sd(x_focal, x_ref)
  if (is.na(pd) || pd < 1e-10) return(list(d = 0, lo = 0, hi = 0))

  d_pt <- (mean(x_focal) - mean(x_ref)) / pd

  boot_d <- numeric(B)
  for (i in seq_len(B)) {
    fb   <- sample(x_focal, n1, replace = TRUE)
    rb   <- sample(x_ref,   n2, replace = TRUE)
    pd_b <- pooled_sd(fb, rb)
    boot_d[i] <- if (!is.na(pd_b) && pd_b > 1e-10) {
      (mean(fb) - mean(rb)) / pd_b
    } else {
      NA_real_
    }
  }
  boot_d <- boot_d[!is.na(boot_d)]
  if (length(boot_d) < 10) return(list(d = round(d_pt, 2), lo = NA_real_, hi = NA_real_))

  ci <- unname(quantile(boot_d, c(0.025, 0.975)))
  list(d = round(d_pt, 2), lo = round(ci[1], 2), hi = round(ci[2], 2))
}

bays_use  <- unique(alpha$BAY_ID)
gears_use <- unique(alpha$GEAR)
q_cols    <- c("q0", "q1", "q2")

effect_rows <- list()
for (g in gears_use) {
  alpha_g <- alpha %>% filter(GEAR == g)
  for (b in bays_use) {
    focal <- alpha_g %>% filter(BAY_ID == b)
    ref   <- alpha_g %>% filter(BAY_ID != b)
    if (nrow(focal) < 3) next
    row <- tibble(BAY_ID = b, GEAR = g)
    for (q in q_cols) {
      res <- cohens_d_boot(focal[[q]], ref[[q]], B = B)
      row[[paste0("d_",  q)]] <- res$d
      row[[paste0("lo_", q)]] <- res$lo
      row[[paste0("hi_", q)]] <- res$hi
    }
    effect_rows[[length(effect_rows) + 1]] <- row
  }
}

effect_tbl <- bind_rows(effect_rows) %>%
  mutate(
    `Cohen's d q=0 [95% CI]` = sprintf("% .2f [%.2f, %.2f]", d_q0, lo_q0, hi_q0),
    `Cohen's d q=1 [95% CI]` = sprintf("% .2f [%.2f, %.2f]", d_q1, lo_q1, hi_q1),
    `Cohen's d q=2 [95% CI]` = sprintf("% .2f [%.2f, %.2f]", d_q2, lo_q2, hi_q2)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 8) Join everything together
# ─────────────────────────────────────────────────────────────────────────────
all_tab <- design_sum %>%
  left_join(
    div_summary %>% select(BAY_ID, GEAR, starts_with("Median")),
    by = c("MAJOR_AREA" = "BAY_ID", "GEAR_clean" = "GEAR")
  ) %>%
  left_join(
    effect_tbl %>% select(BAY_ID, GEAR, starts_with("Cohen")),
    by = c("MAJOR_AREA" = "BAY_ID", "GEAR_clean" = "GEAR")
  ) %>%
  select(
    Estuary, Gear, Years,
    `Sampling months`,
    `Total station visits`,
    `Mean stations per month`,
    `Median q=0 (IQR)`,
    `Median q=1 (IQR)`,
    `Median q=2 (IQR)`,
    `Cohen's d q=0 [95% CI]`,
    `Cohen's d q=1 [95% CI]`,
    `Cohen's d q=2 [95% CI]`
  )

cat("\n--- Table preview ---\n")
print(as.data.frame(all_tab), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 9) Build flextable
# ─────────────────────────────────────────────────────────────────────────────
col_widths <- c(1.30, 0.80, 0.70, 0.65, 0.80, 0.80,
                0.90, 0.90, 0.90, 1.15, 1.15, 1.15)
n_col <- ncol(all_tab)

ft <- flextable(as.data.frame(all_tab)) %>%
  # ── Spanning header groups ────────────────────────────────────────────────
  add_header_row(
    values    = c("", "", "", "Sampling design", "", "",
                  "Median diversity (IQR)", "", "",
                  "Cohen's d vs. coastwide pool [95% CI]", "", ""),
    colwidths = rep(1, n_col)
  ) %>%
  # ── Column widths ─────────────────────────────────────────────────────────
  width(j = seq_len(n_col), width = col_widths) %>%
  # ── Font and size ─────────────────────────────────────────────────────────
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 8, part = "body") %>%
  fontsize(size = 9, part = "header") %>%
  bold(part = "header") %>%
  # ── Alignment ─────────────────────────────────────────────────────────────
  align(j = 1:2,     align = "left",   part = "all") %>%
  align(j = 3:n_col, align = "center", part = "all") %>%
  # ── Merge spanning header cells ───────────────────────────────────────────
  merge_at(i = 1, j = 1:3,   part = "header") %>%
  merge_at(i = 1, j = 4:6,   part = "header") %>%
  merge_at(i = 1, j = 7:9,   part = "header") %>%
  merge_at(i = 1, j = 10:12, part = "header") %>%
  # ── Merge Estuary column in body ──────────────────────────────────────────
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top", part = "body") %>%
  # ── Borders ───────────────────────────────────────────────────────────────
  border_remove() %>%
  hline_top(border    = fp_border(width = 1.5), part = "header") %>%
  hline(i = 1, border = fp_border(width = 0.5, style = "dotted"), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.0), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  vline(j = 3,  border = fp_border(width = 0.5, style = "dotted"), part = "all") %>%
  vline(j = 6,  border = fp_border(width = 0.75), part = "all") %>%
  vline(j = 9,  border = fp_border(width = 0.75), part = "all") %>%
  # ── Shading: alternate by estuary (rows 1-2 shaded, 3-4 plain, etc.) ──────
  bg(i = seq(1, nrow(all_tab), by = 4), bg = "#F0F4F8", part = "body") %>%
  bg(i = seq(2, nrow(all_tab), by = 4), bg = "#F0F4F8", part = "body") %>%
  # ── Row height and padding ────────────────────────────────────────────────
  height_all(height = 0.20) %>%
  hrule(rule = "atleast") %>%
  padding(padding.top = 1, padding.bottom = 1,
          padding.left = 3, padding.right = 3, part = "all") %>%
  line_spacing(space = 1, part = "all")

# ─────────────────────────────────────────────────────────────────────────────
# 10) Caption and footnote text
# ─────────────────────────────────────────────────────────────────────────────
caption_txt <- paste0(
  "Table S1. Sampling design and community diversity summary for the Texas ",
  "Parks and Wildlife Department Coastal Fisheries Marine Resource Monitoring ",
  "Program. Sampling months is the number of calendar months with at least one ",
  "station visit. Total station visits is the cumulative number of individual ",
  "station sampling events across all years. Mean stations per month is the ",
  "average number of stations sampled per monthly survey. Median diversity ",
  "values (with interquartile range, IQR) are coverage-standardized Hill ",
  "numbers computed per monthly time cell (BAY \u00d7 GEAR \u00d7 YEAR \u00d7 MONTH); ",
  "cells with zero catch or failed coverage standardization are excluded. ",
  "q = 0: species richness; q = 1: Shannon diversity (effective species); ",
  "q = 2: dominance-weighted diversity (inverse Simpson effective species). ",
  "Cohen's d compares each estuary's distribution to the pooled distribution ",
  "of all other estuaries within the same gear type; 95% confidence intervals ",
  "were obtained by bootstrap resampling (B = 2,000). Positive d indicates ",
  "higher diversity than the coastwide pool."
)

footnote_txt <- paste0(
  "Abbreviations: BS = Bag Seine; BT = Bay Trawl; IQR = interquartile range. ",
  "Notes: Sabine Lake (MAJOR_AREA = 1) was excluded from the forecasting ",
  "analysis and is not shown."
)

# ─────────────────────────────────────────────────────────────────────────────
# 11) Assemble Word document (landscape)
# ─────────────────────────────────────────────────────────────────────────────
landscape_section <- block_section(
  prop_section(
    page_size = page_size(orient = "landscape"),
    type = "continuous"
  )
)

doc <- read_docx() %>%
  # body_add_fpar(
  #   fpar(
  #     ftext(
  #       "Appendix S17",
  #       fp_text(font.size = 12, bold = TRUE, font.family = "Times New Roman")
  #     )
  #   ),
  #   style = "Normal"
  # ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(
    fpar(
      ftext(
        caption_txt,
        fp_text(font.size = 12, font.family = "Times New Roman")
      ),
      fp_p = fp_par(line_spacing = 2)
    ),
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(
    fpar(
      ftext(
        footnote_txt,
        fp_text(font.size = 10, font.family = "Times New Roman")
      )
    ),
    style = "Normal"
  ) %>%
  body_end_block_section(value = landscape_section)

out_docx <- file.path(tab_dir, "Table_S1.docx")
print(doc, target = out_docx)
cat("\u2705 Table S1 (Appendix S17) saved:", out_docx, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 12) Diagnostic CSV (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
diag_out <- file.path(PATHS$results_dir, "03_Diagnostics", "tabS1_sampling_design")
dir.create(diag_out, recursive = TRUE, showWarnings = FALSE)

div_effect_export <- div_summary %>%
  left_join(
    effect_tbl %>% select(BAY_ID, GEAR,
                          d_q0, lo_q0, hi_q0,
                          d_q1, lo_q1, hi_q1,
                          d_q2, lo_q2, hi_q2),
    by = c("BAY_ID", "GEAR")
  ) %>%
  left_join(
    tibble(BAY_ID = names(bay_labels), Estuary = unname(bay_labels)),
    by = "BAY_ID"
  ) %>%
  arrange(BAY_ID, GEAR)

readr::write_csv(div_effect_export,
                 file.path(diag_out, "diversity_median_cohens_d_by_bay_gear.csv"))
cat("\u2705 Diagnostic CSV saved.\n")
cat("\nOK: tabS1_sampling_design.R completed successfully.\n")
