# =============================================================================
# Script: tabS7_baseline_rmse_by_bay.R
#
# Appendix S23 (Table S7) -- Alpha-Scale Abundance Baseline RMSE by Estuary
#
# Output: Word document (.docx) formatted for PDF conversion
#   - Appendix S23 as section title (bold 12pt)
#   - Table S7 as table caption (12pt double-spaced)
#   - Footnote below table (10pt)
#   - Landscape orientation via body_end_section_landscape()
#   - Table body 8pt, header 9pt (consistent with Tables S3-S5)
#   - Bag Seine panel followed by Bay Trawl panel in one document
#
# Requires: dplyr, tidyr, officer, flextable
# Output: 01_Analysis/tables/Table_S7_Baseline_RMSE.docx
# =============================================================================

rm(list = ls())

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

need <- c("dplyr", "tidyr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(officer); library(flextable)
})

tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Lookups
# =============================================================================
sp_lut <- data.frame(
  SPECIES_CODE = c("601","602","604","608","633","9600","9605","9618"),
  common_name  = c("Bay anchovy","Atlantic croaker","Gulf menhaden","Spot",
                   "Pinfish","White shrimp","Blue crab","Brown shrimp"),
  Group        = c("Fish","Fish","Fish","Fish","Fish",
                   "Crustacean","Crustacean","Crustacean"),
  stringsAsFactors = FALSE)

bay_abbrev <- c(
  "2" = "Galv.",  "3" = "Mata.",  "4" = "S.Ant.",
  "5" = "Aran.",  "6" = "C.Chr.",
  "7" = "ULM",    "8" = "LLM")

# =============================================================================
# 2) Load and prepare data
# =============================================================================
abd_raw <- read.csv(file.path(PATHS$analysis_dir, "Data",
                              "e9_compare_predictability_abundance",
                              "metrics_dictionary_abundance.csv"),
                    stringsAsFactors = FALSE)

base_dat <- abd_raw %>%
  filter(model == "seasonal_naive",
         series == "alpha",
         BAY != "ALL",
         horizon == 1) %>%
  mutate(
    SPECIES_CODE = as.character(SPECIES_CODE),
    bay_abbrev   = recode(as.character(BAY), !!!bay_abbrev)
  ) %>%
  left_join(sp_lut, by = "SPECIES_CODE") %>%
  mutate(
    rmse_fmt    = sprintf("%.3f", rmse),
    common_name = factor(common_name, levels = sp_lut$common_name)
  )

# =============================================================================
# 3) Build wide table for each gear
# =============================================================================
make_wide <- function(gear_val) {
  base_dat %>%
    filter(GEAR == gear_val) %>%
    select(common_name, Group, bay_abbrev, rmse_fmt) %>%
    pivot_wider(names_from = bay_abbrev, values_from = rmse_fmt) %>%
    arrange(Group, common_name) %>%
    select(common_name, Group, all_of(unname(bay_abbrev)))
}

bs_wide <- make_wide("BS")
bt_wide <- make_wide("BT")

cat("BS dims:", nrow(bs_wide), "x", ncol(bs_wide), "\n")

# =============================================================================
# 4) Build flextable function
# =============================================================================
make_ft <- function(df, gear_label) {
  n_tot <- ncol(df)
  n_bay <- n_tot - 2

  flextable(as.data.frame(df)) %>%
    # Spanning gear header row
    add_header_row(
      values    = gear_label,
      colwidths = n_tot) %>%
    set_header_labels(
      common_name = "Species",
      Group       = "Group") %>%
    # ── Font and size ───────────────────────────────────────────────────────
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 8, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    italic(j = 1, part = "body") %>%
    # ── Alignment ───────────────────────────────────────────────────────────
    align(j = 1:2,     align = "left",   part = "all") %>%
    align(j = 3:n_tot, align = "center", part = "all") %>%
    # ── Spanning header shading ─────────────────────────────────────────────
    bg(i = 1, bg = "#2F5496", part = "header") %>%
    color(i = 1, color = "white", part = "header") %>%
    # ── Row shading: fish rows grey, crustacean rows white ──────────────────
    bg(i = which(df$Group == "Fish"),
       bg = "#F0F0F0", part = "body") %>%
    bg(i = which(df$Group == "Crustacean"),
       bg = "#FFFFFF", part = "body") %>%
    # ── Dashed separator between fish and crustacean ────────────────────────
    hline(i = sum(df$Group == "Fish"),
          border = fp_border(width = 0.75, style = "dashed"),
          part = "body") %>%
    # ── Borders ─────────────────────────────────────────────────────────────
    border_remove() %>%
    hline_top(border    = fp_border(width = 1.5), part = "header") %>%
    hline_bottom(border = fp_border(width = 1.0), part = "header") %>%
    hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
    # ── Column widths ────────────────────────────────────────────────────────
    width(j = 1,        width = 1.30) %>%
    width(j = 2,        width = 0.85) %>%
    width(j = 3:n_tot,  width = 0.75) %>%
    # ── Row height and padding ───────────────────────────────────────────────
    height_all(height = 0.20) %>%
    hrule(rule = "atleast") %>%
    padding(padding.top = 1, padding.bottom = 1,
            padding.left = 3, padding.right = 3, part = "all") %>%
    line_spacing(space = 1, part = "all")
}

ft_bs <- make_ft(bs_wide, "Bag Seine")
ft_bt <- make_ft(bt_wide, "Bay Trawl")

# =============================================================================
# 5) Caption and footnote text
# =============================================================================
caption_txt <- paste0(
  "Table S7. Seasonal naive baseline root mean square error (RMSE) for ",
  "population abundance forecasts at the alpha (estuary) scale, h = 1 month. ",
  "Values represent the RMSE of the seasonal naive model (preceding-year same ",
  "calendar month) for each species-estuary-gear combination, computed in units ",
  "of log-transformed catch per unit effort (log CPUE). Higher values indicate ",
  "greater interannual variability and a more challenging forecasting target. ",
  "Species: common name of the target species. ",
  "Group: taxonomic group (Fish or Crustacean). ",
  "Bay columns show RMSE values for each of the seven alpha-scale estuaries. ",
  "Results are shown separately for bag seine (top panel) and bay trawl ",
  "(bottom panel) sampling gears. ",
  "Fish rows are shaded grey; crustacean rows are white. ",
  "RMSE values are rounded to 3 decimal places."
)

footnote_txt <- paste0(
  "Estuary abbreviations: Galv. = Galveston Bay; Mata. = Matagorda Bay; ",
  "S.Ant. = San Antonio Bay; Aran. = Aransas Bay; C.Chr. = Corpus Christi Bay; ",
  "ULM = Upper Laguna Madre; LLM = Lower Laguna Madre. ",
  "Abbreviations: RMSE = Root Mean Square Error; ",
  "log CPUE = log-transformed catch per unit effort."
)

# =============================================================================
# 6) Assemble Word document (landscape)
# =============================================================================
landscape_section <- block_section(
  prop_section(
    page_size = page_size(orient = "landscape"),
    type = "continuous"
  )
)

doc <- read_docx() %>%
  # Appendix title
  # body_add_fpar(
  #   fpar(
  #     ftext(
  #       "Appendix S23",
  #       fp_text(
  #         font.size = 12,
  #         bold = TRUE,
  #         font.family = "Times New Roman"
  #       )
  #     )
  #   ),
  #   style = "Normal"
  # ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Caption
  body_add_fpar(
    fpar(
      ftext(
        caption_txt,
        fp_text(
          font.size = 12,
          font.family = "Times New Roman"
        )
      ),
      fp_p = fp_par(line_spacing = 2)
    ),
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Bag Seine table
  body_add_flextable(ft_bs) %>%
  body_add_par("", style = "Normal") %>%
  
  # Bay Trawl table
  body_add_flextable(ft_bt) %>%
  body_add_par("", style = "Normal") %>%
  
  # Footnote
  body_add_fpar(
    fpar(
      ftext(
        footnote_txt,
        fp_text(
          font.size = 10,
          font.family = "Times New Roman"
        )
      )
    ),
    style = "Normal"
  ) %>%
  
  # End block as one continuous landscape section
  body_end_block_section(value = landscape_section)

# =============================================================================
# 7) Save
# =============================================================================
out_docx <- file.path(tab_dir, "Table_S7_Baseline_RMSE.docx")
print(doc, target = out_docx)
cat("\u2705 Table S7 (Appendix S23) saved:", out_docx, "\n")
cat("BS rows:", nrow(bs_wide), "| BT rows:", nrow(bt_wide), "\n")
