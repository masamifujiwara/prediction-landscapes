# =============================================================================
# Script: tab1_bs_bt_combined.R
#
# Table 1 — SARIMA Forecasting Skill: Bag Seine AND Bay Trawl
#
# Layout (revised):
#   Single table, portrait orientation, 5 columns:
#     Target | Scale | h=1 [95% CI] | h=12 [95% CI] | Ret.
#
#   Section 1: Bag Seine
#     Community diversity: q=0 Alpha, q=0 Gamma, q=1 Alpha, ... q=2 Gamma
#     Population abundance: 8 species (Gamma), ordered by BS h=1 skill desc
#   Section 2: Bay Trawl
#     Community diversity: q=0 Alpha, q=0 Gamma, ... q=2 Gamma
#     Population abundance: 8 species (Gamma), ordered by BT h=1 skill desc
#
# Formatting (ESA Ecological Applications):
#   - Portrait orientation
#   - Table body: 10pt Times New Roman, single-spaced (minimum allowed)
#   - Table caption: 12pt Times New Roman, double-spaced
#   - No landscape; table may extend vertically across pages
#
# ─── RNG fidelity ────────────────────────────────────────────────────────────
# BS: set.seed(42)  — matches original Table 1 exactly
# BT: set.seed(123) — independent stream
#
# Outputs:
#   tables/Table_1_BS_BT_combined_CI.csv   — machine-readable CI lookup
#   tables/Table_1_BS_BT_combined.docx     — formatted Word table
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 0) Locate project root and load PATHS
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
script_dir  <- get_script_dir()
setwd(script_dir)

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

# ─────────────────────────────────────────────────────────────────────────────
# 1) Packages
# ─────────────────────────────────────────────────────────────────────────────
need <- c("dplyr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(officer); library(flextable)
})

# ─────────────────────────────────────────────────────────────────────────────
# 2) Output directory
# ─────────────────────────────────────────────────────────────────────────────
tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Load fold-level data
# ─────────────────────────────────────────────────────────────────────────────
div_folds <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "d1_compare_model_metrics_dictionary_v3", "fc_joined_diversity.rds"))

abd_folds <- readRDS(file.path(PATHS$analysis_dir, "Data",
  "e9_compare_predictability_abundance", "pred_joined_abundance.rds"))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Bootstrap CI function  (identical to tab1_sarima_skill_summary.R §4)
# ─────────────────────────────────────────────────────────────────────────────
boot_skill_ci <- function(error_adv, error_base, B = 2000, alpha = 0.05) {
  ok     <- !is.na(error_adv) & !is.na(error_base)
  e_adv  <- error_adv[ok]
  e_base <- error_base[ok]
  n      <- length(e_adv)
  if (n < 5) return(c(skill = NA, lo = NA, hi = NA))

  rmse_adv  <- sqrt(mean(e_adv^2))
  rmse_base <- sqrt(mean(e_base^2))
  skill_pt  <- 100 * (1 - rmse_adv / rmse_base)

  boot_skills <- vapply(seq_len(B), function(i) {
    idx <- sample(n, n, replace = TRUE)
    ra  <- sqrt(mean(e_adv[idx]^2))
    rb  <- sqrt(mean(e_base[idx]^2))
    100 * (1 - ra / rb)
  }, numeric(1))

  c(skill = skill_pt,
    lo    = unname(quantile(boot_skills, alpha / 2)),
    hi    = unname(quantile(boot_skills, 1 - alpha / 2)))
}

# ─────────────────────────────────────────────────────────────────────────────
# 5) Lookups
# ─────────────────────────────────────────────────────────────────────────────
sp_lut <- data.frame(
  SPECIES_CODE = c("601", "602", "604", "608", "633", "9600", "9605", "9618"),
  common_name  = c("Bay anchovy", "Atlantic croaker", "Gulf menhaden", "Spot",
                   "Pinfish", "White shrimp", "Blue crab", "Brown shrimp"),
  stringsAsFactors = FALSE
)

bays <- as.character(c(2, 3, 4, 5, 6, 7, 8))

# ─────────────────────────────────────────────────────────────────────────────
# 6) BAG SEINE — set.seed(42), exact original order
#    Diversity first, then abundance — same call sequence as original script
# ─────────────────────────────────────────────────────────────────────────────
cat("=== BAG SEINE (seed=42) ===\n")
set.seed(42)

# ── 6a) Diversity BS ──────────────────────────────────────────────────────────
cat("Computing BS diversity CIs...\n")

div_bs_data <- div_folds %>%
  filter(model_clean == "sarima", GEAR == "BS", !is.na(error_base))

div_rows_bs <- list()

for (q_val in c(0, 1, 2)) {

  # Alpha: per-bay point estimates -> median; CI from pooled folds
  for (h_val in c(1, 12)) {
    d_alpha <- div_bs_data %>%
      filter(q == q_val, BAY %in% bays, horizon == h_val)

    bay_skills <- sapply(bays, function(b) {
      d <- div_bs_data %>% filter(q == q_val, BAY == b, horizon == h_val)
      if (nrow(d) == 0) return(NA)
      100 * (1 - sqrt(mean(d$error^2, na.rm = TRUE)) /
                 sqrt(mean(d$error_base^2, na.rm = TRUE)))
    })
    skill_med <- median(bay_skills, na.rm = TRUE)

    ci <- boot_skill_ci(d_alpha$error, d_alpha$error_base)  # RNG draw

    key <- paste0("q", q_val, "_alpha_h", h_val)
    div_rows_bs[[key]] <- data.frame(
      gear = "BS", target = paste0("q = ", q_val), scale = "Alpha",
      q = q_val, h = h_val,
      skill = skill_med, lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BS] q=%d Alpha h=%2d  skill=%6.2f  [%6.2f, %6.2f]\n",
                q_val, h_val, skill_med, ci["lo"], ci["hi"]))
  }

  # Gamma
  for (h_val in c(1, 12)) {
    d_gamma <- div_bs_data %>%
      filter(q == q_val, BAY == "ALL", horizon == h_val)

    ci <- boot_skill_ci(d_gamma$error, d_gamma$error_base)  # RNG draw

    key <- paste0("q", q_val, "_gamma_h", h_val)
    div_rows_bs[[key]] <- data.frame(
      gear = "BS", target = paste0("q = ", q_val), scale = "Gamma",
      q = q_val, h = h_val,
      skill = ci["skill"], lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BS] q=%d Gamma h=%2d  skill=%6.2f  [%6.2f, %6.2f]\n",
                q_val, h_val, ci["skill"], ci["lo"], ci["hi"]))
  }
}

div_long_bs <- bind_rows(div_rows_bs)

# ── 6b) Abundance BS ──────────────────────────────────────────────────────────
# RNG state here is identical to tab1_sarima_skill_summary.R after its
# diversity loops, so abundance BS CIs will match Table 1 exactly.
cat("\nComputing BS abundance CIs...\n")

abd_bs_data <- abd_folds %>%
  filter(model == "sarima", GEAR == "BS", series == "gamma",
         !is.na(error_base)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE")

abd_rows_bs <- list()
for (sp in sp_lut$common_name) {
  for (h_val in c(1, 12)) {
    d  <- abd_bs_data %>% filter(common_name == sp, horizon == h_val)
    ci <- boot_skill_ci(d$error, d$error_base)
    key <- paste0(sp, "_h", h_val)
    abd_rows_bs[[key]] <- data.frame(
      gear = "BS", target = sp, scale = "Gamma", q = NA, h = h_val,
      n_folds = sum(!is.na(d$error) & !is.na(d$error_base)),
      skill = ci["skill"], lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BS] %-20s h=%2d  n=%3d  skill=%6.2f  [%6.2f, %6.2f]\n",
                sp, h_val, abd_rows_bs[[key]]$n_folds,
                ci["skill"], ci["lo"], ci["hi"]))
  }
}
abd_long_bs <- bind_rows(abd_rows_bs)

# ─────────────────────────────────────────────────────────────────────────────
# 7) BAY TRAWL — set.seed(123), fully independent stream
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== BAY TRAWL (seed=123) ===\n")
set.seed(123)

# ── 7a) Diversity BT ──────────────────────────────────────────────────────────
cat("Computing BT diversity CIs...\n")

div_bt_data <- div_folds %>%
  filter(model_clean == "sarima", GEAR == "BT", !is.na(error_base))

div_rows_bt <- list()

for (q_val in c(0, 1, 2)) {

  for (h_val in c(1, 12)) {
    d_alpha <- div_bt_data %>%
      filter(q == q_val, BAY %in% bays, horizon == h_val)

    bay_skills <- sapply(bays, function(b) {
      d <- div_bt_data %>% filter(q == q_val, BAY == b, horizon == h_val)
      if (nrow(d) == 0) return(NA)
      100 * (1 - sqrt(mean(d$error^2, na.rm = TRUE)) /
                 sqrt(mean(d$error_base^2, na.rm = TRUE)))
    })
    skill_med <- median(bay_skills, na.rm = TRUE)

    ci <- boot_skill_ci(d_alpha$error, d_alpha$error_base)

    key <- paste0("q", q_val, "_alpha_h", h_val)
    div_rows_bt[[key]] <- data.frame(
      gear = "BT", target = paste0("q = ", q_val), scale = "Alpha",
      q = q_val, h = h_val,
      skill = skill_med, lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BT] q=%d Alpha h=%2d  skill=%6.2f  [%6.2f, %6.2f]\n",
                q_val, h_val, skill_med, ci["lo"], ci["hi"]))
  }

  for (h_val in c(1, 12)) {
    d_gamma <- div_bt_data %>%
      filter(q == q_val, BAY == "ALL", horizon == h_val)

    ci <- boot_skill_ci(d_gamma$error, d_gamma$error_base)

    key <- paste0("q", q_val, "_gamma_h", h_val)
    div_rows_bt[[key]] <- data.frame(
      gear = "BT", target = paste0("q = ", q_val), scale = "Gamma",
      q = q_val, h = h_val,
      skill = ci["skill"], lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BT] q=%d Gamma h=%2d  skill=%6.2f  [%6.2f, %6.2f]\n",
                q_val, h_val, ci["skill"], ci["lo"], ci["hi"]))
  }
}

div_long_bt <- bind_rows(div_rows_bt)

# ── 7b) Abundance BT ──────────────────────────────────────────────────────────
cat("\nComputing BT abundance CIs...\n")

abd_bt_data <- abd_folds %>%
  filter(model == "sarima", GEAR == "BT", series == "gamma",
         !is.na(error_base)) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(sp_lut, by = "SPECIES_CODE")

missing_sp <- setdiff(sp_lut$common_name, unique(abd_bt_data$common_name))
if (length(missing_sp) > 0)
  warning("Species missing from BT gamma SARIMA data: ",
          paste(missing_sp, collapse = ", "))

abd_rows_bt <- list()
for (sp in sp_lut$common_name) {
  for (h_val in c(1, 12)) {
    d  <- abd_bt_data %>% filter(common_name == sp, horizon == h_val)
    ci <- boot_skill_ci(d$error, d$error_base)
    key <- paste0(sp, "_h", h_val)
    abd_rows_bt[[key]] <- data.frame(
      gear = "BT", target = sp, scale = "Gamma", q = NA, h = h_val,
      n_folds = sum(!is.na(d$error) & !is.na(d$error_base)),
      skill = ci["skill"], lo = ci["lo"], hi = ci["hi"],
      stringsAsFactors = FALSE
    )
    cat(sprintf("  [BT] %-20s h=%2d  n=%3d  skill=%6.2f  [%6.2f, %6.2f]\n",
                sp, h_val, abd_rows_bt[[key]]$n_folds,
                ci["skill"], ci["lo"], ci["hi"]))
  }
}
abd_long_bt <- bind_rows(abd_rows_bt)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Pivot wide per gear
# ─────────────────────────────────────────────────────────────────────────────
pivot_wide_gear <- function(df) {
  has_n <- "n_folds" %in% names(df)
  h1 <- df %>% filter(h == 1) %>%
    select(gear, target, scale, q,
           any_of("n_folds"),
           skill1 = skill, lo1 = lo, hi1 = hi)
  h12 <- df %>% filter(h == 12) %>%
    select(gear, target, scale, q,
           any_of("n_folds"),
           skill12 = skill, lo12 = lo, hi12 = hi)
  if (has_n) {
    h1  <- h1  %>% rename(n_h1  = n_folds)
    h12 <- h12 %>% rename(n_h12 = n_folds)
  }
  left_join(h1, h12, by = c("gear", "target", "scale", "q")) %>%
    mutate(retention = round(100 * skill12 / skill1, 1))
}

wide_div_bs <- pivot_wide_gear(div_long_bs)
wide_div_bt <- pivot_wide_gear(div_long_bt)
wide_abd_bs <- pivot_wide_gear(abd_long_bs)
wide_abd_bt <- pivot_wide_gear(abd_long_bt)

# Order diversity rows: q ascending, Alpha before Gamma
wide_div_bs <- wide_div_bs %>%
  mutate(scale_ord = ifelse(scale == "Alpha", 1, 2)) %>%
  arrange(q, scale_ord) %>%
  select(-scale_ord)

wide_div_bt <- wide_div_bt %>%
  mutate(scale_ord = ifelse(scale == "Alpha", 1, 2)) %>%
  arrange(q, scale_ord) %>%
  select(-scale_ord)

# Order abundance rows: by h=1 skill descending, separately per gear
wide_abd_bs <- wide_abd_bs %>% arrange(desc(skill1))
wide_abd_bt <- wide_abd_bt %>% arrange(desc(skill1))

# ─────────────────────────────────────────────────────────────────────────────
# 9) Verification: BS values vs Table_1_CI_lookup.csv
# ─────────────────────────────────────────────────────────────────────────────
# Join BS and BT side by side for verification only
join_gears <- function(wide_bs, wide_bt, join_cols) {
  wide_bs %>%
    rename_with(~ paste0("bs_", .), .cols = -all_of(join_cols)) %>%
    left_join(
      wide_bt %>% rename_with(~ paste0("bt_", .), .cols = -all_of(join_cols)),
      by = join_cols
    )
}

div_combined <- join_gears(wide_div_bs, wide_div_bt, c("target", "scale", "q"))
abd_combined <- join_gears(wide_abd_bs, wide_abd_bt, c("target", "scale", "q"))

ref_csv <- file.path(tab_dir, "Table_1_CI_lookup.csv")
if (file.exists(ref_csv)) {
  cat("\n--- Verifying BS values against Table_1_CI_lookup.csv ---\n")
  ref <- read.csv(ref_csv, stringsAsFactors = FALSE) %>%
    select(type, target, scale, ref_skill_h1 = skill_h1,
           ref_lo_h1 = ci_lo_h1, ref_hi_h1 = ci_hi_h1,
           ref_skill_h12 = skill_h12,
           ref_lo_h12 = ci_lo_h12, ref_hi_h12 = ci_hi_h12)

  check_div <- div_combined %>%
    mutate(type = "Diversity") %>%
    select(type, target, scale,
           bs_skill1, bs_lo1, bs_hi1, bs_skill12, bs_lo12, bs_hi12) %>%
    left_join(ref, by = c("type", "target", "scale"))

  check_abd <- abd_combined %>%
    mutate(type = "Abundance") %>%
    select(type, target, scale,
           bs_skill1, bs_lo1, bs_hi1, bs_skill12, bs_lo12, bs_hi12) %>%
    left_join(ref, by = c("type", "target", "scale"))

  check_all <- bind_rows(check_div, check_abd) %>%
    mutate(
      d_skill_h1  = round(bs_skill1  - ref_skill_h1,  4),
      d_lo_h1     = round(bs_lo1     - ref_lo_h1,     4),
      d_hi_h1     = round(bs_hi1     - ref_hi_h1,     4),
      d_skill_h12 = round(bs_skill12 - ref_skill_h12, 4),
      d_lo_h12    = round(bs_lo12    - ref_lo_h12,    4),
      d_hi_h12    = round(bs_hi12    - ref_hi_h12,    4),
      ok = d_skill_h1 == 0 & d_lo_h1 == 0 & d_hi_h1 == 0 &
           d_skill_h12 == 0 & d_lo_h12 == 0 & d_hi_h12 == 0
    )

  if (all(check_all$ok, na.rm = TRUE)) {
    cat("  \u2705 All BS values match Table 1 exactly.\n")
  } else {
    cat("  \u26a0 MISMATCH detected:\n")
    print(check_all %>% filter(!ok) %>%
            select(type, target, scale, starts_with("d_")),
          row.names = FALSE)
  }
} else {
  cat("\n  (Table_1_CI_lookup.csv not found -- skipping verification.)\n")
  cat("  Run tab1_sarima_skill_summary.R first to generate the reference.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 10) Export CI lookup CSV  (unchanged from original)
# ─────────────────────────────────────────────────────────────────────────────
ci_lookup_div <- div_combined %>%
  mutate(type = "Diversity") %>%
  select(type, target, scale, q,
         bs_skill_h1 = bs_skill1, bs_ci_lo_h1 = bs_lo1, bs_ci_hi_h1 = bs_hi1,
         bs_skill_h12 = bs_skill12, bs_ci_lo_h12 = bs_lo12, bs_ci_hi_h12 = bs_hi12,
         bs_retention_pct = bs_retention,
         bt_skill_h1 = bt_skill1, bt_ci_lo_h1 = bt_lo1, bt_ci_hi_h1 = bt_hi1,
         bt_skill_h12 = bt_skill12, bt_ci_lo_h12 = bt_lo12, bt_ci_hi_h12 = bt_hi12,
         bt_retention_pct = bt_retention)

ci_lookup_abd <- abd_combined %>%
  mutate(type = "Abundance") %>%
  select(type, target, scale, q,
         bs_n_h1 = bs_n_h1, bs_skill_h1 = bs_skill1,
         bs_ci_lo_h1 = bs_lo1, bs_ci_hi_h1 = bs_hi1,
         bs_n_h12 = bs_n_h12, bs_skill_h12 = bs_skill12,
         bs_ci_lo_h12 = bs_lo12, bs_ci_hi_h12 = bs_hi12,
         bs_retention_pct = bs_retention,
         bt_n_h1 = bt_n_h1, bt_skill_h1 = bt_skill1,
         bt_ci_lo_h1 = bt_lo1, bt_ci_hi_h1 = bt_hi1,
         bt_n_h12 = bt_n_h12, bt_skill_h12 = bt_skill12,
         bt_ci_lo_h12 = bt_lo12, bt_ci_hi_h12 = bt_hi12,
         bt_retention_pct = bt_retention)

ci_lookup_all <- bind_rows(ci_lookup_div, ci_lookup_abd) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

out_csv <- file.path(tab_dir, "Table_1_BS_BT_combined_CI.csv")
write.csv(ci_lookup_all, out_csv, row.names = FALSE)
cat("\n\u2705 CI lookup CSV saved:", out_csv, "\n")

cat("\n--- Full combined skill summary ---\n")
print(as.data.frame(ci_lookup_all), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 11) Build formatted Word table
#
# Layout: 5 columns, portrait, BS stacked above BT
#   Col 1: Target
#   Col 2: Scale
#   Col 3: h = 1 [95% CI]
#   Col 4: h = 12 [95% CI]
#   Col 5: Ret.
#
# Structure:
#   [Bag Seine]          <- gear header, spans all 5 cols
#     Target | Scale | h=1 | h=12 | Ret.   <- column headers
#     [Community diversity]                <- subheader
#       q=0 Alpha ... q=2 Gamma
#     [Population abundance]               <- subheader
#       8 species (ordered by h=1 skill desc)
#   [Bay Trawl]          <- gear header, spans all 5 cols
#     Target | Scale | h=1 | h=12 | Ret.   <- column headers
#     [Community diversity]
#       q=0 Alpha ... q=2 Gamma
#     [Population abundance]
#       8 species (ordered by h=1 skill desc)
#
# Formatting (ESA Ecological Applications):
#   - Body: 10pt Times New Roman, single-spaced
#   - Caption: 12pt Times New Roman, double-spaced
#   - Portrait; table may extend vertically
# ─────────────────────────────────────────────────────────────────────────────

fmt_cell <- function(skill, lo, hi) {
  ifelse(is.na(skill), "\u2014",
         sprintf("%.1f [%.1f, %.1f]", round(skill, 1), round(lo, 1), round(hi, 1)))
}
fmt_ret <- function(r) ifelse(is.na(r), "\u2014", sprintf("%.1f%%", r))

# Build a 5-column display data frame for one gear block
build_gear_block <- function(wide_div, wide_abd) {
  div_disp <- wide_div %>% transmute(
    Target = target,
    Scale  = scale,
    `h = 1 [95% CI]`  = fmt_cell(skill1,  lo1,  hi1),
    `h = 12 [95% CI]` = fmt_cell(skill12, lo12, hi12),
    `Ret.`             = fmt_ret(retention)
  )

  abd_disp <- wide_abd %>% transmute(
    Target = target,
    Scale  = scale,
    `h = 1 [95% CI]`  = fmt_cell(skill1,  lo1,  hi1),
    `h = 12 [95% CI]` = fmt_cell(skill12, lo12, hi12),
    `Ret.`             = fmt_ret(retention)
  )

  list(div = div_disp, abd = abd_disp)
}

bs_block <- build_gear_block(wide_div_bs, wide_abd_bs)
bt_block <- build_gear_block(wide_div_bt, wide_abd_bt)

# Row-index bookkeeping (needed for flextable formatting)
# Structure per gear block:
#   1 row  : gear header      (spans 5 cols, no data col)
#   1 row  : column headers
#   1 row  : "Community diversity" subheader
#   6 rows : diversity data
#   1 row  : "Population abundance" subheader
#   8 rows : abundance data
# Total per block = 18 rows; grand total = 36 rows

# We build a flat data frame with a special "row_type" attribute column
# that we use for formatting, then drop before passing to flextable.

make_block_df <- function(gear_label, block) {
  col_names <- c("Target", "Scale", "h = 1 [95% CI]", "h = 12 [95% CI]", "Ret.")

  gear_hdr  <- setNames(as.data.frame(matrix(
    c(gear_label, "", "", "", ""), nrow = 1)), col_names)
  col_hdr   <- setNames(as.data.frame(matrix(
    col_names, nrow = 1)), col_names)
  div_hdr   <- setNames(as.data.frame(matrix(
    c("Community diversity", "", "", "", ""), nrow = 1)), col_names)
  abd_hdr   <- setNames(as.data.frame(matrix(
    c("Population abundance", "", "", "", ""), nrow = 1)), col_names)

  df <- bind_rows(gear_hdr, col_hdr, div_hdr, block$div, abd_hdr, block$abd)
  df$.row_type <- c(
    "gear_hdr",
    "col_hdr",
    "sub_hdr",
    rep("div_data", nrow(block$div)),
    "sub_hdr",
    rep("abd_data", nrow(block$abd))
  )
  df
}

bs_df <- make_block_df("Bag Seine", bs_block)
bt_df <- make_block_df("Bay Trawl", bt_block)
full_df <- bind_rows(bs_df, bt_df)

# Extract row-type index before dropping helper column
row_type  <- full_df$.row_type
full_df   <- full_df %>% select(-.row_type)
n_rows    <- nrow(full_df)

# Row indices by type
idx_gear_hdr  <- which(row_type == "gear_hdr")   # 2 rows
idx_col_hdr   <- which(row_type == "col_hdr")    # 2 rows
idx_sub_hdr   <- which(row_type == "sub_hdr")    # 4 rows
idx_div_data  <- which(row_type == "div_data")   # 12 rows
idx_abd_data  <- which(row_type == "abd_data")   # 16 rows
last_row      <- n_rows

# Column widths (inches): must fit within 6.5" text width (8.5" - 2x1" margins)
# Target=1.60, Scale=0.65, h1=1.55, h12=1.55, Ret.=0.65 -> total=6.00"
col_w <- c(Target = 1.60, Scale = 0.65,
           `h = 1 [95% CI]` = 1.55, `h = 12 [95% CI]` = 1.55, `Ret.` = 0.65)

ft <- flextable(full_df) %>%
  # ── Fonts and sizes ──────────────────────────────────────────────────────
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  # ── Column widths ─────────────────────────────────────────────────────────
  width(j = 1, width = col_w["Target"]) %>%
  width(j = 2, width = col_w["Scale"]) %>%
  width(j = 3, width = col_w["h = 1 [95% CI]"]) %>%
  width(j = 4, width = col_w["h = 12 [95% CI]"]) %>%
  width(j = 5, width = col_w["Ret."]) %>%
  # ── Remove default header (we baked headers into body rows) ───────────────
  delete_part(part = "header") %>%
  # ── Alignment ─────────────────────────────────────────────────────────────
  align(j = 1:2, align = "left",   part = "body") %>%
  align(j = 3:5, align = "center", part = "body") %>%
  # ── Bold rows ─────────────────────────────────────────────────────────────
  bold(i = idx_gear_hdr, part = "body") %>%
  bold(i = idx_col_hdr,  part = "body") %>%
  bold(i = idx_sub_hdr,  part = "body") %>%
  # ── Borders: remove all, then add selectively ─────────────────────────────
  border_remove() %>%
  # Top rule above first gear header
  hline_top(border = fp_border(width = 1.5), part = "body") %>%
  # Bottom rule of each column header row
  hline(i = idx_col_hdr,
        border = fp_border(width = 0.75), part = "body") %>%
  # Medium rule above each gear header (except the very first)
  hline(i = idx_gear_hdr[-1],
        border = fp_border(width = 1.0), part = "body") %>%
  # Light rule above each subheader
  hline(i = idx_sub_hdr,
        border = fp_border(width = 0.5, style = "dashed"), part = "body") %>%
  # Bottom rule at end of table
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  # ── Merge gear header and subheader rows across all 5 columns ─────────────
  merge_h(i = idx_gear_hdr, part = "body") %>%
  merge_h(i = idx_sub_hdr,  part = "body") %>%
  # ── Shading: diversity data rows light grey ────────────────────────────────
  bg(i = idx_div_data, bg = "#F5F5F5", part = "body") %>%
  # ── Row heights ───────────────────────────────────────────────────────────
  height_all(height = 0.22) %>%
  height(i = idx_gear_hdr, height = 0.26) %>%
  hrule(rule = "atleast") %>%
  # ── Padding ───────────────────────────────────────────────────────────────
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  # ── Line spacing: single for body ─────────────────────────────────────────
  line_spacing(space = 1, part = "body")

# ── Caption (12pt, double-spaced per ESA guidelines) ─────────────────────────
title_text <- paste0(
  "Table 1. SARIMA forecasting skill for community diversity and population ",
  "abundance targets. Values are percentage improvement in RMSE over the ",
  "seasonal naive baseline with bootstrap 95% confidence intervals (B = 2,000 ",
  "resamples, percentile method, fold-level resampling). Diversity rows show ",
  "the median point estimate across seven estuaries (Alpha) and coastwide ",
  "results (Gamma). Abundance rows show gamma-scale (coastwide) results only, ",
  "ordered by h = 1 skill (descending). ",
  "Ret. = skill retention (h = 12 skill as a percentage of h = 1 skill)."
)

doc <- read_docx() %>%
  body_add_fpar(
    fpar(
      ftext(title_text,
            fp_text(font.size = 12, font.family = "Times New Roman")),
      fp_p = fp_par(line_spacing = 2)
    ),
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft)

out_docx <- file.path(tab_dir, "Table_1_BS_BT_combined.docx")
print(doc, target = out_docx)
cat("\u2705 Word table saved:", out_docx, "\n")
