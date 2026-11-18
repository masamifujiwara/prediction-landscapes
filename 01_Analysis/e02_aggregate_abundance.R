# =============================================================================
# Script:  e02_aggregate_abundance.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds the monthly abundance targets for the 8 target taxa: mean catch per
#   haul or tow (CPUE) and log(CPUE + 1), per estuary (alpha) and coastwide
#   (gamma), for both sampling methods (bag seine, bay trawl). Every sample is
#   paired with every target, so a target not caught in a sampled haul counts
#   as a catch of 0 in the mean. April and May 2020 (no sampling) are filled
#   once here and marked, so e03 and every model read one complete series.
#
# Run after: a00, b03, b09
#
# Inputs:
#   01_Analysis/Data.RData                                          (a00)
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   01_Analysis/Data/b9_forecast_settings/settings.rds              (b09)
#
# Outputs:
#   01_Analysis/Data/e2_aggregate_abundance/
#     alpha_abundance_targets.rds  BAY_ID, GEAR, YEAR, MONTH, date,
#                                  SPECIES_CODE, Common_Name, filled,
#                                  n_samples, mean_cpue, log_cpue
#     gamma_abundance_targets.rds  the same without BAY_ID, plus n_bays
#   02_Results/03_Diagnostics/e2_aggregate_abundance/
#     abundance_summary.csv        per sampling method and target: sampled
#                                  months, zero months, mean CPUE (estuary
#                                  and coastwide)
#     e2_coverage_by_series.csv    first and last month, months sampled,
#                                  filled, months before the series start
#     e2_filled_values.csv         the filled months with March and June 2020
#     e2_target_codes.csv          individuals per target and species code,
#                                  1982 to 2012 and 2013 to 2024
#   01_Analysis/Logs/e2_aggregate_abundance/run_<timestamp>.txt
#   Earlier outputs are moved to ARCHIVE_<timestamp> subfolders first.
#
# Method notes:
#   1. Estuaries come from settings.rds (b09, taken from b05); the script stops
#      if station_core has other estuaries.
#   2. Coastwide CPUE pools by sampling effort: the mean over every haul or
#      tow of the sampling method in the month. A coastwide month is kept only
#      when every estuary was sampled. Only leading months may be dropped; an
#      estuary missing inside the record stops the script. Column n_bays.
#   3. Menhaden is a genus-level target, "Menhaden" (Brevoortia spp.), ID 161:
#      codes 604 Gulf menhaden, 658 finescale menhaden and 161 Genus
#      Brevoortia, summed per haul. Since 2013 part of the menhaden catch has
#      been recorded only to genus, mostly small fish in the southern
#      estuaries; with 604 alone, those series would drop because of a change
#      in recording, not in abundance. Family Clupeidae (458) is not included.
#   4. Transform log(CPUE + 1). CPUE is the mean catch per haul or tow, many
#      estuary-months lie between 0 and 1 and some are 0, so log(CPUE) cannot
#      be used.
#   5. April and May 2020 are filled once, here, with
#      fill_gap_seasonal_anomaly() from lib/forecast_rules.R on the
#      log(CPUE + 1) scale (transform log1p, inverse expm1), with monthly means
#      from sampled months up to settings$target_fill$climatology_end. Filled
#      months are marked filled = TRUE, n_samples = 0. This fill was chosen
#      over linear interpolation and the monthly mean alone after hiding April
#      and May in other years and filling them again. A fill below 0 is set to
#      0, because CPUE cannot be negative.
#   6. Each series runs from its first to its last sampled month; any missing
#      month other than settings$gap_months stops the script. "Sampled
#      months" (filled == FALSE) include months with zero catch of a target.
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
source(file.path(PATHS$analysis_dir, "lib", "forecast_rules.R"))

# -----------------------------------------------------------------------------
# 1) Settings
# -----------------------------------------------------------------------------
script_stem <- "e2_aggregate_abundance"   # stem of the output folders

# One row per target; SPECIES_CODE is the target's ID (method note 3).
targets <- tibble::tribble(
  ~SPECIES_CODE, ~Common_Name,
  "602",  "Atlantic croaker",
  "608",  "Spot",
  "633",  "Pinfish",
  "9605", "Blue crab",
  "9618", "Brown shrimp",
  "601",  "Bay anchovy",
  "9600", "White shrimp",
  "161",  "Menhaden"
)
# Species codes counted in each target.
target_codes <- tibble::tribble(
  ~SPECIES_CODE, ~member_code,
  "602",  "602",
  "608",  "608",
  "633",  "633",
  "9605", "9605",
  "9618", "9618",
  "601",  "601",
  "9600", "9600",
  "161",  "604",    # Gulf menhaden
  "161",  "658",    # finescale menhaden
  "161",  "161"     # Genus Brevoortia
)

in_station  <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")

options(scipen = 999, warn = 1, width = 200)

# -----------------------------------------------------------------------------
# 2) Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# 3) Folders, archive, log
# -----------------------------------------------------------------------------
t_start <- Sys.time()
dirs <- script_dirs(script_stem, c("data", "logs", "diag"))
archive_previous(c(dirs$data, dirs$diag))
start_log(dirs$logs, script_stem)
cat("Project folder:", PATHS$project_dir, "\n")

out_alpha  <- file.path(dirs$data, "alpha_abundance_targets.rds")
out_gamma  <- file.path(dirs$data, "gamma_abundance_targets.rds")
out_summ   <- file.path(dirs$diag, "abundance_summary.csv")
out_cover  <- file.path(dirs$diag, "e2_coverage_by_series.csv")
out_filled <- file.path(dirs$diag, "e2_filled_values.csv")
out_codes  <- file.path(dirs$diag, "e2_target_codes.csv")

# -----------------------------------------------------------------------------
# 4) Inputs and checks
# -----------------------------------------------------------------------------
stop_if_missing(c(PATHS$clean_data_rdata, in_station, in_settings),
                c("a00_getdata.R", "b03_harmonize_station_metadata.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
load(PATHS$clean_data_rdata)
for (obj in c("bagseine_organisms", "baytrawl_organisms")) {
  if (!exists(obj)) stop("Data.RData has no object ", obj, call. = FALSE)
  expect_cols(get(obj), c("SAMPLE_ID", "SPECIES_CODE", "catch"), obj)
}
station_core <- readRDS(in_station)
for (el in c("station_bs_core", "station_bt_core")) {
  if (is.null(station_core[[el]])) stop("station_core.rds has no element ", el, call. = FALSE)
  expect_cols(station_core[[el]], c("SAMPLE_ID", "MAJOR_AREA", "YEAR_i", "MONTH_i"), el)
  bays_found <- sort(unique(as.integer(station_core[[el]]$MAJOR_AREA)))
  if (!setequal(bays_found, settings$bays)) {
    stop(el, " has bays ", paste(bays_found, collapse = ", "), "; settings$bays is ",
         paste(settings$bays, collapse = ", "), call. = FALSE)
  }
}
if (!setequal(settings$gears, c("BS", "BT"))) stop("settings$gears is not BS, BT", call. = FALSE)
if (anyDuplicated(target_codes$member_code)) stop("A species code is assigned to two targets", call. = FALSE)
if (!setequal(unique(target_codes$SPECIES_CODE), targets$SPECIES_CODE)) stop("targets and target_codes disagree", call. = FALSE)
gap_months <- as.Date(settings$gap_months)
clim_end   <- as.Date(settings$target_fill$climatology_end)
n_bays_all <- length(settings$bays)
cat("Targets:", paste(targets$Common_Name, collapse = ", "), "\n")
cat("Gap months to fill:", paste(format(gap_months, "%Y-%m"), collapse = ", "),
    "| climatology up to:", format(clim_end, "%Y-%m"), "\n\n")

# -----------------------------------------------------------------------------
# 5) Mean catch per haul or tow, per estuary and coastwide
# -----------------------------------------------------------------------------
aggregate_cpue <- function(station_df, org_df, gear_label) {
  station_base <- station_df %>%
    transmute(SAMPLE_ID, BAY_ID = as.integer(MAJOR_AREA), YEAR = as.integer(YEAR_i), MONTH = as.integer(MONTH_i))
  if (anyDuplicated(station_base$SAMPLE_ID)) stop(gear_label, ": duplicated SAMPLE_ID in station_core", call. = FALSE)

  org_targets <- org_df %>%
    transmute(SAMPLE_ID, member_code = as.character(SPECIES_CODE), catch = as.numeric(catch)) %>%
    inner_join(target_codes, by = "member_code") %>%
    group_by(SAMPLE_ID, SPECIES_CODE) %>%
    summarise(catch = sum(catch), .groups = "drop")
  if (anyNA(org_targets$catch)) stop(gear_label, ": missing catch in a target record", call. = FALSE)
  unmatched <- setdiff(org_targets$SAMPLE_ID, station_base$SAMPLE_ID)
  if (length(unmatched)) stop(gear_label, ": ", length(unmatched), " target records have no station row", call. = FALSE)

  # Every sample x every target; a target not caught counts as 0.
  catch_grid <- tidyr::expand_grid(SAMPLE_ID = station_base$SAMPLE_ID, SPECIES_CODE = targets$SPECIES_CODE) %>%
    left_join(station_base, by = "SAMPLE_ID") %>%
    left_join(org_targets, by = c("SAMPLE_ID", "SPECIES_CODE")) %>%
    mutate(catch = coalesce(catch, 0))

  alpha_df <- catch_grid %>%
    group_by(BAY_ID, YEAR, MONTH, SPECIES_CODE) %>%
    summarise(n_samples = n(), mean_cpue = mean(catch), .groups = "drop") %>%
    mutate(GEAR = gear_label)
  # Coastwide: every haul or tow counts once (method note 2).
  gamma_df <- catch_grid %>%
    group_by(YEAR, MONTH, SPECIES_CODE) %>%
    summarise(n_samples = n(), n_bays = n_distinct(BAY_ID), mean_cpue = mean(catch), .groups = "drop") %>%
    mutate(GEAR = gear_label)
  list(alpha = alpha_df, gamma = gamma_df)
}

cat("Bag seine (BS)...\n")
bs <- aggregate_cpue(station_core$station_bs_core, bagseine_organisms, "BS")
cat("Bay trawl (BT)...\n")
bt <- aggregate_cpue(station_core$station_bt_core, baytrawl_organisms, "BT")
alpha_obs <- bind_rows(bs$alpha, bt$alpha)
gamma_all <- bind_rows(bs$gamma, bt$gamma)

# -----------------------------------------------------------------------------
# 6) Coastwide months need every estuary (method note 2)
# -----------------------------------------------------------------------------
gamma_months <- gamma_all %>% distinct(GEAR, YEAR, MONTH, n_bays) %>%
  mutate(t = month_index(YEAR, MONTH)) %>% arrange(GEAR, t)
dropped <- list()
for (g in unique(gamma_months$GEAR)) {
  gm <- gamma_months %>% filter(GEAR == g)
  first_full <- min(gm$t[gm$n_bays == n_bays_all])
  inside <- gm %>% filter(t > first_full, n_bays < n_bays_all)
  if (nrow(inside) > 0) {
    stop("Coastwide ", g, ": not every bay sampled in ",
         paste(sprintf("%04d-%02d", inside$YEAR, inside$MONTH), collapse = ", "),
         " (inside the record)", call. = FALSE)
  }
  dropped[[g]] <- gm %>% filter(t < first_full)
}
dropped <- bind_rows(dropped)
if (nrow(dropped)) {
  cat("Coastwide months dropped (not every bay sampled):\n")
  print(as.data.frame(dropped %>% transmute(GEAR, month = sprintf("%04d-%02d", YEAR, MONTH), n_bays)), row.names = FALSE)
}
gamma_obs <- gamma_all %>% anti_join(dropped, by = c("GEAR", "YEAR", "MONTH"))

# -----------------------------------------------------------------------------
# 7) Fill April and May 2020 in every series (method note 5)
# -----------------------------------------------------------------------------
n_clamped <- 0L
fill_series <- function(df, key_cols) {
  t <- month_index(df$YEAR, df$MONTH)
  grid <- tibble::tibble(t = seq(min(t), max(t))) %>%
    mutate(YEAR = as.integer(t %/% 12L), MONTH = as.integer(t %% 12L + 1L))
  out <- grid %>% left_join(df, by = c("YEAR", "MONTH"))
  for (k in key_cols) out[[k]] <- df[[k]][1]
  gap <- is.na(out$mean_cpue)
  gap_dates <- month_date(out$YEAR[gap], out$MONTH[gap])
  if (!setequal(gap_dates, gap_months) || anyDuplicated(gap_dates)) {
    stop("Series ", paste(unlist(df[1, key_cols]), collapse = " "), ": missing months are ",
         paste(format(gap_dates, "%Y-%m"), collapse = ", "), "; expected exactly ",
         paste(format(gap_months, "%Y-%m"), collapse = ", "), call. = FALSE)
  }
  v <- fill_gap_seasonal_anomaly(out$mean_cpue, out$YEAR, out$MONTH, fill = gap, clim_end = clim_end,
                                 transform = log1p, inverse = expm1)
  neg <- gap & v < 0
  if (any(neg)) {
    cat("Fill below 0 set to 0:", paste(unlist(df[1, key_cols]), collapse = " "),
        paste(format(month_date(out$YEAR[neg], out$MONTH[neg]), "%Y-%m"), collapse = ", "),
        sprintf("(%.4f)", v[neg]), "\n")
    n_clamped <<- n_clamped + sum(neg)
    v[neg] <- 0
  }
  out$mean_cpue <- v
  out$n_samples[gap] <- 0L
  if ("n_bays" %in% names(out)) out$n_bays[gap] <- 0L
  out$filled <- gap
  out %>% select(-t)
}

alpha_filled <- alpha_obs %>%
  group_by(BAY_ID, GEAR, SPECIES_CODE) %>% group_split() %>%
  lapply(fill_series, key_cols = c("BAY_ID", "GEAR", "SPECIES_CODE")) %>% bind_rows()
gamma_filled <- gamma_obs %>%
  group_by(GEAR, SPECIES_CODE) %>% group_split() %>%
  lapply(fill_series, key_cols = c("GEAR", "SPECIES_CODE")) %>% bind_rows()
cat("Fills set to 0:", n_clamped, "\n")

finish <- function(df) {
  df %>%
    left_join(targets, by = "SPECIES_CODE") %>%
    mutate(date = month_date(YEAR, MONTH), log_cpue = log1p(mean_cpue))
}
alpha_abundance <- finish(alpha_filled) %>%
  select(BAY_ID, GEAR, YEAR, MONTH, date, SPECIES_CODE, Common_Name, filled, n_samples, mean_cpue, log_cpue) %>%
  arrange(BAY_ID, GEAR, SPECIES_CODE, date)
gamma_abundance <- finish(gamma_filled) %>%
  select(GEAR, YEAR, MONTH, date, SPECIES_CODE, Common_Name, filled, n_samples, n_bays, mean_cpue, log_cpue) %>%
  arrange(GEAR, SPECIES_CODE, date)

# -----------------------------------------------------------------------------
# 8) Output checks
# -----------------------------------------------------------------------------
n_alpha_series <- nrow(distinct(alpha_abundance, BAY_ID, GEAR, SPECIES_CODE))
n_gamma_series <- nrow(distinct(gamma_abundance, GEAR, SPECIES_CODE))
if (n_alpha_series != n_bays_all * 2L * nrow(targets)) stop("alpha has ", n_alpha_series, " series", call. = FALSE)
if (n_gamma_series != 2L * nrow(targets)) stop("gamma has ", n_gamma_series, " series", call. = FALSE)
for (d in list(alpha = alpha_abundance, gamma = gamma_abundance)) {
  if (anyNA(d %>% select(-any_of("BAY_ID")))) stop("NA in the output", call. = FALSE)
  if (any(d$mean_cpue < 0)) stop("Negative CPUE in the output", call. = FALSE)
  if (any(d$n_samples[!d$filled] < 1)) stop("A sampled month has no samples", call. = FALSE)
}
if (anyDuplicated(alpha_abundance[, c("BAY_ID", "GEAR", "SPECIES_CODE", "date")])) stop("Duplicated alpha keys", call. = FALSE)
if (anyDuplicated(gamma_abundance[, c("GEAR", "SPECIES_CODE", "date")])) stop("Duplicated gamma keys", call. = FALSE)
if (sum(alpha_abundance$filled) != n_alpha_series * length(gap_months) ||
    sum(gamma_abundance$filled) != n_gamma_series * length(gap_months)) stop("Wrong number of filled months", call. = FALSE)
if (any(gamma_abundance$n_bays[!gamma_abundance$filled] != n_bays_all)) stop("A coastwide month lacks a bay", call. = FALSE)
cat("\nAlpha:", nrow(alpha_abundance), "rows,", n_alpha_series, "series,", sum(alpha_abundance$filled), "filled\n")
cat("Gamma:", nrow(gamma_abundance), "rows,", n_gamma_series, "series,", sum(gamma_abundance$filled), "filled\n")

# -----------------------------------------------------------------------------
# 9) Diagnostics
# -----------------------------------------------------------------------------
summ_g <- gamma_abundance %>% filter(!filled) %>%
  group_by(GEAR, Common_Name, SPECIES_CODE) %>%
  summarise(coast_months = n(), coast_zero_months = sum(mean_cpue == 0),
            coast_mean_cpue = mean(mean_cpue), .groups = "drop")
summ_a <- alpha_abundance %>% filter(!filled) %>%
  group_by(GEAR, SPECIES_CODE) %>%
  summarise(bay_months = n(), bay_zero_months = sum(mean_cpue == 0),
            bay_pct_zero = round(100 * mean(mean_cpue == 0), 1), .groups = "drop")
diag_summary <- summ_g %>% left_join(summ_a, by = c("GEAR", "SPECIES_CODE")) %>%
  arrange(GEAR, desc(coast_mean_cpue))

grid_start <- min(month_index(alpha_abundance$YEAR, alpha_abundance$MONTH))
cover <- bind_rows(
  alpha_abundance %>% mutate(scale = "bay"),
  gamma_abundance %>% mutate(scale = "coastwide", BAY_ID = NA_integer_)
) %>%
  group_by(scale, BAY_ID, GEAR, SPECIES_CODE, Common_Name) %>%
  summarise(first = format(min(date), "%Y-%m"), last = format(max(date), "%Y-%m"),
            months_sampled = sum(!filled), months_filled = sum(filled),
            months_before_start = min(month_index(YEAR, MONTH)) - grid_start, .groups = "drop") %>%
  arrange(scale, BAY_ID, GEAR, SPECIES_CODE)

near <- c(month_date(2020, 3), gap_months, month_date(2020, 6))
filled_values <- bind_rows(
  alpha_abundance %>% mutate(scale = "bay"),
  gamma_abundance %>% mutate(scale = "coastwide", BAY_ID = NA_integer_)
) %>%
  filter(date %in% near) %>%
  transmute(scale, BAY_ID, GEAR, SPECIES_CODE, Common_Name, month = format(date, "%Y-%m"), filled,
            n_samples, mean_cpue = round(mean_cpue, 4), log_cpue = round(log_cpue, 4)) %>%
  arrange(scale, BAY_ID, GEAR, SPECIES_CODE, month)

code_rows <- function(org, st, gear) {
  st <- st %>% transmute(SAMPLE_ID, YEAR = as.integer(YEAR_i))
  org %>% transmute(SAMPLE_ID, member_code = as.character(SPECIES_CODE), catch = as.numeric(catch)) %>%
    inner_join(target_codes, by = "member_code") %>%
    inner_join(st, by = "SAMPLE_ID") %>%
    mutate(period = ifelse(YEAR <= 2012, "ind_1982_2012", "ind_2013_2024")) %>%
    group_by(SPECIES_CODE, member_code, period) %>% summarise(n = sum(catch), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0) %>%
    mutate(GEAR = gear)
}
codes_tab <- bind_rows(code_rows(bagseine_organisms, station_core$station_bs_core, "BS"),
                       code_rows(baytrawl_organisms, station_core$station_bt_core, "BT")) %>%
  left_join(targets, by = "SPECIES_CODE") %>%
  select(GEAR, target = SPECIES_CODE, Common_Name, member_code, any_of(c("ind_1982_2012", "ind_2013_2024"))) %>%
  arrange(GEAR, target, member_code)

# -----------------------------------------------------------------------------
# 10) Write
# -----------------------------------------------------------------------------
saveRDS(alpha_abundance, out_alpha)
saveRDS(gamma_abundance, out_gamma)
readr::write_csv(diag_summary, out_summ)
readr::write_csv(cover, out_cover, na = "")
readr::write_csv(filled_values, out_filled, na = "")
readr::write_csv(codes_tab, out_codes)

cat("\nCoastwide summary (sampled months):\n")
print(as.data.frame(diag_summary %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)
cat("\nMenhaden codes (individuals):\n")
print(as.data.frame(codes_tab %>% filter(target == "161")), row.names = FALSE)

cat("\nSaved:\n")
for (f in c(out_alpha, out_gamma, out_summ, out_cover, out_filled, out_codes)) cat("  ", f, "\n", sep = "")
cat("Run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
end_log("finished")
cat("e02_aggregate_abundance.R completed\n")
