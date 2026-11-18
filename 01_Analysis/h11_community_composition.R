# =============================================================================
# Script:  h11_community_composition.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S3 (Online Resource 2): composition of the sampled community
#   for each sampling method, in three panels:
#     (a) taxonomic composition: number of taxa and percentage of catch;
#     (b) the ten most abundant taxa;
#     (c) the three most abundant taxa in each estuary.
#   It also writes the panels and the full taxon list as CSV files and prints
#   the figures used in the Methods paragraph on the community.
#
# Run after: a00, b04, b09
#
# Inputs:
#   01_Analysis/Data.RData (a00: species, removed_taxa)
#   01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds, bt_events.rds (b04)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09: estuaries)
#
# Outputs:
#   01_Analysis/tables/Table_S3_community_composition.docx
#   01_Analysis/Data/h11_community_composition/panel_a_groups.csv,
#     panel_b_top_taxa.csv, panel_c_top_by_estuary.csv, taxa_full_list.csv
#   01_Analysis/Logs/h11_community_composition/run_<timestamp>.txt
#   Earlier CSV files and the earlier Word file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   Taxa are reported as recorded by the survey, without merges; records
#   identified only to genus or higher are their own rows. The Menidia merge
#   made in a00 stands. Menhaden appears as recorded (three codes); the
#   abundance target "Menhaden" (e02) combines them, which is stated in
#   Online Resource 1, not here.
#   The data are those used in the analysis: the seven estuaries in the
#   settings (Sabine Lake is dropped in a00), the taxa kept by a00 and the
#   catch rows kept by b04 (rows with zero catch are treated as entry errors).
#   The codes removed by a00 are reported in a note, grouped by the reason
#   a00 records in removed_taxa.
#   The lookup's group code 0 (identified above species level) spans fishes
#   and invertebrates and gets its own row; totals are checked so that no
#   taxon leaves the table. Codes still missing from the species list get no
#   group row but a note; they stay in "All taxa recorded" and in the
#   diversity estimates as separate taxa.
#   There is no seasonal breakdown, because every series is deseasonalized
#   before any model is fitted.
#   This script tabulates catch records only; it fits no model.
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

### -----------------------------------------------------------------------------
### 1) Packages
### -----------------------------------------------------------------------------
need <- c("dplyr", "tidyr", "readr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(officer); library(flextable)
})

### -----------------------------------------------------------------------------
### 2) Settings
### -----------------------------------------------------------------------------
N_TOP_OVERALL  <- 10L   # panel (b)
N_TOP_BY_BAY   <- 3L    # panel (c)

TABLE_LABEL    <- "Table S3"   # Online Resource 2 number
IN_SETTINGS <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
stop_if_missing(IN_SETTINGS, "b09_forecast_settings.R")
settings <- load_settings(IN_SETTINGS)
### The seven estuaries of the analysis (Sabine Lake is dropped in a00). Every
### panel covers the same estuaries; the coastwide series pool exactly these.
AREAS_SURVEY <- as.character(settings$bays)
AREAS_ALPHA  <- AREAS_SURVEY
GEARS <- c(BS = "Bag seine", BT = "Bay trawl")

AREA_NAMES <- c(
  "2" = "Galveston Bay",
  "3" = "Matagorda Bay",     "4" = "San Antonio Bay",
  "5" = "Aransas Bay",       "6" = "Corpus Christi Bay",
  "7" = "Upper Laguna Madre","8" = "Lower Laguna Madre")

### The lookup's group codes. 0 is not in the TPWD list; it marks records
### not identified to species (see the Method notes in the header).
GROUP_LABELS <- c(
  "0"  = "Not identified to species",
  "1"  = "Fishes",
  "2"  = "Birds",
  "3"  = "Mammals",
  "4"  = "Amphibians",
  "5"  = "Reptiles",
  "6"  = "Invertebrates",
  "10" = "Vegetation")

### Groups collapsed into one display row, because each is a rounding error on
### its own and a table of three zero rows tells a reader nothing. Full detail
### stays in panel_a_groups.csv.
GROUP_INCIDENTAL <- c("Birds", "Mammals", "Amphibians", "Reptiles", "Vegetation")
LABEL_INCIDENTAL <- "Other vertebrates and vegetation (incidental)"
LABEL_NO_LOOKUP  <- "No entry in the species reference"

### Display order for panel (a).
GROUP_ORDER <- c("Fishes", "Invertebrates", "Not identified to species",
                 LABEL_INCIDENTAL)   # unlisted codes: a note, not a row

### -----------------------------------------------------------------------------
### 3) Directories, archiving, logging
### -----------------------------------------------------------------------------
script_stem <- "h11_community_composition"
dirs <- script_dirs(script_stem, c("data", "logs"))
dir_data_script <- dirs$data
tab_dir <- file.path(PATHS$analysis_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
old <- list.files(dir_data_script, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
say <- function(...) cat(paste0(...), "\n", sep = "")   # the run log captures the console
say("Taxa reported as identified by the survey; no consolidation.")
say("Estuaries: ", paste(AREAS_SURVEY, collapse = ", "))
say("")
### -----------------------------------------------------------------------------
### 4) Inputs
### -----------------------------------------------------------------------------
IN_RDATA <- file.path(PATHS$analysis_dir, "Data.RData")
IN_EVENT <- file.path(PATHS$analysis_dir, "Data", "b4_join_organisms_to_station",
                      c(BS = "bs_events.rds", BT = "bt_events.rds"))
names(IN_EVENT) <- c("BS", "BT")
stopifnot(file.exists(IN_RDATA), all(file.exists(IN_EVENT)))

env_rdata <- new.env()
load(IN_RDATA, envir = env_rdata)
stopifnot(exists("species", envir = env_rdata))
spp_raw <- as.data.frame(env_rdata$species, stringsAsFactors = FALSE)

### The lookup's column names carry the group legend inside the header string, so
### they are matched by position after asserting the shape rather than by name.
stopifnot(ncol(spp_raw) >= 4L)
spp <- data.frame(
  code   = as.character(spp_raw[[1]]),
  latin  = as.character(spp_raw[[2]]),
  common = as.character(spp_raw[[3]]),
  taxa   = as.character(spp_raw[[4]]),
  stringsAsFactors = FALSE)
spp$code <- trimws(spp$code)

### Codes missing from the species list in Data.RData, as identified by
### TPWD. Group codes as in
### the list: 0 not identified to species, 1 fishes, 6 invertebrates.
SPP_SUPPLEMENT <- data.frame(
  code   = c("128", "437", "9273", "9292", "9357"),
  latin  = c("Opsopoeodus emiliae", "Family Hemiramphidae", "Cystobranchus vividus",
             "Angulus texanus", "Armases cinereum"),
  common = c("Pugnose minnow", "Family halfbeaks", "(Marine leech)",
             "Texas tellin", "Squareback marsh crab"),
  taxa   = c("1", "0", "6", "6", "6"),
  stringsAsFactors = FALSE)
spp <- rbind(spp, SPP_SUPPLEMENT[!SPP_SUPPLEMENT$code %in% spp$code, ])

if (anyDuplicated(spp$code)) {
  dups <- unique(spp$code[duplicated(spp$code)])
  stop("Duplicate species codes in the reference table: ",
       paste(head(dups, 10), collapse = ", "),
       ". A join on these would multiply catch rows.")
}

unknown_grp <- setdiff(unique(spp$taxa[!is.na(spp$taxa)]), names(GROUP_LABELS))
if (length(unknown_grp)) {
  stop("Group code(s) absent from GROUP_LABELS: ", paste(unknown_grp, collapse = ", "),
       ". Add them rather than letting the taxa silently leave the table.")
}

say("Species reference: ", nrow(spp), " entries")
say("Group codes present: ",
    paste(sort(unique(spp$taxa)), collapse = ", "),
    "   (0 added by the user for records not identified to species)")
say("")

### -----------------------------------------------------------------------------
### 5) Catch by taxon and gear
### -----------------------------------------------------------------------------
### The species reference wraps many common names in parentheses, and it is its own
### styling rather than information: the convention covers genus- and family-level
### entries such as "(Grass shrimp - unidentified)", Genus Palaemonetes, but equally
### species-level invertebrates with proper binomials, "(Sea squirt)" for Molgula
### manhattensis, "(Striped sea star)" for Luidia clathrata. So the parentheses carry
### no signal a reader of this table could use, and they are stripped for display.
###
### The " - unidentified" suffix does carry signal, since it says the survey could not
### resolve the record below genus, so it is kept and re-set as a parenthetical
### qualifier where it reads as one. The unmodified reference value stays in the CSV
### exports under `label`, so the source is always recoverable.
tidy_label <- function(x) {
  y <- trimws(x)
  outer <- grepl("^\\(.*\\)$", y)
  y[outer] <- sub("^\\((.*)\\)$", "\\1", y[outer])
  y <- sub("\\s*-\\s*unidentified$", " (unidentified)", y)
  trimws(y)
}

label_group <- function(taxa_code) {
  out <- unname(GROUP_LABELS[taxa_code])
  out[is.na(taxa_code)] <- LABEL_NO_LOOKUP
  out[is.na(out)]       <- LABEL_NO_LOOKUP
  out
}

read_gear <- function(gear) {
  d <- readRDS(IN_EVENT[[gear]])
  stopifnot(all(c("SPECIES_CODE", "catch", "MAJOR_AREA") %in% names(d)))
  d <- d %>%
    mutate(SPECIES_CODE = trimws(as.character(SPECIES_CODE)),
           MAJOR_AREA   = trimws(as.character(MAJOR_AREA)),
           catch        = as.numeric(catch))
  if (!setequal(unique(d$MAJOR_AREA), AREAS_SURVEY)) {
    stop("Estuaries in ", gear, " (", paste(sort(unique(d$MAJOR_AREA)), collapse = ", "),
         ") differ from the settings (", paste(AREAS_SURVEY, collapse = ", "), ")", call. = FALSE)
  }
  d %>% filter(!is.na(catch), catch > 0)
}

events <- lapply(names(GEARS), read_gear)
names(events) <- names(GEARS)

### One row per taxon per gear, with its share of that gear's total catch.
by_taxon <- bind_rows(lapply(names(events), function(g) {
  events[[g]] %>%
    group_by(SPECIES_CODE) %>%
    summarise(catch = sum(catch), .groups = "drop") %>%
    mutate(gear = g)
})) %>%
  left_join(spp, by = c("SPECIES_CODE" = "code")) %>%
  mutate(group = label_group(taxa),
         label = ifelse(is.na(common) | !nzchar(common),
                        paste0("Unreferenced code ", SPECIES_CODE), common),
         label_display = tidy_label(label),
         latin = ifelse(is.na(latin), NA_character_, latin)) %>%
  group_by(gear) %>%
  mutate(pct = 100 * catch / sum(catch)) %>%
  ungroup()

for (g in names(GEARS)) {
  sub <- by_taxon %>% filter(gear == g)
  say(GEARS[[g]], ": ", nrow(sub), " taxa with non-zero catch, total catch ",
      format(round(sum(sub$catch)), big.mark = ","))
  stopifnot(abs(sum(sub$pct) - 100) < 1e-6)
}
say("")

### -----------------------------------------------------------------------------
### 6) Panel (a): taxonomic composition by gear
### -----------------------------------------------------------------------------
panel_a_full <- by_taxon %>%
  group_by(gear, group) %>%
  summarise(n_taxa = n(), catch = sum(catch), .groups = "drop") %>%
  group_by(gear) %>%
  mutate(pct = 100 * catch / sum(catch)) %>%
  ungroup()

### Assert nothing left the table between section 5 and here.
chk <- panel_a_full %>% group_by(gear) %>%
  summarise(n = sum(n_taxa), p = sum(pct), .groups = "drop")
ref <- by_taxon %>% group_by(gear) %>% summarise(n = n(), .groups = "drop")
stopifnot(identical(sort(chk$gear), sort(ref$gear)))
stopifnot(all(chk$n[order(chk$gear)] == ref$n[order(ref$gear)]))
stopifnot(all(abs(chk$p - 100) < 1e-6))

panel_a <- panel_a_full %>%
  mutate(display = ifelse(group %in% GROUP_INCIDENTAL, LABEL_INCIDENTAL, group)) %>%
  group_by(gear, display) %>%
  summarise(n_taxa = sum(n_taxa), pct = sum(pct), .groups = "drop") %>%
  filter(display != LABEL_NO_LOOKUP) %>%
  mutate(display = factor(display, levels = GROUP_ORDER)) %>%
  arrange(gear, display)

### Codes still missing from the species list: reported in the note.
unlisted <- by_taxon %>% filter(group == LABEL_NO_LOOKUP)
n_unl_code <- n_distinct(unlisted$SPECIES_CODE); n_unl_ind <- round(sum(unlisted$catch))
n_unl_gear <- function(g) sum(unlisted$gear == g)
n_unl_both <- sum(table(unlisted$SPECIES_CODE) > 1)
num_word <- function(n) if (n >= 1 && n <= 10) c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")[n] else as.character(n)
cap1 <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))
say("Codes not in the species list: ", n_unl_code, " (", n_unl_ind, " individuals): ",
    paste(sort(unique(unlisted$SPECIES_CODE)), collapse = ", "))
if (anyNA(panel_a$display)) stop("A group in panel (a) is missing from GROUP_ORDER", call. = FALSE)

### Values are formatted as strings so one column can carry both the main groups
### and the incidental ones. The incidental groups are small fractions of one per
### cent, where three decimals would print "0.000" and read as a measured zero
### rather than as a trace, so they get a threshold label instead.
fmt_pct <- function(p) {
  ifelse(p >= 0.1,  sprintf("%.1f", p),
  ifelse(p >= 0.01, sprintf("%.2f", p),
  ifelse(p > 0,     "<0.01", "0")))
}

panel_a_wide <- panel_a %>%
  mutate(cell = paste0(format(n_taxa, big.mark = ","), " (", fmt_pct(pct), ")")) %>%
  select(display, gear, cell) %>%
  pivot_wider(names_from = gear, values_from = cell) %>%
  arrange(display) %>%
  mutate(display = as.character(display))

totals <- by_taxon %>% group_by(gear) %>%
  summarise(cell = paste0(format(n(), big.mark = ","), " (100.0)"), .groups = "drop") %>%
  pivot_wider(names_from = gear, values_from = cell) %>%
  mutate(display = "All taxa recorded")
panel_a_wide <- bind_rows(panel_a_wide, totals)

names(panel_a_wide)[1] <- "Taxonomic group"
for (g in names(GEARS)) {
  j <- which(names(panel_a_wide) == g)
  if (length(j)) names(panel_a_wide)[j] <- paste0(GEARS[[g]], ": taxa (% of catch)")
}

say("Panel (a), taxa (% of catch):")
say(paste(capture.output(print(as.data.frame(panel_a_wide), row.names = FALSE)),
          collapse = "\n"))
say("")

### -----------------------------------------------------------------------------
### 7) Panel (b): the most abundant taxa by gear
### -----------------------------------------------------------------------------
top_overall <- by_taxon %>%
  group_by(gear) %>%
  slice_max(order_by = catch, n = N_TOP_OVERALL, with_ties = FALSE) %>%
  arrange(desc(catch), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(entry = ifelse(is.na(latin) | !nzchar(latin),
                        sprintf("%s (%.1f)", label_display, pct),
                        sprintf("%s, %s (%.1f)", label_display, latin, pct)))

panel_b <- top_overall %>%
  select(rank, gear, entry) %>%
  pivot_wider(names_from = gear, values_from = entry) %>%
  arrange(rank)
names(panel_b)[1] <- "Rank"
for (g in names(GEARS)) {
  j <- which(names(panel_b) == g)
  if (length(j)) names(panel_b)[j] <- paste0(GEARS[[g]], ": taxon (% of catch)")
}

say("Panel (b), ", N_TOP_OVERALL, " most abundant taxa per sampling method:")
say(paste(capture.output(print(as.data.frame(panel_b), row.names = FALSE)),
          collapse = "\n"))
say("")

### -----------------------------------------------------------------------------
### 8) Panel (c): the most abundant taxa by estuary and gear
###
### Shares here are within estuary and gear, so each row's percentages are of that
### estuary's catch under that gear, not of the survey total.
### -----------------------------------------------------------------------------
by_area <- bind_rows(lapply(names(events), function(g) {
  events[[g]] %>%
    filter(MAJOR_AREA %in% AREAS_ALPHA) %>%
    group_by(MAJOR_AREA, SPECIES_CODE) %>%
    summarise(catch = sum(catch), .groups = "drop") %>%
    mutate(gear = g)
})) %>%
  left_join(spp, by = c("SPECIES_CODE" = "code")) %>%
  mutate(label = ifelse(is.na(common) | !nzchar(common),
                        paste0("Unreferenced code ", SPECIES_CODE), common),
         label_display = tidy_label(label)) %>%
  group_by(gear, MAJOR_AREA) %>%
  mutate(pct = 100 * catch / sum(catch)) %>%
  ungroup()

panel_c_long <- by_area %>%
  group_by(gear, MAJOR_AREA) %>%
  slice_max(order_by = catch, n = N_TOP_BY_BAY, with_ties = FALSE) %>%
  arrange(desc(catch), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(entry = sprintf("%s (%.0f)", label_display, pct))

### Every estuary must have produced the full set, otherwise a blank cell would
### pass silently into the document.
expect_c <- length(AREAS_ALPHA) * length(GEARS) * N_TOP_BY_BAY
if (nrow(panel_c_long) != expect_c) {
  stop("Panel (c) has ", nrow(panel_c_long), " rows where ", expect_c,
       " were expected. An estuary-gear combination is short of ", N_TOP_BY_BAY,
       " taxa.")
}

rank_cols  <- paste0("r", seq_len(N_TOP_BY_BAY))
rank_names <- paste0(c("First", "Second", "Third", "Fourth",
                       "Fifth")[seq_len(N_TOP_BY_BAY)], " (% of catch)")

panel_c <- panel_c_long %>%
  mutate(Estuary = unname(AREA_NAMES[MAJOR_AREA]),
         Gear    = unname(GEARS[gear])) %>%
  select(Estuary, Gear, rank, entry) %>%
  pivot_wider(names_from = rank, values_from = entry, names_prefix = "r") %>%
  ### Named, not positional: pivot_wider orders new columns by the values it
  ### encounters, so a change upstream could silently transpose the ranks.
  select(Estuary, Gear, all_of(rank_cols)) %>%
  arrange(match(Estuary, unname(AREA_NAMES[AREAS_ALPHA])), Gear)
names(panel_c)[3:(2 + N_TOP_BY_BAY)] <- rank_names
stopifnot(!any(is.na(panel_c)))

say("Panel (c), ", N_TOP_BY_BAY, " most abundant taxa per estuary and sampling method:")
say(paste(capture.output(print(as.data.frame(panel_c), row.names = FALSE)),
          collapse = "\n"))
say("")

### -----------------------------------------------------------------------------
### 8b) Taxa removed upstream (a00), for the note
### -----------------------------------------------------------------------------
removed <- as.data.frame(env_rdata$removed_taxa, stringsAsFactors = FALSE)
if (!all(c("SPECIES_CODE", "gear", "reason") %in% names(removed))) stop("`removed_taxa` lacks SPECIES_CODE, sampling method or reason", call. = FALSE)
# Kinds follow the reasons a00 writes into removed_taxa; a
# reason that matches none of them stops the script, so a new reason in a00
# cannot be reported under the wrong kind.
removed <- removed %>% mutate(kind = dplyr::case_when(
  grepl("presence only", reason)             ~ "presence-only colonial or gelatinous",
  grepl("^sea squirt", reason)               ~ "sea squirt",
  grepl("^colonial cnidarian", reason)       ~ "colonial cnidarian",
  grepl("not in the species list", reason)   ~ "not in the species list",
  grepl("not a target group", reason)        ~ "non-target group",
  TRUE                                       ~ NA_character_))
if (anyNA(removed$kind)) stop("Unknown removal reason in removed_taxa: ",
                              paste(unique(removed$reason[is.na(removed$kind)]), collapse = "; "), call. = FALSE)
removed_by_gear <- removed %>%
  group_by(gear, kind) %>% summarise(codes = n_distinct(SPECIES_CODE), .groups = "drop")
say("Codes removed upstream by a00, by sampling method:")
say(paste(capture.output(print(as.data.frame(removed_by_gear), row.names = FALSE)), collapse = "\n"))
say("")
n_rm <- function(g, k) { v <- removed_by_gear$codes[removed_by_gear$gear == g & removed_by_gear$kind == k]; if (length(v)) v else 0L }
# One clause per kind of removal, with the codes per gear; kinds with no code are left out.
rm_clause <- function(k, what) {
  if (n_rm("BS", k) + n_rm("BT", k) == 0) return(NULL)
  paste0(n_rm("BS", k), " (bag seine) and ", n_rm("BT", k), " (bay trawl) codes of ", what)
}
rm_clauses <- c(
  rm_clause("presence-only colonial or gelatinous",
            "colonial, gelatinous or sessile invertebrates recorded mostly as present without a count"),
  rm_clause("colonial cnidarian", "colonial cnidarians, whose colonies are not counted as individuals"),
  rm_clause("sea squirt", "sea squirts, which may not have been recorded in the early years of the survey"),
  rm_clause("not in the species list", "taxa not in the TPWD species list and not identified"),
  rm_clause("non-target group", "reptiles, amphibians and birds"))
rm_sentence <- if (length(rm_clauses) == 1) rm_clauses else
  paste0(paste(rm_clauses[-length(rm_clauses)], collapse = "; "), "; and ", rm_clauses[length(rm_clauses)])

### -----------------------------------------------------------------------------
### 9) Write the machine-readable outputs
### -----------------------------------------------------------------------------
out_a <- file.path(dir_data_script, "panel_a_groups.csv")
out_b <- file.path(dir_data_script, "panel_b_top_taxa.csv")
out_c <- file.path(dir_data_script, "panel_c_top_by_estuary.csv")
out_f <- file.path(dir_data_script, "taxa_full_list.csv")

write_csv(panel_a_full %>%
            mutate(gear_label = unname(GEARS[gear])) %>%
            arrange(gear, desc(catch)), out_a)
write_csv(top_overall %>%
            mutate(gear_label = unname(GEARS[gear])) %>%
            select(gear, gear_label, rank, SPECIES_CODE, label, label_display,
                   latin, group, catch, pct), out_b)
write_csv(panel_c_long %>%
            mutate(gear_label = unname(GEARS[gear]),
                   estuary    = unname(AREA_NAMES[MAJOR_AREA])) %>%
            select(gear, gear_label, MAJOR_AREA, estuary, rank, SPECIES_CODE,
                   label, label_display, catch, pct), out_c)
write_csv(by_taxon %>%
            mutate(gear_label = unname(GEARS[gear])) %>%
            arrange(gear, desc(catch)) %>%
            select(gear, gear_label, SPECIES_CODE, label, label_display, latin,
                   group, catch, pct), out_f)

say("Wrote:")
for (f in c(out_a, out_b, out_c, out_f)) say("  ", basename(f))
say("")

### -----------------------------------------------------------------------------
### 10) Build the Word table
### -----------------------------------------------------------------------------
TNR <- "Times New Roman"

mk_ft <- function(df, widths) {
  flextable(as.data.frame(df) |> setNames(sub("^Gear$", "Sampling method", names(df)))) %>%
    width(j = seq_len(ncol(df)), width = widths) %>%
    font(fontname = TNR, part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = seq(2, ncol(df)), align = "left", part = "all") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1.5), part = "header") %>%
    hline_bottom(border = fp_border(width = 1), part = "header") %>%
    hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
    padding(padding.top = 2, padding.bottom = 2, part = "all")
}

ft_a <- mk_ft(panel_a_wide, c(2.60, 1.85, 1.85))
ft_b <- mk_ft(panel_b,      c(0.55, 3.00, 3.00))

### Scientific names in italics: genus and species names
### italic; the rank word "Genus" and names of family or higher rank upright.
latin_parts <- function(x) {
  if (is.na(x) || !nzchar(x)) return(list(text = character(0), ital = logical(0)))
  m <- regmatches(x, regexec("^(Genus|Subgenus) (.+)$", x))[[1]]
  if (length(m)) return(list(text = c(paste0(m[2], " "), m[3]), ital = c(FALSE, TRUE)))
  if (grepl("^(Family|Subfamily|Superfamily|Order|Suborder|Infraorder|Class|Subclass|Phylum|Subphylum)\\b", x))
    return(list(text = x, ital = FALSE))
  list(text = x, ital = TRUE)
}
cell_par <- function(label, latin, pct) {
  lp <- latin_parts(latin)
  chunks <- c(list(as_chunk(if (length(lp$text)) paste0(label, ", ") else label)),
              Map(function(t, i) if (i) as_i(t) else as_chunk(t), lp$text, lp$ital),
              list(as_chunk(sprintf(" (%.1f)", pct))))
  do.call(as_paragraph, unname(chunks))
}
for (g in names(GEARS)) {
  j <- which(names(panel_b) == paste0(GEARS[[g]], ": taxon (% of catch)"))
  d <- top_overall %>% filter(gear == g) %>% arrange(rank)
  for (i in seq_len(nrow(d))) {
    if (panel_b$Rank[i] != d$rank[i]) stop("Panel (b) rows out of order", call. = FALSE)
    ft_b <- compose(ft_b, i = i, j = j, part = "body", value = cell_par(d$label_display[i], d$latin[i], d$pct[i]))
  }
}
ft_b <- ft_b %>% font(fontname = TNR, part = "all") %>% fontsize(size = 9, part = "all")
ft_c <- mk_ft(panel_c,      c(1.55, 1.05, 1.60, 1.60, 1.60)) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top", part = "body")

cap_main <- paste0(
  "Composition of the sampled community. Catch records from the Texas ",
  "Parks and Wildlife Department Coastal Fisheries Marine Resource Monitoring ",
  "Program, 1982 to 2024, in the seven estuaries analyzed, summarized separately for the two sampling methods. ",
  "Percentages are of the total catch for that sampling method. Taxa are reported as recorded by the survey, ",
  "except that silverside records were combined at genus level (")
cap_main_italic <- "Menidia"; cap_main_end <- ")."

cap_a <- paste0(
  "(a) Taxonomic composition. Counts are taxa with non-zero catch; percentages in ",
  "parentheses are percentages of the total catch.")
cap_b <- paste0(
  "(b) The ", N_TOP_OVERALL, " most abundant taxa for each sampling method, with the percentage of ",
  "total catch in parentheses.")
cap_c <- paste0(
  "(c) The ", N_TOP_BY_BAY, " most abundant taxa in each estuary, with the percentage of ",
  "that estuary's catch under that sampling method in parentheses.")

foot <- paste0(
  "Notes. \"Not identified to species\" holds mostly records identified only to genus or a higher rank, of ",
  "fishes and invertebrates alike, so it is reported separately rather than assigned to either. Removed before analysis ",
  "and not in this table: ", rm_sentence, ". ",
  if (n_unl_code > 0) paste0(cap1(num_word(n_unl_code)), " codes (", num_word(n_unl_gear("BS")), " recorded by the bag seine and ",
    num_word(n_unl_gear("BT")), " by the bay trawl", if (n_unl_both > 0) paste0(", ", num_word(n_unl_both), " by both") else "",
    "; ", n_unl_ind, " individuals in all) are not in the TPWD species list available to us; they are ",
    "counted in \"All taxa recorded\" and kept as separate taxa in the biodiversity estimates but not assigned to a group. ") else "",
  "All panels cover the seven estuaries analyzed (Sabine Lake is not included), which the coastwide series pool.")

par_txt <- function(txt, size = 10, bold = FALSE, spacing = 1) {
  fpar(ftext(txt, fp_text(font.size = size, bold = bold, font.family = TNR)),
       fp_p = fp_par(line_spacing = spacing))
}

### `body_end_block_section()` wants a block_section, not the prop_section that
### describes the page. `prop_section()` alone fails with
### "inherits(value, \"block_section\") is not TRUE".
landscape_section <- block_section(
  prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1),
    type = "continuous"
  )
)

doc <- read_docx() %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 12, bold = TRUE, font.family = TNR)),
                     ftext(cap_main, fp_text(font.size = 12, font.family = TNR)),
                     ftext(cap_main_italic, fp_text(font.size = 12, italic = TRUE, font.family = TNR)),
                     ftext(cap_main_end, fp_text(font.size = 12, font.family = TNR)),
                     fp_p = fp_par(line_spacing = 2)), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(par_txt(cap_a, bold = TRUE), style = "Normal") %>%
  body_add_flextable(ft_a) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(par_txt(cap_b, bold = TRUE), style = "Normal") %>%
  body_add_flextable(ft_b) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(par_txt(cap_c, bold = TRUE), style = "Normal") %>%
  body_add_flextable(ft_c) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(par_txt(foot, size = 9), style = "Normal") %>%
  body_end_block_section(value = landscape_section)

out_docx <- file.path(tab_dir, "Table_S3_community_composition.docx")
if (file.exists(out_docx)) archive_files(out_docx)
print(doc, target = out_docx)
say("Wrote: ", out_docx)
say("")

### -----------------------------------------------------------------------------
### 11) Figures for the Methods paragraph on the community
###
### Printed so the sentence can be typed from the log rather than from a CSV.
### -----------------------------------------------------------------------------
say("For the Methods paragraph (action 3.6):")
for (g in names(GEARS)) {
  sub  <- by_taxon %>% filter(gear == g)
  grps <- panel_a_full %>% filter(gear == g)
  pick <- function(nm) {
    v <- grps$pct[grps$group == nm]
    if (length(v)) v else 0
  }
  top3 <- top_overall %>% filter(gear == g) %>% arrange(rank) %>% slice_head(n = 3)
  say("  ", GEARS[[g]], ": ", nrow(sub), " taxa; fishes ",
      sprintf("%.0f%%", pick("Fishes")), " of catch, invertebrates ",
      sprintf("%.0f%%", pick("Invertebrates")), "; leading taxa ",
      paste(top3$label_display, collapse = ", "))
}
say("")

### -----------------------------------------------------------------------------
### 12) Provenance
### -----------------------------------------------------------------------------
say("Inputs:")
for (f in c(IN_RDATA, unname(IN_EVENT), IN_SETTINGS)) say("  ", f, "   mtime ", format(file.mtime(f)))
end_log("finished")
cat("h11_community_composition.R completed\n")
