# =============================================================================
# Script:  a00_getdata.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Reads the raw TPWD bag seine (BS) and bay trawl (BT) organism and station
#   tables and the species list, checks them, and merges species codes 615 and
#   212 into 333. Keeps estuaries 2 to 8 (Sabine Lake, bay 1, is dropped here,
#   once, for every later analysis) and removes colonial and gelatinous
#   invertebrates, sea squirts, unidentified codes and taxa outside the target
#   taxon groups (section 5b). Data.RData is the input of the whole pipeline;
#   run_pipeline.R skips this script, which is rerun only when the raw data or
#   these steps change.
#
# Run after: none (first script). Data.RData is read by b03, b04, e02, f07, h11.
#
# Inputs:
#   Raw TPWD data in the folder "0001 TPWD Data Updated" beside the project
#   folder (not made by a script; see RAW_DATA_REL):
#     02 Bag Seine/TPWD BS 05-24.xlsx, TPWD BS 82-04.xlsx
#     03 Bay Trawl/TPWD BT 03-24.xlsx, TPWD BT 82-03.xlsx
#     04 Station Data/TPWD BS 82-24 Corrected 24 hour format.xlsx,
#                     TPWD BT 82-24 Hydro with time.xlsx
#     06 Species Data/TPWD_Spp_codes_MF.xlsx
#
# Outputs:
#   01_Analysis/Data.RData with
#     bagseine_organisms, baytrawl_organisms  SAMPLE_ID, SPECIES_CODE (text), catch
#     STATION_BS, STATION_BT                  all station columns; bay column as
#                                             integer; HYDRO_ID renamed SAMPLE_ID
#     species                                 the TPWD species list, plus codes
#                                             identified by TPWD that the list
#                                             lacks (SPP_SUPPLEMENT)
#     sabine_summary                          Sabine Lake sampling and salinity,
#                                             kept for the Methods text after the
#                                             bay is dropped
#     removed_taxa                            what the taxon exclusions removed,
#                                             per code and gear (bays 2 to 8)
#   An existing Data.RData is moved to 01_Analysis/ARCHIVE_<timestamp>/ first.
#   02_Results/03_Diagnostics/a00_getdata/: columns_<table>.csv,
#     columns_compare_<table>.csv, catch_missing_by_species_<gear>.csv,
#     sabine_lake_summary.csv, nonfish_taxa_review.csv, removed_taxa.csv,
#     species_codes_not_in_list.csv, irregular_values.csv
#   Run log: 01_Analysis/Logs/a00_getdata/run_<timestamp>.txt
#
# Method notes:
#   - Gill net data are not used; only bag seine and bay trawl.
#   - Codes 615, 212 and 333 are merged because TPWD coded the same taxon
#     differently over the survey period; only catch is summed in the merge.
#   - Colonial and gelatinous invertebrates recorded mostly as presence only
#     are removed, because coverage-based Hill numbers need counts of
#     individuals (reasons in section 1, next to EXCLUDE_TAXA). Counted
#     jellyfish are kept. Sessile animals that are not colonial (oysters,
#     mussels, anemones) are kept because they are counted as individuals.
#   - The raw tables are checked (column names, types, missing and irregular
#     values, columns renamed between file versions); nothing is changed by the
#     checks, they are only written to the diagnostics folder.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Find the project folder and build PATHS
# -----------------------------------------------------------------------------
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

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

script_stem <- "a00_getdata"
dirs <- script_dirs(script_stem, c("logs", "diag"))
start_log(dirs$logs, script_stem)

# -----------------------------------------------------------------------------
# 1. Settings and raw data files
# -----------------------------------------------------------------------------
# Bays kept for every analysis. Bay 1 (Sabine Lake) is dropped: it was sampled
# from a later start and is fresher than the other bays (see the summary saved
# in section 5).
KEEP_BAYS <- 2:8
DROP_BAYS <- 1L

# Taxa removed from all analyses (section 5b). Every row of these codes is
# removed, whether or not a count was recorded. They are colonial or gelatinous
# invertebrates: the data checks showed that most of their rows have no catch
# (presence recorded only) and the rest have a count, so whether they entered a sample's species list depended on how the record was
# kept, not on the community. Reasons for removing them:
#   1. Coverage-based Hill numbers need counts of individuals. A bryozoan colony
#      or a mass of comb jellies is not a count of individuals, and a
#      presence-only record cannot enter the coverage estimate.
#   2. Counting was inconsistent between samples, which could create trends in
#      diversity that reflect recording practice rather than the community.
#   3. Bag seines and trawls catch these taxa only incidentally, and their
#      numbers swing with blooms.
#   4. The study's targets are fishes and mobile invertebrates; sessile and
#      benthic taxa are caught incidentally and are not analysis targets.
# Codes are listed one by one, not by taxon group: group 6 (Invertebrates) also
# holds shrimp, crabs and other mobile invertebrates that ARE targets.
# To remove more codes, add rows here; nothing downstream needs to change.
# Also removed, each with its reason in the table:
#   - sea squirts (9055, 9266, 9234): according to TPWD, they may not have been
#     recorded in the early years of the survey.
#   - colonial cnidarians (9356 sea pansy, 9342 by-the-wind sailor, 9124
#     Portuguese man o' war): colonies are not counted as individuals.
#   - the seven codes not in the TPWD species list and still unidentified
#     (3000, 9114, 9276, 9277, 9278, 9279, 9306): they cannot be checked
#     against the rules above.
# Sessile animals that are not colonial (oysters, mussels, jingle shells, pen
# shells, anemones) are kept: they are counted as individuals.
EXCLUDE_TAXA <- data.frame(
  SPECIES_CODE = c("9319", "9343", "9287", "9312", "9249", "9102", "9095", "9105", "9208",
                   "9051", "9068", "9200", "9116", "9314", "9206", "9285", "9079", "9040",
                   "9055", "9266", "9234", "9356", "9342", "9124",
                   "3000", "9114", "9276", "9277", "9278", "9279", "9306"),
  taxon = c("Mnemiopsis mccradyi (comb jelly)", "Beroe ovata (comb jelly)",
            "Phylum Ctenophora (comb jellies)", "Chrysaora quinquecirrha (sea nettle)",
            "Nemopsis bachei (hydromedusa)", "Genus Bougainvillia (hydroid)",
            "Class Hydrozoa (hydroids)", "Hydractinia americana (hydroid)",
            "Order Hydroidea (hydroids)", "Phylum Bryozoa (moss animals)",
            "Zoobotryon verticillatum (bryozoan)", "Bugula neritina (bryozoan)",
            "Amathia alternata (bryozoan)", "Class Ascidiacea (sea squirts)",
            "Phylum Porifera (sponges)", "Order Alcyonacea (soft corals)",
            "Order Pennatulacea (sea pens)", "Oculina diffusa (coral)",
            "Molgula manhattensis (sea squirt)", "Styela plicata (pleated sea squirt)",
            "Ciona intestinalis (sea squirt)", "Renilla muelleri (sea pansy)",
            "Velella velella (by-the-wind sailor)", "Physalia physalis (Portuguese man o' war)",
            paste("code", c("3000", "9114", "9276", "9277", "9278", "9279", "9306"),
                  "(not in the species list)")),
  reason = c(rep("gelatinous; mostly presence only", 5), rep("colonial hydroid; mostly presence only", 4),
             rep("colonial bryozoan; mostly presence only", 4), "colonial or solitary sessile tunicate; mostly presence only",
             rep("sessile; mostly presence only", 4),
             rep("sea squirt; may not have been recorded in the early years", 3),
             rep("colonial cnidarian; colonies not counted as individuals", 3),
             rep("not in the species list; unidentified", 7)),
  stringsAsFactors = FALSE
)

# Taxon groups kept, using the group column of the TPWD species list
# (1 fish, 2 bird, 3 mammal, 4 amphibian, 5 reptile, 6 invertebrates,
# 10 vegetation; 0 is not in the list's legend and marks taxa identified only
# above species level, such as Genus Palaemonetes or Family Mugilidae, which
# hold both fishes and invertebrates). The study's targets are fishes and
# mobile invertebrates, so groups 0, 1 and 6 are kept and every other group is
# removed (for example turtles, alligators, frogs and birds). Codes missing
# from the species list have no group: those identified by TPWD are added to
# the list below (SPP_SUPPLEMENT) and kept; those still unidentified are in
# EXCLUDE_TAXA. Code 9256 is also missing from the list, but its rows have a
# catch of 0 and b04 drops them.
# Counted jellyfish (cannonball, moon jelly) are kept: they are counted as
# individuals and are part of the community the survey records.
# Sea squirts are removed through EXCLUDE_TAXA, not through this rule.
KEEP_TAXA_GROUPS <- c(0L, 1L, 6L)

# Codes missing from the TPWD species list, identified by TPWD; all first
# recorded in recent years of the survey. They are added to the species list
# read below, so that every later script sees their names and groups. Groups as
# in the list: 0 not identified to species, 1 fishes, 6 invertebrates.
SPP_SUPPLEMENT <- data.frame(
  code   = c("128", "437", "9273", "9292", "9357"),
  latin  = c("Opsopoeodus emiliae", "Family Hemiramphidae", "Cystobranchus vividus",
             "Angulus texanus", "Armases cinereum"),
  common = c("Pugnose minnow", "Family halfbeaks", "(Marine leech)",
             "Texas tellin", "Squareback marsh crab"),
  taxa   = c("1", "0", "6", "6", "6"),
  stringsAsFactors = FALSE)

# Raw TPWD data: the folder "0001 TPWD Data Updated" beside the project folder
# (both inside 00 TASKS). Built from the project folder found at start-up, so
# no user name or drive letter is stored. Change RAW_DATA_REL if the layout differs.
RAW_DATA_REL  <- file.path("..", "0001 TPWD Data Updated")
base_data_dir <- normalizePath(file.path(PATHS$project_dir, RAW_DATA_REL),
                               winslash = "/", mustWork = FALSE)

file_bs_1 <- file.path(base_data_dir, "02 Bag Seine", "TPWD BS 05-24.xlsx")
file_bs_2 <- file.path(base_data_dir, "02 Bag Seine", "TPWD BS 82-04.xlsx")
file_bt_1 <- file.path(base_data_dir, "03 Bay Trawl", "TPWD BT 03-24.xlsx")
file_bt_2 <- file.path(base_data_dir, "03 Bay Trawl", "TPWD BT 82-03.xlsx")
file_station_bs <- file.path(base_data_dir, "04 Station Data",
                             "TPWD BS 82-24 Corrected 24 hour format.xlsx")
file_station_bt <- file.path(base_data_dir, "04 Station Data",
                             "TPWD BT 82-24 Hydro with time.xlsx")
file_species <- file.path(base_data_dir, "06 Species Data", "TPWD_Spp_codes_MF.xlsx")

cat("Raw data folder: ", base_data_dir, "\n", sep = "")
stop_if_missing(c(file_bs_1, file_bs_2, file_bt_1, file_bt_2,
                  file_station_bs, file_station_bt, file_species),
                producer = "TPWD raw data (not made by a script)")

output_file <- PATHS$clean_data_rdata   # 01_Analysis/Data.RData

# The species list is read first, so the data checks below can show the name
# and taxon group of each code they report.
species <- readxl::read_excel(file_species)
sp_code_col  <- grep("code", names(species), ignore.case = TRUE, value = TRUE)[1]
sp_latin_col <- grep("latin", names(species), ignore.case = TRUE, value = TRUE)[1]
sp_common_col <- grep("common", names(species), ignore.case = TRUE, value = TRUE)[1]
sp_taxa_col  <- grep("^taxa", names(species), ignore.case = TRUE, value = TRUE)[1]
# Add the codes identified by TPWD (SPP_SUPPLEMENT) that the list lacks. Each
# column keeps the type it has in the list; other columns are left empty.
if (any(is.na(c(sp_code_col, sp_latin_col, sp_common_col, sp_taxa_col)))) {
  stop("Species list: code, latin, common or taxa column not found; cannot add SPP_SUPPLEMENT", call. = FALSE)
}
as_type <- function(x, like) if (is.numeric(like)) as.numeric(x) else as.character(x)
supp_new <- SPP_SUPPLEMENT[!SPP_SUPPLEMENT$code %in% trimws(as.character(species[[sp_code_col]])), ]
if (nrow(supp_new)) {
  add <- species[rep(NA_integer_, nrow(supp_new)), ]
  add[[sp_code_col]]   <- as_type(supp_new$code,   species[[sp_code_col]])
  add[[sp_latin_col]]  <- as_type(supp_new$latin,  species[[sp_latin_col]])
  add[[sp_common_col]] <- as_type(supp_new$common, species[[sp_common_col]])
  add[[sp_taxa_col]]   <- as_type(supp_new$taxa,   species[[sp_taxa_col]])
  species <- dplyr::bind_rows(species, add)
}
cat("Species list: ", nrow(species) - nrow(supp_new), " codes read, ", nrow(supp_new),
    " added from SPP_SUPPLEMENT (", paste(supp_new$code, collapse = ", "), ")\n", sep = "")
# Taxon groups as given in the species list's own column heading.
TAXA_GROUPS <- c("1" = "Fish", "2" = "Bird", "3" = "Mammal", "4" = "Amphibian", "5" = "Reptile",
                 "6" = "Invertebrates", "10" = "Vegetation")
species_names <- function(codes) {
  codes <- trimws(as.character(codes))
  if (is.na(sp_code_col)) return(data.frame(latin = NA_character_, common = NA_character_,
                                            taxon_group = NA_character_)[rep(1, length(codes)), ])
  i <- match(codes, trimws(as.character(species[[sp_code_col]])))
  grp <- if (is.na(sp_taxa_col)) NA_character_ else as.character(species[[sp_taxa_col]][i])
  data.frame(
    latin  = if (is.na(sp_latin_col))  NA_character_ else as.character(species[[sp_latin_col]][i]),
    common = if (is.na(sp_common_col)) NA_character_ else as.character(species[[sp_common_col]][i]),
    taxon_group = ifelse(is.na(i), "not in species list", ifelse(is.na(grp), NA_character_,
                         ifelse(grp %in% names(TAXA_GROUPS), paste(grp, TAXA_GROUPS[grp]), grp))),
    stringsAsFactors = FALSE, row.names = NULL)
}

# -----------------------------------------------------------------------------
# Data-check helpers (used in sections 2 to 4)
# -----------------------------------------------------------------------------
# column_summary(): one row per column with its type, missing values, number of
# distinct values, range (numbers and dates) and a few example values (text).
column_summary <- function(df) {
  data.frame(
    column    = names(df),
    type      = vapply(df, function(x) class(x)[1], character(1)),
    n_rows    = nrow(df),
    n_missing = vapply(df, function(x) sum(is.na(x)), numeric(1)),
    n_unique  = vapply(df, function(x) length(unique(x[!is.na(x)])), numeric(1)),
    min       = vapply(df, function(x) if (is.numeric(x) || inherits(x, c("Date", "POSIXct")))
                         format(suppressWarnings(min(x, na.rm = TRUE))) else "", character(1)),
    max       = vapply(df, function(x) if (is.numeric(x) || inherits(x, c("Date", "POSIXct")))
                         format(suppressWarnings(max(x, na.rm = TRUE))) else "", character(1)),
    examples  = vapply(df, function(x) if (is.character(x) || is.factor(x))
                         paste(utils::head(unique(as.character(x[!is.na(x)])), 3), collapse = " | ") else "",
                       character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Irregular values are collected in one table: table, check, number of rows
# affected, and a few examples. Nothing is removed or changed here; the
# summary at the end says which checks found something.
flags <- list()
flag <- function(table, check, hit, examples = NULL) {
  n <- sum(hit, na.rm = TRUE)
  ex <- if (n > 0 && length(examples)) paste(utils::head(unique(examples[which(hit)]), 5), collapse = " | ") else ""
  flags[[length(flags) + 1L]] <<- data.frame(table = table, check = check, n_rows = n,
                                             examples = ex, stringsAsFactors = FALSE)
}

check_raw_table <- function(df, name) {
  cs <- column_summary(df)
  utils::write.csv(cs, file.path(dirs$diag, paste0("columns_", name, ".csv")), row.names = FALSE)
  cat("\n", name, ": ", nrow(df), " rows, ", ncol(df), " columns\n", sep = "")
  print(cs[, c("column", "type", "n_missing", "n_unique", "min", "max")], row.names = FALSE)
  flag(name, "exact duplicate rows", duplicated(df))
  invisible(cs)
}

# compare_columns(): TPWD has renamed columns between file versions. When two
# tables are stacked with bind_rows(), a renamed column becomes two columns,
# each missing for the other file's rows, without any warning. This lists every
# column of the two tables side by side and flags columns found in only one,
# names that differ only in capitals, spaces or underscores (likely renames),
# and columns whose type differs.
compare_columns <- function(a, b, name_a, name_b, label) {
  cols <- union(names(a), names(b))
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  type_of <- function(df, cl) if (cl %in% names(df)) class(df[[cl]])[1] else ""
  cmp <- data.frame(
    column  = cols,
    in_a    = cols %in% names(a),
    in_b    = cols %in% names(b),
    type_a  = vapply(cols, function(cl) type_of(a, cl), character(1)),
    type_b  = vapply(cols, function(cl) type_of(b, cl), character(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  names(cmp)[2:5] <- c(paste0("in_", name_a), paste0("in_", name_b),
                       paste0("type_", name_a), paste0("type_", name_b))
  only_one <- xor(cols %in% names(a), cols %in% names(b))
  # A column found in only one table whose simplified name matches a column of
  # the other table is probably the same column renamed.
  near <- only_one & vapply(cols, function(cl) {
    other <- if (cl %in% names(a)) names(b) else names(a)
    norm(cl) %in% norm(setdiff(other, cl))
  }, logical(1))
  type_diff <- cols %in% intersect(names(a), names(b)) &
    vapply(cols, function(cl) type_of(a, cl) != type_of(b, cl), logical(1))
  cmp$note <- ifelse(near, "likely renamed (same name apart from capitals or spacing)",
              ifelse(only_one, "only in one table",
              ifelse(type_diff, "different type", "")))
  utils::write.csv(cmp, file.path(dirs$diag, paste0("columns_compare_", label, ".csv")), row.names = FALSE)
  flag(label, paste("columns only in", name_a, "or only in", name_b), only_one, cols)
  flag(label, "column names that differ only in capitals or spacing (likely renamed)", near, cols)
  flag(label, "columns with a different type in the two tables", type_diff,
       paste0(cols, " (", cmp[[4]], " vs ", cmp[[5]], ")"))
  if (any(only_one | type_diff)) {
    cat("  Column differences, ", label, ":\n", sep = "")
    d <- cmp[cmp$note != "", ]
    cat(sprintf("    %-20s %-5s %-5s %s", d$column, ifelse(d[[2]], "yes", "-"),
                ifelse(d[[3]], "yes", "-"), d$note), sep = "\n")
  } else {
    cat("  Columns of ", name_a, " and ", name_b, ": same names and types\n", sep = "")
  }
  invisible(cmp)
}

# Columns that are entirely missing within one raw file of a stacked table
# (another sign of a rename): checked for the columns used later.
flag_empty_in_file <- function(part, file_label, table, cols) {
  for (cl in intersect(cols, names(part)))
    flag(table, paste0(cl, " entirely missing in ", file_label), all(is.na(part[[cl]])))
}

# -----------------------------------------------------------------------------
# 2. Organism data: read and stack the two files per gear
# -----------------------------------------------------------------------------
cat("\nReading organism data (BS and BT)\n")
bs_1 <- read_excel(file_bs_1); bs_2 <- read_excel(file_bs_2)
bt_1 <- read_excel(file_bt_1); bt_2 <- read_excel(file_bt_2)

cat("\nColumn names of the two raw files per gear\n")
compare_columns(bs_1, bs_2, "BS_05_24", "BS_82_04", "organisms_BS_files")
compare_columns(bt_1, bt_2, "BT_03_24", "BT_82_03", "organisms_BT_files")
key_org_cols <- c("SAMPLE_ID", "SPECIES_CODE", "catch")
# Each raw file must have the key columns under these exact names; a renamed
# key column would otherwise lose its values for that file without an error.
expect_cols(bs_1, key_org_cols, "TPWD BS 05-24.xlsx")
expect_cols(bs_2, key_org_cols, "TPWD BS 82-04.xlsx")
expect_cols(bt_1, key_org_cols, "TPWD BT 03-24.xlsx")
expect_cols(bt_2, key_org_cols, "TPWD BT 82-03.xlsx")
flag_empty_in_file(bs_1, "BS 05-24", "organisms_BS", key_org_cols)
flag_empty_in_file(bs_2, "BS 82-04", "organisms_BS", key_org_cols)
flag_empty_in_file(bt_1, "BT 03-24", "organisms_BT", key_org_cols)
flag_empty_in_file(bt_2, "BT 82-03", "organisms_BT", key_org_cols)

bagseine_organisms <- bind_rows(bs_1, bs_2)
baytrawl_organisms <- bind_rows(bt_1, bt_2)

# -----------------------------------------------------------------------------
# 2a. Check the organism data
# -----------------------------------------------------------------------------
check_organisms <- function(df, name, part_1, part_2) {
  expect_cols(df, c("SAMPLE_ID", "SPECIES_CODE", "catch"), name)
  check_raw_table(df, name)
  code_txt <- trimws(as.character(df$SPECIES_CODE))
  flag(name, "SAMPLE_ID missing", is.na(df$SAMPLE_ID))
  flag(name, "SPECIES_CODE missing", is.na(df$SPECIES_CODE))
  flag(name, "SPECIES_CODE not a whole number", !is.na(code_txt) & !grepl("^[0-9]+$", code_txt), code_txt)
  flag(name, "catch missing", is.na(df$catch))
  flag(name, "catch negative", !is.na(df$catch) & df$catch < 0, df$catch)
  flag(name, "catch not a whole number", !is.na(df$catch) & abs(df$catch - round(df$catch)) > 1e-9, df$catch)
  flag(name, "catch zero (rows recorded with no fish)", !is.na(df$catch) & df$catch == 0)
  # Which species codes the rows without a catch value carry, and whether those
  # rows hold other values (for example a length), to see what they record.
  miss <- is.na(df$catch)
  if (any(miss)) {
    has_len <- if ("length_num" %in% names(df)) !is.na(df$length_num) else rep(NA, nrow(df))
    by_code <- data.frame(SPECIES_CODE = code_txt, miss = miss, has_len = has_len) |>
      dplyr::group_by(SPECIES_CODE) |>
      dplyr::summarise(rows = dplyr::n(), rows_catch_missing = sum(miss),
                       catch_missing_with_length = sum(miss & has_len, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(rows_catch_missing > 0) |>
      dplyr::arrange(dplyr::desc(rows_catch_missing))
    by_code <- cbind(by_code, species_names(by_code$SPECIES_CODE))
    utils::write.csv(by_code, file.path(dirs$diag, paste0("catch_missing_by_species_", name, ".csv")),
                     row.names = FALSE)
    cat("  Rows without catch: ", sum(miss), " (", sum(miss & has_len, na.rm = TRUE),
        " of them have a length), in ", nrow(by_code), " codes; taxon groups (rows): ",
        paste(names(tapply(by_code$rows_catch_missing, by_code$taxon_group, sum)),
              tapply(by_code$rows_catch_missing, by_code$taxon_group, sum), sep = " = ", collapse = "; "),
        "\n  Most common:\n", sep = "")
    top <- utils::head(by_code, 8)
    cat(sprintf("    %-6s %7d  %-15s %s (%s)", top$SPECIES_CODE, top$rows_catch_missing,
                top$taxon_group, top$latin, top$common), sep = "\n")
  }
  for (v in intersect(c("catch_per_hour", "length_num"), names(df)))
    flag(name, paste(v, "negative"), !is.na(df[[v]]) & df[[v]] < 0, df[[v]])
  flag(name, "same SAMPLE_ID and SPECIES_CODE on more than one row",
       duplicated(df[, c("SAMPLE_ID", "SPECIES_CODE")]),
       paste(df$SAMPLE_ID, df$SPECIES_CODE))
  # The two files of a gear cover different years; a sample in both would be
  # counted twice after stacking.
  both <- intersect(unique(part_1$SAMPLE_ID), unique(part_2$SAMPLE_ID))
  flag(name, "SAMPLE_ID present in both raw files", df$SAMPLE_ID %in% both, df$SAMPLE_ID)
  if ("YEAR" %in% names(part_1) && "YEAR" %in% names(part_2)) {
    cat("  Years by file: ", paste(range(part_1$YEAR, na.rm = TRUE), collapse = " to "), " and ",
        paste(range(part_2$YEAR, na.rm = TRUE), collapse = " to "), "\n", sep = "")
  }
}
check_organisms(bagseine_organisms, "organisms_BS", bs_1, bs_2)
check_organisms(baytrawl_organisms, "organisms_BT", bt_1, bt_2)
rm(bs_1, bs_2, bt_1, bt_2)

# -----------------------------------------------------------------------------
# 3. Merge species codes 615 and 212 into 333
# -----------------------------------------------------------------------------
# TPWD changed how these taxa are coded over the survey period, so the same
# fish appear under 615, 212 or 333 depending on the year. To keep one
# consistent taxon through time, all three are combined under 333: within each
# sample their catches are summed into one row coded 333. Only SAMPLE_ID,
# SPECIES_CODE and catch are kept afterwards, so only catch is summarised.
# Works whether SPECIES_CODE is read as text or as a number.
combine_species_615_212_333 <- function(df) {
  df <- df |> dplyr::mutate(SPECIES_CODE = as.character(SPECIES_CODE))
  target_codes <- c("615", "212", "333")

  target_samples <- df |>
    dplyr::filter(SPECIES_CODE %in% target_codes) |>
    dplyr::pull(SAMPLE_ID) |>
    unique()

  combined <- df |>
    dplyr::filter(SAMPLE_ID %in% target_samples, SPECIES_CODE %in% target_codes) |>
    dplyr::group_by(SAMPLE_ID) |>
    dplyr::summarize(catch = sum(catch, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(SPECIES_CODE = "333")

  df |>
    dplyr::filter(!(SAMPLE_ID %in% target_samples & SPECIES_CODE %in% target_codes)) |>
    dplyr::bind_rows(combined) |>
    dplyr::arrange(SAMPLE_ID)
}

n_merge <- function(df) sum(as.character(df$SPECIES_CODE) %in% c("615", "212", "333"))
cat("\nRows coded 615, 212 or 333 before the merge: BS ", n_merge(bagseine_organisms),
    ", BT ", n_merge(baytrawl_organisms), "\n", sep = "")

bagseine_organisms <- combine_species_615_212_333(bagseine_organisms) |>
  dplyr::select(SAMPLE_ID, SPECIES_CODE, catch)
baytrawl_organisms <- combine_species_615_212_333(baytrawl_organisms) |>
  dplyr::select(SAMPLE_ID, SPECIES_CODE, catch)

cat("Rows coded 333 after the merge:                BS ", n_merge(bagseine_organisms),
    ", BT ", n_merge(baytrawl_organisms), "\n", sep = "")

# -----------------------------------------------------------------------------
# 4. Station data: read and check
# -----------------------------------------------------------------------------
cat("\nReading station data (BS and BT)\n")
STATION_BS <- read_excel(file_station_bs)
STATION_BT <- read_excel(file_station_bt)

# The bay column is MAJOR_AREA_CODE in the BS table and MAJOR_AREA in the BT table.
check_stations <- function(df, name, bay_col) {
  expect_cols(df, c("HYDRO_ID", bay_col, "YEAR", "MONTH", "COMPLETION_DTTM"), name)
  check_raw_table(df, name)
  bay <- suppressWarnings(as.integer(df[[bay_col]]))
  flag(name, "HYDRO_ID missing", is.na(df$HYDRO_ID))
  flag(name, "HYDRO_ID on more than one row", duplicated(df$HYDRO_ID) & !is.na(df$HYDRO_ID), df$HYDRO_ID)
  flag(name, paste(bay_col, "missing or not a number"), is.na(bay), df[[bay_col]])
  flag(name, paste(bay_col, "outside bays 1 to 8 (dropped)"), !is.na(bay) & !(bay %in% 1:8), bay)
  # MONTH and YEAR are stored as text ("01" to "12"), so they are compared as numbers.
  month_num <- suppressWarnings(as.integer(df$MONTH))
  flag(name, "MONTH not 1 to 12", !is.na(df$MONTH) & !(month_num %in% 1:12), df$MONTH)
  flag(name, "YEAR missing", is.na(df$YEAR))
  flag(name, "COMPLETION_DTTM missing", is.na(df$COMPLETION_DTTM))
  # b03 reads COMPLETION_DTTM as UTC and converts it to Central time to get the
  # sampling year and month. These two checks show which reading agrees with
  # the YEAR and MONTH columns TPWD recorded.
  dt <- df$COMPLETION_DTTM
  if (inherits(dt, "POSIXct")) {
    ym <- sprintf("%04d-%02d", as.integer(df$YEAR), as.integer(df$MONTH))
    ym_as_stored  <- format(dt, "%Y-%m", tz = "UTC")
    ym_as_central <- format(dt, "%Y-%m", tz = "America/Chicago")
    flag(name, "YEAR/MONTH differ from COMPLETION_DTTM read as stored",
         !is.na(dt) & ym_as_stored != ym, paste(ym, format(dt, tz = "UTC")))
    flag(name, "YEAR/MONTH differ from COMPLETION_DTTM converted UTC to Central (as b03 does)",
         !is.na(dt) & ym_as_central != ym, paste(ym, format(dt, tz = "UTC")))
  }
  # Plausible ranges for estuarine water; values outside are reported, not removed.
  ranges <- list(TEMPERATURE = c(0, 40), SALINITY = c(0, 80), DEPTH = c(0, 30),
                 DISS_OXYGEN = c(0, 25), TURBIDITY = c(0, 1000))
  for (v in intersect(names(ranges), names(df))) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    flag(name, sprintf("%s outside %g to %g", v, ranges[[v]][1], ranges[[v]][2]),
         !is.na(x) & (x < ranges[[v]][1] | x > ranges[[v]][2]), x)
  }
  # Values that may be codes for "not measured" rather than measurements.
  if ("TURBIDITY" %in% names(df))
    flag(name, "TURBIDITY exactly 999 (possible code for not measured)",
         !is.na(df$TURBIDITY) & df$TURBIDITY == 999)
  if ("TEMPERATURE" %in% names(df))
    flag(name, "TEMPERATURE exactly 0", !is.na(df$TEMPERATURE) & df$TEMPERATURE == 0)
  cat("  Rows by bay: ", paste(names(table(bay, useNA = "ifany")), table(bay, useNA = "ifany"),
                              sep = "=", collapse = "  "), "\n", sep = "")
}
check_stations(STATION_BS, "stations_BS", "MAJOR_AREA_CODE")
check_stations(STATION_BT, "stations_BT", "MAJOR_AREA")
cat("\nColumn names of the BS and BT station tables\n")
compare_columns(STATION_BS, STATION_BT, "BS", "BT", "stations_BS_vs_BT")

# -----------------------------------------------------------------------------
# 5. Sabine Lake summary, then keep bays 2 to 8
# -----------------------------------------------------------------------------
# Kept so the Methods can describe Sabine Lake (when sampling began, how many
# months were sampled, its salinity) after the bay is dropped from the data.
# Mean salinity is given two ways: over all samples, and as the mean of
# monthly means (each year-month weighted equally).
sabine_one <- function(df, gear, bay_col) {
  s <- df[suppressWarnings(as.integer(df[[bay_col]])) %in% DROP_BAYS, , drop = FALSE]
  sal <- if ("SALINITY" %in% names(s)) suppressWarnings(as.numeric(s$SALINITY)) else rep(NA_real_, nrow(s))
  ym  <- sprintf("%04d-%02d", as.integer(s$YEAR), as.integer(s$MONTH))
  data.frame(
    gear = gear, bay = DROP_BAYS, n_samples = nrow(s),
    first_year = if (nrow(s)) min(s$YEAR, na.rm = TRUE) else NA,
    last_year  = if (nrow(s)) max(s$YEAR, na.rm = TRUE) else NA,
    n_months_sampled = length(unique(ym[!is.na(s$YEAR)])),
    salinity_mean_all_samples = mean(sal, na.rm = TRUE),
    salinity_mean_of_monthly_means = mean(tapply(sal, ym, mean, na.rm = TRUE), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
sabine_summary <- rbind(sabine_one(STATION_BS, "BS", "MAJOR_AREA_CODE"),
                        sabine_one(STATION_BT, "BT", "MAJOR_AREA"))
utils::write.csv(sabine_summary, file.path(dirs$diag, "sabine_lake_summary.csv"), row.names = FALSE)
cat("\nSabine Lake (bay 1), before it is dropped:\n")
for (i in seq_len(nrow(sabine_summary))) with(sabine_summary[i, ], cat(sprintf(
  "  %s: %d samples, %s to %s, %d months sampled; salinity mean %.1f (all samples), %.1f (mean of monthly means)\n",
  gear, n_samples, first_year, last_year, n_months_sampled,
  salinity_mean_all_samples, salinity_mean_of_monthly_means)))

STATION_BS <- STATION_BS |>
  dplyr::mutate(MAJOR_AREA_CODE = as.integer(MAJOR_AREA_CODE)) |>
  dplyr::filter(MAJOR_AREA_CODE %in% KEEP_BAYS) |>
  dplyr::rename(SAMPLE_ID = HYDRO_ID)

STATION_BT <- STATION_BT |>
  dplyr::mutate(MAJOR_AREA = as.integer(MAJOR_AREA)) |>
  dplyr::filter(MAJOR_AREA %in% KEEP_BAYS) |>
  dplyr::rename(SAMPLE_ID = HYDRO_ID)

# Organisms are kept only for samples of the kept bays.
n_bs_before <- nrow(bagseine_organisms); n_bt_before <- nrow(baytrawl_organisms)
bagseine_organisms <- bagseine_organisms |> dplyr::filter(SAMPLE_ID %in% STATION_BS$SAMPLE_ID)
baytrawl_organisms <- baytrawl_organisms |> dplyr::filter(SAMPLE_ID %in% STATION_BT$SAMPLE_ID)
cat("\nOrganism rows kept for bays ", min(KEEP_BAYS), " to ", max(KEEP_BAYS), ": BS ",
    nrow(bagseine_organisms), " of ", n_bs_before, ", BT ", nrow(baytrawl_organisms),
    " of ", n_bt_before, "\n", sep = "")

# -----------------------------------------------------------------------------
# 5b. Remove presence-only colonial and gelatinous invertebrates (EXCLUDE_TAXA)
#     and taxa outside the target taxon groups (KEEP_TAXA_GROUPS); section 1
# -----------------------------------------------------------------------------
# See the reasons next to EXCLUDE_TAXA. A sample whose only records were such
# taxa keeps its station row and becomes a sample with no target catch (b04
# gives it a zero-catch row), which is what it was for the target taxa.
# Before removing, every non-fish code in the data (bays kept) is listed for
# review, with how often it occurs and how often its catch is missing.
nonfish <- dplyr::bind_rows(
  dplyr::mutate(bagseine_organisms, gear = "BS"),
  dplyr::mutate(baytrawl_organisms, gear = "BT")) |>
  dplyr::group_by(SPECIES_CODE, gear) |>
  dplyr::summarise(rows = dplyr::n(), samples = dplyr::n_distinct(SAMPLE_ID),
                   rows_catch_missing = sum(is.na(catch)),
                   total_catch = sum(catch, na.rm = TRUE), .groups = "drop")
# Codes removed by the taxon-group rule: in the species list, with a group
# outside KEEP_TAXA_GROUPS.
grp_num <- if (is.na(sp_taxa_col)) rep(NA_integer_, nrow(species)) else
  suppressWarnings(as.integer(species[[sp_taxa_col]]))
group_excluded_codes <- trimws(as.character(species[[sp_code_col]]))[!is.na(grp_num) &
                                                                       !(grp_num %in% KEEP_TAXA_GROUPS)]
excluded_codes <- union(EXCLUDE_TAXA$SPECIES_CODE, group_excluded_codes)

nonfish <- cbind(nonfish, species_names(nonfish$SPECIES_CODE)) |>
  dplyr::filter(is.na(taxon_group) | taxon_group != "1 Fish") |>
  dplyr::mutate(share_catch_missing = round(rows_catch_missing / rows, 3),
                excluded = SPECIES_CODE %in% excluded_codes) |>
  dplyr::arrange(dplyr::desc(excluded), taxon_group, dplyr::desc(rows))
utils::write.csv(nonfish, file.path(dirs$diag, "nonfish_taxa_review.csv"), row.names = FALSE)
cat("\nNon-fish codes in the data (bays kept): ", length(unique(nonfish$SPECIES_CODE)),
    ", of which excluded: ", length(unique(nonfish$SPECIES_CODE[nonfish$excluded])),
    " (list for review: 03_Diagnostics/a00_getdata/nonfish_taxa_review.csv)\n", sep = "")

removed_taxa <- nonfish |>
  dplyr::filter(excluded) |>
  dplyr::transmute(SPECIES_CODE, gear, latin, common, taxon_group, rows,
                   rows_with_count = rows - rows_catch_missing, rows_catch_missing,
                   samples, catch_removed = total_catch) |>
  dplyr::left_join(EXCLUDE_TAXA[, c("SPECIES_CODE", "reason")], by = "SPECIES_CODE") |>
  dplyr::mutate(reason = ifelse(is.na(reason),
                                paste0("taxon group ", taxon_group, " (not a target group)"), reason))
utils::write.csv(removed_taxa, file.path(dirs$diag, "removed_taxa.csv"), row.names = FALSE)

n_bs_before <- nrow(bagseine_organisms); n_bt_before <- nrow(baytrawl_organisms)
bagseine_organisms <- bagseine_organisms |> dplyr::filter(!(SPECIES_CODE %in% excluded_codes))
baytrawl_organisms <- baytrawl_organisms |> dplyr::filter(!(SPECIES_CODE %in% excluded_codes))
cat("Rows removed (EXCLUDE_TAXA and non-target taxon groups): BS ",
    n_bs_before - nrow(bagseine_organisms),
    ", BT ", n_bt_before - nrow(baytrawl_organisms), "\n", sep = "")
for (g in c("BS", "BT")) {
  r <- removed_taxa[removed_taxa$gear == g, ]
  cat("  ", g, ": ", nrow(r), " codes, ", sum(r$rows), " rows (", sum(r$rows_with_count),
      " with a count, catch ", sum(r$catch_removed), ")\n", sep = "")
}
cat("  Removed by taxon group:\n")
rg <- removed_taxa[!(removed_taxa$SPECIES_CODE %in% EXCLUDE_TAXA$SPECIES_CODE), ]
if (nrow(rg)) cat(sprintf("    %-6s %s  %-16s %s (%s), rows %d", rg$SPECIES_CODE, rg$gear, rg$taxon_group,
                          rg$latin, rg$common, rg$rows), sep = "\n") else cat("    none\n")

# -----------------------------------------------------------------------------
# 6. Species list (read in section 1): checks
# -----------------------------------------------------------------------------
check_raw_table(species, "species")
code_col <- sp_code_col
if (!is.na(code_col)) {
  sp_codes <- trimws(as.character(species[[code_col]]))
  flag("species", paste(code_col, "on more than one row"), duplicated(sp_codes) & !is.na(sp_codes), sp_codes)
  org_codes <- unique(c(bagseine_organisms$SPECIES_CODE, baytrawl_organisms$SPECIES_CODE))
  flag("species", "organism SPECIES_CODE not in the species list",
       !(org_codes %in% sp_codes), org_codes)
  # All organism codes missing from the species list, with how often they occur.
  unknown <- dplyr::bind_rows(
    dplyr::mutate(bagseine_organisms, gear = "BS"),
    dplyr::mutate(baytrawl_organisms, gear = "BT")) |>
    dplyr::filter(!(SPECIES_CODE %in% sp_codes)) |>
    dplyr::group_by(SPECIES_CODE, gear) |>
    dplyr::summarise(rows = dplyr::n(), samples = dplyr::n_distinct(SAMPLE_ID),
                     total_catch = sum(catch, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(rows))
  utils::write.csv(unknown, file.path(dirs$diag, "species_codes_not_in_list.csv"), row.names = FALSE)
  if (nrow(unknown)) {
    cat("Organism codes not in the species list (bays kept):\n")
    cat(sprintf("  %-8s %s  rows %6d  samples %6d  total catch %8.0f", unknown$SPECIES_CODE, unknown$gear,
                unknown$rows, unknown$samples, unknown$total_catch), sep = "\n")
  }
} else {
  cat("No column with 'code' in its name in the species list; code checks skipped.\n")
}

# -----------------------------------------------------------------------------
# 7. Data-check summary and checks before saving
# -----------------------------------------------------------------------------
flags_df <- do.call(rbind, flags)
utils::write.csv(flags_df, file.path(dirs$diag, "irregular_values.csv"), row.names = FALSE)
cat("\nIrregular values (all checks; details in 03_Diagnostics/a00_getdata/irregular_values.csv):\n")
cat(sprintf("  %8d  %-13s %s", flags_df$n_rows, flags_df$table, flags_df$check), sep = "\n")
cat("  Checks that found something: ", sum(flags_df$n_rows > 0), " of ", nrow(flags_df), "\n", sep = "")

cat("\nRows saved\n")
cat("  Station BS:          ", nrow(STATION_BS), "\n")
cat("  Station BT:          ", nrow(STATION_BT), "\n")
cat("  Bag seine organisms: ", nrow(bagseine_organisms), "\n")
cat("  Bay trawl organisms: ", nrow(baytrawl_organisms), "\n")
cat("  Species list:        ", nrow(species), "\n")
cat("  Bays in BS stations: ", paste(sort(unique(STATION_BS$MAJOR_AREA_CODE)), collapse = " "), "\n")
cat("  Bays in BT stations: ", paste(sort(unique(STATION_BT$MAJOR_AREA)), collapse = " "), "\n")

stopifnot(
  "a table is empty" = all(c(nrow(STATION_BS), nrow(STATION_BT), nrow(bagseine_organisms),
                             nrow(baytrawl_organisms), nrow(species)) > 0),
  "SAMPLE_ID missing in a station table" = !anyNA(STATION_BS$SAMPLE_ID) && !anyNA(STATION_BT$SAMPLE_ID),
  "a dropped bay is still present" = !any(c(STATION_BS$MAJOR_AREA_CODE, STATION_BT$MAJOR_AREA) %in% DROP_BAYS),
  "species 615 or 212 left after the merge" =
    !any(c(bagseine_organisms$SPECIES_CODE, baytrawl_organisms$SPECIES_CODE) %in% c("615", "212")),
  "an excluded taxon is still present" =
    !any(c(bagseine_organisms$SPECIES_CODE, baytrawl_organisms$SPECIES_CODE) %in% excluded_codes)
)

# -----------------------------------------------------------------------------
# 8. Save (the previous Data.RData is archived first)
# -----------------------------------------------------------------------------
archive_files(output_file)
save(bagseine_organisms, baytrawl_organisms, STATION_BS, STATION_BT, species, sabine_summary,
     removed_taxa, file = output_file)
cat("\nSaved: ", output_file, "\n", sep = "")

end_log()
