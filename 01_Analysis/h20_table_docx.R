# =============================================================================
# Script:  h20_table_docx.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Joins Tables S1 to S9 into one Word file, Online Resource 2, in order,
#   each with its caption and note, each starting on a new page and keeping
#   its own page orientation. The first page carries the Online Resource
#   title, article title, journal, author, affiliation and e-mail of the
#   corresponding author. Table 1 of the main text stays as the file h01
#   writes. Run last, after the other h scripts.
#
# Run after: h02, h11, h12, h13, h14, h15, h16, h17
#
# Inputs:
#   01_Analysis/tables/Table_S1_sampling_design.docx (h12)
#   01_Analysis/tables/Table_S2_missing_months.docx (h13)
#   01_Analysis/tables/Table_S3_community_composition.docx (h11)
#   01_Analysis/tables/Table_S4_seasonal_fits.docx (h14)
#   01_Analysis/tables/Table_S5_second_measure.docx (h16)
#   01_Analysis/tables/Table_S6_summary.docx (h02)
#   01_Analysis/tables/Table_S7_all_models_diversity.docx (h15)
#   01_Analysis/tables/Table_S8_all_models_abundance.docx (h15)
#   01_Analysis/tables/Table_S9_environment_trends.docx (h17)
#
# Outputs:
#   01_Analysis/tables/Online_Resource_Tables.docx
#   01_Analysis/Logs/h20_table_docx/run_<timestamp>.txt
#   An earlier Online_Resource_Tables.docx is moved to ARCHIVE_<timestamp> first.
#
# Method notes:
#   The tables are not rebuilt here: each h script writes its own table file
#   with its caption, and this script joins those files, so a table is edited
#   in one place only. The list TABLES is the only place that ties a table
#   number to a file.
#   The join works on the document content (officer and xml2): the body of
#   each file is appended, with a section break that carries the previous
#   file's page settings, so landscape and portrait pages survive. Word's own
#   document insertion (body_add_docx) is not used, because it needs Word to
#   finish the merge and cannot be checked elsewhere.
#   Every page is set to US letter, keeping its orientation; the empty
#   portrait section that the landscape files of h11 and h12 end with is
#   dropped.
#   The script stops if a table file is missing, contains links or images
#   (the join copies text and tables only), or lacks its label, and checks
#   that the labels appear in order with one section per table file.
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
need <- c("officer", "xml2")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(officer); library(xml2)})

OR_TABLES_TITLE <- "Online Resource 2: Tables"   # Online Resources: 1 Supplementary methods and results, 2 Tables, 3 Figures
ARTICLE_HEADER <- c(   # first page of each Online Resource file (journal requirement: title, journal, author, affiliation, e-mail)
  "Forecasting estuarine biodiversity and population abundance for conservation: spatial scale, forecast horizon, and choice of measure",
  "Biodiversity and Conservation",
  "Masami Fujiwara",
  "Department of Ecology and Conservation Biology, Texas A&M University, College Station, TX 77843-2258, USA",
  "Corresponding author: Masami Fujiwara, fujiwara@tamu.edu")
script_stem <- "h20_table_docx"
tab_dir <- file.path(PATHS$analysis_dir, "tables")
TABLES <- c("Table S1" = "Table_S1_sampling_design.docx", "Table S2" = "Table_S2_missing_months.docx",
            "Table S3" = "Table_S3_community_composition.docx", "Table S4" = "Table_S4_seasonal_fits.docx",
            "Table S5" = "Table_S5_second_measure.docx", "Table S6" = "Table_S6_summary.docx",
            "Table S7" = "Table_S7_all_models_diversity.docx", "Table S8" = "Table_S8_all_models_abundance.docx",
            "Table S9" = "Table_S9_environment_trends.docx")
dirs <- script_dirs(script_stem, c("logs"))
start_log(dirs$logs, script_stem)
paths <- file.path(tab_dir, TABLES)
stop_if_missing(paths, c("h12_sampling_design.R", "h13_missing_months.R", "h11_community_composition.R", "h14_seasonal_fits.R",
                         "h16_published_measure.R", "h02_summary.R", "h15_all_models.R", "h15_all_models.R",
                         "h17_environment_trends.R"))
out <- file.path(tab_dir, "Online_Resource_Tables.docx")
if (file.exists(out)) archive_files(out)
NS <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")

# body children of a table file: content nodes and the section settings.
# Files written with a landscape block (h11, h12) end with that block's section
# inside the last paragraph, followed by an empty default section: the empty
# section is dropped and the block's settings are used.
has_text <- function(nodes) any(vapply(nodes, function(n) xml_name(n) == "tbl" || nzchar(trimws(xml_text(n))), logical(1)))
parts_of <- function(path, label) {
  body <- xml_find_first(docx_body_xml(read_docx(path)), "//w:body", NS)
  if (length(xml_find_all(body, ".//@r:id", c(NS, r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"))))
    stop(basename(path), " contains links or images; the join copies text and tables only", call. = FALSE)
  if (!grepl(label, xml_text(body), fixed = TRUE)) stop(basename(path), " does not contain its label '", label, "'", call. = FALSE)
  kids <- xml_children(body); is_sect <- xml_name(kids) == "sectPr"
  content <- kids[!is_sect]; sect <- kids[is_sect][[1]]
  inner <- which(vapply(content, function(n) !inherits(xml_find_first(n, "./w:pPr/w:sectPr", NS), "xml_missing"), logical(1)))
  if (length(inner) > 1) stop(basename(path), " has more than one inner section; join not supported", call. = FALSE)
  if (length(inner) == 1) {
    after <- if (inner[1] < length(content)) content[(inner[1] + 1):length(content)] else list()
    if (length(after) && has_text(after)) stop(basename(path), " has content after its inner section; join not supported", call. = FALSE)
    sect <- xml_find_first(content[[inner[1]]], "./w:pPr/w:sectPr", NS)
    keep <- seq_len(inner[1] - 1)
    content <- content[keep]
  }
  list(content = content, sect = sect)
}
# US letter for every page, keeping each section's orientation
to_letter <- function(sect) {
  pg <- xml_find_first(sect, "./w:pgSz", NS)
  if (inherits(pg, "xml_missing")) stop("A section has no page size", call. = FALSE)
  land <- identical(xml_attr(pg, "orient"), "landscape") || as.numeric(xml_attr(pg, "w")) > as.numeric(xml_attr(pg, "h"))
  xml_set_attr(pg, "w:w", if (land) "15840" else "12240"); xml_set_attr(pg, "w:h", if (land) "12240" else "15840")
  if (land) xml_set_attr(pg, "w:orient", "landscape") else xml_set_attr(pg, "w:orient", NULL)
  invisible(sect)
}

# start: title paragraph; its page takes the settings of Table S1
doc <- read_docx() |> body_add_fpar(fpar(ftext(OR_TABLES_TITLE, fp_text(font.size = 14, bold = TRUE, font.family = "Times New Roman"))))
for (h in ARTICLE_HEADER) doc <- body_add_fpar(doc, fpar(ftext(h, fp_text(font.size = 11, font.family = "Times New Roman"))))
body <- xml_find_first(docx_body_xml(doc), "//w:body", NS)
final_sect <- function() xml_find_first(body, "./w:sectPr", NS)
for (k in seq_along(TABLES)) {
  p <- parts_of(paths[k], names(TABLES)[k])
  if (k > 1) {
    # close the previous table's section with its own page settings (next table starts a new page)
    brk <- read_xml(paste0('<w:p xmlns:w="', NS[["w"]], '"><w:pPr>', as.character(final_sect()), '</w:pPr></w:p>'))
    xml_add_sibling(final_sect(), brk, .where = "before")
  }
  for (n in p$content) xml_add_sibling(final_sect(), n, .where = "before")
  xml_replace(final_sect(), p$sect)
}
# every section starts on a new page, on US letter
for (s in xml_find_all(body, ".//w:sectPr", NS)) {
  to_letter(s)
  t <- xml_find_first(s, "./w:type", NS)   # officer writes "continuous" for its landscape blocks
  if (inherits(t, "xml_missing")) {
    pg <- xml_find_first(s, "./w:pgSz", NS)
    xml_add_sibling(pg, read_xml(paste0('<w:type xmlns:w="', NS[["w"]], '" w:val="nextPage"/>')), .where = "before")
  } else xml_set_attr(t, "w:val", "nextPage")
}
print(doc, target = out)

# check: every label present once in order, number of tables and sections
# (a caption paragraph starts with its label followed by a space)
chk <- read_docx(out); paras <- docx_summary(chk)$text
pos <- vapply(names(TABLES), function(l) { k <- which(startsWith(paras, paste0(l, " "))); if (length(k)) k[1] else -1 }, numeric(1))
if (any(pos < 0) || is.unsorted(pos)) stop("Table labels missing or out of order in the joined file", call. = FALSE)
n_tab <- length(xml_find_all(docx_body_xml(chk), "//w:tbl", NS)); n_sec <- length(xml_find_all(docx_body_xml(chk), "//w:sectPr", NS))
cat("Joined", length(TABLES), "files:", n_tab, "tables,", n_sec, "sections\n")
if (n_sec != length(TABLES)) stop("Expected one section per table file", call. = FALSE)
cat("Saved:", out, "\n")
end_log("finished")
cat("h20_table_docx.R completed\n")
