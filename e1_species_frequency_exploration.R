# ─────────────────────────────────────────────────────────────────────────────
# 2) Load Packages
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) Directories
# ─────────────────────────────────────────────────────────────────────────────
script_stem <- "e1_species_frequency_exploration"
diag_out_dir <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir.create(diag_out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Load Data
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(PATHS$clean_data_rdata))
load(PATHS$clean_data_rdata)

# Extract total samples per gear to calculate percentages
total_bs_samples <- n_distinct(STATION_BS$SAMPLE_ID)
total_bt_samples <- n_distinct(STATION_BT$SAMPLE_ID)

cat("Total Bag Seine samples:", total_bs_samples, "\n")
cat("Total Bay Trawl samples:", total_bt_samples, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Calculate Frequency of Occurrence
# ─────────────────────────────────────────────────────────────────────────────

# Helper function to summarize occurrences
calc_freq <- function(org_df, total_samples, gear_label) {
  org_df %>%
    filter(!is.na(SPECIES_CODE) & catch > 0) %>%
    group_by(SPECIES_CODE) %>%
    summarise(
      n_samples_present = n_distinct(SAMPLE_ID),
      .groups = "drop"
    ) %>%
    mutate(
      GEAR = gear_label,
      freq_occ_pct = (n_samples_present / total_samples) * 100
    )
}

bs_freq <- calc_freq(bagseine_organisms, total_bs_samples, "BS")
bt_freq <- calc_freq(baytrawl_organisms, total_bt_samples, "BT")

# Combine and pivot so we can compare gears side-by-side
freq_combined <- bind_rows(bs_freq, bt_freq) %>%
  pivot_wider(
    id_cols = SPECIES_CODE, 
    names_from = GEAR, 
    values_from = c(n_samples_present, freq_occ_pct),
    values_fill = 0
  ) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE))

# ─────────────────────────────────────────────────────────────────────────────
# 6) Join with Species Taxonomy and Highlight Targets
# ─────────────────────────────────────────────────────────────────────────────
# Clean up species lookup column names based on TPWD_Spp_codes_MF.csv
species_lookup <- species %>%
  select(
    SPECIES_CODE = species_code,
    Scientific_Name = `Species (latin name)`,
    Common_Name = `Species (common name)`
  ) %>%
  mutate(SPECIES_CODE = as.character(SPECIES_CODE))

freq_final <- freq_combined %>%
  left_join(species_lookup, by = "SPECIES_CODE") %>%
  select(
    SPECIES_CODE, Common_Name, Scientific_Name, 
    freq_occ_pct_BS, freq_occ_pct_BT, 
    n_samples_present_BS, n_samples_present_BT
  ) %>%
  arrange(desc(freq_occ_pct_BS + freq_occ_pct_BT))

# Mark the species of interest for easy review
target_codes <- c("9618", "9600", "616", "601", "602", "604")

freq_final <- freq_final %>%
  mutate(
    Target_Species = if_else(SPECIES_CODE %in% target_codes, "YES", "")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7) Output to Diagnostics
# ─────────────────────────────────────────────────────────────────────────────
out_file <- file.path(diag_out_dir, "species_frequency_summary.csv")
write_csv(freq_final, out_file)
cat("✅ Frequency summary saved to:", out_file, "\n")

# Print the stats for your currently selected targets to the console
cat("\n--- Currently Identified Target Species ---\n")
print(freq_final %>% filter(Target_Species == "YES") %>% as.data.frame())