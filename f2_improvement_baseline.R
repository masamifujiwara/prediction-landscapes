# ==============================================================================
# Extract exact SARIMA % RMSE improvement for manuscript placeholders
# Targets: Bay Anchovy (601) and Atlantic Croaker (602)
# Conditions: Bag Seine (BS), horizon = 1, model = sarima
# ==============================================================================

library(dplyr)
library(readr)

# 1. Resolve path to the abundance metrics dictionary
resolve_path <- function(subpath, filename) {
  paths <- c(
    file.path("01_Analysis/Data", subpath, filename),
    file.path("Data",             subpath, filename),
    filename
  )
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  stop(paste("File not found:", filename, "in", subpath))
}

abd_path <- resolve_path("e9_compare_predictability_abundance", "metrics_dictionary_abundance.csv")

# 2. Load data
df_abd <- read_csv(abd_path, show_col_types = FALSE)

# 3. Filter for conditions and calculate mean skill
skill_summary <- df_abd %>%
  filter(
    model == "sarima",
    horizon == 1,
    GEAR == "BS",
    SPECIES_CODE %in% c("601", "602") 
  ) %>%
  mutate(
    # Map species codes to common names
    Species = ifelse(SPECIES_CODE == "601", "Bay Anchovy", "Atlantic Croaker"),
    # Standardize scale name
    Scale = ifelse(tolower(series) == "gamma", "Gamma", "Alpha")
  ) %>%
  group_by(Scale, Species) %>%
  # For Alpha, this averages the skill across all individual bays. 
  # For Gamma, it simply returns the single coastwide value.
  summarize(
    Avg_SARIMA_Skill_Pct = mean(skill_rmse_improve_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Species, Scale)

# 4. Print results formatted for the manuscript
cat("\n--- Values for Manuscript Placeholders ---\n")
print(skill_summary)

# Print as a ready-to-copy sentence
bay_anchovy_alpha <- skill_summary$Avg_SARIMA_Skill_Pct[skill_summary$Species == "Bay Anchovy" & skill_summary$Scale == "Alpha"]
croaker_alpha <- skill_summary$Avg_SARIMA_Skill_Pct[skill_summary$Species == "Atlantic Croaker" & skill_summary$Scale == "Alpha"]

bay_anchovy_gamma <- skill_summary$Avg_SARIMA_Skill_Pct[skill_summary$Species == "Bay Anchovy" & skill_summary$Scale == "Gamma"]
croaker_gamma <- skill_summary$Avg_SARIMA_Skill_Pct[skill_summary$Species == "Atlantic Croaker" & skill_summary$Scale == "Gamma"]

cat(sprintf("\nDrafted text with values inserted:\n"))
cat(sprintf("'For example, predicting the local abundance of highly dominant species such as Bay Anchovy and Atlantic Croaker yielded average Alpha-scale SARIMA skill improvements of %.2f%% and %.2f%%, respectively—closely tracking the 22.79%% improvement baseline established by local species richness. Upon spatial aggregation to the coastwide Gamma scale, relative predictive skill for these dominant populations amplified further to %.2f%% and %.2f%%, directly mirroring the portfolio-driven jump to 28.11%% seen in the regional diversity aggregate.'\n", 
            bay_anchovy_alpha, croaker_alpha, bay_anchovy_gamma, croaker_gamma))