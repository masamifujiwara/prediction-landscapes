# prediction-landscapes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24348198.svg)](https://doi.org/10.5281/zenodo.24348198)

**Author:** Masami Fujiwara
**Affiliation:** Department of Ecology and Conservation Biology, Texas A&M University, College Station, TX 77843-2258, USA

R code and prepared data for the paper "Forecasting estuarine biodiversity and population abundance for conservation: spatial scale, forecast horizon, and choice of measure" (revised version, under review at *Biodiversity and Conservation*). The analysis uses the Texas Parks and Wildlife Department (TPWD) Coastal Fisheries monitoring data from seven Texas estuaries (1982 to 2024). It forecasts biodiversity (Hill numbers of order q = 0, 1, and 2) and the abundance of eight common taxa, for each estuary and for the coast as a whole, 1 to 24 months ahead, and evaluates the forecasts by their accuracy and by their skill against two references.

## Citation

Until the paper is published, please cite the preprint of the submitted version:

Fujiwara, M. (2026). A Predictability Framework for Conservation Strategy: Empirical Evidence from Estuarine Biodiversity Forecasting. *bioRxiv*. https://doi.org/10.64898/2026.04.20.719659

The citation will be updated when the revised paper is published. The code in this repository is the version used for the revised paper.

The code and prepared data are archived on Zenodo. To cite the version used for the revised paper (2.0.0), use https://doi.org/10.5281/zenodo.24348198. To cite all versions, use https://doi.org/10.5281/zenodo.19665944.

## Repository contents

```
prediction-landscapes/
├── README.md
├── LICENSE
├── .project_root                 marker file: scripts look upward for it to find the project folder
└── 01_Analysis/
    ├── Data.RData                prepared TPWD data (starting point of the analysis)
    ├── a00_getdata.R ... h20_table_docx.R   59 analysis scripts, run in the order below
    ├── lib/                      functions shared by the scripts
    ├── drivers/                  run_pipeline.R (runs everything in order), check_order.R
    ├── Data/reference_submitted/ values from the tables of the submitted version (read by h16)
    └── figures/source/           drawn source of Figure 1 (read by g01)
```

Each script begins with a header giving its purpose, inputs, outputs, and notes on the method.

## Running the analysis

### Setup

1. Download or clone the repository. The top folder can have any name and sit anywhere; keep `.project_root` at its top.
2. Install R and the packages listed under Software. Each script checks for its packages and stops with a message if one is missing.
3. Run `01_Analysis/b01_setup_paths_packages.R` once. It checks the setup and creates the output folders. Every script finds the project folder by itself, from RStudio (Source), `Rscript`, or `nohup`, on Windows, macOS, or Linux.
4. Run the scripts in the order below, or use `01_Analysis/drivers/run_pipeline.R`.

### Running everything with the driver

`drivers/run_pipeline.R` runs each script in its own R process, in run order, and stops at the first failure. Its settings are near the top:

- `SKIP` holds `a00_getdata.R` by default (see Data below).
- `SKIP_MODELS` is `TRUE` by default and skips the forecast-model scripts (c07 to c12 and e05a to e08), which take days. **For a run from scratch set `SKIP_MODELS <- FALSE`**, because the later scripts need the forecasts these scripts write.
- `START_AT` and `STOP_AFTER` run part of the pipeline.
- `CHECK_ORDER` runs `drivers/check_order.R` first, which checks that every script reads only outputs of scripts that run before it.

Logs of every run are written to `01_Analysis/Logs/`.

### Run order

The run order is file-name order, with one exception: `f02_slow_pattern_diagnostics.R` runs before `f01_forecast_scoring.R`, because f01 reads f02's output. The driver applies this exception (`lib/run_order.R`). Gaps in the numbering are scripts that are not part of the final analysis.

**Data preparation**

| Script | Purpose |
|---|---|
| `a00_getdata.R` | Documents how `Data.RData` was built from the raw TPWD files (not included; see Data) |
| `b01_setup_paths_packages.R` | Setup check and output folders |
| `b03_harmonize_station_metadata.R` | Station table, one row per sample |
| `b04_join_organisms_to_station.R` | Catch records joined to samples |
| `b05_define_spatiotemporal_units.R` | Estuary-month and coast-month units |
| `b06_compute_hill_alpha_q012.R` | Hill numbers for each estuary and month (coverage-standardized, iNEXT) |
| `b07_compute_hill_gamma_q012.R` | Hill numbers for the coast as a whole |
| `b08_build_environment_timeseries.R` | Monthly water temperature and salinity series |
| `b09_forecast_settings.R` | Settings shared by all later scripts (forecast dates, lead times, seed) |

**Biodiversity forecasts**

| Script | Purpose |
|---|---|
| `c00_inputs_check.R` | Checks the biodiversity series before the models run |
| `c01_seasonality_characterization.R` | Seasonal curves of the coastwide series |
| `c02_deseasonalize_alpha_gamma_scaled.R` | Seasonal curve refitted at each forecast date; departures from it |
| `c03_forecast_design_backtesting.R` | Forecast dates of every series |
| `c07_alpha_forecasting_parallel_scaled.R` | SARIMA (automatic, without trend, with linear trend) and damped trend forecasts |
| `c08_alpha_state_space_models_forecasts_parallel_scaled.R` | State-space model forecasts, with and without trend |
| `c09_alpha_arimax_from_env_index_parallel_dynamic_scaled.R` | ARIMAX forecasts (lagged water temperature and salinity) |
| `c10_alpha_state_space_with_env_covariates_parallel_scaled.R` | State-space model with covariates |
| `c12_naive_baseline_diversity.R` | Reference forecasts: cumulative average, 12-month average, last month's value |
| `d01_compare_model_metrics_dictionary_v3_scaled.R` | Joins all biodiversity forecasts |

**Abundance forecasts**

| Script | Purpose |
|---|---|
| `e02_aggregate_abundance.R` | Monthly catch per unit effort of the eight taxa, by estuary and coastwide |
| `e03_deseasonalize_abundance.R` | Seasonal curve refitted at each forecast date; departures from it |
| `e04_forecast_design_abundance.R` | Forecast dates of every series |
| `e05a_baseline_forecasts_abundance.R` | SARIMA and damped trend forecasts |
| `e05b_naive_baseline_abundance.R` | Reference forecasts |
| `e06_state_space_forecasts_abundance.R` | State-space model forecasts |
| `e07_arimax_env_forecasts_abundance.R` | ARIMAX forecasts |
| `e08_ssm_xreg_forecasts_abundance.R` | State-space model with covariates |
| `e09_compare_predictability_abundance.R` | Joins all abundance forecasts |

**Accuracy, skill, and characteristics of the data**

| Script | Purpose |
|---|---|
| `f02_slow_pattern_diagnostics.R` | Slow change in the data (runs before f01) |
| `f01_forecast_scoring.R` | Accuracy and skill of every forecast, moving-block bootstrap ranges, horizons, management situations, variance partition |
| `f03_permutation_entropy.R` | Weighted permutation entropy |
| `f04_environment_trends.R` | Trends in water temperature and salinity |
| `f07_survey_precision.R` | Sampling uncertainty behind the tolerance |

**Figures**

| Script | Output |
|---|---|
| `g01_analysis_overview.R` | Fig. 1 |
| `g02_variance_budget.R` | Fig. 2 |
| `g05_misses_by_lead.R` | Fig. 3 |
| `g06_two_references.R` | Fig. 4 |
| `g03_scale_skill.R` | Fig. 5 |
| `g04_shared_surprises.R` | Fig. 6 |
| `g08_quadrants.R` | Fig. 7 |
| `g11_seasonal_fits.R` | Fig. S1 |
| `g14_horizon_tolerance.R` | Fig. S2 |
| `g12_slow_patterns.R` | Fig. S3 |
| `g17_entropy.R` | Fig. S4 |
| `g13_horizons.R` | Fig. S5 |
| `g15_interval_coverage.R` | Fig. S6 |
| `g07_diversity_periods.R` | Fig. S7 |
| `g16_taxon_skill.R` | Fig. S8 |
| `g20_figure_docx.R` | Word files: Figures 1 to 7 with captions, and Online Resource 3 (Figs. S1 to S8) |

**Tables**

| Script | Output |
|---|---|
| `h01_models.R` | Table 1 |
| `h12_sampling_design.R` | Table S1 |
| `h13_missing_months.R` | Table S2 |
| `h11_community_composition.R` | Table S3 |
| `h14_seasonal_fits.R` | Table S4 |
| `h16_published_measure.R` | Table S5 |
| `h02_summary.R` | Table S6 |
| `h15_all_models.R` | Tables S7 and S8 |
| `h17_environment_trends.R` | Table S9 |
| `h20_table_docx.R` | Word file: Online Resource 2 (Tables S1 to S9) |

The figure and table scripts run after all the others; among them the order does not matter, except that `g20` and `h20` run last.

### Outputs

```
<project folder>/
├── 01_Analysis/
│   ├── Data/<script>/        intermediate data, one folder per script
│   ├── Logs/<script>/        run logs, each ending with sessionInfo()
│   ├── tables/               tables (.docx, .csv)
│   └── figures/              figures (.tiff, .pdf) and figure documents (.docx)
└── 02_Results/               diagnostics and model objects, one folder per script
```

Most scripts name their output folder after themselves; a few keep an older folder name because later scripts read it (the header of each script gives its output folder). Scripts move earlier outputs into an `ARCHIVE_<timestamp>` subfolder before writing, so a rerun never silently overwrites earlier results.

### Run times

Times of the runs used for the paper, with 7 parallel workers. Scripts not listed take a few minutes or less.

| Script | Time |
|---|---|
| `c07` | about 3.5 hours |
| `c09` | about 11.5 hours |
| `c08`, `c10` | under 10 minutes each |
| `e05a` | about 10 hours |
| `e07` | about 38 hours |
| `e06`, `e08` | about 20 minutes each |
| `f01` | about 10 minutes |

Parallel scripts use one worker per physical core minus one, at most 7. Set `N_WORKERS_OVERRIDE` near the top of a script to use a different number (`c10` is set to 4). The long model scripts save checkpoints, so a stopped run resumes where it left off as long as the code, settings and inputs are unchanged.

### Notes on reproducing the results

- Random numbers come from a fixed seed (set in `b09`), so a rerun with the same software gives the same results. Small differences in fitted ARIMA models between machines are possible; each run log records `sessionInfo()`.
- `h16` compares the revised values with the tables of the submitted version, which are stored in `Data/reference_submitted/`.

## Software

The results were produced with R 4.6.1 and these packages: iNEXT 3.0.2, forecast 9.0.2, KFAS 1.6.0, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, zoo 1.9-0, future 1.75.0, future.apply 1.20.2, ggplot2 4.0.3, patchwork 1.3.2, ragg 1.5.2, officer 0.7.6, flextable 0.10.1. Some scripts also use readxl (a00 only), xml2, and rstudioapi (optional). The moving-block bootstrap, permutation entropy, Mann-Kendall test, and Sen slope are the project's own code (`lib/scoring.R`, `f02`, `f03`, `f04`).

## Data

The raw data and their initial processing are documented in `a00_getdata.R`, which reads the raw TPWD files. These files are not included, so `a00` cannot be rerun; the prepared data object `Data.RData` is the starting point of the analysis (start from `b01`). These data are provided only to allow validation and reproduction of the results in this study. They are not intended for use in new analyses or publications without prior permission from the Texas Parks and Wildlife Department.

## License

The code is licensed under the [MIT License](LICENSE). The data in this repository are not covered by this license and are subject to the restrictions above. Use of the data beyond validation of the results presented in this study requires prior permission from the Texas Parks and Wildlife Department.
