# 08_sensitivity_analyses

**Role:** Covariate, robustness, meta-analysis and genetic-interaction reanalysis.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 08_sensitivity_analyses/01_covariate_models.R
Rscript run.R 08_sensitivity_analyses/02_raw_biomarker_comparison.R
Rscript run.R 08_sensitivity_analyses/03_correlation_robustness.R
Rscript run.R 08_sensitivity_analyses/06_crosscohort_covariates.R
Rscript run.R 08_sensitivity_analyses/04_meta_analysis.R approved/table1/data/Table1_analysis_sample_SD.csv /path/to/new/meta --allow-provisional
Rscript run.R 08_sensitivity_analyses/05_genetic_interaction_rebuilt.R
Rscript run.R 08_sensitivity_analyses/compute_R05_stratified_2x2_rebuilt.R
Rscript run.R 08_sensitivity_analyses/plot_R05_forest_interaction_rebuilt.R
Rscript run.R 08_sensitivity_analyses/plot_R05_forest_stratified_rebuilt.R
Rscript run.R 08_sensitivity_analyses/plot_R05_heatmap_2x2_rebuilt.R
```

## Inputs

R01/R02/R03/R06 are separate branches with their own controlled inputs. R04 reads one explicitly versioned Table 1 CSV and refuses mixed models, transformations or SD populations; it no longer embeds historical HRs. Recompute its cohort fits with the approved Table 1 entry described in `../approved/README.txt`. R05 uses UKB_PRS_full.csv, dat_disease_28_rebuilt.rds and clinical/proteomic saved DISCO. Rebuilt interaction and stratified computation precede their matching plots. make_R05_table.R uses the old output family and is historical. add_R05_supp_table.R reads rebuilt estimates plus a supplied workbook template and writes an updated copy under the output root.

## Outputs and interpretation

outputs/R01_* through R06_*. Keep rebuilt and historical R05 outputs separate.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
