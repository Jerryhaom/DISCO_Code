# 05_organ_network

**Role:** Organ saved-score analyses with expanded covariates; figure replay.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 05_organ_network/run_saved_object_audit.R
Rscript run.R 05_organ_network/run_organ_expanded_covariates.R
Rscript run.R 05_organ_network/run_organ_cause_specific.R
Rscript run.R 05_organ_network/run_organ_chronic_disease.R
Rscript run.R 05_organ_network/plot_fig4a_expanded.R
```

## Inputs

Saved organ scores/networks under Result3, UKB mortality, covariates and disease outcomes. IDs must join consistently; confirm time/status/cause and organ score columns against each script input call. Network plotting requires the precomputed tissue mapping and network objects, not just biomarker matrices.

## Outputs and interpretation

outputs/organ_expanded_covariates and outputs/figures. run_figures.R replays historical network figures; expanded-covariate figures are a distinct result family.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
