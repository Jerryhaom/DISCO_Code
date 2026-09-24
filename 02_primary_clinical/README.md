# 02_primary_clinical

**Role:** Saved-object analysis / plotting.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 02_primary_clinical/run_saved_object_audit.R
Rscript run.R 02_primary_clinical/run_figures.R
Rscript run.R 02_primary_clinical/run_fig3ab.R
```

## Inputs

`Result1/NHANES.Clinical.DISCO.rds`, `NHANES.Clinical.DISCO2.rds`, UKB clinical and DNAm comparison objects. NHANES plot data use `id`, `Age`, `mc`, `hd`, `en`; survival objects supply time/status and score columns. No raw-to-score clinical pipeline is claimed.

## Outputs and interpretation

outputs/ and outputs/figures/; Fig 1c/d, Fig 3a/b and audit tables.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
