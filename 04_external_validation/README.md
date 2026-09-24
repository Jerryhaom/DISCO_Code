# 04_external_validation

**Role:** Saved-object validation / plotting.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 04_external_validation/run_saved_object_audit.R
Rscript run.R 04_external_validation/clean_clhls_status.R
Rscript run.R 04_external_validation/run_figures.R
```

## Inputs

`Result4/CHARLS.DISCO.rds`, `CLHLS.DISCO.rds`, `RuLAS.DISCO.rds`; cohort identifiers differ (id/ID), age/mc/hd and survival fields are required. HRS figures additionally need DNAm and comparison inputs. clean_clhls_status writes a separate cleaned object and does not update the historical input automatically.

## Outputs and interpretation

outputs/ and outputs/figures/. Check the saved-object audit before interpreting figures.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
