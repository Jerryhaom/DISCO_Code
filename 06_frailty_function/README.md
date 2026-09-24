# 06_frailty_function

**Role:** Saved-object analysis / plotting.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 06_frailty_function/run_saved_object_audit.R
Rscript run.R 06_frailty_function/run_figures.R
```

## Inputs

Saved Result5 CHARLS/CLHLS/RuLAS function and frailty objects. They contain cohort IDs, age/sex, DISCO and cohort-specific FI/ADL/MMSE fields; see script_input_calls.csv and column_usage.csv for exact accesses. No raw frailty-score construction is claimed.

## Outputs and interpretation

outputs/ and outputs/figures (ED8).

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
