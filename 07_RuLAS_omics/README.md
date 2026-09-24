# 07_RuLAS_omics

**Role:** Reanalysis from saved omics features; historical plotting.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 07_RuLAS_omics/run_saved_object_audit.R
Rscript run.R 07_RuLAS_omics/01_age75_unified.R
Rscript run.R 07_RuLAS_omics/02_discrete_reference.R
```

## Inputs

`Result6/Rugao.PCs.met.rds`, `Rugao.selected.microbe.rds`, `Rugao.Met.surv.rds`, `Rugao.Meta.surv.rds`, `Rugao.Omics.surv.rds`, `Rugao.Phe.DISCO.rds`, `Analysis.2019.x.xlsx`, `Group_info.2019.csv`. Components have pre-existing feature/ID layouts; do not infer interchangeability of survival objects.

## Outputs and interpretation

outputs/07_age75_analysis. run_figures.R is historical saved-object replay and does not automatically consume the age75 reanalysis.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
