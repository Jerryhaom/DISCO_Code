# 03_multimodal_UKB

**Role:** Regeneration from prepared biomarkers and saved selection model; outcome analysis.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 03_multimodal_UKB/01_prepare_prot_data.R
Rscript run.R 03_multimodal_UKB/02_regenerate_prot_disco.R
Rscript run.R 03_multimodal_UKB/03_prot_analysis.R
Rscript run.R 03_multimodal_UKB/05_disease_dual_model.R
Rscript run.R 03_multimodal_UKB/reproduce_ED7_KM_corrected.R
```

## Inputs

`UKB_Data/UKB.Mortality.csv`: id, Age, Gender, time, status, cause. `20240321_prot_instance0.RData` must define prot.dat with eid and numeric proteins. `Result2/Lasso.model.prot.rds` supplies the saved selection model. Disease analyses require their UKB disease objects, plus modality-specific saved scores. See INPUTS.md. Preparation does not retrain LASSO.

## Outputs and interpretation

outputs/03_prot_regeneration, outputs/03_prot_analysis, outputs/03_disease_analysis, outputs/ED7_KM. run_fig2a.R uses the regenerated proteomics object; run_fig2b.R reads historical combined-modality data.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
