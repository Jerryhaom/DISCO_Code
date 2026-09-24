# 01_method_definition

**Role:** Simulations / independent method checks.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 01_method_definition/run_synthetic_demo.R
Rscript run.R 01_method_definition/run_method_smoke_test.R
```

## Inputs

No input data for the two quick-start scripts. `run_all.R` additionally calls simulation_lock and submission_simulation; the latter reads `Simulation/simdata/` including `Simu.Res.rds`. Sensitivity variants may need prior outputs or UKB inputs; consult script_input_calls.csv.

## Outputs and interpretation

outputs/synthetic_demo, outputs/method_smoke_test_summary.csv; each simulation writes its named subdirectory.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
