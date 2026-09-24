# 09_figures_tables

**Role:** Historical figure replay / export tools.

## Recommended entries

Run from the repository root after configuring data and output roots. Where entries represent separate branches, only run the branch needed for your analysis.

```sh
Rscript run.R 09_figures_tables/replay_original_blocks.R
```

## Inputs

Bundled replay_blocks/*.R removes the external source-file requirement. All actual data and saved intermediates remain external. Read provenance.csv, script_input_calls.csv and column_usage.csv. `DISCO_REPLAY_STRICT=1` requires one tag and enables runtime input hashing. Example: `Rscript run.R 09_figures_tables/replay_original_blocks.R Fig2a`.

## Outputs and interpretation

outputs/original_block_replay. Execution status is not visual approval. Python render_replay.py renders PDFs; audit_render.py compares a supplied manuscript PDF; reorganize_to_modules.py targets dated historical run layouts; run_audit.py needs DISCO_DATA_ROOT/manifests/figure_table_manifest.csv with relative final_file paths, status and identifying columns.

The repository-wide [input call inventory](../script_input_calls.csv) gives source locations and expressions; [column_usage.csv](../column_usage.csv) lists statically observed column accesses (not a complete machine-validated schema). All scripts and their preferred/alternative/historical roles are listed in [script_manifest.csv](../script_manifest.csv). Full controlled-data execution is not part of the synthetic validation.
