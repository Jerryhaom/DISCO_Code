Current aggregate entry and historical-source boundaries

run_current_table1_meta.R uses the actual bundled October 7 five-cohort CSV,
not hardcoded HR vectors and not the superseded October 6 approved-source hash.
It checks the exact input contract, shared cohort N/events, per-SD estimands,
log(DM), full core adjustment and warning metadata. It preserves every source
Table1 cell and computes HR meta-analysis plus leave-one-cohort-out estimates.
Crude C-index is reported from the CSV but is not pooled. DerSimonian-Laird with
normal Wald intervals is retained; only stale endpoint metadata are overridden.

The earlier 08_sensitivity_analyses/04_meta_analysis.R is numerically generic,
not itself hardcoded, but is locked to October 6 approved/table1/data/contract.json
and rejects the updated input. Do not weaken that historical contract or call
its old output the current Table1. The generic approved/meta_functions.R is
reused with explicit October 7 aggregate metadata. Other four raw endpoint
chains and UKB provider coverage are not certified by this public replay.

plot_ED4e_current.R reuses the actually identified Test_Distance.R drawing
layout with 11 accepted aggregate rows; official CovTools 0.5.7 SS is used.
Older package paths that retain the 0.5.4 SS value are historical and deprecated
for new SS comparisons. This public drawing entry does not run private NHANES
distance models. Default precision and +/- one-model-SE semantics are explicit.

Optional arguments to either script: INPUT_CSV NEW_OUTPUT_DIRECTORY.
Existing output directories are refused. The default input is entirely bundled;
there is no author's-workstation scientific input fallback. Required R libraries:
digest/jsonlite for Table1/Meta; ggplot2 for the official-SS aggregate drawing.

Provenance: source_paths_relative_to_analysis_workspace and SHA bindings are in
data/current_aggregate_contract.json. License remains the inherited MIT license.
These entries are aggregate replay, not all-paper raw-data reproduction.

CI references:
https://github.com/r-lib/actions/tree/v2/check-r-package
https://github.com/actions/setup-python
https://pybind11.readthedocs.io/en/stable/compiling.html
Prepared remote jobs have not been executed; local tests are separately logged.
