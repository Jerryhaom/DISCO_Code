DISCO 2026-10-08 analysis source integration — local review candidate

Run from this directory with Python 3 and an Rscript installation when required.
Use a new --output directory. Existing directories are rejected. No runner uploads
data, contacts a server, changes a system library, or modifies the supplied inputs.

No-private-input modes:
  python run.py theory-identities --output /path/to/new/theory-check
  python run.py ed4-plot --output /path/to/new/ed4-aggregate-render
  python run.py fig4c-plot --output /path/to/new/fig4c-aggregate-render
  python run.py ph-check --output /path/to/new/frozen-PH-audit
  python run.py ed1-contract --output /path/to/new/ED1-public-source-audit

These consume explicit synthetic calculations or aggregate statistics only.
Plot modes replay supplied complete aggregates; they do not refit cohorts.
PH mode audits already-fitted aggregate HR/PH tables and basis rescaling; it adds
no PH tests and fits no Cox models. Original scope statements accompany each module.

Controlled-input modes: ed4f, ed4h, fig4c, s11b, ed9c, charls, theory-geometry.
ED1 original observation uses ed1-observe with an explicit simulation_source file
matching modules/source_data/original_source_contract.json. Its unchanged original
definitions3–101/function681–788, seed123 and full-logdet entropy are recorded.
Unrelated controlled export callbacks are skipped by an explicit ED1 scope guard.
The full original file and its R environment are caller inputs, not bundled data.
ed1-manifest-check verifies a caller-supplied ED1 synthetic/aggregate export manifest
and its relative output files; neither contract nor manifest checks run the generator.
Provide --config config.json using config.example.json and required_inputs.json.
The old folder layout under input_root/results_root is part of the contract;
it is explicit and must not be guessed. The four S11b endpoint files are named
separately; their original server paths were removed. The S11b original saved
LASSO definitions are retained, not retrained. Its full current target fits need
all four endpoint inputs and the same accepted ID/score bridge. CHARLS needs
actual source wave evidence and icenReg dependencies in the user's configured R
environment. It runs two separate processes for each 200-bootstrap branch.

Source adaptations bind paths and output locations only. The ED4 plotting source
reads byte-preserved aggregate CSVs instead of equivalent aggregate RDSs.
The source manifest links original source hashes and exact mechanical diffs.
No scoring definition, statistical model, seed, fold assignment, endpoint rule,
or covariance convention is replaced by this integration.

The current combined source archive includes the previously public DISCO4HD
aggregate reference registry and a generated synthetic NHANES example. Thus
named disco10_reference() and cal_disco10() references are available after an R
package installation. The registry contains cohort means, standard deviations,
correlations, weights and comparison summaries, not participant records.
Restricted participant tables, fitted participant-level objects and exact
event-time curves are not included or downloaded automatically.

examples/caller_reference_example.R demonstrates an entirely synthetic reference
and independent literal Pearson/DM agreement. After a local R package installation,
run Rscript examples/caller_reference_example.R. Historical current/math_geometry
requires an omitted historical synthetic RDS; use theory-identities for a
self-contained mathematical check. The source_manifest.json records the original
October 8 integration; the enclosing archive's SOURCE_MANIFEST.json records
current file hashes and subsequent edits.

Module names retain the labels at source creation: fig4c refers to the current
Fig. 4d network-degree panel, and s11b refers to the ED12 outcome-power analysis.
Use the latest panel/source mapping when assembling the manuscript.

Original source authorship, provenance and LICENSE files are retained. Scientific
acceptance is based on numeric/structural checks, not word filtering. This analysis source is a
local review candidate. CI results for the separate R package do not validate
restricted cohort analyses or the Python package on other platforms.
