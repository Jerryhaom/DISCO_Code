# October 7 analysis producers

This entry contains the executed day-level endpoint, three-modality disease,
organ-clock and technical simulation corrections. The earlier `approved/` and
`current/` entries are retained as historical pipelines. They must not be used
to label their older outputs as the October 7 results.

Required environment roots:

- `DISCO_PRIVATE_ROOT`: controlled UK Biobank inputs and private run directories.
- `DISCO_RESULTS_ROOT`: analysis workspace with the recorded input/output layout.
- `DISCO_LEGACY_ROOT`: original saved objects and scientific source scripts.
- `DISCO_UKB_SSH`: SSH destination for the explicitly selected private export steps.

The public archive contains source and an aggregate Table 1 input, not cohort
participants, diagnostic dates, fitted participant scores, or credential files.
Source paths are relative; each producer's original and portable hashes and
path-only adaptations are recorded in the accompanying release manifest.

`Rscript run_aggregate_figures.R` rebuilds the two Fig. 3c candidate variants
and two-page CIF sensitivity from the bundled aggregate CSV inputs, then runs
independent numeric checks. It needs R packages data.table, forestploter, gtable,
ggplot2, jsonlite and digest, Python with the QA dependencies, and pdftoppm on PATH.
Set `DISCO_PYTHON` or `DISCO_PDFTOPPM` when those executables are not on PATH.
This figure-only replay does not refit cohort Cox models or certify raw endpoints.

Execution is conditional on the stated input contracts. A source parsing check
does not reproduce a result. The updated endpoint policy starts at recorded
blood collection, uses 31 December 2022 as the administrative cutoff, and keeps
mortality and incident-disease observation rules distinct. That policy does not
establish that every regional provider observed events through the cutoff.

The organ-clock update reuses the frozen train-only 7:3 score models, selected
protein sets, training reference and split. It refits endpoint-dependent Cox
models and age calibration, evaluates concordance in test participants, and
reports PH diagnostics separately. The local algorithms are not the original
authors' pretrained scores or an independent external validation.

The corrected S11b calculation is a transparent analysis of saved aggregates.
It is not recovery of a historical upstream scoring run. The corrected S2–S5
branch uses a signed half-log-determinant change and preserves the distinction
from an absolute-valued theoretical quantity. The theoretical definition still
requires scientific approval.

No publication, repository push, version tag or deployment is performed by this
entry. Third-party sources remain controlled external dependencies unless their
licenses and redistribution terms have been checked.
