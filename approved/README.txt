Approved dated analysis entries, local release candidate (2026-10-06)

These entries supplement the historical analysis tree; they do not silently
replace saved scores or certify final manuscript endpoint definitions.
All source/data/output path adaptations are recorded in source_manifest.json.
No restricted participant input is included. Required source objects are
listed in run.py and are read only from an authorized DISCO data tree.
Generated local participant objects must never be included in a public upload.

Example commands from the DISCO_Code repository:

python3 approved/run.py table1 --data-root /path/to/authorized/DISCO --questionnaire /path/to/wave2.xlsx --output /path/to/new/table1 --allow-provisional
python3 approved/run.py ipaq --data-root /path/to/authorized/DISCO --output /path/to/new/ipaq --allow-provisional
python3 approved/run.py organ_cause --data-root /path/to/authorized/DISCO --output /path/to/new/organ_cause --allow-provisional
python3 approved/run.py organ_cause_ph --data-root /path/to/authorized/DISCO --output /path/to/new/organ_cause_ph --allow-provisional
python3 approved/run.py CIF --data-root /path/to/authorized/DISCO --output /path/to/new/CIF --allow-provisional
python3 approved/run.py nmr --data-root /path/to/authorized/DISCO --output /path/to/new/nmr --allow-provisional
python3 approved/run.py organ_70_30 --data-root /path/to/authorized/DISCO --output /path/to/new/organ70 --preflight

For an explicitly scoped comparison-only replay of the already trained new
70:30 objects (not training), use:
python3 approved/run.py organ_70_30_compare --data-root /path/to/authorized/DISCO --frozen-score-snapshot /path/to/20261006_strict_70_30 --output /path/to/new/comparison --allow-provisional

Replace --preflight with --allow-provisional to train all 18 clocks and run the
new 70:30 branch. This may take substantial time and memory. No old fitted
clock is substituted. Python requires numpy, pandas, scipy, scikit-learn,
matplotlib, joblib and threadpoolctl. R requires the packages used in the
corresponding source. DISCO_R_LIB optionally selects an existing local R
library without recording a workstation path in the released source.

Meta-analysis reads exactly one approved Table1 CSV and refuses changed
source hashes, mixed SD populations, transforms, adjustment formulas, invalid
counts, warned results and nonmatching DISCO/DM participants:

Rscript --vanilla 08_sensitivity_analyses/04_meta_analysis.R approved/table1/data/Table1_analysis_sample_SD.csv /path/to/new/meta --allow-provisional

The two metrics are pooled separately across five independent cohorts. Core
adjustment is age, sex, smoking, alcohol and education; race is additionally
included in NHANES and UKB. Those are similar but not identical adjustments.
HR is per final common cohort analysis-sample SD; DM is log-transformed.
Crude C-index is not the adjusted C-index and is not pooled with HR.

Outputs and scientific limits
- Table1: ten standardized associations and crude C-indices; approved simple
  covariate imputation is explicit. The original clinical score panels remain.
- IPAQ: six UKB associations comparing full Table1, IPAQ-complete common
  samples, and the same-sample additional physical-activity model.
- Organ cause: main age+sex plus separate legacy expanded branches; competing
  deaths remain in the risk set and are censored at their recorded death.
  Per-unit and eligible-sample per-SD HRs are separate. Warned expanded
  Alzheimer's HRs remain NA; no covariates or causes are deleted. The combined
  plot preserves the original disease half, whose raw source dates are not
  certified. The calculation now also runs term/GLOBAL cox.zph diagnostics
  (KM transform, separate BH correction across 18 organs within each branch,
  outcome and term). PH departures are reported, not labelled as passing.
  Diagnostics from the sparse warned models are noninterpretable and stay NA.
- CIF: four sensitivity PDFs plus aggregate AJ/KM curves on identical
  first-event risk sets; zero/ambiguous month ordering rules are explicit in
  PROTOCOL.txt. This does not approve a new final disease-date algorithm.
- NMR: original PCA, pp<.05 selection, reference, weights, plotting source
  blocks and downstream species selection, with the 23-PC dimension repair.
  The source computes participant work files locally. Historical DM remains
  a saved input where its original API cannot be restored; it is not replaced
  by a local BioAge stub. Other original censoring limitations remain.
- Organ 70:30: fixed split, new training-only clocks and DISCO/reference/gap
  fits; matched-protein sensitivity uses these newly selected proteins.
  The requested d display has unequal age information; fair comparisons use
  identical age+sex covariates. Conditional C-index uncertainty omits repeated
  training/imputation uncertainty. TRAIN PH failures must be disclosed.

The 70:30 entry still requires earlier prepared endpoint metadata, the fixed
organ feature sets/mapping, and locally obtained original Oh/Goeminne code.
These third-party files are not re-licensed or bundled under this repository's
MIT license. The local Goeminne README states Academic Non-Commercial terms,
but its LICENSE text is absent from this source snapshot. Obtain and verify
the original terms before redistribution. The Oh snapshot likewise has no
license text in this local copy. The controlled input tree is not a substitute
for recovering and documenting every upstream raw-data construction script.

Historical simulation replay remains conditional on its actual original RDS.
Do not substitute same-named files with different contents (in particular
NEW_R4_blockwise/R4_blockwise.rds). Fig4c's historical protein-network source
chain and all original raw-data ETLs are not recovered by these additions.
Successful numerical reproduction is not final scientific approval, and this
candidate has not been pushed, tagged or deployed.

Explicit day-level endpoint utilities are in endpoint_utilities/. They require
all cutoff, loss, baseline-day and disease/death tie decisions; none is
defaulted. Their 24 synthetic tests do not rebuild actual UKB outcomes.
Outputs containing IDs and true dates must remain PRIVATE.

PRS Joint HR export is fail-closed: f0 main effects cannot be multiplied by
an f1 interaction HR. Without one interaction-fit object containing all three
coefficients and their full covariance, the formatted table leaves Joint HR
and CI unavailable. Optional DISCO_R05_JOINT_MODELS_RDS explicitly supplies
named disease__version Cox f1 objects. The helper reports the joint contrast
for both fitted scores=1 versus both=0, and uses the full covariance. It does
not decide age+sex versus expanded adjustment or approve the raw endpoint.
