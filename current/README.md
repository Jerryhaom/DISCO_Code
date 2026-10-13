# Current correction modules

This folder preserves the earlier correction branches. For the approved 6 October Table 1/IPAQ, organ risk set, CIF/NMR and genuinely retrained 70:30 organ-clock additions, use `../approved/README.txt`. The old 18,648/4,294 comparison and its one-protein Heart limitation are not the new 31,381/13,450 branch.

These modules preserve the executed calculation and plotting segments. Changes
are limited to input/output paths, bundling source dependencies, and removal of
inherited visual-review claims from the numerical HRS audit. Model formulas,
scaling, imputation, selected proteins and plot geometry are unchanged.

Run with Python 3 and R >= 4.5:

```sh
python3 current/run.py organ_d --data-root /path/to/approved/DISCO --output /path/to/new/organ_d
python3 current/run.py organ_abc --data-root /path/to/approved/DISCO --output /path/to/new/organ_abc
python3 current/run.py protein_sensitivity --data-root /path/to/approved/DISCO --output /path/to/new/protein_sensitivity
python3 current/run.py hrs_scores --data-root /path/to/approved/DISCO --output /path/to/new/hrs_scores
python3 current/run.py hrs_audit --data-root /path/to/approved/DISCO --output /path/to/new/hrs_audit
python3 current/run.py math_geometry --data-root /path/to/approved/DISCO --output /path/to/new/math_geometry
```

Every command requires a new output directory and fails on any failed R step.
Participant-level inputs and generated private objects must not be uploaded.
Source files and their SHA-256 identifiers are in `source_manifest.json`.

## Inputs and scope

| Module | Inputs below the approved DISCO root | Outputs and limits |
|---|---|---|
| organ_abc | Same frozen common prediction models, `20261002_local_rebuild/results/common_partial_correlations.csv`, `common_HR.csv`, `common_test_Cindex.csv`, `qa/final_audit.json`, and organ mapping | Independently checks all 18 correlations and 72 HR points, plots a-c and the uniform-age/sex sensitivity d, and reruns 18 joint-model PH diagnostics. Uses rebuilt frozen clock scores, not retraining. |
| organ_d | `Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild/PRIVATE/common_prediction_models.rds`, the `20260928/results/organ_mapping.csv`, and `20261002_rebuilt_figures/numeric` C-index tables | Refits mortality Cox models in the common training sample and evaluates test-sample C-index, paired differences and the d PDF. It does not retrain the protein clocks. |
| protein_sensitivity | Above model inputs; local Oh/Goeminne coefficient CSVs; `UKB_Data/20240321_prot_instance0.RData`; `Result2/Proteins/coding143.tsv`; `Result3/Organ/UKB.prot.imputed.rds` and `UKB.Prot.class.rds`; `Results_final/05_organ_network/corrected_organ/20260926/locked_design.rds` and `pass_1/organ36_exact_scores.rds`; original d numeric tables | Recomputes DISCO from the locally selected protein sets, refits survival models and generates the S9 artwork. Fixed-design sensitivity, not an independently trained full-pipeline validation. One-protein Oh Heart DISCO is undefined. |
| hrs_scores | `Results_final/04_external_validation/HRS_imputed_reference/20261001/numeric/prepared_INTERNAL.rds` and `missingness_inventory.csv`; `Results_final/04_external_validation/HRS_tracker_reanalysis/20261001/numeric/endpoints_INTERNAL.rds` | Recomputes marker-wise median imputation, external/internal DISCO and DM, apparent C-index, missingness sensitivity, and household-disjoint five-fold validation. The prepared input already contains reconciled outcomes; this module does not rebuild Tracker/NDI linkage. |
| hrs_audit | HRS prepared inputs and `median/numeric/models_INTERNAL.rds`; raw `HRS_Data/打包/vbs16v3.sas7bdat`; reconciled endpoint file | Independently checks raw transformations, medians, literal append-and-cor scores, Mahalanobis scores and manually counted concordance. Does not claim a new visual review. |
| math_geometry | Bundled original simulation blocks and synthetic baseline objects only; no cohort data | Recomputes corrected S2 inverse PCA and S3/S4 duplicate-row corrections with original plotting segments. These are corrected analyses, not restoration of an unidentified historical final figure. |

R dependencies: survival, data.table, jsonlite, Rcpp, digest, haven, MASS,
ggplot2, patchwork, ggpubr, grid and scales. Install only those needed for the
selected module. A C++ compiler is required for protein sensitivity.

## Publication boundary

Historical modules outside this folder are retained for provenance. The current
whole-manuscript release remains incomplete: disease definitions, fully adjusted
Table 1, some original plot chains, corrected NMR downstream release, and the
remaining simulation/organ-network modules have open checks. A successful run
here does not approve those other outputs. Use the submission checklist before
selecting final files.

The Original d comparison is DISCO + age versus protein age alone. The Accel
comparison uses training-calibrated age residuals. The corresponding unified
age/sex comparison must be reported alongside these deliberately different
models. Local clock rebuilding is not the authors' original fitted models.
