# Migration and reproducibility scope

This release replaces workstation-specific paths with explicit input/output roots. Original analysis formulae, reference cutoffs, fitting options, seeds and saved-score choices are preserved. The package update does not silently replace scores used in historical manuscript figures.

DISCO4HD 0.2.0 changes default modality decomposition to age-weighted DISCO. Select `weighting="unweighted"` explicitly to study the old structural sum; per-pair means now use correct pair counting. Missing inputs are rejected consistently. Any manuscript analysis that actually used the old decomposition needs a targeted scientific comparison before replacing its results.

The source snapshot before this change is kept outside the public repository in the release workspace. The original 0.1.0 tarball remains available locally. The release directory does not include old mypackage/GOLDBioAge packages, RStudio history or Python caches.

Recommended alternatives: use ED7 `reproduce_ED7_KM_corrected.R`; use R05 `*_rebuilt.R`; prefer `05_disease_dual_model.R` for dual-model disease estimates. These labels identify the corrected analysis family, not certification of every endpoint. Other sensitivity variants are retained as explicitly labelled alternatives in script_manifest.csv.

Clinical/external/frailty modules often start with saved DISCO objects. Access to raw cohort files alone is insufficient for those entries. Proteomics preparation additionally requires the selected LASSO model; its training is not recreated here. `run_all.R` in module 01 also requires the historical simulation dataset, and is not the data-free quick start.

R04 embeds rounded cohort summary estimates, so its output is a reanalysis of those inputs, not a regeneration of individual-level cohort estimates. Figure reconstruction/export utilities referencing dated strict runs remain historical utilities, not preferred starting points.

Validation covers package tests, synthetic examples, syntax, path portability and bundled source provenance. Full clinical outcomes, all alternative sensitivity runs and all figure layouts are not rerun by this release check. Cross-platform CI is provided but requires execution on GitHub; a checked-in workflow is not evidence that those platforms have passed.
