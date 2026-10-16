# DISCO_Code — article analysis and figure code

For the latest disease log(DM) and matched-protein DISCO/clock analyses, start with [submission_20261013/README.md](submission_20261013/README.md). The [October 8 analysis integration](analysis_20261008/README.txt) provides additional sensitivity analyses and aggregate replay. Earlier modules in [current/README.md](current/README.md) and dated directories retain their original scope; see [RELEASE_ENTRIES.md](RELEASE_ENTRIES.md) before selecting an entry.

The 6 October approved Table 1 standardization/imputation, UKB IPAQ, organ cause-specific risk-set correction, CIF sensitivity, NMR correction and new 70:30 organ-clock entries are in [approved/README.txt](approved/README.txt). These dated branches retain explicit provisional endpoint status. They do not make every historical entry a final manuscript pipeline.

Code accompanying **[Human aging reflects increases in entropy across organ networks](https://www.medrxiv.org/content/10.1101/2025.05.06.25327108v2)**, Hao M et al., medRxiv, version 2 (26 April 2026). DOI: [10.1101/2025.05.06.25327108](https://doi.org/10.1101/2025.05.06.25327108). This is a preprint. Earlier scripts may retain terminology from the version 1 title, *Distance of covariance (DISCO), a novel measure of network homeostatic dysregulation, reveals organ system interconnections underlying mortality and disease risk*.

**[R package: DISCO4HD](https://github.com/Jerryhaom/DISCO4HD)** · **[Online calculator](https://jerryhaom.github.io/DISCO-Calculator/)** · **[10-biomarker example](10_biomarker_calculator/README.md)**

## What the article investigates

The study uses Distance of Covariance (DISCO) to measure how an individual's biomarker profile perturbs the correlation structure of a reference population. It investigates aging-related loss of coordination across biological systems using clinical biomarkers, proteomics, metabolomics and microbiome data. The principal cohorts are UK Biobank, NHANES, CHARLS, CLHLS and RuLAS.

The analyses examine associations with mortality, frailty and incident chronic disease, compare DISCO with other aging measures, and study organ/system scores derived from circulating proteins. The article reports that organ scores often predict outcomes beyond diseases of the corresponding organ, with more connected organ scores showing stronger outcome associations. These observational findings motivate a systemic interpretation of aging; the code does not establish causal entropy transfer or a clinical diagnostic threshold. See the linked article for results and limitations.

## How the code maps to the paper

| Scientific question | Main code modules |
|---|---|
| Definition, simulation and sensitivity of DISCO | `01_method_definition` |
| Clinical mortality prediction and aging-measure comparisons | `02_primary_clinical`, `04_external_validation`, `08_sensitivity_analyses` |
| Proteomic and other modality analyses; incident disease | `03_multimodal_UKB`, `07_RuLAS_omics` |
| Organ/system specificity and network relationships | `05_organ_network` |
| Frailty and functional outcomes | `06_frailty_function` |
| Figures, source-block replay and table export | `09_figures_tables` |
| Practical frozen-reference clinical calculator | `10_biomarker_calculator` (companion implementation, not a replacement for the manuscript-wide panels) |

This repository contains analysis, simulation and figure replay scripts. **Restricted participant data and previously computed participant-level objects are not distributed.** Synthetic examples work without cohort data. Some manuscript scripts require approved inputs and saved intermediate objects; a passing demo is not a full reproduction of every manuscript result. Historical and preferred implementations are identified in `script_manifest.csv`.

## Start here

Install the validated **DISCO4HD 0.3.2** commit below (an R source compiler is required for Rcpp), or install the bundled DISCO4HD source with `R CMD INSTALL`. Then clone this analysis repository:

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/DISCO4HD@08db67ce4999ccd3f7e2e7d855f33cc851305ce9")
```

```sh
git clone https://github.com/Jerryhaom/DISCO_Code.git
cd DISCO_Code
```

Run the examples:

```sh
Rscript run.R 01_method_definition/run_synthetic_demo.R
Rscript run.R 01_method_definition/run_method_smoke_test.R
Rscript run.R 10_biomarker_calculator/run_example.R
```

The first checks the package and writes a synthetic score table; the second checks an independent simulation implementation. The third scores fictional clinical profiles with each selectable cohort reference and demonstrates UKB-based z-scores. None reads private participant data. Outputs default to `results/<module>/outputs/`.

For controlled inputs, set two independent roots before running the chosen script:

```sh
export DISCO_DATA_ROOT=/path/to/approved-inputs
export DISCO_OUTPUT_ROOT=/path/to/new-results
Rscript run.R 03_multimodal_UKB/01_prepare_prot_data.R
Rscript run.R 03_multimodal_UKB/02_regenerate_prot_disco.R
Rscript run.R 03_multimodal_UKB/03_prot_analysis.R
```

`DISCO_DATA_ROOT` defaults to `data/` inside this repository. `DISCO_OUTPUT_ROOT` defaults to `results/`. No fallback reads an author's home directory. Relative environment paths are resolved against the initial invocation directory. Use a new output directory; rerunning a script may replace outputs of that script. Do not point it to archival results unless replacement is intended.

Optional overrides: `DISCO_PRS_FILE`, `DISCO_RULAS_WAVE2_FILE`, `DISCO_NHANES4_FILE`. Their defaults and required columns are in [INPUTS.md](INPUTS.md). Direct `Rscript /path/to/module/script.R` invocation is also supported.

## Modules and execution status

| Module | Purpose | Entry guide |
|---|---|---|
| 01_method_definition | Core checks, simulations, R2/R4 sensitivity | [guide](01_method_definition/README.md) |
| 02_primary_clinical | NHANES / UKB saved-object analyses and plots | [guide](02_primary_clinical/README.md) |
| 03_multimodal_UKB | Proteomics regeneration, modalities, incident disease | [guide](03_multimodal_UKB/README.md) |
| 04_external_validation | CHARLS / CLHLS / RuLAS / HRS saved-object validation | [guide](04_external_validation/README.md) |
| 05_organ_network | Organ outcomes, networks, expanded covariates | [guide](05_organ_network/README.md) |
| 06_frailty_function | Saved frailty / function objects and figures | [guide](06_frailty_function/README.md) |
| 07_RuLAS_omics | Omics reanalysis and saved-object figures | [guide](07_RuLAS_omics/README.md) |
| 08_sensitivity_analyses | R01–R06 revision analyses | [guide](08_sensitivity_analyses/README.md) |
| 09_figures_tables | Bundled source-block replay and export utilities | [guide](09_figures_tables/README.md) |
| 10_biomarker_calculator | Frozen 10/8/9-marker reference examples | [guide](10_biomarker_calculator/README.md) |

[script_manifest.csv](script_manifest.csv) assigns a role and release status to each entry. The historical preferred labels identify an older source release; they do not approve a script for the current manuscript. Use the dated approved modules for the corrections listed there. Historical alternatives remain available for provenance and are not all independent final pipelines. R04 meta-analysis now reads the single approved Table 1 CSV and checks its model, transform, SD population and endpoint snapshot contract; the approved Table 1 entry separately recomputes the cohort fits.

## Environment

R >= 4.5 is the tested analysis environment. Install only dependencies needed for your module:

```sh
Rscript install_dependencies.R --check 03_multimodal_UKB
Rscript install_dependencies.R --install 03_multimodal_UKB
```

[dependencies.csv](dependencies.csv) records packages referenced by each script, including bundled replay blocks. `environment_observed.csv` records the local versions available at release validation; it is not a reconstructed historical environment or an all-platform lockfile. Optional Bioconductor dependencies use BiocManager. BioAge appears only as an optional guard against a historical stub and is never installed by this helper. The figure Python utilities need Python 3 with Pillow and the Poppler command-line tools (pdfinfo, pdftoppm, pdftotext) (`python3 -m pip install -r requirements.txt`); core R execution does not need them.

Do not substitute a local BioAge stub for historical calculations. Existing saved DM values remain historical inputs. The proteomics regeneration script contains its own documented formula and computes standard Mahalanobis distance; the synthetic package demo explicitly calls DISCO4HD.

## Figure replay

`09_figures_tables/replay_blocks/` contains the actual expressions previously read from external R scripts. Calls now read whole bundled blocks rather than fixed line numbers in missing legacy files. Most adaptations are path substitutions. The NHANES comparison block also corrects the displayed label from DunedinPACE to DunedinPoAm; its numeric inputs and geometry are unchanged. `provenance.csv` records original file/line ranges and current portable hashes; [RELEASE_CHANGES.md](RELEASE_CHANGES.md) documents further corrections. Restricted data and some previously computed outputs are still necessary; see the module guide. Statistical execution and visual matching remain separate checks.

## Release scope

See [MIGRATION.md](MIGRATION.md) for scientific changes in DISCO4HD, preserved historical analyses and current verification limits. The public release includes portability cleanup and the new calculator example; it does not include a full rerun of every manuscript analysis. Validation scope is recorded in [VALIDATION.md](VALIDATION.md). The analysis code is MIT licensed; this does not grant rights to any external cohort dataset.


## Citation and related software

Please cite the manuscript above when using the method and record the code commit, DISCO4HD version, biomarker panel, transformations and selected reference in your analysis. For historical reproduction, consult [MIGRATION.md](MIGRATION.md) before replacing saved scores with freshly calculated values. The calculator's concise clinical panel should not be described as a statistically established optimal biomarker set.

## Submission update, 13 October 2026

The latest Fig. 3c log-DM and matched-protein mutual-adjustment analyses are in [submission_20261013](submission_20261013/README.md). Their entry point, input requirements and source hashes are documented there.
