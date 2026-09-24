# Input contract

Data are never downloaded by the analysis launcher. Supply authorized local inputs. Original folder labels Result1–Result6 identify input collections under DISCO_DATA_ROOT; they no longer identify author-specific absolute paths.

| Collection / override | Expected input |
|---|---|
| Result1 | NHANES / UKB clinical saved DISCO and survival frames; mc=log-DISCO, hd=historical DM, en=entropy where used. NHANES plotting uses Age/id; UKB commonly uses age/ID. |
| Result2 | Proteomics / combined DISCO, selection model and protein annotation / network objects. Saved scores are not all produced by this release. |
| Result3 | Organ score matrices, tissue mappings, networks and sensitivity results. |
| Result4 | CHARLS, CLHLS, RuLAS saved DISCO / survival objects. |
| Result5 | Frailty / ADL / MMSE analysis objects. |
| Result6 | RuLAS omics feature and survival objects; module 07 lists filenames. |
| UKB_Data | Mortality, clinical/proteomic features, covariates, disease and PRS inputs. Mortality schema is id, Age, Gender, time, status, cause. |
| Simulation/simdata | Historical simulation CSVs with X-prefixed feature columns and survival fields plus Simu.Res.rds. Not required by the synthetic demo. |
| DISCO_PRS_FILE | Defaults to UKB_Data/UKB_PRS_full.csv. eid plus the 28 `262xx-0.0` fields explicitly listed at the R05 fread call. |
| UKB_Data/dat_disease_28_rebuilt.rds | eid with paired time_<disease>/status_<disease> fields for the diseases listed in rebuilt R05; its baseline-exclusion and noncase-censoring conventions are retained. |
| DISCO_RULAS_WAVE2_FILE | Defaults to RuLAS/wave2.xlsx; Sheet2 with R06 RuLAS covariate fields as used by 06_crosscohort_covariates.R. |
| DISCO_NHANES4_FILE | Defaults to NHANES/NHANES4.rda; must load NHANES4 with the R06 covariates. This replaces the hidden dependency on the GOLDBioAge source directory. |
| Manuscript / HRS | External summary CSVs and HRS saved objects for historical replay; see bundled blocks. |

Use module guides for ordering. `script_input_calls.csv` contains every statically found R reader call (including computed path expressions). `column_usage.csv` is an index of explicit `$` accesses to aid schema review, not an assertion that all referenced objects are input data frames. Dynamic columns, matrix layouts and nested R objects must follow the script's scientific contract; the release does not fabricate or auto-convert them.

Generated intermediate files live under DISCO_OUTPUT_ROOT, preserving module/outputs subpaths. If you replay a figure requiring an already computed output, either run its producer or place an approved copy at that documented output location. The launcher does not search the old Results_final tree or mix old results into a new run.

Status is usually 0/1 but follow each cohort script's event definition. Time units must match the original model. Never change age cutoffs, censoring rules, feature ordering or merge keys merely to make a file load.
