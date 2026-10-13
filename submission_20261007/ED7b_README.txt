ED7b saved-result replay

Run from any working directory:
Rscript DISCO_Code/submission_20261007/plot_ED7b_f1_interaction.R /path/to/output

R dependencies: forestploter, gtable and base/recommended grid. If dependencies are in a user library, set DISCO_R_LIBRARY to that library directory. No participant data, private-data path or original DISCO working directory is required by this replay.

The frozen aggregate CSV contains 28 outcomes in each of the clinical and proteomic modalities. The plotted coefficient is exp(beta for zD:zPRS) from the fully adjusted interaction Cox model, not DISCO's main effect and not the (1,1) versus (0,0) Joint contrast. Baseline disease exclusions, actual blood date, death-first ties and the locked 2022-12-31 administrative cutoff were used in the precomputed raw_dates/death_first results. Adjustment: age, sex, smoking, drinking, education and ethnicity; single-sex cancer analyses omit sex. This branch does not include BMI or physical activity. DISCO was standardized across the full modality cohort; each PRS was standardized across its sex-eligible modality cohort before endpoint and complete-covariate exclusions. Per-row covariates are retained in the CSV. The source_result column holds logical attribution names, not required external input paths; only workstation path text was changed, not any statistic.

The proteomic ovarian-cancer row uses the previously saved stricter numerical tolerance for the same model and risk set; it is not a new fitted model chosen by significance. The remaining 55 rows are the existing saved results unchanged.

The visual block is the confirmed original study forest code, with the longer Interaction HR header set to 14 pt; data text and other headers remain 15 pt. The device retains 20 x 12 inches, Helvetica, original point/CI colour and limits. Source algorithms, selection and endpoint models are not refitted by this entry. The old ED7a/full assembled figure is not updated by this command. The primary dataset-provider coverage and original raw upstream certification remain separate tasks.

This local candidate package is not a public-release approval; aggregate/reference redistribution and third-party dependencies still require the study authors' review. No upload is performed by this entry.
