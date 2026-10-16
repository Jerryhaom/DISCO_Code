# Release entry guide

| Analysis | Entry and scope |
|---|---|
| Fig. 3c disease HRs | `submission_20261013/run.py dm`: log(DM), full scored-population SD, age, sex, BMI, smoking, drinking, education, race/ethnicity and physical activity. C-indices remain crude score-only results. Requires authorized inputs. |
| Matched-protein DISCO/clock mortality | `submission_20261013/run.py matched`: age, sex, BMI, smoking, drinking, education and race/ethnicity; separate and mutually adjusted models. Uses frozen training residual SDs and the training analysis population. Requires authorized inputs. |
| CHARLS interval-censoring sensitivities | `analysis_20261008/run.py charls`; aggregate tables distinguish the retained primary Table 1 analysis from interval-censoring sensitivity analyses. Requires authorized inputs. |
| ED4 and Fig. 4d aggregate plots | `analysis_20261008/run.py ed4-plot` and `fig4c-plot`; replay bundled aggregates without refitting cohorts. |
| ED12 outcome-power analysis | `analysis_20261008/run.py s11b`; source-era name retained. Requires authorized inputs. |
| Formula and PH-table checks | `analysis_20261008/run.py theory-identities` and `ph-check`; synthetic identities or frozen aggregate checks. |
| Earlier Table 1 snapshots | `approved/table1` (October 6) and `submission_20261007/release_validation_20261007` (October 7). Preserve their dated contracts; these scripts do not recompute later sensitivity analyses. |

HRS is a benchmark cohort, not a cohort in the five-cohort Table 1 HR meta-analysis.
Do not substitute earlier raw-DM disease results for the current log(DM) entry.
The dated manifests retain source provenance; the archive manifest identifies the
actual release contents. Restricted records and intermediate participant objects
must be supplied by an authorized user and are never downloaded by these runners.

Before distributing an analysis checkout, run `python verify_public_data.py --require-tracked`
and `Rscript verify_release.R`. For a source archive without Git metadata, omit
`--require-tracked`. The public data manifest lists reviewed aggregates explicitly;
adding a data file requires a separate content review and a new manifest entry.
