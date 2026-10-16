# Changes in the 6 October 2026 candidate

This candidate is an analysis-code update, not approval of the full submission.
Original sources and the previous portable release are retained separately.

- Saved-object Cox audit exports use `fit$n` for participants and `fit$nevent`
  for events. `nobs.coxph()` reports the number of events, not participants.
  Seven exports have been corrected. No fitted coefficients are changed.
- Audit entry points attach the packages required by their model expressions
  when invoked in a fresh R process.
- Cox result validation rejects empty output, missing metadata, invalid sample
  and event counts, non-finite results and invalid concordance/standard errors.
  An empty result table can no longer generate a successful audit status.
- The CLHLS audit no longer presents the first numeric column as a participant
  identifier when an ID column is absent. That exported ID is left missing;
  stable identity validation belongs to the original linkage audit.
- The NHANES comparison figure displays DunedinPoAm, consistent with its input
  field and the NHANES documentation. The scores, plot design and C-index inputs
  are unchanged. The corresponding portable source-block hash is updated.
- The combined clinical/protein entropy replay requires positive finite values
  from the plotted `pd$hd` column. The separate legacy HD file is not substituted
  for that column. The saved combined panel has 396 selected features; the
  upstream candidate panel has 425. This does not reproduce selection from raw
  cohort data.
- The analysis smoke workflow now includes the metadata regression test and
  requires the matching DISCO4HD v0.3.1 release. Publish and test the package
  first; the remote workflow has not been run as part of this local review.

The current organ, HRS and corrected simulation modules retain their documented
formulas and input boundaries. Restricted data, local numerical runs and
participant-level outputs are excluded from the public candidate.
# Follow-up audit entry checks on 2026-10-06

- RuLAS audit now reads the endpoint-bearing objects written by Plot_surv.R instead of looking for follow-up in score-only objects. The clinical join and age restriction follow that source and are explicitly marked historical, not validation of loss-to-follow-up or corrected NMR outputs.
- Organ-object finite-value counts now inspect numeric columns independently of character metadata, and total cells use rows times columns. Repeated ages are no longer mislabeled as duplicate participant records; ID duplication remains unavailable when no participant ID is supplied.
- The RuLAS historical aggregate audit executes Plot_subsystems2.R's actual ten-component mean. Its saved `var` does not match that formula, so it is retained only as a separate comparator. Missingness and numerical discrepancies are exported; this change does not replace corrected NMR analyses or modify any legacy input.
