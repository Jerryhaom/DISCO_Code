# Current release review: 16 October 2026

The clean source review restored six reviewed aggregate/contract files previously
excluded by the broad `data/` ignore rule. `verify_public_data.py` now checks the
47 aggregate and metadata inputs by hash; CI also requires them to be tracked.
The metadata, meta-analysis, joint-contrast and launcher checks pass locally.
Five no-private-input analysis modes and the dated Table 1/ED4e aggregate entries
ran successfully. The latest DM and matched-protein input checks found all 8 and
6 required local files; the current R model sources retain their provenance hashes.
This review checked their saved coefficients and confidence intervals without
refitting participant-level models.

The separate R package passed `R CMD check --no-manual` with 99 expectations.
Python 1.0.2 passed 11 tests with its compiled C++ backend on macOS, including R
comparison and real-input guards. Its pure NumPy wheel also passed, with the
5 optional native-only tests explicitly skipped. Python cross-platform execution
is not established by the R package CI. The calculator passed 15 tests and its
embedded evidence check. New local changes are not yet a hosted release.

## Earlier validation records

# Local checklist-action candidate validation (2026-10-06)

This is a separate local candidate, not the earlier upload tree or a remote release. The approved entries genuinely re-executed Table 1 standardization/imputation, UKB IPAQ, CIF, corrected organ cause-specific risk sets, and NMR/core downstream calculations. A comparison-only replay of the newly trained 70:30 clock scores re-executed all 162 C-index models and TRAIN HR/PH; the 18 clocks were not trained again in this software check. Independent aggregate checks match the dated analysis outputs. The full training entry has checked inputs but still requires controlled endpoint intermediates and locally obtained third-party sources.

DISCO4HD 0.3.1 was freshly built and installed into a separate local library. Its 99 expectations and clean-session synthetic demo passed with no skipped tests; the local R CMD check returned `Status: OK`. The first default installation failed because this host's R references absent gfortran libraries. A host-only `FLIBS =` Makevars override was used outside the released package; the package's standard compiler flags were not changed. OpenMP was unavailable and the serial C++ plus PSOCK R backends were tested. These are macOS local checks, not three-platform CI or proof that all hosts have configured toolchains.

The calculator's 15 tests, model-metadata tests, strict Meta input checks, explicit-policy day-endpoint synthetic tests, same-f1 joint-contrast rejection tests, R parsing and source/controlled-file checks were rerun on this candidate. The calculator's raw-DM, fixed-panel, adjusted-C-index evidence is explicitly separate from Table 1's log(DM), core-adjusted HR and crude C-index.

Corrected organ PH diagnostics detect departures; 18 sparse expanded Alzheimer's fits cannot support interpretable PH diagnostics. Survival endpoints remain provisionally retained, not newly finalized. Fig4c's historical network chain, some simulation inputs, missing raw-data ETLs and third-party licensing text remain open. No push, tag, deployment or external CI was performed. Local test success is not full manuscript scientific approval.

# Earlier candidate validation (2026-10-06)

- 197 analysis R files parse; all script-manifest entries resolve. SHA-256 verification passes for the 63 portable replay blocks. Changes beyond input-path substitutions are recorded in RELEASE_CHANGES.md and provenance.csv.
- Nine model-metadata and invalid-output rejection checks pass. Ten controlled-input Cox fits were executed locally: six NHANES/UKB fits and four CLHLS fits. This is not a full rerun of every analysis entry point.
- The web calculator's 15 tests pass. The DISCO4HD 0.3.1 source is unchanged from the checked 2026-10-05 candidate. Its 99 expectations were rerun on 6 October against the locally compiled 0.3.1 package and all passed, including PSOCK after permission for local worker sockets. No fresh package installation, compilation or R CMD check is claimed; those records remain dated 5 October. The default system library still contains 0.1.0 and must not be substituted for the candidate.
- Four additional saved-object audit entry points were executed. CHARLS/RuLAS validation produced six Cox fits, frailty produced two, and the corrected RuLAS omics audit produced eighteen historical/comparator fits. CLHLS validation correctly retains an invalid-status blocker; frailty 2015 lacks endpoints in its saved score object. None is a full raw follow-up validation.
- The RuLAS aggregate audit rebuilds the actual Plot_subsystems2.R score (mean of eight clinical subsystems, Met and Meta). Its crude C-index is 0.712944. The saved `var` column is not that score: 36 missingness-mask differences were found, with maximum overlapping-score difference 5.965618. Both are explicitly distinguished; corrected NMR results are not inferred from this historical replay.
- Organ metadata counts were rerun and independently checked: 51,922 rows, 42 columns, 2,180,724 total cells and 2,076,880 finite numeric cells. A missing ID column is reported as unavailable, not inferred from repeated ages.
- The analysis workflow now requires DISCO4HD v0.3.1. The matching package tag must be published before the workflow can use it; remote three-platform CI has not run.
- Restricted participant inputs are excluded from the public source archive. The archive manifest and per-file checksums have been verified locally. No GitHub upload is claimed.
- Full raw-data ETL, disease-date construction, historical network/simulation inputs, and every manuscript figure are not validated by these checks. Scientific submission approval remains separate from syntax, tests and archive integrity.

The dated sections below are historical records, with their original counts and scope.

# Historical release validation (2026-09-24)

This release adds the clinical-reference calculator example and an article overview to the prior portable analysis release.

- DISCO4HD **0.3.0**: R CMD check OK, 83 passed expectations, no failures/warnings/skips, including PSOCK workers.
- **182 R files parse**, all script-manifest entries resolve, and all 63 bundled replay blocks exist. Five Python files pass AST parsing.
- Both the independent method smoke test and the installed-package synthetic demo pass. The fictional clinical example runs with all five frozen references and verifies a web-calculator score.
- No restricted participant data is added to this repository. Credential-pattern scanning found no matching tokens or private keys in the release text files.
- No full rerun of every controlled-cohort model, or visual reproduction of every manuscript figure, is claimed. Historical validation notes below retain their original scope and counts.

# Local release validation — 2026-09-24

- R 4.5.3 on macOS ARM64.
- DISCO4HD 0.2.0: R CMD check --no-manual, Status: OK; no errors, warnings or notes. Tests: 28 assertions passed, none skipped, including two PSOCK workers.
- 181 R files parse, including 63 bundled replay blocks and the new configuration/verification helpers. Python files parse with ast.
- Package synthetic demo and independent method smoke test passed after copying the repository to a separate directory whose path contains spaces, with no cohort input directory supplied. Explicit relative input/output roots also passed.
- R04 meta-analysis produced five files byte-identical to the original script in separate output directories.
- All 63 bundled replay blocks match the original file slices except documented path substitutions; source and portable hashes are recorded in provenance.csv.
- Full clinical cohorts, alternative sensitivity analyses and visual figure fidelity have not been revalidated. The GitHub workflows have not yet run.
- Local compiler limitation: this host's R Makeconf references a missing Fortran toolchain directory. Local compilation used a temporary user Makevars containing `FLIBS =` with R's existing BLAS/LAPACK. The released Makevars retains standard `$(FLIBS)` linkage for correctly configured toolchains. This override is not a general cross-platform installation instruction. OpenMP was unavailable in this local build; C++ serial and R process-parallel backends were tested.

Full logs and the original-source backup are retained in the local release workspace, outside the public source repository. Existing paper results were not overwritten.
