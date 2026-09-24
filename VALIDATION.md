# Current release validation (2026-09-24)

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
