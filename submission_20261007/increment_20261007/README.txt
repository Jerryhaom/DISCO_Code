New standalone entrypoints: ED7a and complete survival simulation

This increment is the recommended entry for new users. Historical simulation
and sensitivity (SS) outputs in earlier folders are retained for traceability,
not recommended as current corrected results. Earlier ZIP members are immutable;
this supplement neither rewrites their bytes nor silently relabels old results.

Run from any working directory with the required R dependencies installed:
  Rscript --vanilla PATH_TO_INCREMENT/run.R ed7a EMPTY_OUTPUT_DIRECTORY
  Rscript --vanilla PATH_TO_INCREMENT/run.R smoke EMPTY_OUTPUT_DIRECTORY
  Rscript --vanilla PATH_TO_INCREMENT/run.R reference_figures EMPTY_OUTPUT_DIRECTORY
  Rscript --vanilla PATH_TO_INCREMENT/run.R generate_full EMPTY_OUTPUT_DIRECTORY
  Rscript --vanilla PATH_TO_INCREMENT/run.R full EMPTY_OUTPUT_DIRECTORY
Without a mode the default is smoke. Existing outputs are not overwritten.
Simulation modes require no historical dataset or original DISCO installation.
The ED7a mode requires the local curve input described below. DISCO_INCREMENT_RLIB optionally selects a local package library.
install_dependencies.R installs official CRAN sources only to an explicitly
chosen local library. It does not modify a global R installation. All required
tested versions are in package_versions.csv and the CRAN lock description.

ED7a
ED7a redraw requires an authorized local input/ED7a_event_time_curves.csv.
The exact event-time table is not distributed in this public branch.
input/ED7a_C_annotations.csv provides the8 disease annotations and16 score-only
Harrell C values. The age axis is attained age; annotations retain the original
definition, time since blood baseline, NOT an attained-age Cox C.
The curves are1−KM, not competing-risk CIF. Prior death is censored at death.
The baseline date/cutoff/endpoint rules were verified in restricted-data work;
this public entry faithfully redraws the aggregates and does not revalidate raw
UKB clinical records. No private fit data/x/y object is bundled.
source/ED7a_visual_original.R is the actual plotting-only source segment,
not a script chosen by its filename. Font, facets, axes, colours, sizes and
device dimensions are retained by the portable producer.

Simulation
The original configuration/functions/score/model/plotting segments are preserved
in source/, with original paths represented only as repository-relative
attribution and SHA256. They are article-author source under the inherited MIT
notice, not copied third-party package internals. No third-party package source
or binary is embedded; dependency packages retain their own licences.
smoke generates four NEW datasets, Standard p10 rho0.1, and fits16 NEW original
models. It draws a smoke-test figure only. It does NOT reproduce all96 datasets
or full S6/S7. Comparison summaries are used after fitting, not as score/time
or event inputs. generate_full creates all96 NEW datasets but does not fit.
full creates all96 NEW datasets, fits384 models and draws S6/S7 and transparent
ED5; allow several minutes. reference_figures draws all384 S6/S7 and96 ED5 cells
from the bundled verified model summaries. It is full DRAWING replay, not a fresh
96-dataset/384-fit test. Every output records this scope explicitly.

Definitions
Each dataset has1,000 observations; p10/50/100; rho0.1/0.3/0.5/0.7; noise SD0.2/1.
RNG L'Ecuyer-CMRG/Inversion/Rejection matches future.seed=TRUE with original
independent seeds. Predictor scale(sum(X²))+Gaussian noise; coefficient1.
Parameter lambda0.1; Gompertz gamma0.1; Weibull gamma1; exponential no gamma;
parametric administrative censoring at10. Weibull gamma1 equals exponential.
Cox sim.survdata T1000, splineFALSE, censor probability0.3, time divided by100.
Original Cox boundary-duration warnings are retained rather than changing the
generator. Bars are apparent in-sample Harrell C±ONE model SE, not95%CI or
Monte Carlo uncertainty. There is one dataset per scenario cell. S6's historical
N(0,0.2) title is retained; its0.2 parameter is SD, not variance.
The final historical ED5 drawing script is unrecovered. Its authorized transparent
two-metric design now has a complete original-generator/model replacement chain.

Provenance
The original local sources did not have an available Git commit identity.
source_manifest.json records that fact and exact original/extracted SHA256,
line ranges and attribution. It does not invent source commits or release tags.
Official dependency URLs, tested versions and known archive hashes are recorded.
Runtime provenance/sessionInfo come from actual execution, not static examples.
The public candidate has not been uploaded or assigned a release tag.
