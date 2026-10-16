#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Local aggregate-only entry: figures, rasterization, numeric QA and hashes.
root<-.disco_path("results","03_multimodal_UKB/targeted_endpoint_diagnostics_20261007/figure_candidates")
python<-Sys.getenv("DISCO_PYTHON",unset=Sys.which("python3"))
poppler<-Sys.getenv("DISCO_PDFTOPPM",unset=Sys.which("pdftoppm"))
stopifnot(file.exists(python),file.exists(poppler))
stopifnot(system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(root,'code/produce_fig3c_and_cif.R')))==0)
for(n in c('Fig3c_raw_unit_expanded_BMI_NEAR_SOURCE','Fig3c_raw_unit_crude_NEAR_SOURCE','ClinicalProt_CIF_sensitivity'))
  stopifnot(system2(poppler,c('-r','200','-png',shQuote(file.path(root,'output/pdf',paste0(n,'.pdf'))),
                            shQuote(file.path(root,'output/png',n))))==0)
qa<-system2(python,shQuote(file.path(root,'code/independent_numeric_qa.py')),stdout=TRUE)
stopifnot(is.null(attr(qa,'status')))
writeLines(qa,file.path(root,'qa/independent_numeric_QA.json'))
stopifnot(system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(root,'code/finalize_manifest.R')))==0)
cat('Numerically verified candidates rebuilt. Human visual review after changed inputs/layout is still required.\n')
