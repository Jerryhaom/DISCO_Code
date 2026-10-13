#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Independently check saved scores, analytic covariance updates, and two runs.
options(stringsAsFactors=FALSE)
base <- .disco_path("results","01_method_definition/corrected_simulation/20261007_signed_half")
a <- readRDS(file.path(base,'pass_1/numeric_objects.rds'))
b <- readRDS(file.path(base,'pass_2/numeric_objects.rds'))
refs <- readRDS(file.path(base,'pass_1/reference_data.rds'))
stopifnot(identical(a,b))
analytic_entropy <- function(x,ref) {
  n <- nrow(ref); p <- ncol(ref); s <- cov(ref); delta <- x-colMeans(ref)
  # Independent determinant-lemma form of the unbiased covariance update.
  q <- drop(crossprod(delta,solve(s,delta)))
  .5*(p*log((n-1)/n)+log1p(n*q/((n+1)*(n-1))))
}
tests <- do.call(rbind,lapply(grep('^S[1-4]_',names(a),value=TRUE),function(nm) {
  z <- a[[nm]]; ref <- as.matrix(refs[[nm]])
  if(grepl('^S[34]_',nm)) {
    pred <- vapply(seq_len(nrow(ref)),function(i) analytic_entropy(ref[i,],ref),numeric(1))
    error <- max(abs(pred-z$Entropy))
    stopifnot(error<1e-12,nrow(z)==1000L)
  } else error <- NA_real_
  data.frame(panel=nm,n_rows=nrow(z),all_finite=all(is.finite(as.matrix(z))),
    analytic_entropy_max_error=error,two_runs_identical=TRUE)
}))
old <- readRDS(.disco_path("results","09_figures_tables/outputs/strict_extensions_simulation_geometry/20260917T104833_original_S1_S5/pass_1/numeric_objects.rds"))
original_noise <- old$S5_0.7
if(!is.data.frame(original_noise)) stop('Historical S5 object schema changed; verify manually')
cols <- intersect(c('noise','Mahalanobis','DISCO','Entropy','KL'),names(original_noise))
stopifnot(length(cols)==5L,nrow(original_noise)==11000L)
original_noise$Entropy <- .5*original_noise$Entropy
noise_error <- max(abs(as.matrix(a$S5[,cols])-as.matrix(original_noise[,cols])))
stopifnot(noise_error<1e-12)
write.csv(tests,file.path(base,'independent_verification.csv'),row.names=FALSE)
writeLines(c('PASS: pass_1 and pass_2 numeric objects are exactly identical.',
  paste('Independent determinant-lemma entropy maximum error:',max(tests$analytic_entropy_max_error,na.rm=TRUE)),
  paste('S5 maximum score change vs historical source, except Entropy factor 0.5:',noise_error),
  'S5: 1000 underlying subjects x 11 noise SD levels; rho is descriptive, no independent-row P shown.',
  'This verifies this signed entropy implementation, not the competing abs entropy definition or theoretical proof.'),
  file.path(base,'verification.txt'))
print(tests,row.names=FALSE)
cat('S5 score error:',noise_error,'\nCOMPLETE: independent verification.\n')
