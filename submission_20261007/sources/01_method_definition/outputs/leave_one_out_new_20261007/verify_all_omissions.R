# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages(library(Rcpp))
base <- .disco_path("results","01_method_definition/outputs/leave_one_out_new_20261007")
sourceCpp(file.path(base,'pearson_leave_one.cpp'),cacheDir=file.path(base,'compile_cache'))
set.seed(20261007)
ref <- matrix(rnorm(300*58),300,58)
q <- matrix(rnorm(12*58),12,58)
v <- runif(58); w <- abs(v%o%v); diag(w) <- 0; w <- w/sum(w)
actual <- pearson_all_omissions(ref,q,w)
literal <- function(keep,subject) {
  r0 <- cor(ref[,keep,drop=FALSE])
  r1 <- cor(rbind(ref[,keep,drop=FALSE],q[subject,keep]))
  log(nrow(ref)^2*sum((r1-r0)^2*w[keep,keep,drop=FALSE])/sum(w[keep,keep,drop=FALSE]))
}
expected <- vapply(0:58,function(omit) {
  keep <- if(omit==0) 1:58 else setdiff(1:58,omit)
  vapply(1:12,function(i)literal(keep,i),numeric(1))
},numeric(12))
err <- max(abs(actual-expected)); stopifnot(err<1e-10)
constant <- ref; constant[,1] <- 0
checks <- c(constant_reference=inherits(try(pearson_all_omissions(constant,q,w),silent=TRUE),'try-error'),
  zero_weight=inherits(try(pearson_all_omissions(ref,q,w*0),silent=TRUE),'try-error'),
  nonfinite_query=inherits(try(pearson_all_omissions(ref,q*NA_real_,w),silent=TRUE),'try-error'))
stopifnot(all(checks))
writeLines(c(paste('PASS: all 59 columns x12 queries, maximum literal error',format(err,digits=17)),
  paste(names(checks),checks),
  'Synthetic verification only; not UKB evaluation or final ED4i graph.'),file.path(base,'synthetic_all_omissions_verification.txt'))
cat('PASS',format(err,digits=17),'all59 columns; invalid-domain checks PASS\n')
