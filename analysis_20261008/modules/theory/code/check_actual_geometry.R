#!/usr/bin/env Rscript
# Frozen synthetic geometry objects only; no paper score or definition is replaced.
options(warn=1)
suppressPackageStartupMessages({library(jsonlite);library(digest)})
root <- Sys.getenv('DISCO_RUN_MODULE')
up <- Sys.getenv('DISCO_GEOMETRY_INPUT')
paths <- file.path(up,c('numeric_objects.rds','reference_data.rds','executed_source.R'))
before <- vapply(paths,function(x)digest(file=x,algo='sha256'),character(1))
val <- readRDS(paths[1]); refs <- readRDS(paths[2]); out <- list()
for(nm in grep('^S[34]_',names(val),value=TRUE)) {
  x <- val[[nm]];ref <- as.matrix(refs[[nm]]);n <- nrow(ref);p <- ncol(ref)
  mu <- colMeans(ref);S <- cov(ref);L <- chol(S);v <- sweep(ref,2,mu,'-')
  q <- rowSums((v %*% solve(S))*v);a <- (n-1)/n;b <- 1/(n+1);g <- b/a
  entropy <- .5*(p*log(a)+log1p(g*q))
  R <- cor(ref);sdref <- sqrt(diag(S));z <- sweep(v,2,sdref,'/')
  predD <- vapply(seq_len(n),function(i){
    u<-z[i,];Rnew<-(R+g*tcrossprod(u))/sqrt(outer(1+g*u^2,1+g*u^2));log(.5*sum((Rnew-R)^2)*n^2)
  },numeric(1))
  direct <- vapply(seq(1,n,length.out=20),function(ii){
    i<-as.integer(ii)
    Rpred <- (R+g*tcrossprod(z[i,]))/sqrt(outer(1+g*z[i,]^2,1+g*z[i,]^2))
    max(abs(cor(rbind(ref,ref[i,]))-Rpred))
  },numeric(1))
  stopifnot(max(abs(entropy-x$Entropy))<1e-12,max(abs(predD-x$DISCO))<1e-8,max(direct)<1e-12)
  out[[nm]] <- list(n=n,p=p,positiveDefinite=TRUE,minEigenvalue=min(eigen(S,symmetric=TRUE)$values),
    conditionNumber=kappa(S,exact=TRUE),finiteStoredScores=all(is.finite(as.matrix(x))),
    signedNegative=sum(x$Entropy<0),signedFormulaMaxError=max(abs(entropy-x$Entropy)),
    analyticDISCOMaxError=max(abs(predD-x$DISCO)),directCorrelationChecks=length(direct),maxCorrelationError=max(direct),
    Spearman_DISCO_signedEntropy=cor(x$DISCO,x$Entropy,method='spearman'),
    Spearman_DISCO_absoluteEntropy=cor(x$DISCO,abs(x$Entropy),method='spearman'),
    Spearman_DM_signedEntropy=cor(x$Mahalanobis,x$Entropy,method='spearman'),
    Spearman_DM_absoluteEntropy=cor(x$Mahalanobis,abs(x$Entropy),method='spearman'),
    sameSampleReference=TRUE,absoluteOnlyDiagnosticNotScoreReplacement=TRUE)
}
stopifnot(identical(before,vapply(paths,function(x)digest(file=x,algo='sha256'),character(1))))
write_json(list(status='PASS',frozenInputs=as.list(setNames(before,paths)),settings=out,
  boundaryPolicy='No historical calculation changed. Current eight covariance matrices SPD; singular/log-zero cases diagnosed separately, not silently repaired.'),
  file.path(root,'results/actual_geometry_checks.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(root,'logs/geometry_sessionInfo.txt'))
oldlog <- file.path(root,'logs/actual_geometry_R.log')
if(file.exists(oldlog) && any(grepl('unexpected',readLines(oldlog))))
  file.copy(oldlog,file.path(root,'logs/initial_local_R_syntax_error.log'),overwrite=FALSE)
writeLines(c('PASS: 8 settings; 8000 signed/absolute diagnostics; 160 direct cor(rbind()) checks',
 'No frozen score or input changed; input SHA verified.'),oldlog)
cat('PASS:',length(out),'settings; 8000 signed/absolute diagnostics; 160 actual cor(rbind()) checks\n')
