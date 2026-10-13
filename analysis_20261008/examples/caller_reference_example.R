#!/usr/bin/env Rscript
# Entirely fictional seeded example: no cohort, participant or clinical calibration.
suppressPackageStartupMessages(library(DISCO4HD))
set.seed(20261008)
markers <- c('CRP','GLU','HDL')
reference_measurements <- as.data.frame(exp(matrix(rnorm(90),30,3)))
names(reference_measurements) <- markers
target <- as.data.frame(exp(matrix(rnorm(171),57,3)))
names(target) <- markers
transformed <- reference_measurements
transformed$CRP <- log1p(transformed$CRP)
scaling_mean <- colMeans(transformed)
scaling_sd <- vapply(transformed,sd,0)
X <- scale(as.matrix(transformed),center=scaling_mean,scale=scaling_sd)
n <- nrow(X)
weight <- matrix(1,3,3);diag(weight) <- 0;weight <- weight/sum(weight)
ref <- structure(list(cohort='FICTIONAL_EXAMPLE',markers=markers,
 ref_mean=scaling_mean,ref_sd=scaling_sd,mu=colMeans(X),
 ss=crossprod(sweep(X,2,colMeans(X),'-')),Rref=cor(X),weight=weight,
 S=cov(X)+diag(1e-8,3),n_ref=n,units=c('mg/L','mmol/L','mmol/L'),
 source='Seeded fictional lognormal example; uniform illustrative weights'),
 class='disco10_reference')
score <- cal_disco10(target,reference=ref)
y <- target;y$CRP <- log1p(y$CRP)
y <- scale(as.matrix(y),center=scaling_mean,scale=scaling_sd)
literal <- vapply(seq_len(nrow(y)),function(i)
 log(n^2*sum(weight*(cor(X)-cor(rbind(X,y[i,])))^2)),0)
dm <- sqrt(mahalanobis(y,ref$mu,ref$S))
stopifnot(max(abs(literal-score$DISCO))<1e-10,max(abs(dm-score$DM))<1e-10)
cat('CALLER_REFERENCE_PASS fictional_queries=57 markers=3 max_DISCO_error=',
 max(abs(literal-score$DISCO)),' max_DM_error=',max(abs(dm-score$DM)),'\n')
