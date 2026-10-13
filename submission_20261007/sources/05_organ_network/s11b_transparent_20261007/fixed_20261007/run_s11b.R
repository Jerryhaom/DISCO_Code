#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Rebuild the documented transparent weight-power analysis from aggregate inputs.
# The historical plot and the previous mis-keyed long table are preserved.
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(ggplot2); library(readxl)})
args <- commandArgs(TRUE)
src <- if(length(args)>0) args[1] else .disco_path("legacy","Result3/Organ")
out <- if(length(args)>1) args[2] else .disco_path("results","05_organ_network/s11b_transparent_20261007/fixed_20261007")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
inputs <- file.path(src,c('Tissue.pcorr.noverlap2.rds','Tissue.cs1.rds','UKB.tissue.mc.rds','Disease_cindex.xlsx'))
stopifnot(all(file.exists(inputs)))
mat <- readRDS(inputs[1]); cs <- readRDS(inputs[2]); u <- readRDS(inputs[3])
org <- colnames(u)[colnames(u) %in% as.character(unique(cs$tissue))]
stopifnot(length(org)==18L, !anyDuplicated(org), identical(dim(mat),c(18L,18L)),
          all(is.finite(mat)), all(mat>=0), isTRUE(all.equal(mat,t(mat))))
# Matrix order follows the original documented 18-organ cl subset. Check named
# matrices explicitly; do not relabel a named matrix by positional assumption.
if(!is.null(rownames(mat))) stopifnot(identical(rownames(mat),org))
if(!is.null(colnames(mat))) stopifnot(identical(colnames(mat),org))
diag(mat) <- 0
get_cs <- function(variable) {
  d <- cs[cs$variable==variable,]
  stopifnot(!anyDuplicated(d$tissue),setequal(as.character(d$tissue),org),all(is.finite(d$value)))
  d$value[match(org,as.character(d$tissue))]
}
d <- suppressMessages(read_xlsx(inputs[4]))
blocks <- seq.int(2,nrow(d),by=9)
stopifnot(length(blocks)==18L, all(blocks+8L<=nrow(d)))
dn <- as.character(d[[1]][blocks])
dv <- vapply(blocks,function(i) {
  x <- as.numeric(d[[8]][(i+1L):(i+8L)])
  stopifnot(length(x)==8L,all(is.finite(x)))
  mean(x)
},numeric(1))
stopifnot(!anyDuplicated(dn),setequal(dn,org))
targets <- list('Correlation with age'=get_cs('Corr'),
                'C-index for mortality'=get_cs('C-index'),
                'C-index for disease'=dv[match(org,dn)])
powers <- seq(.5,4,by=.5)
# One row per explicit key, never flatten a matrix separately from its labels.
tab <- do.call(rbind,lapply(names(targets),function(n) {
  data.frame(power=powers,
    r=vapply(powers,function(p) cor(rowSums(mat^p),targets[[n]],method='pearson'),numeric(1)),
    outcome=n,n_organs=length(org))
}))
key <- function(x) paste(x$power,x$outcome,sep='|')
stopifnot(nrow(tab)==24L,!anyDuplicated(key(tab)),all(is.finite(tab$r)))
# Independent arithmetic and traversal: power-major order, manual centered dot
# products rather than cor(), then compare by the semantic key.
independent <- do.call(rbind,lapply(powers,function(p) {
  degree <- apply(mat^p,1,sum)
  dc <- degree-mean(degree)
  do.call(rbind,lapply(names(targets),function(n) {
    yc <- targets[[n]]-mean(targets[[n]])
    data.frame(power=p,outcome=n,
      independent_r=sum(dc*yc)/sqrt(sum(dc^2)*sum(yc^2)))
  }))
}))
qa <- tab
qa$independent_r <- independent$independent_r[match(key(tab),key(independent))]
qa$abs_error <- abs(qa$r-qa$independent_r)
stopifnot(!anyNA(qa$independent_r),max(qa$abs_error)<1e-12)
old <- file.path(dirname(out),'S11b_weighted_power_transparent_data.csv')
if(file.exists(old)) {
  legacy <- read.csv(old)
  stopifnot(!anyDuplicated(key(legacy)),setequal(key(legacy),key(tab)))
  qa$legacy_r <- legacy$r[match(key(tab),key(legacy))]
  qa$legacy_abs_error <- abs(qa$r-qa$legacy_r)
}
write.csv(tab,file.path(out,'S11b_weighted_power_data.csv'),row.names=FALSE)
write.csv(qa,file.path(out,'S11b_key_verification.csv'),row.names=FALSE)
organ_values <- data.frame(organ=org,do.call(cbind,targets),check.names=FALSE)
write.csv(organ_values,file.path(out,'S11b_organ_outcomes.csv'),row.names=FALSE)
write.csv(data.frame(input=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_hashes.csv'),row.names=FALSE)
tab$outcome <- factor(tab$outcome,levels=names(targets))
p <- ggplot(tab,aes(power,r,colour=outcome,group=outcome))+
  geom_line(linewidth=.7)+geom_point(size=2)+
  scale_x_continuous(breaks=powers)+
  scale_colour_manual(values=c('#E64B35FF','#4DBBD5FF','#00A087FF'),
    labels=c('Age','Mortality C-index','Disease C-index'),name=NULL)+
  labs(x='Weighted power',y='Pearson correlation coefficient')+
  theme_classic(base_size=12)+
  theme(legend.position='bottom',legend.text=element_text(size=10))
ggsave(file.path(out,'S11b_weighted_power_CORRECTED.pdf'),p,width=6,height=4)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('PASS: all 24 semantic keys agree with independent arithmetic to <1e-12.',
  paste('Maximum arithmetic error:',format(max(qa$abs_error),digits=16)),
  if('legacy_abs_error' %in% names(qa)) paste('Legacy mis-keyed rows:',sum(qa$legacy_abs_error>1e-12),'of',nrow(qa)) else '',
  'This is a transparent new analysis of saved organ aggregates, not a historical full-upstream reproduction.',
  'Disease outcome = mean of the eight stored disease C-indices; mortality/age inputs remain historical saved aggregates.'),
  file.path(out,'verification.txt'))
print(tab,row.names=FALSE)
cat('COMPLETE: S11b keyed correction; original files preserved.\n')
