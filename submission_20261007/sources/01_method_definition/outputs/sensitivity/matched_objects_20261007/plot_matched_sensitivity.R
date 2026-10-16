#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Match the diagnostic plots to their actual saved calculations. This does not
# claim that the old person-level inputs or 10/7 endpoints have been regenerated.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
src <- .disco_path("results","01_method_definition/outputs/sensitivity")
out <- file.path(src,'matched_objects_20261007')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs <- file.path(src,c('R43_exp1_disco_dm.rds','R27_rep50_nested.rds','R27_dm_nested.rds'))
stopifnot(all(file.exists(inputs)))
f <- readRDS(inputs[1])
stopifnot(nrow(f)==8L,!anyDuplicated(f$p),all(is.finite(f$cor_err)),all(is.finite(f$cor_err_se)),
  all(is.finite(f$condnum[f$p<500])),all(is.na(f$condnum[f$p>=500])))
f$correlation_error_SD <- f$cor_err_se
f$condition_status <- ifelse(f$p<500,'estimated','undefined sample-covariance inverse: p >= n_ref=500')
pf1 <- ggplot(f,aes(p,cor_err))+
  geom_errorbar(aes(ymin=cor_err-correlation_error_SD,ymax=cor_err+correlation_error_SD),width=0)+
  geom_line(colour='#E64B35FF',linewidth=.7)+geom_point(colour='#E64B35FF',size=2)+
  scale_x_log10(breaks=f$p,labels=scales::label_number(big.mark=','))+
  labs(x='Number of features',y='Correlation RMS error')+
  theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1))
pf2 <- ggplot(f[f$p<500,],aes(p,condnum))+
  geom_line(colour='#4DBBD5FF',linewidth=.7)+geom_point(colour='#4DBBD5FF',size=2)+
  scale_x_log10(breaks=f$p,labels=scales::label_number(big.mark=','),limits=range(f$p))+
  scale_y_log10(labels=scales::label_number(big.mark=','))+
  labs(x='Number of features',y='Covariance condition number')+
  theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(out,'ED4f_correlation_covariance_DIAGNOSTIC.pdf'),
  pf1+pf2+plot_annotation(tag_levels='a'),width=9,height=3.8)
write.csv(f,file.path(out,'ED4f_matched_values.csv'),row.names=FALSE)
d <- readRDS(inputs[2]); m <- readRDS(inputs[3])
stopifnot(nrow(d)==2100L,nrow(m)==700L,all(is.finite(d$ci)),all(is.finite(m$ci)))
dk <- paste(d$weight,d$model,d$refsize,d$rep,sep='|')
mk <- paste(m$model,m$refsize,m$rep,sep='|')
stopifnot(!anyDuplicated(dk),!anyDuplicated(mk),all(d$rep %in% 1:50),all(m$rep %in% 1:50))
shared <- paste(d$model,d$refsize,d$rep,sep='|')
stopifnot(setequal(unique(shared[d$model!='accel']),mk))
ds <- d %>% group_by(weight,model,refsize) %>% summarise(mean=mean(ci),SD=sd(ci),n_rep=n(),.groups='drop')
ms <- m %>% group_by(model,refsize) %>% summarise(mean=mean(ci),SD=sd(ci),n_rep=n(),.groups='drop')
stopifnot(all(ds$n_rep==50L),all(ms$n_rep==50L))
models <- c('crude','age+sex','accel')
cols <- c('crude'='#E64B35FF','age+sex'='#4DBBD5FF','accel'='#00A087FF')
refsizes <- sort(unique(d$refsize))
ph1 <- ggplot(ds,aes(refsize,mean,colour=model,fill=model,linetype=weight,group=interaction(model,weight)))+
  geom_ribbon(aes(ymin=mean-SD,ymax=mean+SD),alpha=.10,colour=NA)+
  geom_line(linewidth=.7)+geom_point(size=1.4)+
  scale_x_log10(breaks=refsizes,labels=scales::label_number(big.mark=','))+
  scale_colour_manual(values=cols,breaks=models,name=NULL,labels=c('Crude','Age + sex','Age residual'))+
  scale_fill_manual(values=cols,guide='none')+
  scale_linetype_manual(values=c('age-weighted'='solid','uniform'='dashed'),name=NULL,
    labels=c('Age-weighted','Unweighted'))+
  labs(title='DISCO',x='Reference size',y='Mortality C-index')+
  theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1),
    plot.title=element_text(face='bold'),legend.position='bottom',legend.text=element_text(size=10))
ph2 <- ggplot(ms,aes(refsize,mean,colour=model,fill=model,group=model))+
  geom_ribbon(aes(ymin=mean-SD,ymax=mean+SD),alpha=.10,colour=NA)+
  geom_line(linewidth=.7)+geom_point(size=1.4)+
  scale_x_log10(breaks=refsizes,labels=scales::label_number(big.mark=','))+
  scale_colour_manual(values=cols,breaks=c('crude','age+sex'),labels=c('Crude','Age + sex'),name=NULL)+
  scale_fill_manual(values=cols,guide='none')+
  labs(title='DM',x='Reference size',y='Mortality C-index')+
  theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1),
    plot.title=element_text(face='bold'),legend.position='bottom',legend.text=element_text(size=10))
# Explicit breaks/labels prevent factor or alphabetical order from relabelling a
# model, the same key discipline used for the S11b correction.
ggsave(file.path(out,'ED4h_NESTED50_DISCO_DM.pdf'),ph1+ph2+plot_annotation(tag_levels='a'),width=10,height=4.6)
write.csv(ds,file.path(out,'ED4h_DISCO_nested50_summary.csv'),row.names=FALSE)
write.csv(ms,file.path(out,'ED4h_DM_nested50_summary.csv'),row.names=FALSE)
write.csv(data.frame(input=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_hashes.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('PASS: ED4f graph uses R43_exp1_disco_dm.rds, not unrelated 4-panel AUC/C-index object.',
  'ED4f bars = mean +/- SD over five feature subsets; cor_err_se was a historical misnomer.',
  'ED4f is a correlation-estimation/covariance-inversion diagnostic, not direct proof of individual-score or prediction superiority.',
  'For p >= n_ref=500, inverse sample covariance is unavailable by rank; missing condition numbers are not empirical measurements.',
  'PASS: ED4h reads nested 50-replicate DISCO and DM objects, with complete shared (model,refsize,rep) keys.',
  'ED4h ribbons = mean +/- SD across 50 nested reference subsets; they are not 95% prediction CIs.',
  'Saved nested-run endpoints/weights come from the historical 9/21 inputs; no new 10/7 endpoint or input-ID rerun is claimed.',
  'C-indices are apparent, from the fitted evaluation sample; age residual is computed in that sample.',
  'ED4i is outside this plotting entry. A separate 58-marker production entry and small smoke validation are in ../leave_one_out_new_20261007; full new-endpoint leave-one-out calculation is not certified here.'),file.path(out,'verification.txt'))
cat('COMPLETE: matched ED4f diagnostic and paired nested ED4h two-panel plot.\n')
