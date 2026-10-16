#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(jsonlite);library(digest)})
args<-commandArgs(TRUE); module<-if(length(args))normalizePath(args[1])else normalizePath('01_method_definition/ED4fh_full_recompute_20261008')
root<-Sys.getenv('DISCO_INPUT_ROOT');path<-file.path(root,'UKB_Data/UKB.Prot.clean.rds');oldpath<-file.path(root,'Results_final/01_method_definition/outputs/sensitivity/R43_exp1_disco_dm.rds')
sha<-function(p)digest(file=p,algo='sha256');hash<-sha(path);d<-readRDS(path);X<-d$X
stopifnot(is.matrix(X),all(is.finite(X)),ncol(X)>=2920L);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20260918);ri<-sample(nrow(X),500)
rt<-cor(X);pg<-c(50,100,200,400,500,1000,2000,2920);res<-list();detail<-list();sels<-list();tic<-proc.time()[3]
for(p in pg) {
  ce<-cn<-numeric(5)
  for(rep in 1:5) {
    set.seed(1000+p+rep);sel<-sample(ncol(X),p);xr<-X[ri,sel,drop=FALSE];rh<-cor(xr);rr<-rt[sel,sel,drop=FALSE]
    ce[rep]<-sqrt(mean((rh[lower.tri(rh)]-rr[lower.tri(rr)])^2))
    if(p<500){ev<-eigen(cov(xr),symmetric=TRUE,only.values=TRUE)$values;cn[rep]<-max(ev)/max(min(ev),1e-12)}else cn[rep]<-NA_real_
    detail[[length(detail)+1]]<-data.frame(p=p,rep=rep,seed=1000+p+rep,correlation_RMS_error=ce[rep],condition_number=cn[rep]);sels[[paste(p,rep,sep='_')]]<-sel
  }
  res[[length(res)+1]]<-data.frame(p=p,cor_err=mean(ce),cor_err_se=sd(ce),condnum=mean(cn,na.rm=TRUE));cat('ED4f p',p,'done\n');flush.console()
}
res<-do.call(rbind,res);detail<-do.call(rbind,detail);saveRDS(res,file.path(module,'results/ED4f_summary.rds'));write.csv(res,file.path(module,'tables/ED4f_summary.csv'),row.names=FALSE);write.csv(detail,file.path(module,'tables/ED4f_replicates_40.csv'),row.names=FALSE)
saveRDS(list(reference_row=ri,reference_ids=if(!is.null(d$eid))d$eid[ri]else NULL,feature_subsets=sels,feature_names=colnames(X),RNGkind=RNGkind()),file.path(module,'input/private/ED4f_sampling.rds'));Sys.chmod(file.path(module,'input/private/ED4f_sampling.rds'),'0600')
old<-readRDS(oldpath);comparison<-merge(res,old,by='p',suffixes=c('_new','_old'));write.csv(comparison,file.path(module,'results/ED4f_history_comparison.csv'),row.names=FALSE)
delta<-max(abs(comparison$cor_err_new-comparison$cor_err_old),abs(comparison$cor_err_se_new-comparison$cor_err_se_old),abs(comparison$condnum_new-comparison$condnum_old),na.rm=TRUE)
qa<-list(status='PASS',scope='Complete 8 dimensions x 5 independent feature samples using original declared seeds',n_input=nrow(X),p_input=ncol(X),reference_n=500,replicates=5,dimension_count=8,cor_err_se_definition='Sample SD across five feature subsets; not standard error',condition_number_defined_only_p_less_than_500=TRUE,historical_max_abs_difference=delta,input_sha256=hash,input_unchanged=identical(hash,sha(path)),elapsed_sec=unname(proc.time()[3]-tic),R_version=R.version.string)
stopifnot(qa$input_unchanged);write_json(qa,file.path(module,'results/ED4f_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17);print(qa)
