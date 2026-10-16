#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(survival);library(Rcpp);library(jsonlite);library(digest)})
args <- commandArgs(TRUE); module <- if(length(args)) normalizePath(args[1]) else normalizePath('01_method_definition/ED4fh_full_recompute_20261008')
disco_root <- Sys.getenv('DISCO_INPUT_ROOT')
paths <- c(clinical=file.path(disco_root,'Result1/Clean.UKB.Clinical.rds'), metadata=file.path(disco_root,'Result1/UKB.Clinical.DISCO.rds'), historical_DISCO=file.path(disco_root,'Results_final/01_method_definition/outputs/sensitivity/R27_rep50_nested.rds'), historical_DM=file.path(disco_root,'Results_final/01_method_definition/outputs/sensitivity/R27_dm_nested.rds'))
sha <- function(p) digest(file=p,algo='sha256'); hashes <- vapply(paths,sha,'')
sourceCpp(file.path(module,'code/score_kernel.cpp'),cacheDir=file.path(module,'compile_cache'))
clin <- readRDS(paths['clinical']); meta <- readRDS(paths['metadata'])
stopifnot(!anyDuplicated(clin$sampleID), !anyDuplicated(meta$ID), identical(as.character(clin$sampleID),as.character(meta$ID)))
X <- scale(as.matrix(clin[,-1])); ids <- clin$sampleID; rm(clin)
age <- meta$age; sex <- meta$sex; tm <- meta$time; status <- meta$status
stopifnot(ncol(X)==58L,all(is.finite(X)),all(is.finite(age)),all(is.finite(tm)),all(status %in% 0:1),!anyNA(sex))
cc <- apply(X,2,function(x)cor(age,x,use='complete.obs')); wa <- abs(cc%o%cc);diag(wa)<-0;wa<-wa/sum(wa)
wu <- matrix(1,ncol(X),ncol(X));diag(wu)<-0;wu<-wu/sum(wu)
pool <- which(age<=40);stopifnot(length(pool)==5282L)
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20260914); evaluation <- sample(setdiff(seq_len(nrow(X)),pool),50000)
Xa<-X[evaluation,,drop=FALSE]; aa<-age[evaluation];sa<-sex[evaluation];ta<-tm[evaluation];ev<-status[evaluation]
cat('Historical endpoint nonpositive time: full input',sum(tm<=0),'evaluation',sum(ta<=0),'(no endpoint estimand silently changed)\n')
refsizes <- c(50,100,200,500,1000,2000,5282); nrep<-50L
set.seed(123);xr<-X[sample(pool,500),,drop=FALSE]; mu<-colMeans(xr); S<-cov(xr);rr<-cor(xr)
fast<-score_two_weights(Xa[1:12,,drop=FALSE],mu,S,rr,wa,wu,500L)
literal<-t(vapply(1:12,function(i){d<-(cor(rbind(xr,Xa[i,]))-rr)^2;c(log(500^2*sum(wa*d)),log(500^2*sum(wu*d)))},numeric(2)))
literal_error<-max(abs(fast-literal));stopifnot(literal_error<1e-9)
ds<-list();ms<-list(); sampling<-list(); warnings<-list(); score_hashes<-list(); tic<-proc.time()[3]
fit_ci<-function(formula,key) {
  ww<-character();f<-withCallingHandlers(coxph(formula),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
  if(length(ww))warnings[[length(warnings)+1L]]<<-data.frame(key=key,message=unique(ww))
  unname(summary(f)$concordance[1])
}
for(rep in seq_len(nrep)) {
  sampling[[rep]]<-list()
  for(nref in refsizes) {
    set.seed(123L+rep-1L);ri<-sample(pool,nref);sampling[[rep]][[as.character(nref)]]<-ids[ri]
    xr<-X[ri,,drop=FALSE];mu<-colMeans(xr);S<-cov(xr);Rref<-cor(xr)
    scores<-score_two_weights(Xa,mu,S,Rref,wa,wu,as.integer(nref));stopifnot(all(is.finite(scores)))
    score_hashes[[length(score_hashes)+1L]]<-data.frame(rep=rep,refsize=nref,score_sha256=digest(scores,algo='sha256'))
    for(w in seq_len(2)) {
      sc<-scores[,w]; accel<-residuals(lm(sc~aa));wname<-c('age-weighted','uniform')[w]
      vals<-c(fit_ci(Surv(ta,ev)~sc,paste(rep,nref,w,'crude')),fit_ci(Surv(ta,ev)~sc+aa+sa,paste(rep,nref,w,'age+sex')),fit_ci(Surv(ta,ev)~accel,paste(rep,nref,w,'accel')))
      ds[[length(ds)+1L]]<-data.frame(weight=wname,refsize=nref,model=c('crude','age+sex','accel'),rep=rep,ci=vals)
    }
    hd<-sqrt(mahalanobis(Xa,mu,S+diag(1e-8,ncol(Xa))))
    ms[[length(ms)+1L]]<-data.frame(model=c('crude','age+sex'),refsize=nref,rep=rep,ci=c(fit_ci(Surv(ta,ev)~hd,paste(rep,nref,'DM crude')),fit_ci(Surv(ta,ev)~hd+aa+sa,paste(rep,nref,'DM age+sex'))))
  }
  cat(sprintf('ED4h rep %02d/%02d elapsed %.1fs\n',rep,nrep,proc.time()[3]-tic));flush.console()
}
d<-do.call(rbind,ds);m<-do.call(rbind,ms);stopifnot(nrow(d)==2100,nrow(m)==700)
saveRDS(list(evaluation_ids=ids[evaluation],reference_ids=sampling,seeds=list(evaluation=20260914,reference=123:172),RNGkind=RNGkind()),file.path(module,'input/private/ED4h_sampling_IDs.rds'));Sys.chmod(file.path(module,'input/private/ED4h_sampling_IDs.rds'),'0600')
saveRDS(d,file.path(module,'results/ED4h_DISCO_2100.rds'));saveRDS(m,file.path(module,'results/ED4h_DM_700.rds'))
write.csv(d,file.path(module,'tables/ED4h_DISCO_2100.csv'),row.names=FALSE);write.csv(m,file.path(module,'tables/ED4h_DM_700.csv'),row.names=FALSE);write.csv(do.call(rbind,score_hashes),file.path(module,'results/ED4h_score_hashes.csv'),row.names=FALSE)
old_d<-readRDS(paths['historical_DISCO']);old_m<-readRDS(paths['historical_DM']);cmp_d<-merge(d,old_d,by=c('weight','refsize','model','rep'),suffixes=c('_new','_old'));cmp_m<-merge(m,old_m,by=c('refsize','model','rep'),suffixes=c('_new','_old'))
stopifnot(nrow(cmp_d)==2100,nrow(cmp_m)==700);cmp_d$delta<-cmp_d$ci_new-cmp_d$ci_old;cmp_m$delta<-cmp_m$ci_new-cmp_m$ci_old
write.csv(cmp_d,file.path(module,'results/ED4h_DISCO_history_comparison.csv'),row.names=FALSE);write.csv(cmp_m,file.path(module,'results/ED4h_DM_history_comparison.csv'),row.names=FALSE)
qa<-list(status='PASS',scope='Full 50-replicate original seeded statistics; historical endpoint estimand preserved',n_input=nrow(X),p=ncol(X),n_reference=length(pool),n_eval=length(evaluation),deaths=sum(ev),DISCO_rows=nrow(d),DM_rows=nrow(m),model_warning_count=length(warnings),literal_score_max_abs_error=literal_error,historical_DISCO_max_abs_C_difference=max(abs(cmp_d$delta)),historical_DM_max_abs_C_difference=max(abs(cmp_m$delta)),R_version=R.version.string,RNGkind=RNGkind(),elapsed_sec=unname(proc.time()[3]-tic),input_sha256=as.list(hashes),input_unchanged=identical(hashes,vapply(paths,sha,'')))
stopifnot(qa$input_unchanged);if(length(warnings))write.csv(do.call(rbind,warnings),file.path(module,'results/ED4h_model_warnings.csv'),row.names=FALSE)
write_json(qa,file.path(module,'results/ED4h_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17);print(qa)
