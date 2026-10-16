#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE)))
base<-file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final');.libPaths(c(file.path(base,'09_figures_tables/Rlibs'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(Rcpp);library(digest)})
setDTthreads(1);options(warn=1);Sys.umask('0077');sink(file.path(out,'logs/03_DISCO_scores.log'),split=TRUE)
d<-readRDS(file.path(out,'PRIVATE/prepared.rds'));tr<-d$train;te<-d$test;ids<-d$meta$eid
fs<-fromJSON(file.path(out,'numeric/feature_sets.json'));need<-sort(unique(c(unlist(fs$DISCO),unlist(fs$Oh),unlist(fs$Goeminne))))
X<-d$X[,need,drop=FALSE];med<-apply(X[tr,,drop=FALSE],2,median,na.rm=TRUE);stopifnot(all(is.finite(med)))
for(j in seq_len(ncol(X)))X[is.na(X[,j]),j]<-med[j]
stopifnot(all(is.finite(X)));ref<-tr[d$meta$age[tr]<=40];stopifnot(length(ref)>=3,!length(intersect(ids[ref],ids[te])))
rho<-as.vector(cor(X[tr,,drop=FALSE],d$meta$age[tr]));names(rho)<-need;mu<-colMeans(X[ref,,drop=FALSE]);ss<-apply(X[ref,,drop=FALSE],2,sd)
stopifnot(all(is.finite(rho)),all(is.finite(ss)),all(ss>0));
sourceCpp(file.path(Sys.getenv('DISCO_APPROVED_CODE_ROOT'),'resources/exact_score.cpp'),cacheDir=file.path(out,'PRIVATE/cpp_cache'))
probes<-unique(c(1L,ref[c(1,2,length(ref))],te[c(1,2,length(te))],tr[c(100,500,1000)]));checks<-list();weights<-list()
calc<-function(pp,key){
 if(length(pp)<2)return(NULL)
 stopifnot(all(pp%in%need),!anyDuplicated(pp));cc<-cor(X[ref,pp,drop=FALSE]);w<-abs(outer(rho[pp],rho[pp]));diag(w)<-0;stopifnot(sum(w)>0);w<-w/sum(w)
 s<-organ_score(X[,pp,drop=FALSE],mu[pp],ss[pp],cc,w,length(ref),FALSE)
 direct<-vapply(probes,function(i)log(length(ref)^2*sum(w*(cor(rbind(X[ref,pp,drop=FALSE],X[i,pp,drop=FALSE]))-cc)^2)),0.)
 err<-max(abs(s[probes]-direct));stopifnot(all(is.finite(s)),err<1e-9)
 checks[[key]]<<-data.frame(score=key,p=length(pp),direct_probes=length(probes),max_abs_error=err,train_reference_N=length(ref),reference_test_overlap=0)
 weights[[key]]<<-list(proteins=pp,correlation=cc,weights=w)
 cat('Exact train-reference DISCO',key,'P',length(pp),'direct max error',err,'\n');flush.console();s
}
scores<-lapply(names(fs$DISCO),function(o)calc(fs$DISCO[[o]],paste0(o,'_Original')));names(scores)<-names(fs$DISCO)
saveRDS(list(scores=do.call(cbind,scores),ids=ids,X=X,train=tr,test=te,reference=ref,meta=d$meta,medians=med,agecor=rho,mu=mu,sd=ss,weights=weights),file.path(out,'PRIVATE/DISCO_original.rds'),compress=FALSE)
fwrite(rbindlist(checks),file.path(out,'qa/DISCO_original_direct_checks.csv'));fwrite(data.frame(protein=need,training_median=unname(med),training_age_r=unname(rho[need]),reference_mean=unname(mu[need]),reference_sd=unname(ss[need]),train_N=length(tr),reference_N=length(ref)),file.path(out,'numeric/DISCO_training_parameters.csv'))
write_json(list(all_train_only=TRUE,train_N=length(tr),test_N=length(te),reference_N=length(ref),reference_test_overlap=0,score_kernel='Original exact correlation-difference kernel; same no-self-exclusion definition for training references',kernel_sha256=digest(file.path(Sys.getenv('DISCO_APPROVED_CODE_ROOT'),'resources/exact_score.cpp'),algo='sha256',file=TRUE),all_direct_checks_pass=TRUE),file.path(out,'qa/DISCO_original_complete.json'),pretty=TRUE,auto_unbox=TRUE)
pf<-list.files(file.path(out,'PRIVATE'),full.names=TRUE,recursive=FALSE);Sys.chmod(pf[!file.info(pf)$isdir],'0600');Sys.chmod(pf[file.info(pf)$isdir],'0700');sink()
