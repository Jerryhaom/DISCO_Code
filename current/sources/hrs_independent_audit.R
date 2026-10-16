#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(haven);library(survival)})
dir.create(file.path(Sys.getenv('DISCO_RUN_OUTPUT'),'numeric'),recursive=TRUE,showWarnings=FALSE)
base<-Sys.getenv('DISCO_RUN_OUTPUT')
p<-readRDS(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/04_external_validation/HRS_imputed_reference/20261001/numeric/prepared_INTERNAL.rds'))
o<-readRDS(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/04_external_validation/HRS_imputed_reference/20261001/median/numeric/models_INTERNAL.rds'));d<-o$data
bl<-as.data.frame(read_sas(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'HRS_Data/打包/vbs16v3.sas7bdat')))
id<-sprintf('%06d%03d',as.integer(bl$HHID),as.integer(bl$PN));j<-match(p$id,id)
raw<-as.matrix(bl[j,p$features]);raw[,'PGLUFF']<-raw[,'PGLUFF']/18.0182
raw[,'PCRP']<-log1p(raw[,'PCRP']/10)
logs<-c('PALT','PAST','PCMVGE','PCYSC','PDHEASE','PFERTN','PNTBNPE','PTGF');raw[,logs]<-log1p(raw[,logs])
stopifnot(identical(is.na(raw),is.na(as.matrix(p$X[,p$features]))),max(abs(raw-as.matrix(p$X[,p$features])),na.rm=TRUE)<1e-12)
train<-!p$id%in%p$cid
med<-apply(raw[train,],2,median,na.rm=TRUE)
stopifnot(max(abs(med-o$medians[p$features]))<1e-12)
filled<-raw;for(k in seq_len(ncol(filled)))filled[is.na(filled[,k]),k]<-med[k]
ref<-p$drawage<60;refmat<-filled[ref,]
ac<-apply(filled[train,],2,function(v)cor(v,p$X$AGE[train]));w<-abs(ac%o%ac);diag(w)<-0;w<-w/sum(w)
# Literal append-and-cor calculation, on reversed reference order, including
# leave-one-out reference members and nonmembers. No optimized scoring helper.
probe<-unique(c(head(which(ref[p$eval_index]),8),head(which(!ref[p$eval_index]),8)))
errd<-errm<-numeric()
for(i in probe){rr<-refmat
 if(ref[p$eval_index[i]])rr<-rr[-match(p$eval_index[i],which(ref)),,drop=FALSE]
 rr<-rr[nrow(rr):1,];v<-filled[p$eval_index[i],]
 disco<-log(sum((cor(rr)-cor(rbind(rr,v)))^2*w)*nrow(rr)^2)
 z<-(v-colMeans(rr))/apply(rr,2,sd);dm<-log(drop(t(z)%*%solve(cor(rr),z)))
 errd<-c(errd,abs(disco-d$HRS35_DISCO[i]));errm<-c(errm,abs(dm-d$HRS35_DM[i]))
}
manual<-function(time,status,lp){num<-den<-0
 for(i in which(status==1)){ii<-which(time>time[i]|(time==time[i]&status==0));den<-den+length(ii);num<-num+sum(lp[i]>lp[ii])+.5*sum(lp[i]==lp[ii])}
 num/den}
cerr<-numeric()
for(m in c('NHANES18_DISCO','NHANES18_DM','HRS35_DISCO','HRS35_DM')){
 dd<-d[nrow(d):1,];dd$s<-dd[[m]]
 f<-coxph(Surv(time,status)~s+age,data=dd,ties='efron')
 c<-manual(dd$time,dd$status,predict(f,type='lp'))
 expected<-o$table$C[o$table$metric==m&o$table$model=='adj'];cerr<-c(cerr,abs(c-expected))
}
stopifnot(max(errd)<1e-7,max(errm)<1e-7,max(cerr)<1e-10)
old<-readRDS(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/04_external_validation/HRS_tracker_reanalysis/20261001/numeric/endpoints_INTERNAL.rds'))
oi<-match(old$data$person_id,d$person_id)
stopifnot(length(oi)==3923,!anyNA(oi),max(abs(d$time[oi]-old$data$time))<1e-12,all(d$status[oi]==old$data$status))
out<-list(primary_N=nrow(d),deaths=sum(d$status),reference_N=sum(ref),median_pool_N=sum(train),
 raw_source_transform_error=max(abs(raw-as.matrix(p$X[,p$features])),na.rm=TRUE),
 independent_literal_DISCO_max_error=max(errd),independent_DM_max_error=max(errm),manual_concordance_max_error=max(cerr),
 old_3923_endpoints_identical=TRUE,old_3923_NHANES_scores_identical=TRUE,
 no_MICE_results_used=TRUE,simple_imputation='Training-only per-marker median; no stochastic draws',
 sources=c(blood=file.path(Sys.getenv('DISCO_DATA_ROOT'), 'HRS_Data/打包/vbs16v3.sas7bdat'),
  original_plot=file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Result1/DNAm/HRS.cmp.v2.R')))
jsonlite::write_json(out,file.path(base,'numeric/independent_audit.json'),pretty=TRUE,auto_unbox=TRUE,digits=14)
writeLines(capture.output(sessionInfo()),file.path(base,'numeric/sessionInfo.txt'))
print(out)
