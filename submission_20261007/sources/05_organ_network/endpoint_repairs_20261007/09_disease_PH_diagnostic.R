#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Refit only primary 288 disease models for PH diagnostics; verify identity.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1);Sys.umask('0077')
root<-.disco_path("private","organ_endpoint_repairs_20261007_run02")
out<-file.path(root,'disease_PH');stopifnot(!dir.exists(out));dir.create(out,mode='0700')
input<-.disco_path("private","month_organ_execution_20261007/input_organ.csv")
death<-.disco_path("private","month_organ_execution_20261007/private_organ_endpoint.csv")
source8<-.disco_path("private","month_disease_prs_execution_20261007/private_endpoint_8.csv")
patch8<-.disco_path("private","targeted_endpoint_diagnostics_20261007.JaNZnJ/private_organ_52_endpoint_8.csv")
m<-fread(input,colClasses=c(eid='character'));ep<-fread(death,colClasses=c(eid='character'))
ii<-match(m$eid,ep$eid);stopifnot(!anyNA(ii));m[,time:=ep$time[ii]];m[,event:=as.integer(ep$event[ii])]
m[,sex:=factor(sex)];m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
branches<-list(age_sex=c('age','sex'),expanded=c('age','sex','bmi','race','edu','smoke','drink'))
organs<-c('tongue','salivary.gland','stomach','adipose.tissue','choroid.plexus','retina','lung','pancreas',
 'kidney','heart.muscle','skin','esophagus','skeletal.muscle','intestine','bone.marrow','liver','brain','lymphoid.tissue')
dis_names<-c(obesity='Obesity',hypertension='Hypertension',diabetes='Diabetes',cancer='Cancer',mi='MI',stroke='Stroke',copd='COPD',dementia='Dementia')
e<-rbindlist(list(fread(source8,colClasses=c(eid='character')),fread(patch8,colClasses=c(eid='character'))),use.names=TRUE)
stopifnot(!anyDuplicated(e[,.(eid,disease)]));e<-e[eid%in%m$eid]
expected<-fread(file.path(root,'organ/organ_disease_models.csv'))[policy=='death_first']
checks<-PH<-errors<-list()
for(dis in names(dis_names)){
  ee<-e[disease==dis];idx<-match(m$eid,ee$eid);stopifnot(!anyNA(idx))
  for(branch in names(branches))for(org in organs){
    z<-m[,branches[[branch]],with=FALSE];z[,score:=m[[org]]]
    z[,time:=ee$time[idx]];z[,event:=as.integer(ee$event[idx]==1)];z[,reason:=ee$reason[idx]]
    z<-z[is.na(reason)|reason==''];need<-c('score','time','event',branches[[branch]])
    z<-z[complete.cases(z[,..need])&time>0]
    ff<-reformulate(c('score',branches[[branch]]),'Surv(time,event)');ws<-character()
    f<-withCallingHandlers(coxph(ff,data=z,ties='efron',x=TRUE,y=TRUE,model=TRUE),
      warning=function(w){ws<<-c(ws,conditionMessage(w));invokeRestart('muffleWarning')})
    b<-unname(coef(f)['score']);se<-sqrt(vcov(f)['score','score'])
    wanted_branch<-branch
    target<-expected[organ==org&branch==wanted_branch&outcome==dis_names[[dis]]]
    stopifnot(nrow(target)==1,f$n==target$N,f$nevent==target$events,abs(b-target$beta)<1e-9,abs(se-target$SE)<1e-9)
    key<-paste('incident_disease',branch,org,dis_names[[dis]],'competing_death_censor','death_first',sep='|')
    err<-'';ph<-tryCatch(cox.zph(f,transform='km')$table,error=function(e){err<<-conditionMessage(e);NULL})
    checks[[length(checks)+1]]<-data.table(model_key=key,N=f$n,events=f$nevent,beta_abs_error=abs(b-target$beta),
      SE_abs_error=abs(se-target$SE),warning_count=length(ws),warning=paste(ws,collapse=' || '),PH_available=!is.null(ph),PH_error=err)
    if(!is.null(ph))PH[[length(PH)+1]]<-data.table(model_key=key,term=rownames(ph),chisq=ph[,'chisq'],df=ph[,'df'],P=ph[,'p'])
    else errors[[length(errors)+1]]<-data.table(model_key=key,error=err)
  }
  cat('PH completed',dis_names[[dis]],'\n');flush.console()
}
ck<-rbindlist(checks);ph<-rbindlist(PH);stopifnot(nrow(ck)==288)
fwrite(ck,file.path(out,'model_identity_checks.csv'));fwrite(ph,file.path(out,'disease_PH.csv'))
if(length(errors))fwrite(rbindlist(errors),file.path(out,'PH_errors.csv'))
write_json(list(status='EXECUTED_DIAGNOSTIC',models=288,PH_available=sum(ck$PH_available),
  warning_count=sum(ck$warning_count),model_beta_max_error=max(ck$beta_abs_error),model_SE_max_error=max(ck$SE_abs_error),
  PH_score_significant=sum(ph[term=='score']$P<.05),PH_GLOBAL_significant=sum(ph[term=='GLOBAL']$P<.05),
  source_sha256=digest(file=file.path(root,'09_disease_PH_diagnostic.R'),algo='sha256'),
  caveat='Diagnostic does not change primary model or prove a constant HR'),
  file.path(out,'execution_QA.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));cat('DONE disease PH diagnostics\n')
