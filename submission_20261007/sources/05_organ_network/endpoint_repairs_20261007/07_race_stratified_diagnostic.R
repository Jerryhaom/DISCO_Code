#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Post-hoc separation diagnostic; not a replacement primary model.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1);Sys.umask('0077')
root<-.disco_path("private","organ_endpoint_repairs_20261007_run01")
out<-file.path(root,'race_diagnostic');stopifnot(!dir.exists(out));dir.create(out,mode='0700')
m<-fread(.disco_path("private","month_organ_execution_20261007/input_organ.csv"),colClasses=c(eid='character'))
ep<-fread(.disco_path("private","month_organ_execution_20261007/private_organ_endpoint.csv"),colClasses=c(eid='character'))
ii<-match(m$eid,ep$eid);stopifnot(!anyNA(ii));m[,time:=ep$time[ii]];m[,event:=as.integer(ep$event[ii])]
m[,sex:=factor(sex)];m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
main<-fread(file.path(root,'organ/organ_mortality_models.csv'))[branch=='expanded'&outcome=="Alzheimer's disease"&riskset=='competing_death_censor']
R<-W<-list()
for(org in main$organ){
  z<-m[,.(time,event,cause,age,sex,bmi,race,edu,smoke,drink)];z[,score:=m[[org]]]
  z<-z[complete.cases(z[,.(time,event,score,age,sex,bmi,race,edu,smoke,drink)])&time>0]
  z[,event:=as.integer(event==1&!is.na(cause)&cause=="Alzheimer's disease")]
  ws<-character();f<-withCallingHandlers(coxph(Surv(time,event)~score+age+sex+bmi+edu+smoke+drink+strata(race),
    data=z,ties='efron',x=TRUE),warning=function(w){ws<<-c(ws,conditionMessage(w));invokeRestart('muffleWarning')})
  beta<-coef(f)['score'];SE<-sqrt(vcov(f)['score','score']);SD<-main$SD[match(org,main$organ)]
  stopifnot(all(is.finite(c(beta,SE))),SE>0,f$n==50572,f$nevent==145)
  R[[length(R)+1]]<-data.table(organ=org,outcome="Alzheimer's disease",N=f$n,events=f$nevent,
    HR_unit=exp(beta),HR_SD=exp(beta*SD),lower_SD=exp((beta-qnorm(.975)*SE)*SD),upper_SD=exp((beta+qnorm(.975)*SE)*SD),
    P=2*pnorm(-abs(beta/SE)),SD=SD,unstratified_HR_SD=main$HR_SD[match(org,main$organ)],
    warning_count=length(ws),warning=paste(ws,collapse=' || '),
    formula=paste(deparse(formula(f)),collapse=''),status='post-hoc race-stratified separation diagnostic; not main model replacement')
  if(length(ws))W[[length(W)+1]]<-data.table(organ=org,warning=ws)
}
r<-rbindlist(R);r[,HR_SD_change:=HR_SD-unstratified_HR_SD]
fwrite(r,file.path(out,'race_stratified_Alzheimer_diagnostic.csv'));if(length(W))fwrite(rbindlist(W),file.path(out,'warnings.csv'))
write_json(list(status='EXECUTED_DIAGNOSTIC_ONLY',models=nrow(r),warning_count=sum(r$warning_count),
  N=50572,events=145,NonWhite_N=3054,NonWhite_events=0,max_abs_score_HR_SD_change=max(abs(r$HR_SD_change)),
  policy='race strata retains all eligible participants; no race coefficient estimated; original flagged estimates retained separately',
  caveats=c('Does not establish a constant score HR or cure PH violations','No proof of robustness from finite coefficients alone'),
  source_sha256=digest(file=file.path(root,'07_race_stratified_diagnostic.R'),algo='sha256')),
  file.path(out,'execution_QA.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));cat('DONE race-stratified diagnostic\n')
