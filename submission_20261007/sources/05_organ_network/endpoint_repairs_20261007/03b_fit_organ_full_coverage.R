#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1);Sys.umask('0077')
root<-.disco_path("private","organ_endpoint_repairs_20261007_run02");out<-file.path(root,'organ')
stopifnot(!dir.exists(out));dir.create(out,mode='0700')
oldroot<-.disco_path("private","month_organ_execution_20261007")
input<-file.path(oldroot,'input_organ.csv');death<-file.path(oldroot,'private_organ_endpoint.csv')
disease<-.disco_path("private","month_disease_prs_execution_20261007/private_endpoint_8.csv")
disease_patch<-.disco_path("private","targeted_endpoint_diagnostics_20261007.JaNZnJ/private_organ_52_endpoint_8.csv")
m<-fread(input,colClasses=c(eid='character'));ep<-fread(death,colClasses=c(eid='character'))
stopifnot(!anyDuplicated(m$eid),!anyDuplicated(ep$eid));ii<-match(m$eid,ep$eid);stopifnot(!anyNA(ii))
m[,time:=ep$time[ii]];m[,event:=as.integer(ep$event[ii])]
m[,sex:=factor(sex)];m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
organs<-c('tongue','salivary.gland','stomach','adipose.tissue','choroid.plexus','retina','lung','pancreas',
  'kidney','heart.muscle','skin','esophagus','skeletal.muscle','intestine','bone.marrow','liver','brain','lymphoid.tissue')
causes<-setNames(c('Malignant Neoplasms','Heart Diseases','Cerebrovascular diseases','Respiratory system',"Alzheimer's disease",'Others'),
  c('Cancer','Heart disease','Cerebrovascular disease','Respiratory disease',"Alzheimer's disease",'Others'))
branches<-list(age_sex=c('age','sex'),expanded=c('age','sex','bmi','race','edu','smoke','drink'))
R<-W<-coefrows<-cross<-flows<-PH<-disease_flows<-list()
fit_one<-function(dat,org,branch,outcome,riskset,sd_score,domain,policy='not_applicable'){
  key<-paste(domain,branch,org,outcome,riskset,policy,sep='|');warnings_here<-character();error_here<-''
  ff<-reformulate(c('score',branches[[branch]]),'Surv(time,event)')
  f<-tryCatch(withCallingHandlers(coxph(ff,data=dat,ties='efron',x=TRUE,y=TRUE,model=TRUE),
    warning=function(w){warnings_here<<-c(warnings_here,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error_here<<-conditionMessage(e);NULL})
  r<-data.table(domain=domain,branch=branch,organ=org,outcome=outcome,riskset=riskset,policy=policy,N=nrow(dat),events=sum(dat$event),
    SD=sd_score,SD_population='branch common complete mortality population, prior to disease-specific prevalence exclusions',
    beta=NA_real_,SE=NA_real_,HR_unit=NA_real_,HR_SD=NA_real_,lower_unit=NA_real_,upper_unit=NA_real_,lower_SD=NA_real_,upper_SD=NA_real_,P=NA_real_,C=NA_real_,C_SE=NA_real_,
    status=if(is.null(f))'FIT_ERROR'else'FITTED_NO_WARNING',warning_count=length(warnings_here),warning=paste(warnings_here,collapse=' || '),
    error=error_here,score_reliable_numeric=FALSE,needs_model_review=length(warnings_here)>0||is.null(f),formula=paste(deparse(ff),collapse=''),endpoint='actual blood day + cutoff2022-12-31')
  if(!is.null(f)){
    b<-coef(f);V<-vcov(f);ss<-summary(f);se<-sqrt(V['score','score']);bb<-b['score']
    r[,`:=`(N=f$n,events=f$nevent,beta=bb,SE=se,HR_unit=exp(bb),HR_SD=exp(bb*sd_score),
      lower_unit=exp(bb-qnorm(.975)*se),upper_unit=exp(bb+qnorm(.975)*se),lower_SD=exp((bb-qnorm(.975)*se)*sd_score),upper_SD=exp((bb+qnorm(.975)*se)*sd_score),
      P=2*pnorm(-abs(bb/se)),C=unname(ss$concordance[1]),C_SE=unname(ss$concordance[2]))]
    finite_score<-all(is.finite(c(bb,se)))&&se>0
    score_named_warning<-any(grepl('variable\\s+1\\s*;',warnings_here))
    r[,score_reliable_numeric:=finite_score&&!score_named_warning]
    r[,status:=if(!finite_score||score_named_warning)'UNRELIABLE_SCORE'else if(length(warnings_here))'WARNING_REVIEW_MODEL'else'FITTED_NO_WARNING']
    # Preserve finite score estimates when warning concerns a nuisance variable.
    # No scientific reliability claim is inferred merely from finite coefficients.
    coefrows[[length(coefrows)+1]]<<-data.table(model_key=key,term=names(b),beta=unname(b),SE=sqrt(diag(V)),finite=is.finite(b)&is.finite(sqrt(diag(V))))
    for(w in warnings_here){
      hit<-regexec('variable\\s+([0-9 ,]+)\\s*;',w);z<-regmatches(w,hit)[[1]]
      indices<-if(length(z)>1)as.integer(strsplit(gsub(' ','',z[2]),',',fixed=TRUE)[[1]])else integer()
      terms<-names(b)[indices];W[[length(W)+1]]<<-data.table(model_key=key,warning=w,implicated_terms=paste(terms,collapse=','),score_implicated='score'%in%terms)
    }
    if(domain=='mortality'&&riskset!='legacy_delete'){
      z<-tryCatch(cox.zph(f,transform='km')$table,error=function(e)NULL)
      if(!is.null(z))PH[[length(PH)+1]]<<-data.table(model_key=key,term=rownames(z),chisq=z[,'chisq'],df=z[,'df'],P=z[,'p'])
    }
  }else W[[length(W)+1]]<<-data.table(model_key=key,warning=error_here,implicated_terms='',score_implicated=NA)
  for(v in intersect(c('race','sex'),names(dat))){
    tt<-dat[,.(N=.N,events=sum(event)),by=v];setnames(tt,v,'level')
    cross[[length(cross)+1]]<<-data.table(model_key=key,variable=v,level=as.character(tt$level),N=tt$N,events=tt$events)
  }
  R[[length(R)+1]]<<-r
}
for(branch in names(branches))for(org in organs){
  z<-m[,c('eid','time','event','cause',branches[[branch]]),with=FALSE];z[,score:=m[[org]]]
  needs<-c('time','event','score',branches[[branch]]);z<-z[complete.cases(z[,..needs])&time>0]
  sd_score<-sd(z$score);stopifnot(is.finite(sd_score),sd_score>0)
  flows[[length(flows)+1]]<-data.table(branch=branch,organ=org,input_N=nrow(m),common_eligible_N=nrow(z),all_cause_deaths=sum(z$event),SD=sd_score)
  fit_one(z,org,branch,'All-cause','all_cause',sd_score,'mortality')
  for(oc in names(causes)){
    target<-z$event==1&!is.na(z$cause)&z$cause==causes[[oc]]
    legacy<-copy(z[target|event==0]);correct<-copy(z);correct[,event:=as.integer(target)]
    fit_one(legacy,org,branch,oc,'legacy_delete',sd_score,'mortality');fit_one(correct,org,branch,oc,'competing_death_censor',sd_score,'mortality')
  }
  cat('Mortality complete',branch,org,'\n');flush.console()
}
# Keep an aggregate-only checkpoint if disease coverage fails later.
fwrite(rbindlist(R),file.path(out,'organ_mortality_models_checkpoint.csv'))
e<-rbindlist(list(fread(disease,colClasses=c(eid='character')),fread(disease_patch,colClasses=c(eid='character'))),use.names=TRUE);stopifnot(!anyDuplicated(e[,.(eid,disease)]))
e<-e[eid%in%m$eid]
dis_names<-c(obesity='Obesity',hypertension='Hypertension',diabetes='Diabetes',cancer='Cancer',mi='MI',stroke='Stroke',copd='COPD',dementia='Dementia')
for(dis in names(dis_names)){
  ee<-e[disease==dis];idx<-match(m$eid,ee$eid);covered<-!is.na(idx)
  stopifnot(all(covered))
  for(branch in names(branches))for(org in organs){
    z<-m[,c('eid',branches[[branch]]),with=FALSE];z[,score:=m[[org]]]
    z[,time:=ee$time[idx]];z[,event:=ee$event[idx]];z[,reason:=ee$reason[idx]];z[,tie:=ee$tie[idx]]
    # A missing endpoint row is not evidence of disease absence/censoring.
    z<-z[covered]
    source_N<-nrow(z);excluded_prevalent<-sum(z$reason=='prevalent',na.rm=TRUE)
    excluded_unknown<-sum(z$reason=='known_disease_unknown_date',na.rm=TRUE)
    z<-z[is.na(reason)|reason==''];needs<-c('time','event','score',branches[[branch]])
    before_complete_N<-nrow(z)
    z<-z[complete.cases(z[,..needs])&time>0]
    disease_flows[[length(disease_flows)+1]]<-data.table(disease=dis_names[[dis]],branch=branch,organ=org,
      input_N=nrow(m),covered_endpoint_N=source_N,not_in_endpoint_source=sum(!covered),
      prevalent_excluded=excluded_prevalent,unknown_diagnosis_date_excluded=excluded_unknown,
      other_covariate_time_score_excluded=before_complete_N-nrow(z),eligible_N=nrow(z),
      primary_events=sum(z$event==1),same_day_ties=sum(z$tie==1,na.rm=TRUE))
    sd_score<-flows[[which(vapply(flows,function(x)x$branch==branch&&x$organ==org,FALSE))]]$SD
    for(pol in c('death_first','disease_first')){
      zz<-copy(z);if(pol=='disease_first')zz[tie==1,event:=1L];zz[,event:=as.integer(event==1)]
      fit_one(zz,org,branch,dis_names[[dis]],'competing_death_censor',sd_score,'incident_disease',pol)
    }
  }
  cat('Disease complete',dis_names[[dis]],'\n');flush.console()
}
r<-rbindlist(R);stopifnot(nrow(r[domain=='mortality'])==468,nrow(r[domain=='incident_disease'])==576)
fwrite(r[domain=='mortality'],file.path(out,'organ_mortality_models.csv'));fwrite(r[domain=='incident_disease'],file.path(out,'organ_disease_models.csv'))
fwrite(rbindlist(flows),file.path(out,'organ_common_sample_flows.csv'));fwrite(rbindlist(coefrows),file.path(out,'full_coefficient_diagnostics.csv'))
fwrite(rbindlist(disease_flows),file.path(out,'organ_disease_sample_flows.csv'))
if(length(W))fwrite(rbindlist(W),file.path(out,'model_warnings.csv'))
fwrite(rbindlist(cross),file.path(out,'race_sex_event_cross_tabs.csv'));fwrite(rbindlist(PH),file.path(out,'mortality_PH.csv'))
paths<-c(input,death,disease,disease_patch,file.path(root,'03b_fit_organ_full_coverage.R'))
write_json(list(status='EXECUTED',mortality_rows=468,disease_rows=576,warning_rows=length(W),
  models_needing_review=sum(r$needs_model_review),unreliable_score_models=sum(!r$score_reliable_numeric),
  policy='primary death-first; same-day disease-first sensitivity; competing death censoring retained',
  formulas=list(age_sex='score+age+sex',expanded_sensitivity='score+age+sex+BMI+race+education+smoking+alcohol'),
  HR_units='per original score unit and per fixed branch/common-population SD both exported',
  warning_policy='finite score estimates retained with explicit review flags; nuisance separation not automatically proof of score reliability',
  endpoint_coverage_policy='participants absent from endpoint source explicitly excluded, not assigned disease-free status; counts exported by disease/organ/branch',
  input_code_sha256=setNames(vapply(paths,function(p)digest(file=p,algo='sha256'),''),basename(paths))),file.path(out,'execution_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));cat('DONE organ mortality/disease repairs\n')
