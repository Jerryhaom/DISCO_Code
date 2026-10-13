#!/usr/bin/env Rscript
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(survival);library(data.table);library(jsonlite);library(digest)})
args<-commandArgs(TRUE);stopifnot(length(args)==1,args[1]%in%c('pass_1','pass_2'))
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);here<-dirname(normalizePath(entry))
out<-file.path(here,args[1]);dir.create(out,recursive=TRUE,showWarnings=FALSE)
root<-Sys.getenv('DISCO_DATA_ROOT');setDTthreads(1)
inputs<-c(scores=file.path(root,'Result2/UKB.ClinicalProt.DISCO.rds'),
 disease=file.path(root,'UKB_Data/dat_disease.rds'),mortality=file.path(root,'UKB_Data/UKB.Mortality.csv'),
 original_plot=file.path(root,'Results_final/03_multimodal_UKB/reproduce_ED7_KM_facet.R'),
 protocol=file.path(here,'PROTOCOL.txt'))
hash<-vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))
p<-readRDS(inputs['scores']);d<-readRDS(inputs['disease']);m<-fread(inputs['mortality'])
stopifnot(!anyDuplicated(p$ID),!anyDuplicated(d$eid),!anyDuplicated(m$id))
dm<-match(p$ID,d$eid);mm<-match(p$ID,m$id);stopifnot(!anyNA(dm),!anyNA(mm))
stopifnot(isTRUE(all.equal(p$time,m$time[mm],check.attributes=FALSE)),
 identical(as.integer(p$status),as.integer(m$status[mm])),all(p$status%in%c(0,1)))
d<-d[dm,];p$mcg<-as.integer(cut(p$mc,quantile(p$mc),include.lowest=TRUE,right=TRUE))
p$hdg<-as.integer(cut(p$hd,quantile(p$hd),include.lowest=TRUE,right=TRUE))
stopifnot(all(is.finite(p$mc)),all(is.finite(p$hd)),!anyNA(p$mcg),!anyNA(p$hdg))
stopifnot(max(abs(p$age*12-round(p$age*12)))<1e-6,
 max(abs(p$time*12-round(p$time*12)))<1e-6)
p$time<-round(p$time*12)/12
diseases<-c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')
labels<-c('Obesity','Hypertension','Diabetes','Cancer','Myocardial infarction','Stroke','COPD','Dementia')
groups<-c('DISCO Q1','DISCO Q4','DM Q1','DM Q4')
curves<-risk<-flows<-checks<-landmark_counts<-list();warnings<-character()
fit_logged<-function(expr)withCallingHandlers(expr,warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
# Independent event-time product updates: no survfit call and no saved fit reuse.
direct<-function(start,stop,event,times){
 S<-1;F<-0;K<-1;D<-0;ans<-matrix(NA_real_,length(times),4)
 for(i in seq_along(times)){
  t<-times[i];at<-start<t & stop>=t;Y<-sum(at)
  n1<-sum(stop==t & event==1);n2<-sum(stop==t & event==2)
  if(n1+n2>0){stopifnot(Y>=n1+n2);F<-F+S*n1/Y;D<-D+S*n2/Y;S<-S*(1-(n1+n2)/Y);K<-K*(1-n1/Y)}
  ans[i,]<-c(F,D,S,1-K)
 };colnames(ans)<-c('CIF','Death','Event_free','KM');ans
}
for(j in seq_along(diseases)){
 key<-diseases[j];dt<-d[[paste0('time_',key)]];ds<-d[[paste0('status_',key)]]
 stopifnot(max(abs(dt*12-round(dt*12)),na.rm=TRUE)<1e-6);dt<-round(dt*12)/12
 valid<-is.finite(dt)&is.finite(p$time)&is.finite(p$age)&ds%in%c(0,1)&p$time>0
 positive<-valid & dt>0
 conflict<-positive & ds==1 & p$status==1 & dt>p$time+1e-8
 same_month<-positive & ds==1 & p$status==1 & abs(dt-p$time)<1e-8
 flows[[key]]<-data.frame(disease=key,N_score=nrow(p),invalid=sum(!valid),negative=sum(valid&dt<0),
  zero=sum(valid&dt==0),positive=sum(positive),diagnosis_after_death=sum(conflict),
  diagnosis_death_same_month=sum(same_month),primary_eligible=sum(positive&!conflict&!same_month),
  observed_diagnosis_after_mortality_admin_end=sum(positive&ds==1&p$status==0&dt>p$time+1e-8),
  disease_nonevent_observation_after_mortality_end=sum(positive&ds==0&dt>p$time+1e-8))
 for(policy in c('exclude_ambiguous','disease_first','death_first')){
  keep<-positive&!conflict;if(policy=='exclude_ambiguous')keep<-keep&!same_month
  target<-ds==1 & if(policy=='death_first')dt<p$time-1e-8 else dt<=p$time+1e-8
  competing<-p$status==1 & p$time<=dt+1e-8 & !target
  evt<-ifelse(target,1L,ifelse(competing,2L,0L));end<-pmin(dt,p$time)
  for(g in seq_along(groups)){
   sel<-keep & if(g==1)p$mcg==1 else if(g==2)p$mcg==4 else if(g==3)p$hdg==1 else p$hdg==4
   z<-data.frame(age=p$age[sel],time=end[sel],event=evt[sel]);z$event_f<-factor(z$event,0:2,c('censor','disease','death'));z$id<-seq_len(nrow(z))
   stopifnot(nrow(z)>0,all(z$time>0))
   risk[[paste(key,policy,g)]]<-data.frame(disease=key,policy,group=groups[g],N=nrow(z),
     disease_events=sum(z$event==1),competing_deaths=sum(z$event==2),censored=sum(z$event==0))
   for(scale in c('followup','attained_age')){
    z$start<-if(scale=='followup')0 else round(z$age*12);z$stop<-z$start+round(z$time*12)
    f<-fit_logged(survfit(Surv(start,stop,event_f)~1,data=z,id=id,conf.type='plain'))
    k<-fit_logged(survfit(Surv(start,stop,event==1)~1,data=z,conf.type='log-log'))
    stopifnot(identical(f$time,k$time))
    ic<-match('disease',f$states);dc<-match('death',f$states)
    a<-data.frame(disease=key,label=labels[j],policy,group=groups[g],time_scale=scale,time=f$time/12,
      CIF=f$pstate[,ic],CIF_low=f$lower[,ic],CIF_high=f$upper[,ic],CIF_death=f$pstate[,dc],
      event_free=f$pstate[,1],KM=1-k$surv,KM_low=1-k$upper,KM_high=1-k$lower,n_risk=f$n.risk[,1],
      N=nrow(z),disease_events=sum(z$event==1),competing_deaths=sum(z$event==2))
    brute<-direct(z$start,z$stop,z$event,f$time)
    err<-max(abs(a$CIF-brute[,'CIF']),abs(a$CIF_death-brute[,'Death']),abs(a$KM-brute[,'KM']))
    if(err>=1e-10)cat('independent discrepancy',key,policy,groups[g],scale,err,
      'CIF',max(abs(a$CIF-brute[,'CIF'])),'death',max(abs(a$CIF_death-brute[,'Death'])),
      'KM',max(abs(a$KM-brute[,'KM'])),'\n')
    stopifnot(err<1e-10,max(abs(rowSums(f$pstate)-1))<1e-10,
      all(diff(a$CIF)>=-1e-12),all(a$CIF<=a$KM+1e-10),all(a$CIF>=-1e-12),all(a$CIF<=1+1e-12))
    checks[[paste(key,policy,g,scale)]]<-data.frame(disease=key,policy,group=groups[g],time_scale=scale,
      independent_max_error=err,state_sum_max_error=max(abs(rowSums(f$pstate)-1)),CIF_below_KM=TRUE)
    curves[[paste(key,policy,g,scale)]]<-a
    ts<-if(scale=='followup')c(5,10,15) else c(60,70,80)
    for(t in ts)landmark_counts[[paste(scale,policy,key,groups[g],t)]]<-
      sum(z$start<t*12 & z$stop>=t*12)
   }
  }
 }
 cat('completed',key,'\n');flush.console()
}
C<-rbindlist(curves);fwrite(C,file.path(out,'curves_aggregate.csv'));fwrite(rbindlist(risk),file.path(out,'risk_sets.csv'))
fwrite(rbindlist(flows),file.path(out,'sample_flow.csv'));fwrite(rbindlist(checks),file.path(out,'independent_checks.csv'))
landmark<-list()
for(scale in unique(C$time_scale))for(policy_name in unique(C$policy))for(dis in diseases)for(group_name in groups){
 a<-C[time_scale==scale & policy==policy_name & disease==dis & group==group_name]
 ts<-if(scale=='followup')c(5,10,15) else c(60,70,80)
 for(t in ts){b<-a[time<=t];if(nrow(b)){
  v<-tail(b,1);v$time<-t;v$n_risk<-landmark_counts[[paste(scale,policy_name,dis,group_name,t)]]
  v$within_observed_support<-t>=min(a$time)&t<=max(a$time)
  landmark[[paste(scale,policy_name,dis,group_name,t)]]<-v
 }}
}
fwrite(rbindlist(landmark),file.path(out,'landmark_risks.csv'))
after<-vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1));stopifnot(identical(hash,after))
write_json(list(pass=args[1],status='COMPLETED_CONDITIONAL_SAVED_ENDPOINT_SENSITIVITY',
 N_score=nrow(p),input_hash=setNames(as.list(hash),inputs),all_inputs_unchanged=TRUE,
 original_dates_verified=FALSE,positive_followup_only=TRUE,primary_policy='exclude_ambiguous',
 warnings=unique(warnings),N_curve_comparisons=length(checks),max_independent_error=max(rbindlist(checks)$independent_max_error),
 mortality_match=TRUE,methods=c('Aalen-Johansen','1-KM, same first-event risk set'),
 adjustment='Unadjusted score quartiles; not age+sex-adjusted',
 cutpoints=list(DISCO=as.list(quantile(p$mc)),DM=as.list(quantile(p$hd)))),
 file.path(out,'verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));cat('completed',args[1],'\n')
