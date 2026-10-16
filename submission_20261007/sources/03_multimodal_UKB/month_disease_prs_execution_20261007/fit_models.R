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
setDTthreads(1)
a<-commandArgs(TRUE);mod<-a[1];task<-a[2]
entry<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE)))
out<-file.path(root,paste(mod,task,sep='_'));dir.create(out,showWarnings=FALSE,mode='0700')
m<-fread(file.path(root,paste0('input_',mod,'.csv')));stopifnot(!anyDuplicated(m$eid))
m[,zD:=as.numeric(scale(mc))];m[,zDM:=as.numeric(scale(log(hd)))];m[,sex:=factor(sex)]
m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
# Missing ethnicity imputed to White (author-approved C02 rule).
m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
m[,qD:=as.integer(cut(mc,quantile(mc,0:4/4,na.rm=TRUE),include.lowest=TRUE))]
m[,qDM:=as.integer(cut(hd,quantile(hd,0:4/4,na.rm=TRUE),include.lowest=TRUE))]
ep<-fread(file.path(root,paste0('private_endpoint_',task,'.csv')))
results<-list();vcovs<-list();betas<-list();flows<-list();curves<-list();checks<-list();errors<-list();scales<-list()
context<-''
logged<-function(expr) tryCatch(withCallingHandlers(expr,warning=function(w){
 errors[[length(errors)+1L]]<<-data.table(context=context,type='warning',message=conditionMessage(w));cat('WARNING',context,conditionMessage(w),'\n');invokeRestart('muffleWarning')
}),error=function(e){errors[[length(errors)+1L]]<<-data.table(context=context,type='error',message=conditionMessage(e));cat('ERROR',context,conditionMessage(e),'\n');NULL})
fit<-function(dd,rhs,disease,policy,version,model,terms,joint=FALSE){
 context<<-paste(disease,policy,version,model)
 need<-all.vars(as.formula(paste('~ time+status+',rhs)));dd<-dd[complete.cases(dd[,..need]) & is.finite(time) & time>=0]
 if(nrow(dd)<500 || sum(dd$status)<5){errors[[length(errors)+1L]]<<-data.table(context=context,type='insufficient',message=paste('N',nrow(dd),'events',sum(dd$status)));return(NULL)}
 f<-logged(coxph(as.formula(paste('Surv(time,status)~',rhs)),data=dd,ties='efron',singular.ok=FALSE))
 if(is.null(f))return(NULL)
 b<-coef(f);V<-vcov(f);s<-summary(f)
 meta<-data.table(modality=mod,task=task,disease=disease,policy=policy,version=version,model=model)
 vv<-as.data.table(as.table(V));setnames(vv,c('term1','term2','covariance'));vcovs[[length(vcovs)+1L]]<<-cbind(meta,vv)
 betas[[length(betas)+1L]]<<-cbind(meta,data.table(term=names(b),beta=as.numeric(b),se=sqrt(diag(V))))
 for(term in c(terms,if(joint)'Joint_11_vs_00')){
  L<-setNames(rep(0,length(b)),names(b))
  if(term=='Joint_11_vs_00')L[c('zD','zPRS','zD:zPRS')]<-1 else L[term]<-1
  estimate<-sum(L*b);se<-sqrt(as.numeric(t(L)%*%V%*%L));z<-estimate/se
  results[[length(results)+1L]]<<-cbind(meta,data.table(term=term,N=f$n,events=f$nevent,
   beta=estimate,se=se,HR=exp(estimate),low=exp(estimate-qnorm(.975)*se),high=exp(estimate+qnorm(.975)*se),P=2*pnorm(-abs(z)),
   C=unname(s$concordance[1]),C_se=unname(s$concordance[2]),covariates=rhs))
 }
 invisible(f)
}
aj<-function(dd,disease,policy,group){
 if(!nrow(dd))return(NULL)
 context<<-paste('CIF',disease,policy,group)
 dd<-copy(dd);dd[,event_f:=factor(event,levels=0:2,labels=c('censor','disease','death'))]
 f<-logged(survfit(Surv(time,event_f)~1,data=dd,conf.type='plain'))
 if(is.null(f))return(NULL)
 ix<-match('disease',f$states);dx<-match('death',f$states)
 # Independent event-table Aalen-Johansen and cause-specific Kaplan-Meier.
 tab<-dd[,.(n=.N,n1=sum(event==1),n2=sum(event==2)),by=time][order(time)]
 tab[,Y:=nrow(dd)-shift(cumsum(n),fill=0)]
 tab[,S:=cumprod(1-(n1+n2)/Y)];tab[,F:=cumsum(shift(S,fill=1)*n1/Y)]
 tab[,KM:=1-cumprod(1-n1/Y)]
 stopifnot(identical(as.numeric(tab$time),as.numeric(f$time)))
 err<-max(abs(tab$F-f$pstate[,ix]));stopifnot(err<1e-10,all(tab$F<=tab$KM+1e-10),max(abs(rowSums(f$pstate)-1))<1e-10)
 checks[[length(checks)+1L]]<<-data.table(modality=mod,disease,policy,group,max_AJ_error=err,N=nrow(dd),events=sum(dd$event==1),deaths=sum(dd$event==2))
 # Monthly display grid only; fitting and event updates use fractional months.
 grid<-sort(unique(c(0,seq(12,180,12),max(f$time))));grid<-grid[grid<=max(f$time)]
 k<-findInterval(grid,f$time);get<-function(x,initial=0)c(initial,x)[k+1]
 curves[[length(curves)+1L]]<<-data.table(modality=mod,disease,policy,group,time_months=grid,
  CIF=get(f$pstate[,ix]),CIF_low=pmax(0,get(f$lower[,ix])),CIF_high=pmin(1,get(f$upper[,ix])),death_CIF=get(f$pstate[,dx]),
  KM=get(tab$KM),N=nrow(dd),at_risk=vapply(grid,function(t)sum(dd$time>=t),integer(1)))
}
if(task=='8'){
 old<-as.data.table(readRDS(.disco_path("private","\u7ed3\u5c40/UKB/\u65b0\u53d1\u75be\u75c5/dat_disease.rds")))
 diseases<-c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')
 for(dis in diseases){
  e<-ep[disease==dis];idx<-match(m$eid,e$eid);stopifnot(!anyNA(idx));d<-copy(m)
  for(v in c('time','event','tie','reason','after_death','after_loss','after_cutoff','unknown_date','prevalent'))set(d,j=v,value=e[[v]][idx])
  flows[[dis]]<-d[,.(N=.N),by=.(reason,event,tie)][,`:=`(modality=mod,disease=dis)]
  for(policy in c('death_first','disease_first')){
   dd<-copy(d[is.na(reason)|reason=='']);if(policy=='disease_first')dd[tie==1,event:=1L]
   dd[,status:=as.integer(event==1)]
   for(metric in c('zD','zDM')){
    fit(dd,metric,dis,policy,'raw_dates','crude',metric)
    fit(dd,paste(metric,'age+sex+smoke+drink+edu+race',sep='+'),dis,policy,'raw_dates','extended',metric)
    fit(dd,paste(metric,'age+sex+bmi+smoke+drink+edu+race',sep='+'),dis,policy,'raw_dates','extended_BMI',metric)
   }
   aj(dd,dis,policy,'All')
   for(metric in c('qD','qDM'))for(q in c(1,4))aj(dd[get(metric)==q],dis,policy,paste(metric,q,sep='_'))
  }
  oo<-copy(m);j<-match(oo$eid,old$eid);oo[,time:=old[[paste0('time_',dis)]][j]*12];oo[,status:=old[[paste0('status_',dis)]][j]]
  for(metric in c('zD','zDM'))for(model in c('crude','extended','extended_BMI')){
   rhs<-switch(model,crude=metric,extended=paste(metric,'age+sex+smoke+drink+edu+race',sep='+'),extended_BMI=paste(metric,'age+sex+bmi+smoke+drink+edu+race',sep='+'))
   fit(oo,rhs,dis,'legacy','old_endpoint_same_covariates',model,metric)
  }
  cat('completed 8',mod,dis,'\n');flush.console()
 }
}else{
 names28<-c('amd','ad','asthma','af','bipolar','bowel_ca','breast_ca','cvd','coeliac','cad','crohn','ovarian_ca','hypertension','ischaemic_stroke','melanoma','ms','osteoporosis','parkinson','glaucoma','prostate_ca','psoriasis','ra','schizophrenia','sle','t1d','t2d','ulcerative_colitis','vte')
 fields<-paste0(c(26204,26206,26210,26212,26214,26218,26220,26223,26225,26227,26229,26232,26244,26248,26252,26254,26258,26260,26265,26267,26269,26273,26275,26278,26283,26285,26287,26289),'-0.0')
 prs<-fread(file.path(root,'PRS_input.csv'),select=c('eid',fields));stopifnot(!anyDuplicated(prs$eid))
 old<-fread(.disco_path("private","bioage_v2/outputs/dat_disease_28_rebuilt.csv"))
 for(jj in seq_along(names28)){
  dis<-names28[jj];e<-ep[disease==dis];d<-copy(m)
  sexrestriction<-if(dis%in%c('breast_ca','ovarian_ca'))'Female' else if(dis=='prostate_ca')'Male' else NULL
  if(!is.null(sexrestriction))d<-d[sex==sexrestriction]
  p<-prs[[fields[jj]]][match(d$eid,prs$eid)];d[,zPRS:=as.numeric(scale(p))]
  scales[[dis]]<-data.table(modality=mod,disease=dis,D_mean=mean(m$mc),D_sd=sd(m$mc),PRS_mean=mean(p,na.rm=TRUE),PRS_sd=sd(p,na.rm=TRUE),scaling_N=nrow(d),prs_available=sum(!is.na(p)))
  oi<-match(d$eid,old$eid);ot<-old[[paste0('time_',dis)]][oi]*12;os<-old[[paste0('status_',dis)]][oi]
  ot[is.na(ot)]<-d$time[is.na(ot)]*12;os[is.na(os)]<-0
  idx<-match(d$eid,e$eid);stopifnot(!anyNA(idx))
  for(v in c('time','event','tie','reason'))set(d,j=v,value=e[[v]][idx])
  flows[[dis]]<-d[,.(N=.N),by=.(reason,event,tie)][,`:=`(modality=mod,disease=dis)]
  cv<-paste(c('age',if(is.null(sexrestriction))'sex','smoke','drink','edu','race'),collapse='+')
  for(policy in c('death_first','disease_first','legacy')){
   dd<-copy(d)
   if(policy=='legacy'){dd[,time:=ot];dd[,status:=os];version<-'old_endpoint_same_covariates'} else {
    dd<-dd[is.na(reason)|reason==''];if(policy=='disease_first')dd[tie==1,event:=1L]
    dd[,status:=as.integer(event==1)];version<-'raw_dates'
   }
   fit(dd,paste('zD+zPRS',cv,sep='+'),dis,policy,version,'f0_additive',c('zD','zPRS'))
   fit(dd,paste('zD*zPRS',cv,sep='+'),dis,policy,version,'f1_interaction',c('zD','zPRS','zD:zPRS'),TRUE)
  }
  cat('completed 28',mod,dis,'\n');flush.console()
 }
}
save<-function(x,n)if(length(x))fwrite(rbindlist(x,fill=TRUE),file.path(out,n))
save(results,'models.csv');save(vcovs,'full_vcov.csv');save(betas,'coefficients.csv');save(flows,'sample_flow.csv');save(curves,'CIF_aggregate.csv');save(checks,'CIF_QA.csv');save(errors,'model_messages.csv');save(scales,'scaling.csv')
r<-rbindlist(results,fill=TRUE);new<-r[version=='raw_dates' & policy=='death_first'];oldr<-r[version=='old_endpoint_same_covariates']
comparison<-merge(new,oldr,by=c('modality','task','disease','model','term'),suffixes=c('_new','_old'),all=TRUE)
fwrite(comparison,file.path(out,'old_new_comparison.csv'))
write_json(list(modality=mod,task=task,status='EXECUTED',model_rows=nrow(r),errors=errors,
 input_N=nrow(m),time_unit='days/(365.25/12)',joint='f1 contrast (zD,zPRS)=(1,1) vs (0,0): betaD+betaPRS+betaDxPRS; full covariance',
 race_missing='missing/blank imputed to White (author-approved C02 rule)',date_tie='death-first analysis convention; disease-first sensitivity',
 code_sha256=digest(file=normalizePath(gsub("~+~", " ", entry, fixed=TRUE)),algo='sha256')),file.path(out,'execution_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('DONE',mod,task,'\n')
