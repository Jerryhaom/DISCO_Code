#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Complete missing ClinicalProt outcomes and preserve historical raw mc/raw hd scale.
# This supplements, rather than replaces, prior standardized/log-DM analyses.
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1)
a <- commandArgs(TRUE);stopifnot(length(a)==1L);out <- normalizePath(a[1])
base <- .disco_path("private","month_disease_prs_execution_20261007")
entry <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
input_paths <- c(clinical=file.path(base,'input_clinical.csv'),
                 proteomics=file.path(base,'input_proteomics.csv'),
                 clinicalprot=file.path(out,'input_clinicalprot.csv'))
sources <- c(input_paths,endpoint=file.path(base,'private_endpoint_8.csv'),
             clinical_models=file.path(base,'clinical_8/models.csv'),
             proteomics_models=file.path(base,'proteomics_8/models.csv'))
hash_before <- vapply(sources,function(p)digest(file=p,algo='sha256'),character(1))
ep <- fread(sources['endpoint'])
raw <- list();standard <- list();scaling <- list();flows <- list();warnings <- list()
coefficients <- list();vcovs <- list();curves <- list();checks <- list()
context <- ''
capture_fit <- function(expr) tryCatch(withCallingHandlers(expr,warning=function(w){
 warnings[[length(warnings)+1L]]<<-data.table(context=context,type='warning',message=conditionMessage(w));
 invokeRestart('muffleWarning')
}),error=function(e){warnings[[length(warnings)+1L]]<<-data.table(context=context,type='error',message=conditionMessage(e));NULL})
aj <- function(dd,mod,dis,policy,group){
 if(!nrow(dd))return(invisible(NULL))
 context<<-paste('ClinicalProt CIF',dis,policy,group)
 dd<-copy(dd);dd[,event_f:=factor(event,levels=0:2,labels=c('censor','disease','death'))]
 f<-capture_fit(survfit(Surv(time,event_f)~1,data=dd,conf.type='plain'))
 if(is.null(f))return(invisible(NULL))
 ix<-match('disease',f$states);dx<-match('death',f$states)
 tab<-dd[,.(n=.N,n1=sum(event==1),n2=sum(event==2)),by=time][order(time)]
 tab[,Y:=nrow(dd)-shift(cumsum(n),fill=0)]
 tab[,S:=cumprod(1-(n1+n2)/Y)];tab[,F:=cumsum(shift(S,fill=1)*n1/Y)]
 tab[,KM:=1-cumprod(1-n1/Y)]
 err<-max(abs(tab$F-f$pstate[,ix]));stopifnot(err<1e-10,all(tab$F<=tab$KM+1e-10),max(abs(rowSums(f$pstate)-1))<1e-10)
 checks[[length(checks)+1L]]<<-data.table(modality=mod,disease=dis,policy=policy,group=group,max_AJ_error=err,N=nrow(dd),events=sum(dd$event==1),deaths=sum(dd$event==2))
 grid<-sort(unique(c(0,seq(12,180,12),max(f$time))));grid<-grid[grid<=max(f$time)]
 k<-findInterval(grid,f$time);get<-function(x,initial=0)c(initial,x)[k+1]
 curves[[length(curves)+1L]]<<-data.table(modality=mod,disease=dis,policy=policy,group=group,time_months=grid,
  CIF=get(f$pstate[,ix]),CIF_low=pmax(0,get(f$lower[,ix])),CIF_high=pmin(1,get(f$upper[,ix])),death_CIF=get(f$pstate[,dx]),
  KM=get(tab$KM),N=nrow(dd),at_risk=vapply(grid,function(t)sum(dd$time>=t),integer(1)))
}
fitted <- 0L;converted <- 0L
for(mod in names(input_paths)){
 m<-fread(input_paths[mod]);stopifnot(!anyDuplicated(m$eid),all(is.finite(m$mc)),all(is.finite(m$hd)),all(m$hd>0))
 m[,loghd:=log(hd)];m[,sex:=factor(sex)]
 m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
 m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
 m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
 m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
 m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
 m[,qD:=as.integer(cut(mc,quantile(mc,0:4/4),include.lowest=TRUE))]
 m[,qDM:=as.integer(cut(hd,quantile(hd,0:4/4),include.lowest=TRUE))]
 sdv<-c(mc=sd(m$mc),hd=sd(m$hd),loghd=sd(m$loghd));means<-c(mc=mean(m$mc),hd=mean(m$hd),loghd=mean(m$loghd))
 scaling[[mod]]<-data.table(modality=mod,input_N=nrow(m),score=c('mc','hd','loghd'),mean=means,SD=sdv,
  SD_population='full saved modality score population, before baseline prevalence/date/covariate exclusions')
 saved<-if(mod!='clinicalprot')fread(file.path(base,paste0(mod,'_8/models.csv')))[version=='raw_dates' & policy=='death_first'] else NULL
 for(dis in c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')){
  e<-ep[disease==dis];stopifnot(!anyDuplicated(e$eid));idx<-match(m$eid,e$eid);stopifnot(!anyNA(idx))
  d<-copy(m);for(v in c('time','event','tie','reason'))set(d,j=v,value=e[[v]][idx])
  flows[[length(flows)+1L]]<-d[,.(N=.N),by=.(reason,event,tie)][,`:=`(modality=mod,disease=dis)]
  dd<-d[is.na(reason)|reason==''];dd[,status:=as.integer(event==1)]
  for(metric in c('mc','hd','loghd'))for(model in c('crude','extended','extended_BMI')){
   rhs<-switch(model,crude=metric,extended=paste(metric,'age+sex+smoke+drink+edu+race',sep='+'),
               extended_BMI=paste(metric,'age+sex+bmi+smoke+drink+edu+race',sep='+'))
   context<-paste(mod,dis,metric,model)
   if(mod!='clinicalprot' & metric%in%c('mc','loghd')){
    oldterm<-if(metric=='mc')'zD' else 'zDM'
    model_name<-model;z<-saved[disease==dis & term==oldterm & model==model_name]
    stopifnot(nrow(z)==1L)
    estimate<-z$beta/sdv[metric];se<-z$se/sdv[metric]
    r<-data.table(N=z$N,events=z$events,beta=estimate,se=se,HR=exp(estimate),low=exp(estimate-qnorm(.975)*se),
                 high=exp(estimate+qnorm(.975)*se),P=z$P,C=z$C,C_se=z$C_se)
    origin<-'exact linear re-expression of saved standardized model';converted<-converted+1L
   }else{
    need<-all.vars(as.formula(paste('~time+status+',rhs)))
    x<-dd[complete.cases(dd[,..need]) & is.finite(time) & time>=0]
    f<-capture_fit(coxph(as.formula(paste('Surv(time,status)~',rhs)),data=x,ties='efron',singular.ok=FALSE,
                        control=coxph.control(eps=1e-11,iter.max=60)))
    if(is.null(f))next
    b<-coef(f);V<-vcov(f);s<-summary(f);estimate<-b[metric];se<-sqrt(V[metric,metric])
    r<-data.table(N=f$n,events=f$nevent,beta=estimate,se=se,HR=exp(estimate),low=exp(estimate-qnorm(.975)*se),
                 high=exp(estimate+qnorm(.975)*se),P=2*pnorm(-abs(estimate/se)),C=unname(s$concordance[1]),C_se=unname(s$concordance[2]))
    origin<-'new targeted fit; saved scores not retrained';fitted<-fitted+1L
    meta<-data.table(modality=mod,disease=dis,metric=metric,model=model)
    coefficients[[length(coefficients)+1L]]<-cbind(meta,data.table(term=names(b),beta=unname(b),se=sqrt(diag(V))))
    vv<-as.data.table(as.table(V));setnames(vv,c('term1','term2','covariance'));vcovs[[length(vcovs)+1L]]<-cbind(meta,vv)
   }
   r[,`:=`(modality=mod,disease=dis,metric=if(metric=='mc')'DISCO' else if(metric=='hd')'DM' else 'logDM',
            score_column=metric,model=model,policy='death_first',version='raw_dates',unit='one saved score unit',
            covariates=rhs,source=origin,SD=unname(sdv[metric]),SD_input_N=nrow(m))]
   raw[[length(raw)+1L]]<-r
   if(metric!='hd'){
    z<-copy(r);z[,`:=`(unit='1 SD of full saved modality score population',beta=beta*SD,se=se*SD)]
    z[,`:=`(HR=exp(beta),low=exp(beta-qnorm(.975)*se),high=exp(beta+qnorm(.975)*se))]
    standard[[length(standard)+1L]]<-z
   }
  }
  if(mod=='clinicalprot')for(policy in c('death_first','disease_first')){
   cdd<-copy(dd);if(policy=='disease_first')cdd[tie==1,event:=1L]
   aj(cdd,mod,dis,policy,'All')
   for(metric in c('qD','qDM'))for(q in c(1,4))aj(cdd[get(metric)==q],mod,dis,policy,paste(metric,q,sep='_'))
  }
  cat('completed',mod,dis,'targeted new fits',fitted,'\n');flush.console()
 }
}
save_table<-function(x,name)if(length(x))fwrite(rbindlist(x,fill=TRUE),file.path(out,name))
save_table(raw,'three_modalities_raw_unit_models.csv');save_table(standard,'three_modalities_perSD_models.csv')
save_table(scaling,'three_modalities_scaling.csv');save_table(flows,'three_modalities_sample_flow.csv')
save_table(coefficients,'three_modalities_new_fit_coefficients.csv');save_table(vcovs,'three_modalities_new_fit_full_vcov.csv')
save_table(curves,'ClinicalProt_CIF_aggregate.csv');save_table(checks,'ClinicalProt_CIF_QA.csv');save_table(warnings,'three_modalities_model_messages.csv')
r<-rbindlist(raw);stopifnot(nrow(r)==216L,fitted==120L,converted==96L)
hr<-r[metric%in%c('DISCO','DM') & model!='crude'];c0<-r[metric%in%c('DISCO','DM') & model=='crude',
 .(modality,disease,metric,crude_N=N,crude_events=events,Crude_C=C,Crude_C_se=C_se,
   Crude_C_low=pmax(0,C-qnorm(.975)*C_se),Crude_C_high=pmin(1,C+qnorm(.975)*C_se))]
panel<-merge(hr,c0,by=c('modality','disease','metric'),all.x=TRUE)
fwrite(panel,file.path(out,'Table19_Fig3c_raw_unit_candidate.csv'))
hash_after<-vapply(sources,function(p)digest(file=p,algo='sha256'),character(1));stopifnot(identical(hash_before,hash_after))
write_json(list(status='EXECUTED',new_fit_count=fitted,exact_scale_conversion_count=converted,raw_unit_model_rows=nrow(r),
 three_modalities=c('Clinical','Prot','ClinicalProt'),combined_N=47173,source_hashes=as.list(hash_before),sources_unchanged=TRUE,
 endpoint='Existing approved raw_dates death_first; no endpoint rules or raw dates changed',
 raw_DM_warning='raw hd is not log(hd). It was newly fitted; never obtained by rescaling logDM coefficients',
 score_SD='full saved modality before disease/covariate exclusions, unlike Table1 final common Cox SD',
 covariates='current expanded_BMI branch retained: age+sex+bmi+smoke+drink+edu_q0 factor+race; missing race White',
 historical_covariate_difference='near-source reproduction uses derived education (edu_high/all qualifications); current new-endpoint branch retains edu_q0 definition. This difference must be adjudicated before final substitution',
 C_index='Apparent Harrell C; crude C from own score-only complete-case sample. Adjusted HR has separately reported denominator',
 figure_status='Numerical candidate only; original plotting segment and final manuscript/table version not modified',
 provider_coverage='Not established by this targeted run; approved 2022-12-31 is not proof of complete provider coverage',
 code_sha256=digest(file=normalizePath(gsub("~+~", " ", entry, fixed=TRUE)),algo='sha256')),file.path(out,'three_modalities_execution_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'three_modalities_sessionInfo.txt'))
cat('DONE: 120 targeted fits and 96 exact scale conversions; no historical output overwritten\n')
