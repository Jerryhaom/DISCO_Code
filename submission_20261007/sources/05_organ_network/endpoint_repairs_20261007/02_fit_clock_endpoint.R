#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Frozen 7:3 scores, actual blood-day endpoint. All nuisance fits use training rows.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1);Sys.umask('0077')
root<-.disco_path("private","organ_endpoint_repairs_20261007_run01")
out<-file.path(root,'clock');stopifnot(!dir.exists(out));dir.create(out,mode='0700')
scores<-.disco_path("private","month_7030_execution_20261007/scores_7030.csv")
endpoint<-.disco_path("private","month_organ_execution_20261007/private_organ_endpoint.csv")
d<-fread(scores,colClasses=c(eid='character'));ep<-fread(endpoint,colClasses=c(eid='character'))
matched<-fread(file.path(root,'matched_scores.csv'),colClasses=c(eid='character'))
stopifnot(!anyDuplicated(d$eid),!anyDuplicated(ep$eid),!anyDuplicated(matched$eid))
ie<-match(d$eid,ep$eid);im<-match(d$eid,matched$eid);stopifnot(!anyNA(ie),!anyNA(im))
d[,time:=ep$time[ie]];d[,event:=as.integer(ep$event[ie])];d[,sex:=factor(sex)]
for(k in setdiff(names(matched),'eid'))d[[k]]<-matched[[k]][im]
stopifnot(all(is.finite(d$time)),all(d$time>0),all(d$event%in%0:1),sum(d$is_train)==31381,sum(d$is_test)==13450,
  all(d$is_train+d$is_test==1),sum(d$event[d$is_train==1])==3318,sum(d$event[d$is_test==1])==1423)
orgs<-c('Adipose','Brain','Heart','Intestine','Kidney','Liver','Lung','Muscle','Pancreas')
CI<-pairs<-HR<-PH<-cal<-corr<-piece<-warns<-flows<-list()
label_now<-''
fit<-function(f,dat,label){label_now<<-label;withCallingHandlers(coxph(f,data=dat,ties='efron',x=TRUE,y=TRUE,model=TRUE),
  warning=function(w){warns[[length(warns)+1]]<<-data.table(model=label,warning=conditionMessage(w));invokeRestart('muffleWarning')})}
termrow<-function(f,t,o,m,model,scale){b<-coef(f)[t];s<-sqrt(vcov(f)[t,t]);data.table(organ=o,method=m,model=model,term=t,
  N=f$n,deaths=f$nevent,beta=unname(b),SE=unname(s),HR=exp(b),lower95=exp(b-qnorm(.975)*s),upper95=exp(b+qnorm(.975)*s),P=2*pnorm(-abs(b/s)),
  formula=paste(deparse(formula(f)),collapse=''),scale=scale,fit_sample='frozen_common_train',endpoint='actual blood day; cutoff2022-12-31; fractional months')}
for(o in orgs){
  a<-d[,.(eid,time,event,age,sex,is_train,is_test)]
  a$D<-d[[paste0('D_',o)]];a$O<-d[[paste0('Oh_',o)]];a$G<-d[[paste0('Goeminne_',o)]]
  a$DO<-d[[paste0('M_Oh_',o)]];a$DG<-d[[paste0('M_Goeminne_',o)]]
  tr<-a[is_train==1];te<-a[is_test==1];stopifnot(all(complete.cases(tr)),all(complete.cases(te)))
  for(k in c('D','O','G','DO','DG')){
    f<-lm(reformulate('age',k),data=tr);den<-sd(resid(f));mu<-mean(tr[[k]]);rawsd<-sd(tr[[k]])
    stopifnot(den>0,rawsd>0)
    tr[[paste0(k,'a')]]<-resid(f)/den;te[[paste0(k,'a')]]<-(te[[k]]-predict(f,newdata=te))/den
    tr[[paste0(k,'z')]]<-(tr[[k]]-mu)/rawsd;te[[paste0(k,'z')]]<-(te[[k]]-mu)/rawsd
    stopifnot(abs(cor(tr[[paste0(k,'a')]],tr$age))<1e-10,abs(sd(tr[[paste0(k,'a')]])-1)<1e-10)
    cal[[length(cal)+1]]<-data.table(organ=o,score=k,train_N=nrow(tr),test_N_used=0,mean=mu,raw_SD=rawsd,
      age_intercept=unname(coef(f)[1]),age_slope=unname(coef(f)[2]),age_residual_SD=den)
  }
  for(m in c('Oh','Goeminne')){
    k<-if(m=='Oh')'O'else'G'
    ft<-lm(D~age+sex,tr);fp<-lm(reformulate(c('age','sex'),k),tr)
    rd<-te$D-predict(ft,newdata=te);rp<-te[[k]]-predict(fp,newdata=te);r<-cor(rd,rp)
    corr[[length(corr)+1]]<-data.table(organ=o,method=m,N=nrow(te),r=r,
      lower95=tanh(atanh(r)-qnorm(.975)/sqrt(nrow(te)-3)),upper95=tanh(atanh(r)+qnorm(.975)/sqrt(nrow(te)-3)),
      P=2*pt(-abs(r*sqrt((nrow(te)-2)/(1-r*r))),nrow(te)-2),adjustment='age+sex nuisance fits learned in train; applied to test')
    for(scale in c('age_residual_SD','raw_score_SD')){
      vD<-if(scale=='age_residual_SD')'Da'else'Dz';vP<-paste0(k,if(scale=='age_residual_SD')'a'else'z')
      fs<-list(D_only=reformulate(c('age','sex',vD),'Surv(time,event)'),P_only=reformulate(c('age','sex',vP),'Surv(time,event)'),
        joint=reformulate(c('age','sex',vD,vP),'Surv(time,event)'))
      for(n in names(fs)){
        f<-fit(fs[[n]],tr,paste(o,m,scale,n,sep='|'));stopifnot(all(is.finite(coef(f))))
        for(t in intersect(c(vD,vP),names(coef(f))))HR[[length(HR)+1]]<-termrow(f,t,o,m,n,scale)
        if(n=='joint'){
          ph<-cox.zph(f,transform='km')$table
          PH[[length(PH)+1]]<-data.table(organ=o,method=m,basis=scale,term=rownames(ph),chisq=ph[,'chisq'],df=ph[,'df'],P=ph[,'p'],N=f$n,deaths=f$nevent)
        }
      }
    }
    # Exploratory 0-5yr / >5yr sensitivity, not a replacement primary model.
    long<-survSplit(Surv(time,event)~.,data=as.data.frame(tr),cut=60,start='tstart',episode='period')
    long$late<-as.integer(long$period==2)
    vP<-paste0(k,'a')
    ff<-reformulate(c('strata(period)','age','sex','Da',vP,'age:late','sex:late','Da:late',paste0(vP,':late')),'Surv(tstart,time,event)')
    pf<-fit(ff,long,paste(o,m,'exploratory_5yr',sep='|'))
    for(t in c('Da',vP))for(per in 1:2){
      b<-coef(pf)[t];V<-vcov(pf);se2<-V[t,t];it<-paste0(t,':late')
      if(per==2){b<-b+coef(pf)[it];se2<-se2+V[it,it]+2*V[t,it]}
      se<-sqrt(se2)
      piece[[length(piece)+1]]<-data.table(organ=o,method=m,term=t,period=if(per==1)'0-60 months'else'>60 months',
        HR=exp(b),lower95=exp(b-qnorm(.975)*se),upper95=exp(b+qnorm(.975)*se),P=2*pnorm(-abs(b/se)),
        interaction_P=2*pnorm(-abs(coef(pf)[it]/sqrt(V[it,it]))),N=nrow(tr),events=sum(long$event[long$period==per]),
        scale='training age-residual SD',status='post-hoc exploratory; not main HR replacement',formula=paste(deparse(ff),collapse=''))
    }
  }
  fs<-list(Fair_Age_sex=Surv(time,event)~age+sex,Fair_DISCO=Surv(time,event)~age+sex+D,
    Fair_Oh=Surv(time,event)~age+sex+O,Fair_Goeminne=Surv(time,event)~age+sex+G,
    Fair_DISCO_Oh=Surv(time,event)~age+sex+D+O,Fair_DISCO_Goeminne=Surv(time,event)~age+sex+D+G,
    Display_DISCO=Surv(time,event)~age+D,Display_Oh=Surv(time,event)~O,Display_Goeminne=Surv(time,event)~G,
    Accel_DISCO=Surv(time,event)~Da,Accel_Oh=Surv(time,event)~Oa,Accel_Goeminne=Surv(time,event)~Ga,
    Matched_Display_DISCO_Oh=Surv(time,event)~age+DO,Matched_Accel_DISCO_Oh=Surv(time,event)~DOa,Matched_Fair_DISCO_Oh=Surv(time,event)~age+sex+DO,
    Matched_Display_DISCO_Goeminne=Surv(time,event)~age+DG,Matched_Accel_DISCO_Goeminne=Surv(time,event)~DGa,Matched_Fair_DISCO_Goeminne=Surv(time,event)~age+sex+DG)
  fits<-lapply(names(fs),function(n)fit(fs[[n]],tr,paste(o,n,sep='|')));names(fits)<-names(fs)
  cc<-do.call(concordance,c(unname(fits),list(newdata=te,timewt='n',influence=2)));cv<-cc$concordance;V<-cc$var
  stopifnot(all(is.finite(V)),length(cv)==18,min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
  for(j in seq_along(fits)){
    lp<-predict(fits[[j]],newdata=te,type='lp');direct<-concordance(Surv(te$time,te$event)~lp,reverse=TRUE)
    stopifnot(abs(cv[j]-direct$concordance)<1e-12,abs(V[j,j]-direct$var)<1e-12)
    se<-sqrt(V[j,j]);CI[[length(CI)+1]]<-data.table(organ=o,model=names(fits)[j],N=nrow(te),deaths=sum(te$event),Cindex=cv[j],SE=se,
      lower95=max(0,cv[j]-qnorm(.975)*se),upper95=min(1,cv[j]+qnorm(.975)*se),formula=paste(deparse(fs[[j]]),collapse=''),
      fit_sample='frozen_common_train',evaluation='frozen_common_test',CI='paired IJ conditional on frozen trained scores/models')
  }
  contrasts<-list(c('Fair_DISCO','Fair_Oh'),c('Fair_DISCO','Fair_Goeminne'),c('Fair_DISCO_Oh','Fair_Oh'),c('Fair_DISCO_Goeminne','Fair_Goeminne'),
    c('Fair_DISCO_Oh','Fair_DISCO'),c('Fair_DISCO_Goeminne','Fair_DISCO'),c('Display_DISCO','Display_Oh'),c('Display_DISCO','Display_Goeminne'),
    c('Accel_DISCO','Accel_Oh'),c('Accel_DISCO','Accel_Goeminne'))
  for(m in c('Oh','Goeminne'))contrasts<-c(contrasts,list(c(paste0('Matched_Fair_DISCO_',m),paste0('Fair_',m)),
    c(paste0('Matched_Fair_DISCO_',m),'Fair_DISCO'),c(paste0('Matched_Display_DISCO_',m),paste0('Display_',m)),c(paste0('Matched_Accel_DISCO_',m),paste0('Accel_',m))))
  for(v in contrasts){
    i<-match(v[1],names(fits));j<-match(v[2],names(fits));delta<-cv[i]-cv[j];sv<-V[i,i]+V[j,j]-2*V[i,j];stopifnot(is.finite(sv),sv>=-1e-12)
    se<-sqrt(max(0,sv));pairs[[length(pairs)+1]]<-data.table(organ=o,contrast=paste(v,collapse=' minus '),N=nrow(te),deaths=sum(te$event),
      delta_Cindex=delta,SE=se,lower95=delta-qnorm(.975)*se,upper95=delta+qnorm(.975)*se,P=if(se>0)2*pnorm(-abs(delta/se))else NA_real_)
  }
  flows[[length(flows)+1]]<-data.table(organ=o,train_N=nrow(tr),train_deaths=sum(tr$event),test_N=nrow(te),test_deaths=sum(te$event),reference_N=344,reference_test_overlap=0)
  cat('COMPLETED',o,'18 C models; both HR/PH bases; matched score new endpoint; 5yr exploratory\n');flush.console()
}
ci<-rbindlist(CI);h<-rbindlist(HR);ph<-rbindlist(PH);pa<-rbindlist(pairs);pa[,q_BH:=p.adjust(P,'BH'),by=contrast]
ph[,q_BH:=p.adjust(P,'BH'),by=.(basis,term)]
fwrite(ci,file.path(out,'all_test_Cindex.csv'));fwrite(h,file.path(out,'training_HR_both_scales.csv'));fwrite(h[scale=='age_residual_SD'],file.path(out,'training_HR.csv'))
fwrite(ph,file.path(out,'training_PH_both_bases.csv'));fwrite(pa,file.path(out,'paired_Cindex_differences.csv'));fwrite(rbindlist(cal),file.path(out,'training_calibration.csv'))
fwrite(rbindlist(corr),file.path(out,'test_train_adjusted_correlations.csv'));fwrite(rbindlist(piece),file.path(out,'exploratory_5yr_piecewise_HR.csv'));fwrite(rbindlist(flows),file.path(out,'cohorts.csv'))
if(length(warns))fwrite(rbindlist(warns),file.path(out,'warnings.csv'))
display<-ci[grepl('^(Display_|Accel_|Matched_Display_|Matched_Accel_)',model)]
display[,group:=ifelse(grepl('Accel_',model),'Accel','Original')];display[,score:=sub('^(Matched_Display_|Matched_Accel_|Display_|Accel_)','',model)]
display[,available:=TRUE];setnames(display,'Cindex','C');fwrite(display[!grepl('^Matched_',model)],file.path(out,'panel_d_Cindex.csv'))
fwrite(display,file.path(out,'panel_d_matched_combined_Cindex.csv'))
stopifnot(nrow(ci)==162,nrow(h)==144,nrow(ph)==180,nrow(pa)==162,nrow(display)==90)
paths<-c(scores,endpoint,file.path(root,'matched_scores.csv'),file.path(root,'02_fit_clock_endpoint.R'))
write_json(list(status='EXECUTED',train_N=31381,train_deaths=3318,test_N=13450,test_deaths=1423,reference_N=344,
  clock_models_retrained=FALSE,frozen_20261006_models_reused=TRUE,endpoint='actual blood day + cutoff2022-12-31 + death priority',
  Cindex_rows=162,HR_rows_both_scales=144,PH_rows_both_bases=180,S9_estimated_cells=90,d_original_cells=54,
  independent_direct_C_SE_checks=324,warning_count=length(warns),piecewise_status='post-hoc exploratory train-only sensitivity, not main replacement',
  input_code_sha256=setNames(vapply(paths,function(p)digest(file=p,algo='sha256'),''),basename(paths))),file.path(out,'execution_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));cat('DONE clock endpoint repair\n')
