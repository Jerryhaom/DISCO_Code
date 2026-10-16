#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1);Sys.umask('0077')
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1);args<-commandArgs(TRUE);stopifnot(length(args)==1);work<-normalizePath(args[1])
score_path<-Sys.getenv('DISCO_ORGAN_INPUT');death_path<-Sys.getenv('DISCO_MORTALITY_ENDPOINT');disease_path<-Sys.getenv('DISCO_DISEASE_ENDPOINT');patch_path<-Sys.getenv('DISCO_ORGAN_DISEASE_PATCH')
paths<-c(score_path,death_path,disease_path,patch_path);hashes<-setNames(vapply(paths,function(p)digest(file=p,algo='sha256'),''),basename(paths));s<-fread(score_path,colClasses=c(eid='character'));ep<-fread(death_path,colClasses=c(eid='character'));stopifnot(!anyDuplicated(s$eid),!anyDuplicated(ep$eid));j<-match(s$eid,ep$eid);stopifnot(!anyNA(j));org<-c('tongue','salivary.gland','stomach','adipose.tissue','choroid.plexus','retina','lung','pancreas','kidney','heart.muscle','skin','esophagus','skeletal.muscle','intestine','bone.marrow','liver','brain','lymphoid.tissue');stopifnot(all(org%in%names(s)));warnings<-list();rows<-list()
one<-function(score,time,event,organ,outcome) {
  ok<-is.finite(score)&is.finite(time)&time>0&!is.na(event);z<-data.frame(score=score[ok],time=time[ok],event=event[ok]);ww<-character();fit<-withCallingHandlers(coxph(Surv(time,event)~score,data=z,ties='efron'),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')});ss<-summary(fit);if(length(ww))warnings[[length(warnings)+1]]<<-data.frame(organ=organ,outcome=outcome,message=unique(ww));data.frame(organ=organ,outcome=outcome,N=fit$n,events=fit$nevent,C=unname(ss$concordance[1]),C_SE=unname(ss$concordance[2]),beta=unname(coef(fit)[1]),P=ss$coefficients[1,5],warning_count=length(ww),estimand='apparent crude score-only Cox Harrell C')
}
for(o in org)rows[[length(rows)+1]]<-one(s[[o]],ep$time[j],ep$event[j],o,'All-cause mortality')
e<-rbindlist(list(fread(disease_path,colClasses=c(eid='character')),fread(patch_path,colClasses=c(eid='character'))),use.names=TRUE);e<-e[eid%in%s$eid];stopifnot(!anyDuplicated(e[,.(eid,disease)]));dis<-c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')
for(d in dis) {
  ee<-e[disease==d];ix<-match(s$eid,ee$eid);stopifnot(!anyNA(ix));reason<-ee$reason[ix];allowed<-is.na(reason)|reason=='';tt<-ee$time[ix];event<-as.integer(ee$event[ix]==1L);tt[!allowed]<-NA_real_
  for(o in org)rows[[length(rows)+1]]<-one(s[[o]],tt,event,o,d)
  cat('Disease target complete:',d,'\n');flush.console()
}
r<-rbindlist(rows);stopifnot(nrow(r)==18*9,all(is.finite(r$C)));fwrite(r,file.path(work,'current_target_models_162.csv'));if(length(warnings))fwrite(rbindlist(warnings),file.path(work,'model_warnings.csv'))
write_json(list(status='PASS',model_rows=nrow(r),organ_n=18,disease_n=8,N_scores=nrow(s),model_warning_count=length(warnings),endpoint='Approved actual blood date / 2022-12-31 cutoff / death-first; prevalent and known-unknown diagnosis excluded; all51922 IDs explicitly covered',input_sha256=as.list(hashes),input_unchanged=identical(hashes,setNames(vapply(paths,function(p)digest(file=p,algo='sha256'),''),basename(paths))),private_individual_records_remain_on_server=TRUE,estimand='Crude apparent C from each original organ definition (all atlas organ genes), not age+sex adjusted predictive C or held-out performance'),file.path(work,'current_targets_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17);writeLines(capture.output(sessionInfo()),file.path(work,'sessionInfo.txt'));cat('COMPLETE162targets\n')
