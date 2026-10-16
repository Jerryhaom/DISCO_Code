#!/usr/bin/env Rscript
# Refit the original age+sex model and the explicitly separate legacy expanded model.
# The only endpoint change is retaining competing deaths as censored observations.
options(stringsAsFactors=FALSE,warn=1,contrasts=c("contr.treatment","contr.poly"))
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(survival);library(data.table);library(jsonlite);library(digest)})
setDTthreads(1)
args<-commandArgs(TRUE);stopifnot(length(args)==1,args%in%c("pass_1","pass_2"))
root<-Sys.getenv("DISCO_DATA_ROOT")
entry<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
base<-dirname(entry);out<-file.path(base,args);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
started<-Sys.time()
paths<-c(organ_scores=file.path(root,"Result3/Organ/UKB.odis.rds"),
  alignment=file.path(root,"Result3/Organ/UKB.prot.imputed.rds"),
  covariates=file.path(root,"UKB_Data/extended_covariates_clinical.rds"),
  race=file.path(root,"UKB_Data/ukb_race_21000.csv"),mortality=file.path(root,"UKB_Data/UKB.Mortality.csv"),
  original_plot=file.path(root,"Result3/Organ/Plot_Fig4E.R"),
  legacy_expanded=file.path(root,"Results_final/05_organ_network/outputs/organ_expanded_covariates/organ_cause_specific_expanded.csv"),
  disease_HR=file.path(root,"Result3/Organ/mortality/res.rds"),
  disease_stars=file.path(root,"Result3/Organ/mortality/pl.rds"),
  organ_order=file.path(root,"Result3/Organ/Tissue.order.rds"),script=entry)
hash<-function(p)digest(file=p,algo="sha256")
before<-vapply(paths,hash,character(1))
odis<-readRDS(paths["organ_scores"]);prot<-readRDS(paths["alignment"]);eid<-prot[,1];rm(prot);gc(FALSE)
cv<-readRDS(paths["covariates"]);race<-fread(paths["race"],colClasses="character");mort<-fread(paths["mortality"])
stopifnot(length(eid)==nrow(odis),!anyDuplicated(eid),!anyDuplicated(cv$ID),!anyDuplicated(race$eid),!anyDuplicated(mort$id))
mi<-match(eid,mort$id);ci<-match(eid,cv$ID);ri<-match(as.character(eid),race$eid)
stopifnot(!anyNA(mi))
for(n in c("age","sex","time","status")){
  z<-c(age="Age",sex="Gender",time="time",status="status")[[n]]
  stopifnot(isTRUE(all.equal(odis[[n]],mort[[z]][mi],check.attributes=FALSE)))
}
ca<-mort$cause[mi];ca[!is.na(ca)&ca==""]<-NA_character_;ca[!is.na(ca)&ca=="Diabetes mellitus"]<-"Others"
stopifnot(identical(as.character(odis$cause),as.character(ca)),all(odis$status%in%c(0,1)),
  !any(odis$status==1&is.na(ca)),!any(odis$status==0&!is.na(ca)),all(is.finite(odis$time)),all(odis$time>=0))
rr<-race$race[ri];race_imputed<-is.na(rr)|rr=="";rr[race_imputed]<-"White"
d0<-data.frame(time=odis$time,event=odis$status,cause=ca,age=odis$age,sex=factor(odis$sex),
  bmi=cv$bmi_raw[ci],race=factor(rr,levels=c("White","Non-White")),edu=factor(cv$edu_high[ci]),
  smoke=factor(cv$smoke_status[ci],levels=c(0,1,2)),alcohol=factor(cv$alcohol_freq[ci],levels=1:6,ordered=TRUE))
organs<-c("tongue","salivary.gland","stomach","adipose.tissue","choroid.plexus","retina","lung","pancreas",
  "kidney","heart.muscle","skin","esophagus","skeletal.muscle","intestine","bone.marrow","liver","brain","lymphoid.tissue")
causes<-c(Cancer="Malignant Neoplasms","Heart disease"="Heart Diseases","Cerebrovascular disease"="Cerebrovascular diseases",
  "Respiratory disease"="Respiratory system","Alzheimer's disease"="Alzheimer's disease",Others="Others")
branches<-list(age_sex=c("age","sex"),expanded_legacy=c("age","sex","bmi","race","edu","smoke","alcohol"))
stopifnot(all(organs%in%names(odis)),setequal(na.omit(unique(ca)),unname(causes)))
rows<-coefs<-warns<-support<-flows<-proofs<-ph_rows<-list();z<-qnorm(.975)
fit_one<-function(d,org,branch,outcome,method,sd_fixed,n_eligible,competing_n){
  ws<-character();err<-NA_character_;form<-reformulate(c("score",branches[[branch]]),"Surv(time,event)")
  f<-tryCatch(withCallingHandlers(coxph(form,data=d,ties="efron",na.action=na.fail,x=TRUE,y=TRUE),
    warning=function(w){ws<<-c(ws,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){err<<-conditionMessage(e);NULL})
  key<-paste(branch,org,outcome,method,sep="|")
  r<-data.frame(key=key,branch=branch,organ=org,outcome=outcome,riskset=method,N=nrow(d),events=sum(d$event),
    eligible_N=n_eligible,competing_deaths=competing_n,SD_analysis_population=sd_fixed,
    beta=NA_real_,SE=NA_real_,HR_per_unit=NA_real_,lower_per_unit=NA_real_,upper_per_unit=NA_real_,
    HR_per_SD=NA_real_,lower_per_SD=NA_real_,upper_per_SD=NA_real_,P=NA_real_,C_index=NA_real_,C_SE=NA_real_,
    warning_count=length(ws),warnings=paste(ws,collapse=" || "),error=err,model_status="FIT_ERROR",
    stable_for_display=FALSE,formula=paste(deparse(form),collapse=""),PH_verified=FALSE)
  if(!is.null(f)){
    s<-summary(f);cf<-coef(f);se<-sqrt(diag(vcov(f)));r$N<-f$n;r$events<-f$nevent
    r$beta<-unname(cf["score"]);r$SE<-unname(se["score"]);r$P<-s$coefficients["score","Pr(>|z|)"]
    r$HR_per_unit<-exp(r$beta);r$lower_per_unit<-exp(r$beta-z*r$SE);r$upper_per_unit<-exp(r$beta+z*r$SE)
    r$HR_per_SD<-exp(r$beta*sd_fixed);r$lower_per_SD<-exp((r$beta-z*r$SE)*sd_fixed);r$upper_per_SD<-exp((r$beta+z*r$SE)*sd_fixed)
    r$C_index<-unname(s$concordance[1]);r$C_SE<-unname(s$concordance[2])
    numeric_ok<-all(is.finite(cf))&&all(is.finite(se))&&!anyNA(cf)&&all(is.finite(unlist(r[c("beta","SE","HR_per_unit","HR_per_SD")])))
    r$model_status<-if(length(ws))"WARNING_REVIEW_REQUIRED" else if(!numeric_ok)"COEFFICIENT_REVIEW_REQUIRED" else "FITTED_NO_WARNING"
    r$stable_for_display<-length(ws)==0&&numeric_ok
    if(method!="legacy_delete"){
      ph_warnings<-character();ph_error<-NA_character_
      ph<-tryCatch(withCallingHandlers(cox.zph(f,transform="km",terms=TRUE,singledf=FALSE,global=TRUE),
        warning=function(w){ph_warnings<<-c(ph_warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
        error=function(e){ph_error<<-conditionMessage(e);NULL})
      t<-if(is.null(ph))matrix(c(NA_real_,NA_real_,NA_real_),nrow=1,
        dimnames=list("GLOBAL",c("chisq","df","p"))) else ph$table
      ph_rows[[length(ph_rows)+1]]<<-data.frame(key=key,branch=branch,organ=org,outcome=outcome,
        riskset=method,N=f$n,events=f$nevent,term=rownames(t),chisq=t[,"chisq"],df=t[,"df"],P=t[,"p"],
        transform="km",terms_test=TRUE,singledf=FALSE,fit_estimable=r$stable_for_display,
        diagnostic_status=if(!r$stable_for_display)"NOT_INTERPRETABLE_WARNED_FIT" else if(is.null(ph))"PH_ERROR" else if(length(ph_warnings))"PH_WARNING_REVIEW" else "DIAGNOSTIC_COMPUTED",
        warning_count=length(ph_warnings),warnings=paste(unique(ph_warnings),collapse=" || "),error=ph_error,
        endpoint_status="PROVISIONAL saved mortality dates; no final raw endpoint approval",stringsAsFactors=FALSE)
    }
    coefs[[length(coefs)+1]]<<-data.frame(key=key,branch=branch,organ=org,outcome=outcome,riskset=method,
      term_index=seq_along(cf),term=names(cf),beta=unname(cf),SE=unname(se),HR=exp(cf),lower=exp(cf-z*se),upper=exp(cf+z*se))
    if(length(ws))for(w in ws){
      ix<-integer();if(grepl("before variable",w,fixed=TRUE)){
        seg<-sub(".*before variable[[:space:]]*","",w);seg<-sub(";.*","",seg);ix<-as.integer(strsplit(trimws(seg),"[ ,]+")[[1]])
      }
      warns[[length(warns)+1]]<<-data.frame(key=key,branch=branch,organ=org,outcome=outcome,riskset=method,
        warning=w,implicated_terms=paste(names(cf)[ix],collapse=","))
    }
    # Algebraic scaling is independently checked by a fresh standardized fit for all-cause and each target in one organ.
    if(org=="brain"&&method!="legacy_delete"&&r$stable_for_display){
      dz<-d;dz$score<-(dz$score-mean(dz$score))/sd_fixed
      fz<-coxph(form,data=dz,ties="efron",na.action=na.fail);sz<-summary(fz)
      er<-max(abs(c(coef(fz)["score"]-r$beta*sd_fixed,sqrt(vcov(fz)["score","score"])-r$SE*sd_fixed,
        sz$concordance[1]-r$C_index)))
      stopifnot(er<1e-9)
      proofs[[length(proofs)+1]]<<-data.frame(branch=branch,organ=org,outcome=outcome,standardized_refit_error=er)
    }
  }
  if(outcome!="All-cause")for(v in intersect(c("sex","race","edu","smoke","alcohol"),branches[[branch]])){
    for(lev in levels(d[[v]])){ix<-as.character(d[[v]])==lev
      support[[length(support)+1]]<<-data.frame(key=key,branch=branch,organ=org,outcome=outcome,riskset=method,
        covariate=v,level=lev,N=sum(ix),target_events=sum(d$event[ix]))
    }
  }
  r
}
for(branch in names(branches))for(org in organs){
  d<-d0;d$score<-odis[[org]];need<-c("time","event","score",branches[[branch]])
  keep<-complete.cases(d[need]);eligible<-droplevels(d[keep,,drop=FALSE])
  stopifnot(all(is.finite(eligible$score)),all(is.finite(eligible$age)))
  sd_fixed<-sd(eligible$score);stopifnot(is.finite(sd_fixed),sd_fixed>0)
  flows[[length(flows)+1]]<-data.frame(branch=branch,organ=org,input_N=nrow(d),eligible_N=nrow(eligible),
    missing_covariate_excluded=sum(!keep),all_cause_deaths=sum(eligible$event),zero_followup=sum(eligible$time==0),SD=sd_fixed,
    score_sd_population="all complete eligible participants for this organ/model branch, shared across all outcomes")
  rows[[length(rows)+1]]<-fit_one(eligible,org,branch,"All-cause","all_cause",sd_fixed,nrow(eligible),0)
  for(outcome in names(causes)){
    target<-!is.na(eligible$cause)&eligible$cause==causes[[outcome]]
    competing<-eligible$event==1&!target
    legacy<-eligible[c(which(target),which(eligible$event==0)),,drop=FALSE]
    correct<-eligible;correct$event<-as.integer(target)
    stopifnot(nrow(correct)-nrow(legacy)==sum(competing),sum(correct$event)==sum(legacy$event),identical(correct$time,eligible$time))
    for(method in c("legacy_delete","competing_death_censor"))rows[[length(rows)+1]]<-
      fit_one(if(method=="legacy_delete")legacy else correct,org,branch,outcome,method,sd_fixed,nrow(eligible),sum(competing))
  }
  cat(branch,org,"fits",length(rows),"\n");flush.console()
}
r<-rbindlist(rows);stopifnot(nrow(r)==468,!anyDuplicated(r$key),!any(r$model_status=="FIT_ERROR"))
pub<-copy(r[riskset!="legacy_delete"]);stopifnot(nrow(pub)==252)
pub[stable_for_display==FALSE,c("HR_per_unit","lower_per_unit","upper_per_unit","HR_per_SD","lower_per_SD","upper_per_SD","P"):=NA_real_]
le<-r[riskset=="legacy_delete"];cr<-r[riskset=="competing_death_censor"]
delta<-merge(le,cr,by=c("branch","organ","outcome"),suffixes=c("_legacy","_correct"))
delta[,`:=`(delta_N=N_correct-N_legacy,delta_events=events_correct-events_legacy,
  HR_percent_change=100*(HR_per_unit_correct/HR_per_unit_legacy-1),
  logHR_difference=beta_correct-beta_legacy,C_index_difference=C_index_correct-C_index_legacy)]
stopifnot(nrow(delta)==216,all(delta$delta_events==0),all(delta$N_correct==delta$eligible_N_correct),
  all(delta$delta_N==delta$competing_deaths_correct))
write<-function(x,n)fwrite(x,file.path(out,n),na="NA",quote="auto")
write(r,"all_model_results_diagnostic.csv");write(pub,"corrected_model_results.csv");write(delta,"legacy_corrected_differences.csv")
write(rbindlist(flows),"analysis_sample_SD.csv");write(rbindlist(coefs),"all_coefficients.csv")
write(if(length(warns))rbindlist(warns) else data.table(key=character(),warning=character()),"model_warnings.csv")
write(rbindlist(support),"categorical_event_support.csv");write(rbindlist(proofs),"independent_SD_refit_checks.csv")
ph_tab<-rbindlist(ph_rows);stopifnot(uniqueN(ph_tab$key)==252L)
ph_tab[,q_BH_18_organs:=p.adjust(P,"BH"),by=.(branch,outcome,term)]
ph_tab[,nominal_PH_departure:=is.finite(P)&P<.05]
ph_tab[,BH_PH_departure:=is.finite(q_BH_18_organs)&q_BH_18_organs<.05]
write(ph_tab,"corrected_model_PH_diagnostics.csv")
write(ph_tab[,.(models=uniqueN(key),estimable_models=uniqueN(key[fit_estimable]),
  computed=sum(is.finite(P)),nominal_departures=sum(nominal_PH_departure),BH_departures=sum(BH_PH_departure),
  noninterpretable=sum(diagnostic_status!="DIAGNOSTIC_COMPUTED")),by=.(branch,outcome,term)],"PH_family_summary.csv")
after<-vapply(paths,hash,character(1));stopifnot(identical(before,after))
write(data.frame(role=names(paths),path=unname(paths),sha256_before=unname(before),sha256_after=unname(after),unchanged=before==after),"input_manifest.csv")
summary<-list(pass=args,original_data_N=nrow(odis),deaths=sum(odis$status),zero_followup=sum(odis$time==0),
  model_count=nrow(r),corrected_public_rows=nrow(pub),corrected_warning_models=sum(!pub$stable_for_display),
  all_warning_models=sum(r$warning_count>0),failed_models=sum(r$model_status=="FIT_ERROR"),
  extended_covariate_unmatched=sum(is.na(ci)),race_unmatched=sum(is.na(ri)),
  race_missing_or_blank_imputed_white=sum(race_imputed),all_input_hashes_unchanged=identical(before,after),
  primary="age+sex only; original Fig4a model retained",expanded="separate legacy expanded branch age+sex+bmi+race+edu+smoke+alcohol",
  riskset_change="competing deaths retained and censored at their unchanged recorded death follow-up",
  SD="within organ/model branch complete eligible population; N-1 sample SD; never target-death-only SD",
  warnings="Affected public HR/CI/P are NA; diagnostic unmasked estimates retained separately; covariates and causes not deleted or merged",
  PH_verified=FALSE,PH_diagnostic_models=uniqueN(ph_tab$key),PH_transform="km",
  PH_nominal_departure_rows=sum(ph_tab$nominal_PH_departure),PH_BH_departure_rows=sum(ph_tab$BH_PH_departure),
  PH_limits="Conditional diagnostics do not certify final endpoints or proportional hazards; warned fits remain noninterpretable",
  incident_dates_verified=FALSE,time_units="unaltered stored mortality time",endpoints_changed=FALSE,
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
write_json(summary,file.path(out,"run_metadata.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("COMPLETE",nrow(r),"fits; corrected warning models",sum(!pub$stable_for_display),"\n")
