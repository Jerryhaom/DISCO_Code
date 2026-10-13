#!/usr/bin/env Rscript
# Verify interval estimand, propagate aggregate-only Meta, and audit original-panel groups.
Sys.umask('0077')
root <- Sys.getenv('DISCO_RESULTS_INPUT')
out <- Sys.getenv('DISCO_RUN_MODULE')
.libPaths(c(file.path(out,'rlib'),.libPaths()))
suppressPackageStartupMessages({library(icenReg);library(survival);library(jsonlite);library(digest)})
primary <- read.csv(file.path(out,'results/primary_interval_PH_results.csv'))
repeat_run <- read.csv(file.path(out,'results/repeat_interval_PH_results.csv'))
numeric_fields <- c('N','deaths','primary_SD','beta','SE','HR','lower95','upper95','P')
dif <- max(abs(as.matrix(primary[,numeric_fields])-as.matrix(repeat_run[,numeric_fields])),na.rm=TRUE)
stopifnot(dif==0)
qa <- list();full <- readRDS(file.path(out,'private/input/wave_interval_endpoint_frame.rds'))
d <- full[!full$no_postbaseline_followup,];d$age<-d$age-mean(full$age)
for(metric in c('DISCO','DM')){
 x<-readRDS(file.path(out,'private/models',paste0('primary_',metric,'_interval_PH.rds')))
 y<-readRDS(file.path(out,'private/models',paste0('repeat_',metric,'_interval_PH.rds')))
 stopifnot(nrow(x$bsMat)==200,all(is.finite(x$bsMat)),max(abs(x$bsMat-y$bsMat))==0,all(eigen(vcov(x),symmetric=TRUE)$values>0))
 raw <- if(metric=='DISCO')d$mc else d$log_hd
 primary_raw <- if(metric=='DISCO')full$mc else full$log_hd
 d$raw_score<-raw;d$score_z<-(raw-mean(primary_raw))/sd(primary_raw)
 rawfit<-ic_sp(Surv(lower_nominal,upper_nominal,type='interval2')~raw_score+age+sex+smoke+drink+education,data=d,model='ph',bs_samples=0)
 # Monotone wave recoding: the semiparametric PH likelihood depends on inspection order,
 # not an invented uniform number of calendar years for each person's actual visits.
 d$lower_order<-match(d$lower_nominal,c(0,2,4,7))-1
 d$upper_order<-match(d$upper_nominal,c(0,2,4,7))-1
 orderfit<-ic_sp(Surv(lower_order,upper_order,type='interval2')~score_z+age+sex+smoke+drink+education,data=d,model='ph',bs_samples=0)
 raw_error<-abs(coef(rawfit)['raw_score']*sd(primary_raw)-coef(x)['score_z'])
 order_error<-max(abs(coef(orderfit)-coef(x)))
 stopifnot(raw_error<1e-6,order_error<1e-9,x$iterations<10000)
 qa[[metric]]<-data.frame(metric,point_and_bootstrap_repeat_max_diff=max(abs(x$bsMat-y$bsMat)),bootstrap_N=nrow(x$bsMat),
  covariance_positive_definite=TRUE,iterations=x$iterations,log_likelihood=x$llk,raw_score_scale_beta_error=raw_error,
  monotone_wave_order_beta_error=order_error,ordinary_Harrell_C_not_computed_for_interval_endpoint=TRUE)
}
write.csv(do.call(rbind,qa),file.path(out,'results/interval_model_numerical_QA.csv'),row.names=FALSE)
old_path<-file.path(root,'04_external_validation/CHARLS_followup_sensitivity_20261007/results/primary_CHARLS_followup_sensitivity.csv')
old<-read.csv(old_path)
summary_rows <- list()
for(metric in c('DISCO','DM')){
 for(state in c('primary_original','exclude_unverified_2018_alive')){
  s<-old[old$metric==metric&old$state==state,]
  summary_rows[[length(summary_rows)+1L]]<-data.frame(metric,analysis=state,N=s$N,deaths=s$deaths,
   HR=s$adjusted_HR_per_primary_SD,lower95=s$adjusted_CI95_lower,upper95=s$adjusted_CI95_upper,P=s$adjusted_P,
   primary_SD=s$fixed_primary_SD,ordinary_Harrell_C=s$crude_C_index,
   endpoint='Legacy wave-assigned 2/4/6-year deaths, 8-year censored; retained only as historical comparison')
 }
 for(method in unique(primary$method)){
  s<-primary[primary$metric==metric&primary$method==method,]
  summary_rows[[length(summary_rows)+1L]]<-data.frame(metric,analysis=method,N=s$N,deaths=s$deaths,HR=s$HR,lower95=s$lower95,
    upper95=s$upper95,P=s$P,primary_SD=s$primary_SD,ordinary_Harrell_C=NA_real_,endpoint=s$endpoint)
 }
}
comparison<-do.call(rbind,summary_rows)
retro<-read.csv(file.path(out,'results/retrospective_primary_interval_PH_results.csv'))
retro_repeat<-read.csv(file.path(out,'results/retrospective_repeat_interval_PH_results.csv'))
retro_diff<-max(abs(as.matrix(retro[,numeric_fields])-as.matrix(retro_repeat[,numeric_fields])),na.rm=TRUE)
stopifnot(retro_diff==0)
for(metric in c('DISCO','DM')){
 x<-readRDS(file.path(out,'private/models',paste0('retrospective_primary_',metric,'_interval_PH.rds')))
 y<-readRDS(file.path(out,'private/models',paste0('retrospective_repeat_',metric,'_interval_PH.rds')))
 stopifnot(max(abs(x$bsMat-y$bsMat))==0,all(is.finite(x$bsMat)),nrow(x$bsMat)==200)
}
for(i in seq_len(nrow(retro))){s<-retro[i,];comparison<-rbind(comparison,data.frame(metric=s$metric,
 analysis=paste0('retrospective_',s$method),N=s$N,deaths=s$deaths,HR=s$HR,lower95=s$lower95,upper95=s$upper95,
 P=s$P,primary_SD=s$primary_SD,ordinary_Harrell_C=NA_real_,
 endpoint='Observed wave interval plus10 explicit solar deaths before2018 from2020 exit; 1 date/alive conflict excluded'))}
write.csv(comparison,file.path(out,'tables/CHARLS_endpoint_sensitivity_comparison.csv'),row.names=FALSE)
# Freeze the other four current cohorts; replace CHARLS only in a separately labelled Meta sensitivity.
meta_fun<-file.path(Sys.getenv('DISCO_RELEASE_ROOT'),'../approved/meta_functions.R')
meta_input<-file.path(root,'03_multimodal_UKB/month_endpoint_execution_20261007/Meta_inputs.csv')
original<-read.csv(meta_input,stringsAsFactors=FALSE);changed<-original
source(meta_fun)
for(metric in c('DISCO','DM')){
 i<-which(changed$cohort=='CHARLS'&changed$metric==metric)
 s<-primary[primary$metric==metric&primary$method=='semiparametric_interval_PH',]
 stopifnot(length(i)==1,nrow(s)==1,abs(changed$score_SD[i]-s$primary_SD)<1e-12)
 changed$n[i]<-s$N;changed$events[i]<-s$deaths;changed$HR[i]<-s$HR
 changed$lower95[i]<-s$lower95;changed$upper95[i]<-s$upper95;changed$beta[i]<-s$beta;changed$SE[i]<-s$SE
 changed$model_formula[i]<-'Interval-censored semiparametric PH ~ score_z_fixed + age + sex + smoke + drink + education'
 changed$SD_population[i]<-'Original common CHARLS Table1 population N=11601; sample SD fixed'
 changed$endpoint_policy_version[i]<-'CHARLS_observed_wave_interval_PH_20261008'
 changed$endpoint_status[i]<-'Actual observed alive/dead wave intervals; no precise individual calendar-time claim'
}
stopifnot(identical(changed[changed$cohort!='CHARLS',],original[original$cohort!='CHARLS',]))
retrospective_changed<-changed
for(metric in c('DISCO','DM')){
 i<-which(retrospective_changed$cohort=='CHARLS'&retrospective_changed$metric==metric)
 s<-retro[retro$metric==metric&retro$method=='semiparametric_interval_PH',]
 retrospective_changed$n[i]<-s$N;retrospective_changed$events[i]<-s$deaths;retrospective_changed$HR[i]<-s$HR
 retrospective_changed$lower95[i]<-s$lower95;retrospective_changed$upper95[i]<-s$upper95
 retrospective_changed$beta[i]<-s$beta;retrospective_changed$SE[i]<-s$SE
 retrospective_changed$endpoint_policy_version[i]<-'CHARLS_strict_solar_retrospective_interval_PH_20261008'
 retrospective_changed$endpoint_status[i]<-'10 unambiguous Gregorian retrospective deaths before2018; 1 date/alive conflict isolated'
}
stopifnot(identical(retrospective_changed[retrospective_changed$cohort!='CHARLS',],original[original$cohort!='CHARLS',]))
pool_rows<-list()
for(metric in c('DISCO','DM'))for(left in c('none',unique(original$cohort)))for(state in c('primary_current','CHARLS_interval_PH_sensitivity','CHARLS_retrospective_interval_PH_sensitivity')){
 v<-switch(state,primary_current=original,CHARLS_interval_PH_sensitivity=changed,CHARLS_retrospective_interval_PH_sensitivity=retrospective_changed)
 v<-v[v$metric==metric&if(left=='none')TRUE else v$cohort!=left,]
 p<-dl_pool(v$beta,v$SE,metric,left);p$state<-state
 p$estimand<-'Association per original cohort SD; CHARLS interval PH sensitivity uses bootstrap SE, other four cohorts unchanged'
 pool_rows[[length(pool_rows)+1L]]<-p
}
pool<-do.call(rbind,pool_rows)
base<-pool[pool$state=='primary_current',];new<-pool[pool$state!='primary_current',]
delta<-merge(base,new,by=c('metric','left_out'),suffixes=c('_primary','_interval'))
delta$HR_difference<-delta$HR_interval-delta$HR_primary
stopifnot(all(delta$HR_difference[delta$left_out=='CHARLS']==0))
original$state<-'primary_current';changed$state<-'CHARLS_interval_PH_sensitivity';retrospective_changed$state<-'CHARLS_retrospective_interval_PH_sensitivity'
write.csv(rbind(original,changed,retrospective_changed),file.path(out,'tables/Meta_inputs_primary_interval.csv'),row.names=FALSE)
write.csv(pool[pool$left_out=='none',],file.path(out,'tables/Meta_pooled_primary_interval.csv'),row.names=FALSE)
write.csv(pool[pool$left_out!='none',],file.path(out,'tables/Meta_LOO_primary_interval.csv'),row.names=FALSE)
write.csv(delta,file.path(out,'tables/Meta_paired_differences.csv'),row.names=FALSE)
# Reproduce actual original group-selection semantics, not guessed figure file names.
ed8source<-file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/CHARLS.2011.disco.rds')
ed3source<-file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result4/CHARLS.DISCO.rds')
ed8<-readRDS(ed8source);ed3<-readRDS(ed3source)
m8<-match(full$id,ed8$id);m3<-match(full$id,ed3$id)
stopifnot(!anyNA(m8),!anyNA(m3),max(abs(full$mc-ed8$mc[m8]))==0,max(abs(full$mc-ed3$mc[m3]))==0,
          max(abs(full$hd-ed8$hd[m8]))==0,max(abs(full$time-ed8$time[m8]))==0,
          max(abs(full$status-ed8$status[m8]))==0,max(abs(full$time-ed3$time[m3]))==0,
          max(abs(full$status-ed3$status[m3]))==0)
full$fl<-ed8$fl[m8]
full$g_ED8<-as.numeric(cut(full$mc,breaks=quantile(ed8$mc))) # original include.lowest default FALSE.
ed8_selected<-full[!is.na(full$fl)&full$g_ED8%in%c(1,4),]
ed8_selected$group<-paste(ed8_selected$fl,ed8_selected$g_ED8)
ed8_selected$group<-factor(ed8_selected$group,levels=names(table(ed8_selected$group))[c(2,1,4,3,6,5)])
ed3_groups<-list()
for(metric in c('DISCO','DM')){
 v<-if(metric=='DISCO')full$mc else full$hd
 source_v<-if(metric=='DISCO')ed3$mc else ed3$hd
 g<-as.numeric(cut(v,breaks=quantile(source_v,probs=seq(0,1,.25)),include.lowest=TRUE,right=TRUE))
 u<-full[g%in%c(1,4),];u$group<-paste(metric,ifelse(g[g%in%c(1,4)]==1,'Low Risk','High Risk'))
 ed3_groups[[metric]]<-u
}
ed3_selected<-do.call(rbind,ed3_groups)
aggregate_groups<-function(data,panel){
 rows<-list()
 for(g in unique(as.character(data$group))){
   u<-data[as.character(data$group)==g,]
   for(state in c('historical_wave_assigned','observed_wave_interval')){
    z<-if(state=='historical_wave_assigned')u else u[!u$no_postbaseline_followup,]
    rows[[length(rows)+1L]]<-data.frame(panel,group=g,state,N=nrow(z),deaths=sum(z$status),sensitivity_excluded_N=nrow(u)-nrow(z),
     confirmed_alive2013=sum(z$alive13),confirmed_alive2015=sum(z$alive15),confirmed_alive2018=sum(z$alive18))
   }
 }
 do.call(rbind,rows)
}
groups<-rbind(aggregate_groups(ed3_selected,'ED3b_CHARLS'),aggregate_groups(ed8_selected,'ED8b_CHARLS'))
retrospective_full<-readRDS(file.path(out,'private/input/retrospective_strict_interval_frame.rds'))
augment_retrospective_groups<-function(data,panel){
 j<-match(data$id,retrospective_full$id);stopifnot(!anyNA(j))
 data$status<-retrospective_full$status[j]
 data$no_postbaseline_followup<-retrospective_full$analysis_excluded[j]
 t<-aggregate_groups(data,panel);t<-t[t$state=='observed_wave_interval',]
 t$state<-'observed_wave_interval_2020_retrospective';t
}
groups<-rbind(groups,augment_retrospective_groups(ed3_selected,'ED3b_CHARLS'),augment_retrospective_groups(ed8_selected,'ED8b_CHARLS'))
write.csv(groups,file.path(out,'tables/ED3b_ED8b_group_endpoint_propagation.csv'),row.names=FALSE)
stopifnot(nrow(ed8_selected)==5472L,sum(ed8_selected$status)==654L)
saveRDS(list(ED3b=ed3_selected,ED8b=ed8_selected),file.path(out,'private/input/panel_interval_group_frames.rds'))
Sys.chmod(file.path(out,'private/input/panel_interval_group_frames.rds'),'0600')
manifest_paths<-c(old_path,meta_fun,meta_input,ed8source,ed3source,
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result4/Validate_Charls.R'),file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Asso_DISCO_FI.R'),
 file.path(root,'00_manifest/parallel_completion_20261007/current_candidates.json'))
write.csv(data.frame(path=manifest_paths,sha256=vapply(manifest_paths,function(p)digest(file=p,algo='sha256'),character(1))),
 file.path(out,'input/propagation_source_manifest.csv'),row.names=FALSE)
write_json(list(status='INTERVAL_PH_REPEAT_AND_PROPAGATION_QA_PASS',independent_R_process_numeric_max_diff=dif,
  bootstrap_repeated_exactly=TRUE,interval_model_count_per_process=6,additional_raw_scale_and_wave_order_check_count=4,
  strict_retrospective_independent_R_process_numeric_max_diff=retro_diff,strict_retrospective_bootstrap_repeated_exactly=TRUE,
  approved_primary_SD_and_covariates_held_fixed=TRUE,other_four_Meta_cohort_inputs_unchanged=TRUE,
  original_ED3b_and_ED8b_objects_same_endpoints_and_scores_verified_by_ID=TRUE,
  ED8b_historical_selected_N=nrow(ed8_selected),ED8b_historical_deaths=sum(ed8_selected$status),
  interval_endpoints_are_wave_order_not_precision_calendar_followup=TRUE,
  original_Table1_Meta_ED3b_ED8b_and_CURRENT_untouched=TRUE,
  no_ordinary_Harrell_C_claim_for_interval_endpoints=TRUE),file.path(out,'results/propagation_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
print(pool[pool$left_out=='none',c('metric','state','HR','lower95','upper95','I2')],row.names=FALSE)
print(groups)
