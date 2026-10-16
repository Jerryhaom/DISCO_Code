#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
options(stringsAsFactors=FALSE)
.libPaths(c(.disco_path("results","09_figures_tables/Rlibs"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
entry<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]));root<-dirname(entry)
old<-file.path(root,'../organ_clock_comparison/20261006_strict_70_30')
c<-fread(file.path(root,'output/clock/all_test_Cindex.csv'))
p<-fread(file.path(root,'output/clock/paired_Cindex_differences.csv'))
h<-fread(file.path(root,'output/clock/training_HR_both_scales.csv'))
ph<-fread(file.path(root,'output/clock/training_PH_both_bases.csv'))
stopifnot(nrow(c)==162,nrow(p)==162,nrow(h)==144,nrow(ph)==180,
  !anyDuplicated(c[,.(organ,model)]),!anyDuplicated(h[,.(organ,method,model,term,scale)]))
pe<-vapply(seq_len(nrow(p)),function(i){
  models<-strsplit(p$contrast[i],' minus ',fixed=TRUE)[[1]]
  v<-c[organ==p$organ[i]&model%in%models];v<-v$Cindex[match(models,v$model)]
  abs((v[1]-v[2])-p$delta_Cindex[i])
},0.)
stopifnot(max(pe)<1e-12,max(abs(p$P-2*pnorm(-abs(p$delta_Cindex/p$SE))))<1e-12,
  max(abs(p$lower95-(p$delta_Cindex-qnorm(.975)*p$SE)))<1e-12,
  max(abs(h$HR-exp(h$beta)))<1e-12,
  max(abs(h$P-2*pnorm(-abs(h$beta/h$SE))))<1e-12)
on<-fread(file.path(old,'numeric/panel_d_matched_combined_Cindex.csv'))
nn<-fread(file.path(root,'output/clock/panel_d_matched_combined_Cindex.csv'))
s<-merge(on[,.(organ,group,score,old_C=C)],nn[,.(organ,group,score,new_C=C)],by=c('organ','group','score'))
stopifnot(nrow(s)==90);s[,change:=new_C-old_C]
fwrite(s,file.path(root,'qa/S9_old_vs_new_endpoint_actual90.csv'))
pb<-ph[term%in%c('Da','Dz'),.(N=.N,nominal_PH_significant=sum(P<.05),BH_significant=sum(q_BH<.05)),by=basis]
global<-ph[term=='GLOBAL',.(N=.N,nominal_PH_significant=sum(P<.05),BH_significant=sum(q_BH<.05)),by=basis]
fair<-p[contrast%in%c('Fair_DISCO minus Fair_Oh','Fair_DISCO minus Fair_Goeminne')]
joint<-p[contrast%in%c('Fair_DISCO_Oh minus Fair_Oh','Fair_DISCO_Goeminne minus Fair_Goeminne')]
fair_summary<-fair[,.(N=.N,positive=sum(delta_Cindex>0),nominal_significant=sum(P<.05),BH_significant=sum(q_BH<.05),mean_delta=mean(delta_Cindex))]
joint_summary<-joint[,.(N=.N,positive=sum(delta_Cindex>0),nominal_significant=sum(P<.05),BH_significant=sum(q_BH<.05),mean_delta=mean(delta_Cindex))]
mc<-readRDS(file.path(old,'PRIVATE/matched_scores.rds'))
fs<-fromJSON(file.path(old,'numeric/feature_sets.json'));params<-fread(file.path(old,'numeric/DISCO_training_parameters.csv'))
selected_counts<-data.table(score_key=names(mc$selected),P=vapply(mc$selected,length,0L),reference_N=344,available=TRUE)
fwrite(selected_counts,file.path(root,'qa/matched_train_selected_protein_counts.csv'))
extra<-setdiff(params$protein,unique(unlist(fs$DISCO)))
selected_extra<-intersect(extra,unique(unlist(mc$selected)))
oldcheck<-fromJSON(file.path(old,'qa_independent_r/independent_R_checks.json'))
md<-oldcheck[grepl('matched|selected|supplement|extra|direct',oldcheck$check),]
stopifnot(all(md$pass))
phase<-fread(file.path(root,'output/clock/exploratory_5yr_piecewise_HR.csv'))
piece_summary<-unique(phase[term=='Da',.(organ,method,interaction_P)])[,.(N=.N,nominal_time_interaction_significant=sum(interaction_P<.05))]
ans<-list(status='AGGREGATE_CHECKS_PASS',C_rows=162,paired_rows=162,HR_rows_both_scales=144,PH_rows_both_bases=180,
  paired_Cindex_arithmetic_max_error=max(pe),S9_actual_recomputed_cells=90,
  S9_max_absolute_endpoint_Cindex_change=max(abs(s$change)),S9_all_cells_numerically_unchanged=all(s$change==0),
  DISCO_term_PH_by_basis=pb,global_PH_by_basis=global,fair_same_age_sex_comparison=fair_summary,
  joint_increment_over_clock=joint_summary,Brain_fair=as.data.frame(fair[organ=='Brain']),
  protein_parameter_union=nrow(params),fixed_atlas_union=length(unique(unlist(fs$DISCO))),
  candidate_union_extra_original_columns=length(extra),nonzero_selected_union_extra_original_columns=length(selected_extra),
  prior_independent_matched_direct_checks=nrow(md),prior_matched_direct_check_failures=sum(!md$pass),
  posthoc_5yr_DISCO_interaction=piece_summary,
  caveats=c('BH corrections exported by contrast family (nine organs), not across every comparison in the study',
    'PH basis counts refer to training joint models, not held-out test predictions',
    '5-year cut is post-hoc exploratory and does not prove proportional hazards',
    'All C-index CI conditional on frozen fitted pipeline; no split/retraining uncertainty',
    'This script validates aggregate arithmetic; direct C-index and paired covariance tests executed on server by02 script'))
if(file.exists(file.path(root,'output/organ/organ_disease_models.csv'))){
  om<-fread(file.path(root,'output/organ/organ_mortality_models.csv'))
  od<-fread(file.path(root,'output/organ/organ_disease_models.csv'))
  df<-fread(file.path(root,'output/organ/organ_disease_sample_flows.csv'))
  stopifnot(nrow(om)==468,nrow(od)==576,!anyDuplicated(om[,.(branch,organ,outcome,riskset)]),
    !anyDuplicated(od[,.(branch,organ,outcome,policy)]))
  stopifnot(max(abs(om$HR_unit-exp(om$beta)),na.rm=TRUE)<1e-12,
    max(abs(om$HR_SD-exp(om$beta*om$SD)),na.rm=TRUE)<1e-12,
    max(abs(od$HR_SD-exp(od$beta*od$SD)),na.rm=TRUE)<1e-12)
  ans$organ_mortality_counts<-om[,.(N=.N,warnings=sum(warning_count>0),unreliable_score=sum(!score_reliable_numeric)),by=.(branch,riskset)]
  ans$organ_disease_counts<-od[,.(N=.N,warnings=sum(warning_count>0),unreliable_score=sum(!score_reliable_numeric)),by=.(branch,policy)]
  primary_mortality<-om[riskset!='legacy_delete'];primary_disease<-od[policy=='death_first']
  stopifnot(nrow(primary_mortality)==252,nrow(primary_disease)==288)
  add_PH<-function(r,ph){
    r<-copy(r);r[,model_key:=paste(domain,branch,organ,outcome,riskset,policy,sep='|')]
    ps<-ph[term=='score',.(model_key,PH_score_P=P)];pg<-ph[term=='GLOBAL',.(model_key,PH_GLOBAL_P=P)]
    stopifnot(!anyDuplicated(ps$model_key),!anyDuplicated(pg$model_key))
    r[,PH_score_P:=ps$PH_score_P[match(model_key,ps$model_key)]]
    r[,PH_GLOBAL_P:=pg$PH_GLOBAL_P[match(model_key,pg$model_key)]]
    r[,PH_available:=is.finite(PH_score_P)&is.finite(PH_GLOBAL_P)]
    r[,PH_review:=ifelse(PH_available,PH_score_P<.05|PH_GLOBAL_P<.05,TRUE)]
    r[,PH_status:=ifelse(PH_available,'computed; nominal diagnostic',
      'not returned by cox.zph; check model warning, not a passed test')]
    r
  }
  primary_mortality<-add_PH(primary_mortality,fread(file.path(root,'output/organ/mortality_PH.csv')))
  fwrite(primary_mortality[PH_available==FALSE],file.path(root,'qa/mortality_PH_unavailable_models.csv'))
  ans$primary_mortality_PH<-primary_mortality[,.(N=.N,PH_available=sum(PH_available),
    score_nominal_significant=sum(PH_score_P<.05,na.rm=TRUE),GLOBAL_nominal_significant=sum(PH_GLOBAL_P<.05,na.rm=TRUE)),by=branch]
  if(file.exists(file.path(root,'output/disease_PH/disease_PH.csv'))){
    primary_disease<-add_PH(primary_disease,fread(file.path(root,'output/disease_PH/disease_PH.csv')))
    ans$primary_disease_PH<-primary_disease[,.(N=.N,PH_available=sum(PH_available),
      score_nominal_significant=sum(PH_score_P<.05,na.rm=TRUE),GLOBAL_nominal_significant=sum(PH_GLOBAL_P<.05,na.rm=TRUE)),by=branch]
  }
  fwrite(primary_mortality,file.path(root,'output/organ/material_primary_mortality_252.csv'))
  fwrite(primary_disease,file.path(root,'output/organ/material_primary_disease_288.csv'))
  ans$disease_age_sex_flows<-unique(df[branch=='age_sex',.(disease,input_N,covered_endpoint_N,not_in_endpoint_source,prevalent_excluded,unknown_diagnosis_date_excluded,eligible_N,primary_events,same_day_ties)])
}
write_json(ans,file.path(root,'qa/aggregate_checks.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
print(ans[c('status','S9_max_absolute_endpoint_Cindex_change','DISCO_term_PH_by_basis','fair_same_age_sex_comparison','joint_increment_over_clock','posthoc_5yr_DISCO_interaction')])
