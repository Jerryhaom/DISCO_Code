#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
.libPaths(c(.disco_path("results","09_figures_tables/Rlibs"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
entry<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]));root<-dirname(entry)
files<-c(list.files(root,pattern='[.]R$',full.names=TRUE),file.path(root,'README.txt'),
  list.files(file.path(root,'output'),full.names=TRUE,recursive=TRUE),
  list.files(file.path(root,'panels'),pattern='[.](pdf|png)$',full.names=TRUE),file.path(root,'Fig4c_remaining.txt'),
  list.files(file.path(root,'qa'),pattern='[.](json|csv|txt|log|png)$',full.names=TRUE))
files<-sort(unique(files[file.exists(files)&!file.info(files)$isdir]))
stopifnot(!any(grepl('PRIVATE|private_',basename(files))),length(list.files(file.path(root,'panels'),pattern='[.]pdf$'))==8)
tab<-data.table(path=sub(paste0('^',root,'/'),'',files),bytes=file.info(files)$size,
  sha256=vapply(files,function(f)digest(file=f,algo='sha256'),''))
fwrite(tab,file.path(root,'manifest_files_sha256.csv'))
clock<-fromJSON(file.path(root,'output/clock/execution_QA.json'))
organ<-fromJSON(file.path(root,'output/organ/execution_QA.json'))
qa<-fromJSON(file.path(root,'qa/aggregate_checks.json'))
input_map<-rbindlist(list(
  data.table(run='clock_run01',role=c('frozen_scores','new_death_endpoint','matched_scores','executed_code'),
    path=c(.disco_path("private","month_7030_execution_20261007/scores_7030.csv"),.disco_path("private","month_organ_execution_20261007/private_organ_endpoint.csv"),
      .disco_path("private","organ_endpoint_repairs_20261007_run01/matched_scores.csv"),'02_fit_clock_endpoint.R'),sha256=unlist(clock$input_code_sha256,use.names=FALSE)),
  data.table(run='final_organ_run02',role=c('frozen_organ_scores','new_death_endpoint','source8_endpoint','52_extra_endpoint','executed_code'),
    path=c(.disco_path("private","month_organ_execution_20261007/input_organ.csv"),.disco_path("private","month_organ_execution_20261007/private_organ_endpoint.csv"),
      .disco_path("private","month_disease_prs_execution_20261007/private_endpoint_8.csv"),
      .disco_path("private","targeted_endpoint_diagnostics_20261007.JaNZnJ/private_organ_52_endpoint_8.csv"),'03b_fit_organ_full_coverage.R'),
    sha256=unlist(organ$input_code_sha256,use.names=FALSE))))
stopifnot(digest(file=file.path(root,'02_fit_clock_endpoint.R'),algo='sha256')==input_map[run=='clock_run01'&role=='executed_code']$sha256,
  digest(file=file.path(root,'03b_fit_organ_full_coverage.R'),algo='sha256')==input_map[run=='final_organ_run02'&role=='executed_code']$sha256)
fwrite(input_map,file.path(root,'source_input_sha256_map.csv'))
raceqa<-fromJSON(file.path(root,'output/race_diagnostic/execution_QA.json'))
diseaseqa<-fromJSON(file.path(root,'output/disease_PH/execution_QA.json'))
stopifnot(digest(file=file.path(root,'07_race_stratified_diagnostic.R'),algo='sha256')==raceqa$source_sha256,
  digest(file=file.path(root,'09_disease_PH_diagnostic.R'),algo='sha256')==diseaseqa$source_sha256,
  digest(file=file.path(root,'05_plot_Fig4a_new_endpoint.R'),algo='sha256')=='b366872b14af61f03ae5fd7a0b1d6486e3810cb6d6041012d9d0912401ab7d22')
files<-sort(unique(c(files,file.path(root,'source_input_sha256_map.csv'),file.path(root,'export_manifest.json'))))
tab<-data.table(path=sub(paste0('^',root,'/'),'',files),bytes=file.info(files)$size,
  sha256=vapply(files,function(f)digest(file=f,algo='sha256'),''))
fwrite(tab,file.path(root,'manifest_files_sha256.csv'))
stopifnot(clock$Cindex_rows==162,organ$mortality_rows==468,organ$disease_rows==576,
  all(qa$disease_age_sex_flows$not_in_endpoint_source==0),file.exists(file.path(root,'qa/visual_review.txt')))
write_json(list(status='EXECUTED_TECHNICAL_REPAIRS; scientific interpretation and manuscript synchronization still required',
  date='2026-10-07',scope='05_organ_network; old source/output preserved; no individual-level records downloaded',
  server_clock_run=.disco_path("private","organ_endpoint_repairs_20261007_run01"),
  server_final_organ_run=.disco_path("private","organ_endpoint_repairs_20261007_run02"),
  final_mortality_table='output/organ/material_primary_mortality_252.csv',
  final_disease_table='output/organ/material_primary_disease_288.csv',
  organ_death_population=51922,organ_allcause_deaths=5615,disease_endpoint_covered_population=51922,
  clock_train_N=31381,clock_train_deaths=3318,clock_test_N=13450,clock_test_deaths=1423,clock_reference_N=344,
  clock_pipeline='Frozen 20261006 train-only scores; new endpoint Cox refits; new training age calibration; held-out C evaluation. Clock algorithms not retrained in this endpoint repair run.',
  primary_plot_HR_basis='ED10 b/c training age-residual SD; score+age+sex, before/after same-organ alternative score',
  raw_HR_alternative='output/clock/training_HR_both_scales.csv raw_score_SD rows; do not mix with residual-SD PH',
  Fig4a_primary='panels/Fig4a_age_sex_NEW_ENDPOINT_per_unit.pdf; original raw-score per-unit age+sex design',
  Fig4a_separate_SD_candidate='panels/Fig4a_age_sex_NEW_ENDPOINT_per_SD.pdf; common mortality population SD before disease exclusions',
  S9_actual_recomputed_cells=90,ED10b_points=36,ED10c_points=36,ED10d_cells=54,
  PH_training_joint_residual_basis_DISCO_significant=18,PH_training_joint_raw_basis_DISCO_significant=16,
  expanded_Alzheimer_warning_models=36,race_stratified_diagnostic_models=18,race_stratified_diagnostic_warnings=0,
  disease_PH_models=288,disease_PH_available=288,disease_PH_score_nominal_significant=180,
  mortality_PH_missing_expanded_Alzheimer_models=18,
  outstanding=c('PH violations remain; do not interpret b/c as time-invariant causal/constant effects',
    '5-year segmented model is post-hoc train-only sensitivity, not a primary replacement',
    'Display Original comparison gives DISCO age but raw ProtAge no explicit age; fairness claims must use same age+sex models',
    'Strict split is internal and previously seen at study level; not external replication of original authors fitted models',
    'Conditional IJ CI excludes repeated split, clock retraining and imputation uncertainty',
    'Old KM draft lacks verified producer; not silently declared reproduced',
    'Submission composite/legends/tables/full text must be synchronized by main task; this run does not edit those files',
    'ED10/S9 PDFs retain original Helvetica PDF device with unembedded core fonts; final submission font embedding still needs checking'),
  input_hash_index='source_input_sha256_map.csv',
  frozen_52_patch_server_backup=.disco_path("private","organ_endpoint_repairs_20261007_run02/private_organ_52_endpoint_8_frozen.csv"),
  file_hash_index='manifest_files_sha256.csv'),file.path(root,'manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('Manifest created:',nrow(tab),'files; eight PDFs; final organ full coverage\n')
