#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages({library(jsonlite);library(digest)})
root <- .disco_path("results","")
out <- file.path(root,'03_multimodal_UKB/month_endpoint_execution_20261007')
sourcefile <- file.path(root,'08_software/20261006_checklist_actions/repos/DISCO_Code/approved/meta_functions.R')
source(sourcefile)
oldfile <- file.path(root,'02_primary_clinical/table1_sd_standardization/20261006/Table1_analysis_sample_SD.csv')
old <- read.csv(oldfile);models <- read.csv(file.path(out,'UKB_models_old_new.csv'))
for(m in c('DISCO','DM')) {
 a <- old[old$cohort=='UKB clinical'&old$metric==m,];b <- models[models$version=='old'&models$metric==m,]
 stopifnot(a$n==b$n,a$events==b$events,abs(a$adjusted_HR_per_analysis_SD-b$adjusted_HR_per_analysis_SD)<1e-9,
 abs(a$crude_C_index-b$crude_C_index)<1e-9)
}
# Retain only current estimates/metadata shared with server output, avoiding stale legacy fields.
cols <- intersect(names(old),names(models))
new <- rbind(old[old$cohort!='UKB clinical',cols],models[models$version=='new',cols])
write.csv(new,file.path(out,'Table1_updated_five_cohorts.csv'),row.names=FALSE)
inputs <- table1_meta_input(new)
inputs$endpoint_policy_version <- ifelse(inputs$cohort=='UKB clinical','UKB_actual_blood_20221231_fractional_months_20261007','other_cohort_saved_endpoint_20261006')
inputs$endpoint_status <- ifelse(inputs$cohort=='UKB clinical','Raw ZIP rebuilt; independent R QA passed','Existing aggregate; endpoint not rebuilt in this task')
write.csv(inputs,file.path(out,'Meta_inputs.csv'),row.names=FALSE)
pool <- function(t,version) {
 inp <- table1_meta_input(t)
 rows <- list()
 for(m in c('DISCO','DM')) {
  d <- inp[inp$metric==m,]
  for(left in c('none',d$cohort)) {
   s <- if(left=='none') d else d[d$cohort!=left,]
   r <- dl_pool(s$beta,s$SE,m,left);r$version <- version
   r$endpoint_policy_version <- if(version=='new') 'UKB_rebuilt_20261007_other_four_existing' else 'existing_20261006'
   rows[[length(rows)+1]] <- r
  }
 }
 do.call(rbind,rows)
}
p <- rbind(pool(old,'old'),pool(new,'new'))
write.csv(p[p$left_out=='none',],file.path(out,'Meta_old_new.csv'),row.names=FALSE)
write.csv(p[p$left_out!='none',],file.path(out,'Meta_LOO_old_new.csv'),row.names=FALSE)
compare <- merge(models[models$version=='old',],models[models$version=='new',],by=c('cohort','metric'),suffixes=c('_old','_new'))
write.csv(compare,file.path(out,'UKB_old_new_comparison.csv'),row.names=FALSE)
write.csv(data.frame(path=c(oldfile,sourcefile),sha256=vapply(c(oldfile,sourcefile),function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'meta_source_hashes.csv'),row.names=FALSE)
write_json(list(old_UKB_aggregate_reproduced=TRUE,independent_cohorts=5,
 other_four_cohorts='Existing Table1_analysis_sample_SD aggregates, unchanged',
 method='DerSimonian-Laird random effects; normal Wald CI; SE from log CI',
 scope='Association meta-analysis, not C-index pooling; no cross-metric difference test',
 limitations='Other four cohort endpoints not rebuilt by this task; crude concordance is apparent in-sample'),file.path(out,'meta_verification.json'),pretty=TRUE,auto_unbox=TRUE)
print(p[p$left_out=='none',c('metric','version','HR','lower95','upper95','I2')],row.names=FALSE)
