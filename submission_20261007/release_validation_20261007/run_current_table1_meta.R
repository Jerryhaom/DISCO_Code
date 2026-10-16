#!/usr/bin/env Rscript
# Current locked Table 1 aggregate entry; no HR literals or participant fits.
entry <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
here <- dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE),mustWork=TRUE))
repo <- dirname(dirname(here))
suppressPackageStartupMessages({library(digest);library(jsonlite)})
source(file.path(repo,'approved','meta_functions.R'))
args <- commandArgs(TRUE)
input <- if(length(args)>=1L)normalizePath(args[1],mustWork=TRUE) else file.path(dirname(here),'aggregate','Table1_updated_five_cohorts.csv')
out <- if(length(args)>=2L)args[2] else file.path(getwd(),'current_Table1_Meta_output')
if(file.exists(out))stop('Choose a new output directory; existing output is not overwritten')
contract <- fromJSON(file.path(here,'data','current_aggregate_contract.json'))
stopifnot(digest(file=input,algo='sha256')==contract$table1_sha256)
raw <- read.csv(input,stringsAsFactors=FALSE,check.names=FALSE)
d <- table1_meta_input(raw)
# The generic old estimator is retained, but its old endpoint metadata is not.
d$endpoint_policy_version <- ifelse(d$cohort=='UKB clinical','UKB_actual_blood_20221231_fractional_months_20261007','other_four_saved_cohort_endpoint_20261006')
d$endpoint_status <- ifelse(d$cohort=='UKB clinical','Locked 20261007 aggregate; raw UKB endpoint reconstruction audited separately','Retained Table1 aggregate; raw endpoint not rebuilt by this public entry')
rows <- list()
for(metric in c('DISCO','DM')) {
  x <- d[d$metric==metric,]
  for(left in c('none',x$cohort)) {
    keep <- if(left=='none')x else x[x$cohort!=left,]
    s <- dl_pool(keep$beta,keep$SE,metric,left)
    s$endpoint_policy_version <- 'UKB_locked_20261007_other_four_retained_20261006'
    rows[[length(rows)+1L]] <- s
  }
}
pooled <- do.call(rbind,rows)
dir.create(out,recursive=TRUE)
write.csv(raw,file.path(out,'Table1_current_full_precision.csv'),row.names=FALSE)
write.csv(d,file.path(out,'Meta_inputs_current.csv'),row.names=FALSE)
write.csv(pooled[pooled$left_out=='none',],file.path(out,'Meta_current.csv'),row.names=FALSE)
write.csv(pooled[pooled$left_out!='none',],file.path(out,'Meta_LOO_current.csv'),row.names=FALSE)
write_json(list(source_sha256=contract$table1_sha256,model_rows=10L,independent_cohorts=5L,
  HR_scale='DISCO or natural-log DM per effective within-cohort analysis-sample SD',
  concordance='Crude apparent score-only C-index; not pooled',
  scope='Portable aggregate replay and DerSimonian-Laird/Wald meta-analysis; no raw cohort fits',
  endpoint_limitation='Other four cohort endpoints are retained; UKB observation-window limitations are inherited'),
  file.path(out,'audit_status.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(pooled[pooled$left_out=='none',c('metric','HR','lower95','upper95','I2')],row.names=FALSE)
