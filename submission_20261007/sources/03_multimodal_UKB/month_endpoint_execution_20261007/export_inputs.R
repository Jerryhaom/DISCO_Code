#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages({library(digest);library(jsonlite)})
root <- .disco_path("results","")
out <- file.path(root,'03_multimodal_UKB/month_endpoint_execution_20261007')
remote <- .disco_path("private","month_endpoint_execution_20261007_0WHBMrmO")
paths <- c(.disco_path("legacy","UKB_Data/extended_covariates_clinical.rds"),
 .disco_path("legacy","UKB_Data/covariates_raw_full.rds"),
 .disco_path("legacy","UKB_Data/ukb_race_21000.csv"),
 file.path(root,'02_primary_clinical/table1_sd_standardization/prepare_approved_cohorts.R'),
 file.path(root,'90_archive/sensitivity_analyses_20260924/run_DISCO_DM_table.R'))
hash <- function() vapply(paths,function(p) digest(file=p,algo='sha256'),character(1))
before <- hash()
source(paths[4]); prep <- prepare_approved_cohorts(root)
original <- readRDS(paths[1]); d <- prep$cohorts[['UKB clinical']]$d
stopifnot(nrow(d)==nrow(original),!anyNA(original$ID),!anyDuplicated(original$ID),
 identical(d$mc,original$mc),identical(d$hd,original$hd),
 identical(d$age,original$age),identical(d$time,original$time),identical(d$status,original$status))
# Explicitly reconstruct the preparation's ID order from its original source.
env <- new.env(parent=globalenv())
for(e in parse(paths[5])) { if(is.call(e)&&identical(e[[1]],as.name('for'))) break; invisible(eval(e,env)) }
stopifnot(identical(original$ID,env$ukb$ID),identical(d$mc,env$ukb$mc),identical(d$hd,env$ukb$hd))
d <- data.frame(eid=as.character(original$ID),d,check.names=FALSE)
# Individual data goes directly through SSH stdin; no local export file.
con <- pipe(sprintf('ssh -o BatchMode=yes "${DISCO_UKB_SSH:?Set DISCO_UKB_SSH}" "umask 077; cat > %s/analysis_input.csv"',remote),'w')
write.table(d,con,sep=',',row.names=FALSE,col.names=TRUE,quote=TRUE,na='NA'); close(con)
after <- hash();stopifnot(identical(before,after))
write.csv(data.frame(path=paths,sha256=before),file.path(out,'local_input_hashes.csv'),row.names=FALSE)
write_json(list(n=nrow(d),unique_ids=length(unique(d$eid)),original_ID_order_verified=TRUE,
 scores_identical_to_original=TRUE,source_hashes_unchanged=TRUE,remote_directory=remote,
 input_columns=names(d),missing_rules=prep$audit),file.path(out,'input_protection.json'),pretty=TRUE,auto_unbox=TRUE)
