#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages({library(data.table);library(digest)})
out <- .disco_path("results","03_multimodal_UKB/month_disease_prs_execution_20261007")
remote <- .disco_path("private","month_disease_prs_execution_20261007")
root <- .disco_path("legacy","UKB_Data")
hashes <- list(); counts <- list()
race <- fread(file.path(root,'ukb_race_21000.csv')); stopifnot(!anyDuplicated(race$eid))
for (mod in c('clinical','proteomics')) {
 p <- file.path(root,paste0('extended_covariates_',mod,'.rds')); h <- digest(file=p,algo='sha256')
 d <- as.data.table(readRDS(p)); stopifnot(!anyDuplicated(d$ID)); setnames(d,'ID','eid')
 d[,race:=race$race[match(eid,race$eid)]]
 con <- pipe(sprintf('ssh -o BatchMode=yes "${DISCO_UKB_SSH:?Set DISCO_UKB_SSH}" "umask 077; cat > %s/input_%s.csv"',remote,mod),'w')
 write.table(d,con,sep=',',row.names=FALSE,na='NA'); close(con)
 stopifnot(h==digest(file=p,algo='sha256'))
 hashes[[mod]] <- data.frame(path=p,sha256=h); counts[[mod]] <- data.frame(modality=mod,N=nrow(d))
}
fwrite(rbindlist(hashes),file.path(out,'local_input_hashes.csv')); fwrite(rbindlist(counts),file.path(out,'input_counts.csv'))
