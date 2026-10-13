#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Join the saved 396-marker combined score to the matching saved UKB covariates.
# Individual rows are streamed to a task-private server file, never to stdout.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
setDTthreads(1)
a <- commandArgs(TRUE);stopifnot(length(a)==2L)
out <- normalizePath(a[1]);remote <- a[2]
base <- .disco_path("legacy","")
paths <- c(combined=file.path(base,'Result2/UKB.ClinicalProt.DISCO.rds'),
           covariates=file.path(base,'UKB_Data/extended_covariates_proteomics.rds'),
           pure_protein=file.path(base,'Result2/UKB.Prot.DISCO.rds'),
           race=file.path(base,'UKB_Data/ukb_race_21000.csv'),
           combined_features=file.path(base,'Result2/Selected.UKB.Clinicalprot.rds'),
           pure_features=file.path(base,'Result2/Selected.UKB.Prot.rds'))
hash_before <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
cp <- as.data.table(readRDS(paths['combined']));p <- as.data.table(readRDS(paths['covariates']))
pure <- as.data.table(readRDS(paths['pure_protein']));race <- fread(paths['race'])
stopifnot(!anyDuplicated(cp$ID),!anyDuplicated(p$ID),!anyDuplicated(pure$ID),!anyDuplicated(race$eid))
j <- match(p$ID,pure$ID);stopifnot(!anyNA(j),all(p$mc==pure$mc[j]),all(p$hd==pure$hd[j]))
ix <- match(cp$ID,p$ID);stopifnot(!anyNA(ix),all(cp$age==p$age[ix]),
                              all(as.character(cp$sex)==as.character(p$sex[ix])))
d <- copy(p[ix]);d[,`:=`(mc=cp$mc,hd=cp$hd)]
setnames(d,'ID','eid');d[,race:=race$race[match(eid,race$eid)]]
stopifnot(nrow(d)==47173L,all(is.finite(d$mc)),all(is.finite(d$hd)),all(d$hd>0))
private_path <- file.path(remote,'input_clinicalprot.csv')
# Exclusive creation prevents replacing a private input after an unknown outcome.
py <- paste0("import sys,os,shutil; os.umask(0o077); f=open(",
             encodeString(private_path,quote='"'),",'xb'); shutil.copyfileobj(sys.stdin.buffer,f); f.close()")
cmd <- paste('ssh -o BatchMode=yes "${DISCO_UKB_SSH:?Set DISCO_UKB_SSH}"',shQuote(paste('python3 -c',shQuote(py))))
con <- pipe(cmd,'w');write.table(d,con,sep=',',row.names=FALSE,na='NA');rc <- close(con)
stopifnot(rc==0L)
hash_after <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
stopifnot(identical(hash_before,hash_after))
q <- list(status='EXECUTED',modality='ClinicalProt',score_N=nrow(cp),covariate_join_N=nrow(d),
          score_source=paths['combined'],covariate_source=paths['covariates'],
          score_columns=list(DISCO='saved mc',DM='saved raw hd'),
          selected_features_combined=length(readRDS(paths['combined_features'])),
          selected_features_pure_protein=ncol(readRDS(paths['pure_features']))-1L,
          pure_proteomics_identity='extended_covariates_proteomics mc/hd equal UKB.Prot.DISCO exactly in all 47859 participants',
          demographic_join='all 47173 age and sex exactly matched by ID',
          weights_or_scores_retrained=FALSE,individual_export='task-private SSH stream, not downloaded',
          private_server_path=private_path,source_hashes=as.list(hash_before),sources_unchanged=TRUE)
write_json(q,file.path(out,'ClinicalProt_input_QA.json'),pretty=TRUE,auto_unbox=TRUE)
print(q[c('status','modality','score_N','covariate_join_N','selected_features_combined','selected_features_pure_protein')])
