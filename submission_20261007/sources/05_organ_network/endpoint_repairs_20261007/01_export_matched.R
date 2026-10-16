#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
setDTthreads(1)
root<-.disco_path("results","05_organ_network")
frozen<-file.path(root,'organ_clock_comparison/20261006_strict_70_30')
dest<-.disco_path("private","organ_endpoint_repairs_20261007_run01")
ds<-readRDS(file.path(frozen,'PRIVATE/DISCO_original.rds'))
mt<-readRDS(file.path(frozen,'PRIVATE/matched_scores.rds'))
d<-data.table(eid=as.character(ds$ids))
for(o in colnames(ds$scores))for(m in c('Oh','Goeminne')){
  k<-paste(m,o,sep='_');stopifnot(k%in%names(mt$matched),length(mt$matched[[k]])==nrow(d),all(is.finite(mt$matched[[k]])))
  d[[paste0('M_',m,'_',o)]]<-mt$matched[[k]]
}
stopifnot(nrow(d)==44831,!anyDuplicated(d$eid),length(ds$train)==31381,length(ds$test)==13450,length(ds$reference)==344)
con<-pipe(sprintf('ssh -S none -o BatchMode=yes "${DISCO_UKB_SSH:?Set DISCO_UKB_SSH}" "umask 077; cat > %s/matched_scores.csv"',dest),'w')
write.table(d,con,sep=',',row.names=FALSE,na='NA');close(con)
paths<-c(file.path(frozen,'PRIVATE/DISCO_original.rds'),file.path(frozen,'PRIVATE/matched_scores.rds'),file.path(frozen,'numeric/new_train_selected_proteins.json'))
write_json(list(N=nrow(d),columns=names(d)[-1],train_N=31381,test_N=13450,reference_N=344,
  source_hashes=setNames(vapply(paths,function(p)digest(file=p,algo='sha256'),''),basename(paths)),
  transfer='existing frozen scores to authorized server only; no participant data downloaded'),
  file.path(root,'endpoint_repairs_20261007/export_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Exported matched scores:',nrow(d),'participants;',ncol(d)-1,'score columns; reference344\n')
