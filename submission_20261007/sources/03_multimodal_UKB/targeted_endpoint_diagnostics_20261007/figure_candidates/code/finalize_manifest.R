#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
root<-.disco_path("results","03_multimodal_UKB/targeted_endpoint_diagnostics_20261007/figure_candidates")
paths<-list.files(root,recursive=TRUE,full.names=TRUE)
paths<-paths[!dir.exists(paths)&!basename(paths)%in%c('artifact_manifest.csv')]
tab<-data.table(path=substring(paths,nchar(root)+2),bytes=file.info(paths)$size,
                mtime=as.character(file.info(paths)$mtime),
                sha256=vapply(paths,function(p)digest(file=p,algo='sha256'),''))
fwrite(tab,file.path(root,'artifact_manifest.csv'))
cat(nrow(tab),'producer/source/input/QA/figure files hashed; includes3PDFs,4PNGs.\n')
