#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
suppressPackageStartupMessages({library(data.table);library(digest)})
out<-.disco_path("results","03_multimodal_UKB/targeted_endpoint_diagnostics_20261007")
files<-list.files(out,full.names=TRUE,recursive=FALSE)
files<-files[!dir.exists(files)&basename(files)!='targeted_artifact_manifest.csv']
d<-data.table(file=basename(files),bytes=file.info(files)$size,mtime=as.character(file.info(files)$mtime),
              sha256=vapply(files,function(p)digest(file=p,algo='sha256'),character(1)))
fwrite(d,file.path(out,'targeted_artifact_manifest.csv'))
cat(nrow(d),'local code/aggregate/documented-evidence files hashed; no private input downloaded\n')
