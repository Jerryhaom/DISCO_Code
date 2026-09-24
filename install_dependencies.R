#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
entry <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
repo <- dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE),mustWork=TRUE))
mode <- if (length(args)) args[1] else "--check"
if (!mode %in% c("--check","--install")) stop("Use --check or --install [module]")
d <- read.csv(file.path(repo,"dependencies.csv"))
if (length(args)>1) d <- d[startsWith(d$script,paste0(args[2],"/")),]
if (!nrow(d)) stop("No scripts match the requested module")
d <- unique(d[c("package","repository")])
d$installed <- vapply(d$package,requireNamespace,logical(1),quietly=TRUE)
if ("DISCO4HD" %in% d$package && requireNamespace("DISCO4HD",quietly=TRUE))
  d$installed[d$package=="DISCO4HD"] <- packageVersion("DISCO4HD") >= "0.3.0"
print(d,row.names=FALSE)
if (mode=="--check" && any(!d$installed & d$repository!="guard-only")) quit(status=1L)
if (mode=="--install") {
  missing <- d[!d$installed,]
  cran <- missing$package[missing$repository=="CRAN"]
  if (length(cran)) install.packages(cran,repos="https://cloud.r-project.org")
  bio <- missing$package[missing$repository=="Bioconductor"]
  if (length(bio)) {
    if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
    BiocManager::install(bio,ask=FALSE,update=FALSE)
  }
  if (any(missing$repository=="local")) message("Install the DISCO4HD >= 0.3.0 from Jerryhaom/DISCO4HD separately.")
}
