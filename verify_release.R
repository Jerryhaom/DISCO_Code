#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
repo <- dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE), mustWork=TRUE))
files <- list.files(repo, "[.]R$", full.names=TRUE, recursive=TRUE)
files <- files[!grepl("/(data|results)/", files)]
for (f in files) {
  parse(f)
  lines <- readLines(f, warn=FALSE)
  if (basename(f) != "verify_release.R" && any(grepl("/Users/|~/Desktop/", lines))) stop("Workstation path remains in ", f)
}
blocks <- read.csv(file.path(repo, "09_figures_tables/replay_blocks/provenance.csv"))
stopifnot(all(file.exists(file.path(repo,"09_figures_tables/replay_blocks",paste0(blocks$block,".R")))))
manifest <- read.csv(file.path(repo,"script_manifest.csv"))
stopifnot(all(file.exists(file.path(repo, manifest$script))))
cat("PASS:",length(files),"R files parse;",nrow(blocks),"bundled blocks exist; manifest entries resolve\n")
