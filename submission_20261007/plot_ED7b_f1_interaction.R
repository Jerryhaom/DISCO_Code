#!/usr/bin/env Rscript
# Reproduce the saved, fully adjusted DISCO-by-PRS interaction panel.
options(device=function(...) grDevices::pdf(file=NULL,...))
entry <- normalizePath(sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
base <- dirname(entry)
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(base,'output_ED7b')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out)
extra_library <- Sys.getenv('DISCO_R_LIBRARY','')
if(nzchar(extra_library)) .libPaths(c(extra_library,.libPaths()))
suppressPackageStartupMessages({library(forestploter);library(grid);library(gtable)})
d <- read.csv(file.path(base,'ED7b_f1_interaction_56.csv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(d)==56L,all(d$model=='f1_interaction'),all(d$term=='zD:zPRS'),
 all(d$endpoint_version=='raw_dates'),all(d$policy=='death_first'))
source(file.path(base,'ED7b_visual_block.R'),local=environment(),encoding='UTF-8')
stopifnot(file.exists(file.path(out,'ED7b_f1_interaction_NEW_ENDPOINT_CANDIDATE.pdf')))
cat('Saved-result interaction panel generated; no participant scoring or Cox refit.\n')
