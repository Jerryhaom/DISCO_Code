#!/usr/bin/env Rscript
# Rebuild figure candidates from bundled, non-participant aggregate inputs.
args <- commandArgs(trailingOnly=FALSE)
file <- sub('^--file=', '', args[grepl('^--file=', args)])
stopifnot(length(file)==1L)
entry <- dirname(normalizePath(file))
if (!nzchar(Sys.getenv('DISCO_RESULTS_ROOT')))
  Sys.setenv(DISCO_RESULTS_ROOT=file.path(entry,'sources'))
producer <- file.path(Sys.getenv('DISCO_RESULTS_ROOT'),
 '03_multimodal_UKB/targeted_endpoint_diagnostics_20261007/figure_candidates/code/run_all.R')
stopifnot(file.exists(producer))
source(producer, chdir=FALSE)
