#!/usr/bin/env Rscript
# Purely synthetic demonstration; no restricted inputs are required.
entry <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(entry)),'../project_config.R'))
if(!requireNamespace('DISCO4HD',quietly=TRUE) || packageVersion('DISCO4HD')<'0.3.1')
  stop('Install the audited DISCO4HD >= 0.3.1 source package first')
Sys.setenv(DISCO_EXAMPLE_OUTPUT=disco_results_path('10_biomarker_calculator/synthetic_nhanes'))
source(system.file('examples','synthetic_nhanes.R',package='DISCO4HD'))
