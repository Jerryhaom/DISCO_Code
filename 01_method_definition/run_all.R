#!/usr/bin/env Rscript
# Portable release bootstrap (source/data/output roots are independent).
local({
  frames <- vapply(sys.frames(), function(e) if (is.null(e$ofile)) "" else as.character(e$ofile)[1], character(1))
  frames <- frames[nzchar(frames)]
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  entry <- if (length(frames)) tail(frames, 1) else if (length(arg)) sub("^--file=", "", arg[1]) else ""
  cfg <- file.path(dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE), mustWork = TRUE)), "..", "project_config.R")
  source(cfg, local = .GlobalEnv)
})
disco_init_module("01_method_definition")

module <- disco_results_path("01_method_definition")
setwd(module)
status <- system2("Rscript", c(shQuote(disco_code_path("01_method_definition", "run_method_smoke_test.R"))))
if (status != 0) stop("method smoke test failed")
status <- system2("Rscript", c(shQuote(disco_code_path("01_method_definition", "run_simulation_lock.R"))))
if (status != 0) stop("simulation lock failed")
status <- system2("Rscript", c(shQuote(disco_code_path("01_method_definition", "run_submission_simulation.R"))))
if (status != 0) stop("submission simulation audit failed")
status
