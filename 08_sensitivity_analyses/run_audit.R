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
disco_init_module("08_sensitivity_analyses")


options(stringsAsFactors = FALSE)
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

status <- read.csv("analysis_status.csv", check.names = FALSE)
status$script_exists <- file.exists(status$script)
status$output_exists <- file.exists(status$output)
status$ready <- status$status == "completed" & status$script_exists & status$output_exists

write.csv(status, file.path(out_dir, "revision_analysis_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PLANNING_ONLY",
  paste0("analysis_rows=", nrow(status)),
  paste0("runnable_scripts=", sum(status$script_exists)),
  paste0("completed_outputs=", sum(status$output_exists)),
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PLANNING_ONLY\n")
