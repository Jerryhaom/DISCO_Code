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
disco_init_module("04_external_validation")


# Build the valid CLHLS survival analysis object from the existing DISCO object.
# The CLHLS codebook defines -9 as lost to follow-up. It is therefore excluded
# from the binary mortality analysis; the legacy object is not modified.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))
root <- disco_project_root()
source_file <- file.path(root, "Result4/CLHLS.DISCO.rds")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

legacy <- readRDS(source_file)
stopifnot(is.data.frame(legacy), all(c("status", "time", "mc", "hd") %in% names(legacy)))

status_before <- as.data.frame(table(legacy$status, useNA = "ifany"), stringsAsFactors = FALSE)
names(status_before) <- c("status", "n")
valid <- !is.na(legacy$status) & legacy$status %in% c(0, 1)
clean <- legacy[valid, , drop = FALSE]
clean$status <- as.integer(clean$status)

status_after <- as.data.frame(table(clean$status, useNA = "ifany"), stringsAsFactors = FALSE)
names(status_after) <- c("status", "n")
status_after$status <- as.character(status_after$status)
status_before$status <- as.character(status_before$status)
status_audit <- rbind(
  data.frame(stage = "before", status_before),
  data.frame(stage = "after", status_after)
)

cox_rows <- lapply(c("Surv(time,status)~mc", "Surv(time,status)~mc+age+sex", "Surv(time,status)~hd", "Surv(time,status)~log(hd)+age+sex"), function(f) {
  fit <- survival::coxph(as.formula(f), data = clean)
  s <- summary(fit)
  data.frame(formula = f, n = fit$n, events = fit$nevent,
             c_index = s$concordance[1], se = s$concordance[2])
})

disco_validate_cox_rows(do.call(rbind, cox_rows))
saveRDS(clean, file.path(out_dir, "CLHLS.DISCO.valid_status01.rds"))
write.csv(status_audit, file.path(out_dir, "CLHLS_status_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, cox_rows), file.path(out_dir, "CLHLS_valid_status_cox.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_VALID_STATUS_OBJECT_CREATED",
  "codebook_definition=-9 means lost to follow-up",
  paste0("legacy_n=", nrow(legacy)),
  paste0("valid_n=", nrow(clean)),
  paste0("excluded_invalid_status_n=", sum(!valid)),
  paste0("events=", sum(clean$status == 1)),
  "legacy_files_modified=false"
), file.path(out_dir, "clhls_clean_status.txt"))
cat("PASS_VALID_STATUS_OBJECT_CREATED\n")
