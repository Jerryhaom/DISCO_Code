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


# Audit external-validation saved objects without modifying legacy files.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))
root <- disco_project_root()
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

specs <- list(
  CHARLS = list(file = file.path(root, "Result4/CHARLS.DISCO.rds"), id = "id", age = "age", expected_n = NA_integer_),
  CLHLS = list(file = file.path(root, "Result4/CLHLS.DISCO.rds"), id = NULL, age = "age", expected_n = 2546L),
  RuLAS = list(file = file.path(root, "Result4/RuLAS.DISCO.rds"), id = "ID", age = "age", expected_n = 1823L)
)

audit <- list(); cox_rows <- list()
for (label in names(specs)) {
  s <- specs[[label]]
  d <- readRDS(s$file)
  id_col <- if (is.null(s$id)) NA_character_ else s$id
  id <- if (is.na(id_col)) NULL else d[[id_col]]
  age_col <- s$age
  audit[[label]] <- data.frame(
    cohort = label, file = basename(s$file), n = nrow(d), expected_n = s$expected_n,
    n_matches_expected = is.na(s$expected_n) || nrow(d) == s$expected_n,
    id_column = id_col, unique_ids = if(is.null(id)) NA_integer_ else length(unique(id)),
    duplicated_ids = if(is.null(id)) NA else anyDuplicated(id) > 0,
    age_column = age_col, na_age = sum(is.na(d[[age_col]])),
    na_mc = sum(is.na(d$mc)), na_hd = sum(is.na(d$hd)),
    has_time = all(c("time", "status") %in% names(d)),
    na_time = if ("time" %in% names(d)) sum(is.na(d$time)) else NA,
    na_status = if ("status" %in% names(d)) sum(is.na(d$status)) else NA,
    invalid_status = if ("status" %in% names(d)) sum(!is.na(d$status) & !d$status %in% c(0, 1)) else NA,
    events = if ("status" %in% names(d)) sum(d$status == 1, na.rm = TRUE) else NA
  )
  if (all(c("time", "status") %in% names(d))) {
    # Never allow Surv() to silently convert undocumented status codes to NA.
    if (any(!is.na(d$status) & !d$status %in% c(0, 1))) next
    for (formula_text in c("Surv(time,status)~mc", "Surv(time,status)~mc+age", "Surv(time,status)~mc+age+sex")) {
      if (all(c("mc", "age", "sex") %in% names(d)) || !grepl("sex", formula_text)) {
        fit <- survival::coxph(as.formula(formula_text), data = d)
        ci <- summary(fit)$concordance
        cox_rows[[length(cox_rows) + 1L]] <- data.frame(
          cohort = label, formula = formula_text, n = fit$n, events = fit$nevent,
          c_index = unname(ci[1]), se = unname(ci[2])
        )
      }
    }
  }
}

# HRS is part of the PNAS submission package and may not be present as a local RDS
# in this repository. Record its documented audit values explicitly.
hrs_doc <- data.frame(
  cohort = "HRS", source = "PNAS/Submission_Rev1/Supplementary_Materials.docx",
  n = 4018L, events = 766L, reference = "NHANES age<30 external; N=6466",
  status = "documented_in_PNAS_submission"
)

disco_validate_cox_rows(do.call(rbind, cox_rows))
write.csv(do.call(rbind, audit), file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, cox_rows), file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
write.csv(hrs_doc, file.path(out_dir, "hrs_documented_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=BLOCKED_INVALID_STATUS_CODE_IN_CLHLS",
  "charls_clhls_rulas_saved_objects_audited=true",
  "hrs_audited_from_submission_documentation=true",
  "clhls_requires_explicit_handling_of_status_minus_9=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("BLOCKED_INVALID_STATUS_CODE_IN_CLHLS\n")
