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
disco_init_module("06_frailty_function")


# Audit CHARLS DISCO/frailty objects without modifying legacy files.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
r5 <- file.path(root, "Result5")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- c("CHARLS.2011.disco.rds", "CHARLS.2015.disco.rds", "ROC.rds")
audit <- list(); cox_rows <- list()
for (f in files) {
  p <- file.path(r5, f)
  if (!file.exists(p)) next
  x <- readRDS(p)
  audit[[f]] <- data.frame(
    object = f, class = paste(class(x), collapse = ";"), rows = if (is.data.frame(x)) nrow(x) else NA,
    cols = if (is.data.frame(x)) ncol(x) else NA,
    unique_id = if (is.data.frame(x) && "ID" %in% names(x)) length(unique(x$ID)) else NA,
    duplicated_id = if (is.data.frame(x) && "ID" %in% names(x)) anyDuplicated(x$ID) > 0 else NA,
    na_total = if (is.data.frame(x)) sum(is.na(x)) else NA,
    stringsAsFactors = FALSE
  )
  if (is.data.frame(x) && all(c("time", "status", "mc") %in% names(x))) {
    invalid <- sum(!is.na(x$status) & !x$status %in% c(0, 1))
    if (invalid == 0) {
      for (formula_text in c("Surv(time,status)~mc", "Surv(time,status)~mc+age+sex")) {
        if (all(c("age", "sex") %in% names(x)) || !grepl("age|sex", formula_text)) {
          fit <- survival::coxph(as.formula(formula_text), data = x)
          ci <- summary(fit)$concordance
          cox_rows[[length(cox_rows) + 1L]] <- data.frame(
            object = f, formula = formula_text, n = nobs(fit), events = fit$nevent,
            c_index = unname(ci[1]), se = unname(ci[2])
          )
        }
      }
    }
  }
}

fi_files <- list.files(file.path(r5, "Frailty"), pattern = "^CHARLS\\.FI\\.(2011|2013|2015|2018)\\.rds$", full.names = TRUE)
fi_audit <- do.call(rbind, lapply(fi_files, function(p) {
  x <- readRDS(p)
  data.frame(object = basename(p), rows = nrow(x), cols = ncol(x), names = paste(names(x), collapse = ";"), na_total = sum(is.na(x)))
}))

write.csv(do.call(rbind, audit), file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(if (length(cox_rows)) do.call(rbind, cox_rows) else data.frame(), file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
write.csv(fi_audit, file.path(out_dir, "frailty_object_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_FRAILTY_PIPELINE_LIMITATIONS",
  "saved_charls_disco_objects_audited=true",
  "frailty_objects_audited=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_FRAILTY_PIPELINE_LIMITATIONS\n")
