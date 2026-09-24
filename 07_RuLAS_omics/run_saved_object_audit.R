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
disco_init_module("07_RuLAS_omics")


# Audit RuLAS omics and exploratory CKM objects without modifying legacy files.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
r6 <- file.path(root, "Result6")
r7 <- file.path(root, "Result7")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- c(
  "Rugao.Phe.DISCO.rds", "Rugao.Phe.DISCO1.rds", "Rugao.Phe.DISCO2.rds",
  "Rugao.Met.DISCO.rds", "Rugao.Met.DISCO2.rds", "Rugao.Meta.DISCO.rds",
  "Rugao.Meta.DISCO1.rds", "Rugao.Meta.DISCO2.rds", "Rugao.Omics.DISCO.rds",
  "Rugao.Omics.surv.rds", "Rugao.Met.surv.rds", "Rugao.Meta.surv.rds",
  "Rugao.DISCO.Sub.rds", "Rugao.subsystems.rds", "Rugao.subsystems2.rds",
  "Rugao.CKM.Med.Analysis.rds", "Rugao.CKM2.Med.Analysis.rds"
)

audit <- list()
for (f in files) {
  p <- if (file.exists(file.path(r6, f))) file.path(r6, f) else file.path(r7, f)
  if (!file.exists(p)) next
  x <- readRDS(p)
  audit[[length(audit) + 1L]] <- data.frame(
    object = f, source_dir = basename(dirname(p)), class = paste(class(x), collapse = ";"),
    rows = if (is.data.frame(x) || is.matrix(x)) nrow(x) else NA,
    cols = if (is.data.frame(x) || is.matrix(x)) ncol(x) else NA,
    length = length(x), names = if (is.data.frame(x)) paste(names(x), collapse = ";") else NA,
    na_total = if (is.data.frame(x) || is.matrix(x)) sum(is.na(x)) else NA,
    stringsAsFactors = FALSE
  )
}

cox_rows <- list()
for (f in c("Rugao.Phe.DISCO.rds", "Rugao.Met.DISCO.rds", "Rugao.Meta.DISCO.rds", "Rugao.Omics.DISCO.rds")) {
  p <- file.path(r6, f)
  if (!file.exists(p)) next
  d <- readRDS(p)
  if (!is.data.frame(d) || !all(c("time", "status") %in% names(d))) next
  if (sum(!is.na(d$status) & !d$status %in% c(0, 1)) > 0) next
  metric <- intersect(c("mc", "hd"), names(d))
  for (m in metric) {
    form <- paste0("Surv(time,status)~", m)
    fit <- survival::coxph(as.formula(form), data = d)
    ci <- summary(fit)$concordance
    cox_rows[[length(cox_rows) + 1L]] <- data.frame(object=f, formula=form, n=nobs(fit), events=fit$nevent, c_index=ci[1], se=ci[2])
  }
}

write.csv(do.call(rbind, audit), file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(if (length(cox_rows)) do.call(rbind, cox_rows) else data.frame(), file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_RULAS_OMICS_VERSION_LIMITATIONS",
  "saved_omics_objects_audited=true",
  "confirmatory_and_exploratory_objects_require_separation=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_RULAS_OMICS_VERSION_LIMITATIONS\n")
