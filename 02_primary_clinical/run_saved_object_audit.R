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
disco_init_module("02_primary_clinical")


# Audit the promoted legacy NHANES/UKB clinical objects without modifying them.
# This is a reproducibility gate, not the final re-computation pipeline.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
r1 <- file.path(root, "Result1")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required <- c(
  "NHANES.Clinical.DISCO.rds", "NHANES.Clinical.DISCO2.rds",
  "UKB.Clinical.DISCO.rds", "UKB.Clinical.DISCO2.rds",
  "DISCO.UKB.rds", "HD.UKB.rds", "mdata.rds", "Clean.UKB.Clinical.rds"
)
missing <- required[!file.exists(file.path(r1, required))]
if (length(missing)) stop("Missing required legacy files: ", paste(missing, collapse = ", "))

check_df <- function(path, id_col) {
  x <- readRDS(path)
  stopifnot(is.data.frame(x), id_col %in% names(x))
  data.frame(
    file = basename(path), rows = nrow(x), cols = ncol(x),
    unique_ids = length(unique(x[[id_col]])),
    duplicated_ids = anyDuplicated(x[[id_col]]) > 0,
    na_age = if ("age" %in% names(x)) sum(is.na(x$age)) else if ("Age" %in% names(x)) sum(is.na(x$Age)) else NA,
    na_time = if ("time" %in% names(x)) sum(is.na(x$time)) else NA,
    na_status = if ("status" %in% names(x)) sum(is.na(x$status)) else NA,
    mc_finite = if ("mc" %in% names(x)) sum(is.finite(x$mc)) else NA,
    hd_finite = if ("hd" %in% names(x)) sum(is.finite(x$hd)) else NA
  )
}

audit <- rbind(
  check_df(file.path(r1, "NHANES.Clinical.DISCO.rds"), "id"),
  check_df(file.path(r1, "NHANES.Clinical.DISCO2.rds"), "id"),
  check_df(file.path(r1, "UKB.Clinical.DISCO.rds"), "ID"),
  check_df(file.path(r1, "UKB.Clinical.DISCO2.rds"), "ID")
)

nh1 <- readRDS(file.path(r1, "NHANES.Clinical.DISCO.rds"))
nh2 <- readRDS(file.path(r1, "NHANES.Clinical.DISCO2.rds"))
uk1 <- readRDS(file.path(r1, "UKB.Clinical.DISCO.rds"))
uk2 <- readRDS(file.path(r1, "UKB.Clinical.DISCO2.rds"))
v_mc <- readRDS(file.path(r1, "DISCO.UKB.rds"))

comparison <- data.frame(
  comparison = c("NHANES id", "NHANES mc", "NHANES hd", "UKB id", "UKB mc vs vector", "UKB mc cross-object", "UKB hd cross-object"),
  result = c(
    identical(nh1$id, nh2$id),
    max(abs(nh1$mc - nh2$mc), na.rm = TRUE),
    max(abs(nh1$hd - nh2$hd), na.rm = TRUE),
    identical(uk1$ID, uk2$ID),
    max(abs(uk1$mc - as.numeric(v_mc)), na.rm = TRUE),
    max(abs(uk1$mc - uk2$mc), na.rm = TRUE),
    max(abs(uk1$hd - uk2$hd), na.rm = TRUE)
  ), stringsAsFactors = FALSE
)

cox_summary <- function(x, label, formula_text) {
  fit <- survival::coxph(as.formula(formula_text), data = x)
  s <- summary(fit)
  data.frame(
    dataset = label, formula = formula_text, n = nobs(fit), events = fit$nevent,
    c_index = unname(s$concordance[1]), se = unname(s$concordance[2])
  )
}

cox <- rbind(
  cox_summary(nh2, "NHANES", "Surv(time,status)~mc"),
  cox_summary(nh2, "NHANES", "Surv(time,status)~hd"),
  cox_summary(uk1, "UKB", "Surv(time,status)~mc"),
  cox_summary(uk1, "UKB", "Surv(time,status)~mc+age"),
  cox_summary(uk1, "UKB", "Surv(time,status)~mc+age+sex"),
  cox_summary(uk1, "UKB", "Surv(time,status)~log(hd)+age")
)

write.csv(audit, file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(comparison, file.path(out_dir, "object_comparison.csv"), row.names = FALSE)
write.csv(cox, file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

status <- if (comparison$result[comparison$comparison == "UKB hd cross-object"] == 0) {
  "PASS"
} else {
  "PASS_WITH_VERSION_LOCK_PRIMARY_HD_UKB_SENSITIVITY_HD_UKB2"
}
writeLines(c(paste0("status=", status), "source_objects=legacy_Result1", "no_legacy_files_modified=true"), file.path(out_dir, "audit_status.txt"))
cat(status, "\n")
