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
suppressPackageStartupMessages(library(survival))
root <- disco_project_root()
r6 <- file.path(root, "Result6")
r7 <- file.path(root, "Result7")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- c(
  "Rugao.Phe.DISCO.rds", "Rugao.Phe.DISCO1.rds", "Rugao.Phe.DISCO2.rds",
  "Rugao.Met.DISCO.rds", "Rugao.Met.DISCO2.rds", "Rugao.Meta.DISCO.rds",
  "Rugao.Meta.DISCO1.rds", "Rugao.Meta.DISCO2.rds", "Rugao.Omics.DISCO.rds",
  "Rugao.Omics.surv.rds", "Rugao.Met.surv.rds", "Rugao.Meta.surv.rds", "Rugao.2019.mortality.rds",
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

# Plot_surv.R uses the .DISCO2 clinical score joined by ID to mortality;
# the other three endpoint-bearing objects are saved by that same script.
# These are historical scores/endpoints, not approval of corrected NMR results.
specs <- c(Clinical="Rugao.Phe.DISCO2.rds", Metabolomics="Rugao.Met.surv.rds",
           Microbiota="Rugao.Meta.surv.rds", StoredOmics="Rugao.Omics.surv.rds",
           Fig2fAggregate="Rugao.DISCO.Sub.rds")
cox_rows <- list(); source_rows <- list(); definition_rows <- list()
for (label in names(specs)) {
  f <- unname(specs[label])
  p <- file.path(r6, f)
  if (!file.exists(p)) next
  d <- readRDS(p)
  endpoint_source <- f
  if (label == "Fig2fAggregate") {
    # Plot_subsystems2.R uses mean(eight clinical subsystem scores, Met, Meta),
    # not Rugao.Omics.surv.rds$mc. Check the actual plotted score definition.
    parts <- c("Heart","Oxygen.Transport","Immunity","Kidney","Liver",
               "Lipid.and.Glucose","Electrolyte","Vascular.health","Met","Meta")
    stopifnot(all(c(parts,"var","Age","Sex","ID") %in% names(d)),
              !anyNA(d$ID), !anyDuplicated(d$ID))
    rebuilt <- rowMeans(d[,parts,drop=FALSE])
    definition_rows[[length(definition_rows)+1L]] <- data.frame(
      object=f,source="Plot_subsystems2.R: mean columns 3:12 without na.rm",
      n=nrow(d),stored_var_missing=sum(is.na(d$var)),rebuilt_var_missing=sum(is.na(rebuilt)),
      missing_mask_differences=sum(is.na(rebuilt)!=is.na(d$var)),
      max_absolute_difference=max(abs(rebuilt-d$var),na.rm=TRUE),
      stored_var_matches_plotted_definition=isTRUE(all.equal(rebuilt,d$var,tolerance=1e-12,check.attributes=FALSE)))
    # The saved var is not assumed to match the plotted formula. Preserve it
    # as a comparator, then rebuild the historical plot's declared score.
    d$stored_var <- d$var; d$var <- rebuilt
    d$age <- d$Age; d$sex <- d$Sex
    endpoint_source <- paste(f,"; Plot_subsystems2.R historical endpoint conventions")
  }
  if (label == "Clinical") {
    death_path <- file.path(r6,"Rugao.2019.mortality.rds")
    if (!file.exists(death_path)) stop("Missing clinical mortality join input")
    mortality <- readRDS(death_path)
    stopifnot(all(c("ID","age") %in% names(d)),
              all(c("x2019ID","time","status") %in% names(mortality)),
              !anyNA(d$ID), !anyNA(mortality$x2019ID),
              !anyDuplicated(d$ID), !anyDuplicated(mortality$x2019ID))
    matched <- mortality[match(d$ID,mortality$x2019ID),,drop=FALSE]
    # Preserve the documented historical convention; loss-to-follow-up is not verified here.
    missing_status <- is.na(matched$status)
    matched$time[missing_status] <- 24
    matched$status[missing_status] <- 0
    d$time <- matched$time; d$status <- matched$status
    d <- d[which(d$age >= 60),,drop=FALSE]
    endpoint_source <- "Rugao.2019.mortality.rds; historical unmatched status -> censored at 24 months"
  }
  if (!is.data.frame(d) || !all(c("time", "status") %in% names(d))) next
  if (sum(!is.na(d$status) & !d$status %in% c(0, 1)) > 0) stop("Invalid status in ",f)
  source_rows[[length(source_rows)+1L]] <- data.frame(cohort=label,score_object=f,
    endpoint_source=endpoint_source,n=nrow(d),missing_time=sum(is.na(d$time)),
    missing_age=if("age" %in% names(d))sum(is.na(d$age)) else NA_integer_)
  metric <- if(label == "Fig2fAggregate") c("var","stored_var") else intersect(c("mc", "hd"), names(d))
  for (m in metric) {
    forms <- paste0("Surv(time,status)~",m)
    if(all(c("age","sex") %in% names(d))) forms <- c(forms,paste0(forms,"+age+sex"))
    for(form in forms) {
      fit <- survival::coxph(as.formula(form), data = d)
      ci <- summary(fit)$concordance
      cox_rows[[length(cox_rows) + 1L]] <- data.frame(cohort=label,object=f, formula=form, n=fit$n, events=fit$nevent, c_index=unname(ci[1]), se=unname(ci[2]))
    }
  }
}

write.csv(do.call(rbind, audit), file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
disco_validate_cox_rows(do.call(rbind, cox_rows))
write.csv(if (length(cox_rows)) do.call(rbind, cox_rows) else data.frame(), file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
write.csv(do.call(rbind,source_rows),file.path(out_dir,"cox_input_sources.csv"),row.names=FALSE)
write.csv(do.call(rbind,definition_rows),file.path(out_dir,"aggregate_definition_check.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_RULAS_OMICS_VERSION_LIMITATIONS",
  "saved_omics_objects_audited=true",
  "confirmatory_and_exploratory_objects_require_separation=true",
  "historical_endpoint_objects_not_raw_followup_validation=true",
  "corrected_nmr_downstream_not_replaced_by_this_audit=true",
  "aggregate_historical_plot_formula_rebuilt_saved_var_not_assumed_equivalent=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_RULAS_OMICS_VERSION_LIMITATIONS\n")
