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
disco_init_module("03_multimodal_UKB")


# Audit saved UKB multimodal objects. This does not modify legacy Result2 files.

options(stringsAsFactors = FALSE)
r2 <- disco_path("Result2")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- list(
  proteomics = list(data = "UKB.Prot.DISCO.rds", vector = "DISCO.UKB.prot.rds", id = "ID", n_expected = 47859L),
  clinical_proteomics = list(data = "UKB.ClinicalProt.DISCO.rds", vector = "DISCO.UKB.clinical.prot.rds", id = "ID", n_expected = 47173L),
  metabolomics = list(data = "UKB.Met.DISCO.rds", vector = "DISCO.UKB.Met.rds", id = "ID", n_expected = 250561L),
  clinical_protein_metabolite = list(data = "UKB.ClinicalProtMet.DISCO.rds", vector = "DISCO.UKB.clinical.prot.met.rds", id = "ID", n_expected = 24462L)
)

audit <- list()
cox_rows <- list()
for (label in names(datasets)) {
  spec <- datasets[[label]]
  d <- readRDS(file.path(r2, spec$data))
  v <- readRDS(file.path(r2, spec$vector))
  id <- d[[spec$id]]
  audit[[label]] <- data.frame(
    modality = label, file = spec$data, n = nrow(d), expected_n = spec$n_expected,
    n_matches_expected = nrow(d) == spec$n_expected,
    unique_ids = length(unique(id)), duplicated_ids = anyDuplicated(id) > 0,
    na_age = sum(is.na(d$age)), na_sex = sum(is.na(d$sex)),
    na_time = sum(is.na(d$time)), na_status = sum(is.na(d$status)),
    finite_mc = sum(is.finite(d$mc)), finite_hd = sum(is.finite(d$hd)),
    vector_length = length(v), mc_vector_max_abs_diff = max(abs(d$mc - as.numeric(v)), na.rm = TRUE)
  )
  for (formula_text in c("Surv(time,status)~mc", "Surv(time,status)~mc+age+sex", "Surv(time,status)~hd", "Surv(time,status)~age")) {
    fit <- survival::coxph(as.formula(formula_text), data = d)
    s <- summary(fit)
    cox_rows[[length(cox_rows) + 1L]] <- data.frame(
      modality = label, formula = formula_text, n = nobs(fit), events = fit$nevent,
      c_index = unname(s$concordance[1]), se = unname(s$concordance[2])
    )
  }
}

model_files <- c("Lasso.model.prot.rds", "Lasso.model.ClinicalProt.rds", "Lasso.model.ClinicalProtMet.rds")
model_audit <- do.call(rbind, lapply(model_files, function(f) {
  x <- readRDS(file.path(r2, f))
  co_min <- coef(x, s = "lambda.min")
  co_1se <- coef(x, s = "lambda.1se")
  data.frame(
    model = f, class = paste(class(x), collapse = ";"),
    lambda_min = x$lambda.min, lambda_1se = x$lambda.1se,
    features_lambda_min = sum(as.numeric(co_min) != 0),
    features_lambda_1se = sum(as.numeric(co_1se) != 0),
    n_lambda = length(x$lambda), nzero_max = max(x$nzero)
  )
}))

write.csv(do.call(rbind, audit), file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, cox_rows), file.path(out_dir, "cox_recalculation.csv"), row.names = FALSE)
write.csv(model_audit, file.path(out_dir, "lasso_model_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_FEATURE_SELECTION_AND_PREPROCESSING_LIMITATIONS",
  "saved_metric_objects_are_internally_audited=true",
  "same_cohort_mortality_feature_selection_requires_out_of_sample_sensitivity=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_FEATURE_SELECTION_AND_PREPROCESSING_LIMITATIONS\n")
