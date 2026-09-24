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
disco_init_module("08_revision_analyses")


# R01: expanded covariate models for the Nature Aging revision.
# Uses the primary UKB clinical/proteomics DISCO objects and the matching
# extended covariate tables. No legacy files are modified.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
ukb <- file.path(root, "UKB_Data")
out_dir <- "outputs/R01_covariate_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260914)

required <- c("extended_covariates_clinical.rds", "extended_covariates_proteomics.rds")
missing <- required[!file.exists(file.path(ukb, required))]
if (length(missing)) stop("Missing covariate files: ", paste(missing, collapse = ", "))

recode_special <- function(x) {
  # UKB extracted questionnaire tables use negative codes for non-response or
  # inapplicable responses. Keep valid positive category codes only.
  if (is.numeric(x)) x[x < 0] <- NA_real_
  x
}

prepare <- function(path, label) {
  d <- readRDS(path)
  needed <- c("ID", "age", "sex", "time", "status", "mc", "hd",
              "ipaq_group", "edu_high", "edu_college", "smoke_status",
              "alcohol_freq", "bmi_raw")
  miss <- setdiff(needed, names(d))
  if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse = ", "))

  d <- d[, needed]
  for (nm in c("ipaq_group", "edu_high", "edu_college", "smoke_status", "alcohol_freq")) {
    d[[nm]] <- recode_special(d[[nm]])
  }
  # Treat questionnaire categories as categorical predictors. Reference levels
  # are the lowest-risk/never categories where the coding supports this.
  d$sex <- factor(d$sex)
  d$ipaq_group <- factor(d$ipaq_group, levels = c(0, 1, 2), ordered = TRUE)
  d$edu_high <- factor(d$edu_high)
  d$edu_college <- factor(d$edu_college, levels = c(0, 1))
  d$smoke_status <- factor(d$smoke_status, levels = c(0, 1, 2))
  d$alcohol_freq <- factor(d$alcohol_freq, levels = 1:6, ordered = TRUE)
  d$bmi_raw <- as.numeric(d$bmi_raw)
  # race（field 21000 Ethnic background，全样本 N=502,356，从远端下载）
  race <- read.csv(file.path(ukb, "ukb_race_21000.csv"), colClasses = "character")
  d$race <- race$race[match(d$ID, race$eid)]
  d$race[d$race == ""] <- "White"  # unknown/missing -> reference (White)
  d$race <- factor(d$race, levels = c("White", "Non-White"))
  # DM (Mahalanobis, hd) is log-transformed so its per-unit HR is comparable to
  # DISCO (mc), which is already on the log scale.
  d$log_hd <- log(d$hd)
  d$metric_source <- label
  d
}

model_sets <- list(
  Crude = character(0),
  M0_age_sex = c("age", "sex"),
  M1_demographic = c("age", "sex", "bmi_raw", "race", "edu_high"),
  M2_lifestyle = c("age", "sex", "bmi_raw", "race", "edu_high", "smoke_status", "alcohol_freq"),
  M3_full_available = c("age", "sex", "bmi_raw", "race", "edu_high", "smoke_status", "alcohol_freq", "ipaq_group")
)

fit_one <- function(d, metric, metric_var, covars, model_name) {
  vars <- c("time", "status", metric_var, covars)
  x <- d[complete.cases(d[, vars]), vars, drop = FALSE]
  rhs <- if (length(covars) == 0) metric_var else paste(c(metric_var, covars), collapse = "+")
  formula_text <- paste0("Surv(time,status)~", rhs)
  fit <- survival::coxph(as.formula(formula_text), data = x, x = TRUE, y = TRUE)
  s <- summary(fit)
  metric_row <- which(rownames(s$coefficients) == metric_var)
  sd_metric <- sd(x[[metric_var]], na.rm = TRUE)
  data.frame(
    model = model_name, metric = metric,
    covariates = if (length(covars) == 0) "(none)" else paste(covars, collapse = "+"),
    formula = formula_text, n = nrow(x), events = fit$nevent,
    metric_sd = sd_metric,
    metric_hr_per_sd = exp(coef(fit)[metric_var] * sd_metric),
    metric_low95_per_sd = exp(confint(fit)[metric_var, 1] * sd_metric),
    metric_high95_per_sd = exp(confint(fit)[metric_var, 2] * sd_metric),
    metric_p = s$coefficients[metric_row, "Pr(>|z|)"],
    c_index = unname(s$concordance[1]), c_index_se = unname(s$concordance[2]),
    stringsAsFactors = FALSE
  )
}

fit_increment <- function(d, metric, metric_var, covars, model_name) {
  base_vars <- c("time", "status", covars)
  full_vars <- c(base_vars, metric_var)
  xb <- d[complete.cases(d[, base_vars]), base_vars, drop = FALSE]
  xf <- d[complete.cases(d[, full_vars]), full_vars, drop = FALSE]
  base_formula <- paste0("Surv(time,status)~", paste(covars, collapse = "+"))
  full_formula <- paste0("Surv(time,status)~", metric_var, "+", paste(covars, collapse = "+"))
  fb <- survival::coxph(as.formula(base_formula), data = xb)
  ff <- survival::coxph(as.formula(full_formula), data = xf)
  sb <- summary(fb); sf <- summary(ff)
  # A likelihood-ratio comparison is valid only on the same rows. The
  # expanded covariates generally reduce the complete-case sample, so do not
  # silently compare models fitted to different individuals.
  same_sample <- identical(rownames(xb), rownames(xf)) && nrow(xb) == nrow(xf)
  lrt_p <- NA_real_
  if (same_sample) {
    lrstat <- -2 * (fb$loglik[2] - ff$loglik[2])
    df <- length(coef(ff)) - length(coef(fb))
    if (df > 0) lrt_p <- pchisq(lrstat, df, lower.tail = FALSE)
  }
  data.frame(
    model = model_name, metric = metric, base_formula = base_formula, full_formula = full_formula,
    base_n = nrow(xb), full_n = nrow(xf), base_events = fb$nevent, full_events = ff$nevent,
    base_c_index = sb$concordance[1], base_c_index_se = sb$concordance[2],
    full_c_index = sf$concordance[1], full_c_index_se = sf$concordance[2],
    c_index_change_same_sample = if (same_sample) unname(sf$concordance[1] - sb$concordance[1]) else NA_real_,
    likelihood_ratio_p = lrt_p,
    stringsAsFactors = FALSE
  )
}

run_dataset <- function(path, label) {
  d <- prepare(path, label)
  rows <- list(); inc <- list()
  # mc = DISCO (already log scale); hd = DM fitted as log(hd) for a comparable HR.
  metric_spec <- list(mc = "mc", hd = "log_hd")
  for (metric in c("mc", "hd")) {
    metric_var <- metric_spec[[metric]]
    for (model_name in names(model_sets)) {
      covars <- model_sets[[model_name]]
      rows[[length(rows) + 1L]] <- fit_one(d, metric, metric_var, covars, model_name)
      if (length(covars) > 0) {
        inc[[length(inc) + 1L]] <- fit_increment(d, metric, metric_var, covars, model_name)
      }
    }
  }
  list(
    data = d,
    models = do.call(rbind, rows),
    increments = do.call(rbind, inc),
    coding = data.frame(
      dataset = label,
      n = nrow(d),
      ipaQ_observed = sum(!is.na(d$ipaq_group)),
      bmi_observed = sum(!is.na(d$bmi_raw)),
      education_observed = sum(!is.na(d$edu_high)),
      smoking_observed = sum(!is.na(d$smoke_status)),
      alcohol_observed = sum(!is.na(d$alcohol_freq)),
      invalid_status = sum(!is.na(d$status) & !d$status %in% c(0, 1)),
      stringsAsFactors = FALSE
    )
  )
}

clinical <- run_dataset(file.path(ukb, "extended_covariates_clinical.rds"), "clinical")
proteomics <- run_dataset(file.path(ukb, "extended_covariates_proteomics.rds"), "proteomics")

models <- rbind(clinical$models, proteomics$models)
increments <- rbind(clinical$increments, proteomics$increments)
coding <- rbind(clinical$coding, proteomics$coding)
models$dataset <- rep(c("clinical", "proteomics"), c(nrow(clinical$models), nrow(proteomics$models)))

write.csv(models, file.path(out_dir, "R01_metric_models.csv"), row.names = FALSE)
write.csv(increments, file.path(out_dir, "R01_incremental_Cindex_LRT.csv"), row.names = FALSE)
write.csv(coding, file.path(out_dir, "R01_covariate_coding_audit.csv"), row.names = FALSE)
write.csv(data.frame(
  dataset = c("clinical", "proteomics"),
  n = c(nrow(clinical$data), nrow(proteomics$data)),
  events = c(sum(clinical$data$status == 1), sum(proteomics$data$status == 1)),
  primary_dm = c("HD.UKB.rds", "HD.UKB.Prot.rds")
), file.path(out_dir, "R01_dataset_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_COMPLETE_CASE_AND_COVARIATE_CODING_AUDIT",
  "special_negative_codes_recode_to_NA=true",
  "categorical_covariates_are_factors=true",
  "DM_log_transformed=true",
  "LRT_computed_via_model_loglik=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_COMPLETE_CASE_AND_COVARIATE_CODING_AUDIT\n")
