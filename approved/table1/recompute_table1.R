#!/usr/bin/env Rscript
# Within-cohort standardization in the effective Table 1 analysis sample.
suppressPackageStartupMessages({library(survival); library(digest); library(jsonlite)})
results_root <- file.path(Sys.getenv("DISCO_DATA_ROOT"),"Results_final")
code <- Sys.getenv("DISCO_APPROVED_CODE_ROOT")
source_script <- file.path(code,"table1/source_prepare.R")
preparation_script <- file.path(code,"table1/prepare_approved_cohorts.R")
reference_csv <- file.path(code,"table1/data/Table1_existing_full_models.csv")
crude_csv <- file.path(code,"table1/data/Table1_crude_C_index_details.csv")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)))
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[[1]] else file.path(dirname(script_path), "20261006")
if (dir.exists(out)) stop("Output exists; use a new output directory.")
raw_inputs <- c(file.path(Sys.getenv("DISCO_DATA_ROOT"),c("Result4/CHARLS.DISCO.rds","Result4/CLHLS.DISCO.rds",
 "Result4/CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds","Result4/RuLAS.DISCO.rds",
 "Result1/NHANES.Clinical.DISCO2.rds","Result1/NHANES.Clinical.DISCO.rds","UKB_Data/NHANES_covariates.rds",
 "UKB_Data/extended_covariates_clinical.rds","UKB_Data/covariates_raw_full.rds","UKB_Data/ukb_race_21000.csv")),
 Sys.getenv("DISCO_RULAS_QUESTIONNAIRE"))
input_paths <- unique(c(script_path,source_script,preparation_script,reference_csv,crude_csv,raw_inputs))
before <- vapply(input_paths, function(p) digest(file=p, algo="sha256"), character(1))
source(preparation_script)
prep <- prepare_approved_cohorts(results_root)
preparation_warnings <- prep$preparation_warnings
reference <- read.csv(reference_csv, stringsAsFactors=FALSE)
crude_reference <- read.csv(crude_csv, stringsAsFactors=FALSE)
rows <- list()
checks <- list()

for (cohort in names(prep$cohorts)) {
  d <- prep$cohorts[[cohort]]$d
  has_race <- prep$cohorts[[cohort]]$race
  for (v in c("mc", "hd", "time", "status", "age", "sex")) d[[v]] <- as.numeric(unclass(d[[v]]))
  invalid_status <- sum(!is.na(d$status) & !d$status %in% c(0,1))
  d$status[!d$status %in% c(0,1)] <- NA
  d$log_hd <- log(d$hd)
  for (v in c("sex", "smoke", "drink", "education", if (has_race) "race")) d[[v]] <- factor(d[[v]])
  covs <- c("age", "sex", "smoke", "drink", "education", if (has_race) "race")
  keep <- complete.cases(d[c("mc", "time", "status", covs)])
  stopifnot(identical(keep, complete.cases(d[c("log_hd", "time", "status", covs)])))
  evaluation <- d[keep, , drop=FALSE]
  stopifnot(all(is.finite(evaluation$mc)), all(is.finite(evaluation$log_hd)),
    all(is.finite(evaluation$time)), all(evaluation$time >= 0))

  for (label in c("DISCO", "DM")) {
    metric <- if (label=="DISCO") "mc" else "log_hd"
    sd_effective <- sd(evaluation[[metric]])
    sd_original <- sd(d[[metric]])
    stopifnot(is.finite(sd_effective), sd_effective > 0)
    evaluation$score_z <- (evaluation[[metric]]-mean(evaluation[[metric]]))/sd_effective
    scaled_fml <- as.formula(paste("Surv(time,status)~", paste(c("score_z", covs), collapse="+")))
    raw_fml <- as.formula(paste("Surv(time,status)~", paste(c(metric, covs), collapse="+")))
    warnings <- character()
    fit_model <- function(fml) withCallingHandlers(coxph(fml, data=evaluation, ties="efron"),
      warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    scaled_fit <- fit_model(scaled_fml)
    raw_fit <- fit_model(raw_fml)
    crude_fit <- fit_model(as.formula(paste0("Surv(time,status)~", metric)))
    sm <- summary(scaled_fit)
    b <- unname(coef(scaled_fit)["score_z"])
    se <- unname(sqrt(diag(vcov(scaled_fit)))["score_z"])
    ci <- unname(confint(scaled_fit)["score_z",])
    p <- unname(sm$coefficients["score_z","Pr(>|z|)"])
    raw_b <- unname(coef(raw_fit)[metric])
    raw_ci <- unname(confint(raw_fit)[metric,])
    old <- reference[reference$cohort==cohort & reference$metric==label, , drop=FALSE]
    old_crude <- crude_reference[crude_reference$cohort==cohort & crude_reference$metric==label, , drop=FALSE]
    stopifnot(nrow(old)==1L, nrow(old_crude)==1L,
      scaled_fit$n==old$n, scaled_fit$nevent==old$events,
      abs(sd(evaluation$score_z)-1)<1e-12)
    ci_crude <- unname(summary(crude_fit)$concordance)
    direct <- concordance(as.formula(paste0("Surv(time,status)~", metric)),
      data=evaluation, reverse=unname(coef(crude_fit)[metric])>0, timewt="n")
    proof <- c(abs(exp(b)-exp(raw_b*sd_effective)),
      abs(exp(ci)-exp(raw_ci*sd_effective)))
    p_difference <- p-old$P
    adjusted_c_difference <- unname(sm$concordance[1])-old$C_index
    crude_c_difference <- ci_crude[1]-old_crude$C_index_crude
    stopifnot(max(proof)<1e-10, abs(p_difference)<1e-12,
      abs(adjusted_c_difference)<1e-12, abs(crude_c_difference)<1e-12,
      abs(ci_crude[1]-direct$concordance)<1e-12,
      abs(ci_crude[2]-sqrt(direct$var))<1e-12,
      abs(exp(raw_b*sd_original)-old$HR)<1e-10)
    key <- paste(cohort, label, sep="__")
    hr <- exp(b)
    lower <- exp(ci[1]); upper <- exp(ci[2])
    rows[[key]] <- data.frame(cohort=cohort, metric=label, n=scaled_fit$n, events=scaled_fit$nevent,
      model=if (has_race) "A" else "B", adjusted_HR_per_analysis_SD=hr,
      adjusted_HR_CI95_lower=lower, adjusted_HR_CI95_upper=upper,
      adjusted_HR_report=sprintf("%.3f (%.3f\u2013%.3f)",hr,lower,upper),
      adjusted_P=p, adjusted_neg_log10_P=-log(2)/log(10)-pnorm(abs(b/se),lower.tail=FALSE,log.p=TRUE)/log(10),
      P_report=if (p==0) "<1\u00d710\u207b\u00b3\u2070\u2070" else sprintf("%.2e",p),
      crude_C_index=ci_crude[1], crude_C_index_SE=ci_crude[2],
      crude_C_index_CI95_lower=max(0,ci_crude[1]-qnorm(.975)*ci_crude[2]),
      crude_C_index_CI95_upper=min(1,ci_crude[1]+qnorm(.975)*ci_crude[2]),
      crude_C_index_report=sprintf("%.3f",ci_crude[1]),
      score_mean_analysis_sample=mean(evaluation[[metric]]), score_SD_analysis_sample=sd_effective,
      score_SD_previous_input=sd_original, input_n=nrow(d), excluded_n=sum(!keep),
      invalid_status_n=invalid_status, previous_HR=old$HR,
      previous_CI95_lower=old$CI95_lower, previous_CI95_upper=old$CI95_upper,
      HR_change=hr-old$HR, raw_score_formula=paste(deparse(raw_fml),collapse=""),
      standardized_formula=paste(deparse(scaled_fml),collapse=""),
      scaling="One within-cohort sample SD in the effective Table 1 Cox analysis population; DM log-transformed before scaling",
      warning_count=length(warnings), warnings=paste(unique(warnings),collapse=" | "),
      stringsAsFactors=FALSE)
    checks[[key]] <- data.frame(cohort=cohort, metric=label,
      HR_scaled_refit_vs_rescaled_raw=proof[1], CI_lower_difference=proof[2], CI_upper_difference=proof[3],
      P_difference=p_difference, adjusted_C_difference=adjusted_c_difference,
      crude_C_difference=crude_c_difference, crude_direct_difference=ci_crude[1]-direct$concordance)
    cat(sprintf("%s N=%d deaths=%d HR %.9f -> %.9f; SD %.9f -> %.9f\n",
      key,scaled_fit$n,scaled_fit$nevent,old$HR,hr,sd_original,sd_effective))
  }
}
tab <- do.call(rbind,rows); rownames(tab) <- NULL
check_tab <- do.call(rbind,checks); rownames(check_tab) <- NULL
stopifnot(nrow(tab)==10L, !anyDuplicated(paste(tab$cohort,tab$metric)), nrow(check_tab)==10L)
after <- vapply(input_paths,function(p) digest(file=p,algo="sha256"),character(1))
stopifnot(identical(before,after))
dir.create(out,recursive=TRUE)
write.csv(tab,file.path(out,"Table1_analysis_sample_SD.csv"),row.names=FALSE)
write.csv(check_tab,file.path(out,"standardization_verification.csv"),row.names=FALSE)
write.csv(tab[c("cohort","metric","n","events","model","adjusted_HR_report","P_report","crude_C_index_report")],
  file.path(out,"Table1_display.csv"),row.names=FALSE)
write.csv(data.frame(path=input_paths,sha256=before),file.path(out,"source_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
write_json(prep$audit,file.path(out,"approved_missing_data_rules.json"),
  pretty=TRUE,auto_unbox=TRUE)
writeLines(c("Table 1 standardization", "",
  "Adjusted HRs are expressed per one SD of DISCO or natural-log-transformed DM, calculated within each cohort's effective Cox analysis sample.",
  "SD uses the sample estimator with denominator N-1. The same participants are used for DISCO and DM within each cohort.",
  "Covariates and source missing-data coding are unchanged. Only the score scaling population is changed.",
  "Crude C-index uses only DISCO or log(DM) on the same participants. HR scaling does not change P values or concordance.",
  "CLHLS excludes 333 records with invalid mortality status before calculating either score's SD. Its effective N is 2,213 with 828 deaths.",
  "UKB missing race is explicitly imputed as White. CHARLS missing education is explicitly filled with the mode of observed binary-coded education in its scored cohort. These two rules were author-approved on 2026-10-06.",
  "Other Table 1 covariate-selection decisions remain open.", "",
  capture.output(print(tab[c("cohort","metric","n","events","adjusted_HR_report","crude_C_index_report")],row.names=FALSE))),
  file.path(out,"README.txt"))
write_json(list(status="UNIFIED_WITHIN_COHORT_EFFECTIVE_ANALYSIS_SAMPLE_SD_VERIFIED",
  run_time=as.character(Sys.time()), input_hashes_unchanged=TRUE,
  HR_scaling="Per sample SD of DISCO or log(DM) in the effective within-cohort analysis sample",
  n_models=nrow(tab), model_warning_count=sum(tab$warning_count),
  preparation_warnings=unique(preparation_warnings),
  max_scaled_refit_HR_CI_difference=max(abs(as.matrix(check_tab[3:5]))),
  max_P_difference=max(abs(check_tab$P_difference)),
  max_crude_C_difference=max(abs(check_tab$crude_C_difference)),
  changed_HR_rows=paste(tab$cohort[abs(tab$HR_change)>1e-10],tab$metric[abs(tab$HR_change)>1e-10]),
  limitations=c("Covariate sets, categorical imputation and UKB race coding are unchanged",
    "Crude C-index is apparent in-sample Harrell concordance, not independently validated")),
  file.path(out,"verification.json"),pretty=TRUE,auto_unbox=TRUE)
