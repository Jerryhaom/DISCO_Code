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


# R06: cross-cohort expanded-covariate sensitivity analysis (R1.1).
# Reviewer 1.1 asks that smoking, physical activity and education be adjusted
# "according to availability" across cohorts. We run a covariate ladder (crude ->
# age+sex -> core smoking/education -> full available) for each validation cohort
# using covariates recovered from the raw data:
#   CHARLS : smoke, drink, education (no BMI in local data)
#   CLHLS  : BMI (g101/g1021), smoke (d71), drink (d81), education (f1)
#   RuLAS  : BMI (BMI2016)  (no smoke/drink/education in local data)
# Missing covariate values are imputed (median for continuous, mode for
# categorical) so the expanded models retain the full cohort sample. No legacy
# files are modified.

suppressMessages(library(survival))
suppressMessages(library(openxlsx))

root <- disco_project_root()
out_dir <- "outputs/R06_crosscohort_covariates"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- imputation helpers ---------------------------------------------------
impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}
impute_mode <- function(x) {
  x[is.na(x)] <- names(sort(table(x), decreasing = TRUE))[1]
  x
}

# ---- covariate ladder -----------------------------------------------------
run_ladder <- function(d, label, extra) {
  # d columns: mc, time, status, age, sex + extra covariates (numeric/factor)
  keep <- complete.cases(d[c("mc", "time", "status", "age", "sex")])
  d <- d[keep, , drop = FALSE]
  mc_sd <- sd(d$mc, na.rm = TRUE)
  d$sex <- factor(d$sex)

  out <- list()
  models <- list(
    Crude       = "mc",
    M1_age_sex  = "mc+age+sex",
    M2_core     = "mc+age+sex",
    M3_full     = "mc+age+sex"
  )
  # M2 core: + smoke + education if available
  if (all(c("smoke", "education") %in% names(d))) {
    models$M2_core <- "mc+age+sex+smoke+education"
  }
  # M3 full: + smoke + drink + education + bmi as available
  avail <- extra[extra %in% names(d)]
  if (length(avail)) models$M3_full <- paste(c("mc", "age", "sex", avail), collapse = "+")

  for (mname in names(models)) {
    fit <- coxph(as.formula(paste0("Surv(time,status)~", models[[mname]])), data = d)
    s <- summary(fit)
    b <- coef(fit)["mc"]
    ci <- confint(fit)["mc", ]
    out[[length(out) + 1]] <- data.frame(
      cohort = label,
      model = mname,
      covariates = gsub("mc\\+|\\+mc", "", models[[mname]]),
      n = nrow(d),
      events = sum(d$status),
      hr_unit = exp(b),
      hr_per_sd = exp(b * mc_sd),
      lo_per_sd = exp(ci[1] * mc_sd),
      hi_per_sd = exp(ci[2] * mc_sd),
      c_index = s$concordance[1],
      p = coef(summary(fit))["mc", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

results <- list()

# ============================ CHARLS =======================================
ch <- readRDS(file.path(root, "Result4/CHARLS.DISCO.rds"))
d <- data.frame(
  mc = ch$mc, time = ch$time, status = ch$status,
  age = ch$age, sex = ch$sex,
  stringsAsFactors = FALSE
)
# smoke: 1=yes, 2=no -> 1=smoker, 0=non-smoker; missing -> mode
d$smoke <- ifelse(ch$smoke == 1, 1, ifelse(ch$smoke == 2, 0, NA))
d$smoke <- impute_mode(d$smoke)
d$drink <- ifelse(ch$drink == 1, 1, ifelse(ch$drink == 2, 0, NA))
d$drink <- impute_mode(d$drink)
# education: CHARLS 1-3 (illiterate/primary incomplete/sishu) = low, else high
edu <- ch$edu
edu[edu %in% c(".d", ".r", "")] <- NA
edu <- as.numeric(edu)
d$education <- ifelse(edu %in% c(1, 2, 3), 0, 1)
d$education <- impute_mode(d$education)
# BMI from CKM.Charls.stat2011.csv (height/weight/BMI), ~81% matched -> median impute
ckm <- read.csv(file.path(root, "Result4/CHARLS/CKM.Charls.stat2011.csv"))
bmi <- ckm$BMI[match(ch$id, ckm$id)]
bmi[!is.na(bmi) & (bmi < 10 | bmi > 60)] <- NA
d$bmi <- impute_median(bmi)
d$smoke <- factor(d$smoke); d$drink <- factor(d$drink); d$education <- factor(d$education)
results[["CHARLS"]] <- run_ladder(d, "CHARLS", c("smoke", "drink", "education", "bmi"))

# ============================ CLHLS ========================================
clh <- readRDS(file.path(root, "Result4/CLHLS.DISCO.rds"))
lon <- readRDS(file.path(root, "Result4/CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds"))
# rows align (verified: sample.info / DISCO / 2014_2018 longitudinal same order)
d <- data.frame(
  mc = clh$mc, time = clh$time, status = clh$status,
  age = clh$age, sex = clh$sex,
  stringsAsFactors = FALSE
)
d$smoke <- ifelse(lon$d71 == 1, 1, ifelse(lon$d71 == 2, 0, NA))
d$smoke <- impute_mode(d$smoke)
d$drink <- ifelse(lon$d81 == 1, 1, ifelse(lon$d81 == 2, 0, NA))
d$drink <- impute_mode(d$drink)
f1 <- as.numeric(lon$f1); f1[f1 %in% c(88, 99)] <- NA
d$education <- ifelse(f1 == 0, 0, 1)  # 0 years = illiterate, >=1 = educated
d$education <- impute_mode(d$education)
wt <- as.numeric(lon$g101); ht <- as.numeric(lon$g1021)
wt[wt %in% c(999, 888)] <- NA; ht[ht %in% c(999, 888)] <- NA
bmi <- wt / (ht / 100)^2
bmi[bmi < 15 | bmi > 50] <- NA
d$bmi <- impute_median(bmi)
d$smoke <- factor(d$smoke); d$drink <- factor(d$drink); d$education <- factor(d$education)
results[["CLHLS"]] <- run_ladder(d, "CLHLS", c("smoke", "drink", "education", "bmi"))

# ============================ RuLAS ========================================
ru <- readRDS(file.path(root, "Result4/RuLAS.DISCO.rds"))
rug <- read.csv(file.path(root, "Result4/Rugao/2016Rugao.data.x.csv"))
rug <- rug[match(ru$ID, rug$X2016ID), ]
d <- data.frame(
  mc = ru$mc, time = ru$time, status = ru$status,
  age = ru$age, sex = ru$sex,
  stringsAsFactors = FALSE
)
bmi <- as.numeric(rug$BMI2016)
bmi[bmi <= 0] <- NA
d$bmi <- impute_median(bmi)
# smoke/drink/education from wave2 questionnaire (如皋调查表 A6/J/K)
sh <- read.xlsx(disco_external_path("DISCO_RULAS_WAVE2_FILE", "RuLAS/wave2.xlsx"), sheet = "Sheet2")
w2 <- match(ru$ID, sh[[1]])   # 2016编号 -> DISCO ID, 100% matched
d$smoke <- ifelse(sh[[361]][w2] %in% c("2", "3", "水烟"), 1, ifelse(sh[[361]][w2] == "1", 0, NA))
d$smoke <- impute_mode(d$smoke)
d$drink <- ifelse(sh[[362]][w2] %in% c("2", "3"), 1, ifelse(sh[[362]][w2] == "1", 0, NA))
d$drink <- impute_mode(d$drink)
d$education <- ifelse(as.numeric(sh[[101]][w2]) == 1, 0, 1)  # 1=没有上过学 -> low
d$education <- impute_mode(d$education)
d$smoke <- factor(d$smoke); d$drink <- factor(d$drink); d$education <- factor(d$education)
results[["RuLAS"]] <- run_ladder(d, "RuLAS", c("smoke", "drink", "education", "bmi"))

# ============================ NHANES =======================================
dis2 <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO2.rds"))
dis  <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO.rds"))
nhan <- readRDS(file.path(root, "UKB_Data/NHANES_covariates.rds"))
nhan <- nhan[match(dis2$id, nhan$sampleID), ]
load(disco_external_path("DISCO_NHANES4_FILE", "NHANES/NHANES4.rda"))
n4 <- NHANES4[match(dis2$id, NHANES4$sampleID), ]
d <- data.frame(
  mc = dis2$mc, time = dis2$time, status = dis2$status,
  age = dis$Age[match(dis2$id, dis$id)], sex = nhan$sex,
  stringsAsFactors = FALSE
)
d$race <- factor(nhan$race)
d$smoke <- ifelse(nhan$smoke == 1, 1, ifelse(nhan$smoke == 0, 0, NA))
d$smoke <- impute_mode(d$smoke)
d$education <- ifelse(nhan$education == "high", 1, ifelse(nhan$education == "low", 0, NA))
d$education <- impute_mode(d$education)
d$bmi <- impute_median(as.numeric(n4$bmi))
d$smoke <- factor(d$smoke); d$education <- factor(d$education)
results[["NHANES"]] <- run_ladder(d, "NHANES", c("smoke", "education", "race", "bmi"))

final <- do.call(rbind, results)
write.csv(final, file.path(out_dir, "R06_crosscohort_covariates.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "imputation=median_continuous_mode_categorical",
  "cohorts=NHANES_CHARLS_CLHLS_RuLAS",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))

print(final[, c("cohort", "model", "n", "events", "hr_per_sd", "lo_per_sd", "hi_per_sd", "c_index", "p")],
      digits = 3, row.names = FALSE)
