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

# UKB 慢病 incidence 双模型（Model1 crude + Model2 扩展协变量 adjusted）
# 对齐 Supplementary Table 17 的协变量：age+sex+BMI+race+education+smoking+drinking
# 不动 Result2 旧文件

options(stringsAsFactors = FALSE)
suppressMessages({library(survival); library(Hmisc)})
root <- disco_project_root()
out <- disco_results_path("03_multimodal_UKB/outputs/03_disease_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

diseases <- c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')

# 三个 modality 的 DISCO/DM
mods <- list(
  ClinicalProt = readRDS(file.path(root,'Result2/UKB.ClinicalProt.DISCO.rds')),
  Prot         = readRDS(disco_results_path('03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds')),
  Clinical     = readRDS(file.path(root,'Result1/UKB.Clinical.DISCO.rds'))
)
dis <- readRDS(file.path(root,'UKB_Data/dat_disease.rds'))

# 协变量：extended_covariates（edu/bmi/smoke/alcohol）+ 全样本 race（ukb_race_21000.csv）
cov_clin <- readRDS(file.path(root,'UKB_Data/extended_covariates_clinical.rds'))
race <- read.csv(file.path(root,'UKB_Data/ukb_race_21000.csv'), colClasses='character')

recode <- function(x) { x[x < 0] <- NA; x }

run_one <- function(disease, d) {
  st <- dis[[paste0('status_',disease)]]; tm <- dis[[paste0('time_',disease)]]
  # 对齐协变量
  cv <- cov_clin[match(d$ID, cov_clin$ID), c('bmi_raw','smoke_status','alcohol_freq','edu_high')]
  rc <- race$race[match(as.character(d$ID), race$eid)]
  dd <- data.frame(time=tm[match(d$ID, dis$eid)], status=st[match(d$ID, dis$eid)],
                   mc=d$mc, loghd=log(d$hd), age=d$age, sex=d$sex,
                   bmi=as.numeric(recode(cv$bmi_raw)), smoke=as.factor(recode(cv$smoke_status)),
                   alcohol=as.factor(recode(cv$alcohol_freq)), edu=as.factor(recode(cv$edu_high)),
                   race=as.factor(as.character(rc)))
  dd$sex <- as.factor(dd$sex)
  # 排除 prevalent 和缺失
  dd <- dd[!is.na(dd$time) & dd$time >= 0 & dd$status %in% c(0,1) & !is.na(dd$mc), ]
  # 动态协变量：只保留 ≥2 水平且无非 NA 单水平的 factor；BMI 连续
  cov_vars <- character()
  for (nm in c('smoke','alcohol','edu','race')) {
    v <- dd[[nm]]
    if (length(unique(v[!is.na(v)])) >= 2) cov_vars <- c(cov_vars, nm)
  }
  if (length(unique(dd$bmi[!is.na(dd$bmi)])) >= 2) cov_vars <- c(cov_vars, 'bmi')
  out <- list()
  for (metric in c('mc','loghd')) {
    f1 <- coxph(as.formula(paste0('Surv(time,status)~',metric)), data=dd)
    f2 <- coxph(as.formula(paste0('Surv(time,status)~', metric, '+age+sex+', paste(cov_vars, collapse='+'))), data=dd)
    s1 <- summary(f1); s2 <- summary(f2)
    out[[metric]] <- data.frame(
      disease=disease, metric=metric, N=f1$n, events=f1$nevent,
      M1_crude_HR=exp(coef(f1)[metric]), M1_C=unname(s1$concordance[1]),
      M2_adj_HR=exp(coef(f2)[metric]),
      M2_HR_low=exp(confint(f2)[metric,1]), M2_HR_high=exp(confint(f2)[metric,2]),
      M2_C=unname(s2$concordance[1]),
      M2_covariates=paste(c('age','sex',cov_vars), collapse='+'))
  }
  do.call(rbind, out)
}

res <- list()
for (modname in names(mods)) {
  for (disease in diseases) {
    r <- run_one(disease, mods[[modname]])
    r$modality <- modname
    res[[length(res)+1]] <- r
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out, 'disease_dual_model.csv'), row.names=FALSE)

cat('=== ClinicalProt DISCO 慢病双模型 ===\n')
sub <- res[res$metric=='mc' & res$modality=='ClinicalProt', c('disease','events','M1_crude_HR','M1_C','M2_adj_HR','M2_C')]
print(sub, row.names=FALSE)
cat('\nmanuscript Table 17 (ClinicalProt DISCO): Obesity HR 1.11 C 0.615, Hypertension 1.41 C 0.642, Diabetes 1.73 C 0.721\n')
writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo_dual.txt'))
writeLines(c('status=PASS_DUAL_MODEL','M1=crude','M2=age+sex+BMI+race+education+smoking+drinking','legacy_files_modified=false'), file.path(out,'dual_status.txt'))
cat('\nDONE\n')
