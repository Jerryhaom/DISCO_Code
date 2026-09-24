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

# UKB 慢病 incidence 的 DISCO/DM HR 重算（对比 Supplementary Table 17）
# 不动 Result2 旧文件

options(stringsAsFactors = FALSE)
suppressMessages(library(survival))
root <- disco_project_root()
out <- disco_results_path("03_multimodal_UKB/outputs/03_disease_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

diseases <- c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')

# 三个 modality 的 DISCO/DM
mods <- list(
  ClinicalProt = readRDS(file.path(root,'Result2/UKB.ClinicalProt.DISCO.rds')),
  Prot         = readRDS(file.path(root,'Result2/UKB.Prot.DISCO.rds')),
  Clinical     = readRDS(file.path(root,'Result1/UKB.Clinical.DISCO.rds'))
)
dis <- readRDS(file.path(root,'UKB_Data/dat_disease.rds'))

run_one <- function(disease, d) {
  st <- dis[[paste0('status_',disease)]]
  tm <- dis[[paste0('time_',disease)]]
  # 对齐
  idx <- match(d$ID, dis$eid)
  dd <- data.frame(time=tm[idx], status=st[idx], mc=d$mc, loghd=log(d$hd), age=d$age, sex=d$sex)
  # 排除 prevalent（time<0）和缺失
  dd <- dd[!is.na(dd$time) & dd$time >= 0 & !is.na(dd$status) & dd$status %in% c(0,1) & !is.na(dd$mc), ]
  out <- list()
  for (metric in c('mc','loghd')) {
    f <- coxph(as.formula(paste0('Surv(time,status)~',metric)), data=dd)
    s <- summary(f)
    out[[metric]] <- data.frame(
      disease=disease, metric=metric, N=f$n, events=f$nevent,
      HR=exp(coef(f)[metric]), HR_low=exp(confint(f)[metric,1]), HR_high=exp(confint(f)[metric,2]),
      C=unname(s$concordance[1]))
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
write.csv(res, file.path(out, 'disease_HR_Cindex_recomputed.csv'), row.names=FALSE)

# 对比 Table 17 的 DISCO HR（ClinicalProt Obesity 等）
cat('=== 重算 DISCO 慢病 HR（crude）===\n')
sub <- res[res$metric=='mc', c('disease','modality','events','HR','C')]
print(sub[sub$modality=='ClinicalProt',], row.names=FALSE)
cat('\nmanuscript Table 17 (ClinicalProt DISCO): Obesity 1.11 C 0.615, Hypertension 1.41 C 0.642, Diabetes 1.73 C 0.721\n')
writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo.txt'))
writeLines(c('status=PASS','model=crude Cox','prevalent_time<0_excluded=true','legacy_files_modified=false'), file.path(out,'status.txt'))
cat('\nDONE\n')
