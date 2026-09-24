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

# RuLAS omics: 正确的层级 reference 版本（方案 B）
# metabolomics -> phe DISCO 25%; microbiome -> metabolomics DISCO 25%; aggregate -> age<75
options(stringsAsFactors = FALSE)
suppressMessages({library(survival); library(Hmisc); library(parallel)})
r6 <- disco_path("Result6")
out <- disco_results_path("07_RuLAS_omics/outputs/07_age75_analysis")
ncores <- max(1, detectCores() - 1)

disco_core <- function(i, d, d2, weight, use_spearman) {
  cc1 <- cor(rbind(d2, d[i,,drop=FALSE]), method = if(use_spearman)"spearman" else "pearson")
  sum((cor(d2, method=if(use_spearman)"spearman" else "pearson") - cc1)^2 * weight)
}
dmc_par <- function(d, d2, weight, use_spearman, ncores=8) {
  n0 <- nrow(d2); idx <- seq_len(nrow(d))
  chunks <- split(idx, cut(idx, breaks=min(ncores,length(idx)), labels=FALSE))
  raw <- unlist(mclapply(chunks, function(ch){
    cc_ref <- cor(d2, method=if(use_spearman)"spearman" else "pearson")
    sapply(ch, function(i){
      cc1 <- cor(rbind(d2, d[i,,drop=FALSE]), method=if(use_spearman)"spearman" else "pearson")
      sum((cc_ref - cc1)^2 * weight)
    })
  }, mc.cores=ncores))
  log(raw[order(unlist(chunks))] * n0 * n0)
}
age_weight <- function(X, age){
  cc <- apply(X,2,function(x) cor(age,x,use="pairwise.complete.obs"))
  w <- abs(cc %*% t(cc)); diag(w)<-0; w<-w/sum(w); w
}

# 数据
pc <- readRDS(file.path(r6,'Rugao.PCs.met.rds'))
mic <- readRDS(file.path(r6,'Rugao.selected.microbe.rds'))
surv_met <- readRDS(file.path(r6,'Rugao.Met.surv.rds'))
surv_mic <- readRDS(file.path(r6,'Rugao.Meta.surv.rds'))
surv_agg <- readRDS(file.path(r6,'Rugao.Omics.surv.rds'))
phe <- readRDS(file.path(r6,'Rugao.Phe.DISCO.rds'))  # clinical DISCO (dis$mc)

X_met <- as.matrix(pc[match(surv_met$ID, pc$ID), -1]); X_met <- apply(X_met,2,as.numeric)
age_met <- surv_met$age
X_mic <- as.matrix(mic[match(surv_mic$ID, mic$ID), -1]); X_mic <- apply(X_mic,2,as.numeric)
age_mic <- surv_mic$age

# --- metabolomics: reference = phe DISCO 25% ---
phe_mc <- phe$mc[match(surv_met$ID, phe$ID)]
ref_met_idx <- which(phe_mc < quantile(phe_mc, na.rm=TRUE)[2])
cat('metabolomics reference n=', length(ref_met_idx), '\n')
w_met <- age_weight(X_met, age_met)
met_disco <- dmc_par(X_met, X_met[ref_met_idx,,drop=FALSE], w_met, FALSE, ncores)

# --- microbiome: reference = metabolomics DISCO 25%（对齐到 microbiome 样本）---
met_disco_mic <- met_disco[match(surv_mic$ID, surv_met$ID)]
ref_mic_idx <- which(met_disco_mic < quantile(met_disco_mic, na.rm=TRUE)[2])
cat('microbiome reference n=', length(ref_mic_idx), '\n')
w_mic <- age_weight(X_mic, age_mic)
mic_disco <- dmc_par(X_mic, X_mic[ref_mic_idx,,drop=FALSE], w_mic, TRUE, ncores)

# --- aggregate: reference = age<75 ---
rugao.phe <- readxl::read_xlsx(file.path(r6,'Analysis.2019.x.xlsx'))
group.phe <- read.csv(file.path(r6,'Group_info.2019.csv'))
rugao.phe <- rugao.phe[,c(colnames(rugao.phe)[1:6], group.phe$label)]
colnames(rugao.phe)[7:ncol(rugao.phe)] <- group.phe$Pheno
rugao.phe1 <- apply(rugao.phe[7:ncol(rugao.phe)],2,function(x) suppressWarnings(as.numeric(x)))
rugao.phe1 <- apply(rugao.phe1,2,function(x) Hmisc::impute(as.numeric(x)))
samples <- surv_agg$ID
clin <- rugao.phe1[match(samples, rugao.phe$X2019ID), 10:ncol(rugao.phe1)]
pc_a <- as.matrix(pc[match(samples, pc$ID), -1]); pc_a <- apply(pc_a,2,as.numeric)
mic_a <- as.matrix(mic[match(samples, mic$ID), -1]); mic_a <- apply(mic_a,2,as.numeric)
X_agg <- cbind(clin, pc_a, mic_a); X_agg <- apply(X_agg,2,as.numeric)
age_agg <- surv_agg$age
ref_agg_idx <- which(age_agg < 75)
cat('aggregate reference n=', length(ref_agg_idx), '\n')
w_agg <- age_weight(X_agg, age_agg)
agg_disco <- dmc_par(X_agg, X_agg[ref_agg_idx,,drop=FALSE], w_agg, TRUE, ncores)

eval_surv <- function(mc, surv, label, ref_desc){
  d <- data.frame(time=surv$time, status=as.numeric(surv$status), age=surv$age, sex=surv$sex, mc=mc)
  d <- d[d$status %in% c(0,1),]
  sd_m <- sd(d$mc, na.rm=TRUE)
  f1 <- coxph(Surv(time,status)~mc, data=d)
  f2 <- coxph(Surv(time,status)~mc+age+sex, data=d)
  b1 <- coef(f1)['mc']; se1 <- sqrt(diag(vcov(f1)))['mc']
  b2 <- coef(f2)['mc']; se2 <- sqrt(diag(vcov(f2)))['mc']
  data.frame(modality=label, reference=ref_desc, N=f1$n, events=f1$nevent, sd_mc=sd_m,
    HR_crude=exp(b1), HR_crude_sd=exp(b1*sd_m), lo_crude_sd=exp((b1-1.96*se1)*sd_m), hi_crude_sd=exp((b1+1.96*se1)*sd_m),
    C_crude=summary(f1)$concordance[1],
    HR_adj=exp(b2), HR_adj_sd=exp(b2*sd_m), lo_adj_sd=exp((b2-1.96*se2)*sd_m), hi_adj_sd=exp((b2+1.96*se2)*sd_m),
    C_adj=summary(f2)$concordance[1])
}
res <- rbind(
  eval_surv(met_disco, surv_met, 'metabolomics', 'phe DISCO <25%'),
  eval_surv(mic_disco, surv_mic, 'microbiome', 'met DISCO <25%'),
  eval_surv(agg_disco, surv_agg, 'aggregate', 'age<75')
)
write.csv(res, file.path(out,'discrete_ref_HR_Cindex.csv'), row.names=FALSE)
print(res[,c('modality','reference','N','events','sd_mc','HR_crude','HR_crude_sd','HR_adj','HR_adj_sd','C_adj')], digits=4, row.names=FALSE)
writeLines(c('status=PASS_DISCRETE_REFERENCE','reference=met_pheDISCO25/mic_metDISCO25/agg_age75','legacy_files_modified=false'), file.path(out,'status_discrete.txt'))
cat('\nDONE\n')
