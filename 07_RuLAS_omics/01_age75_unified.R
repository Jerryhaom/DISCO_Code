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


# RuLAS omics: age<75 统一 reference 版本（重算 DISCO/DM）
# 对比原版（分位数 reference）的 HR/C-index
# 不动 Result6 旧文件

options(stringsAsFactors = FALSE)
suppressMessages({library(survival); library(Hmisc); library(parallel)})
r6 <- disco_path("Result6")
out <- disco_results_path("07_RuLAS_omics/outputs/07_age75_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ncores <- max(1, detectCores() - 1)

# ---- DISCO 函数（mclapply 并行分块）----
disco_core <- function(i, d, d2, weight, use_spearman) {
  cc1 <- cor(rbind(d2, d[i,,drop=FALSE]), method = if (use_spearman) "spearman" else "pearson")
  sum((cor(d2, method = if (use_spearman) "spearman" else "pearson") - cc1)^2 * weight)
}
dmc_par <- function(d, d2, weight, use_spearman, ncores = 8) {
  n0 <- nrow(d2)
  # 分块，减少 fork 开销
  idx <- seq_len(nrow(d))
  chunks <- split(idx, cut(idx, breaks = min(ncores, length(idx)), labels = FALSE))
  raw <- unlist(mclapply(chunks, function(ch) {
    cc_ref <- cor(d2, method = if (use_spearman) "spearman" else "pearson")
    sapply(ch, function(i) {
      cc1 <- cor(rbind(d2, d[i,,drop=FALSE]), method = if (use_spearman) "spearman" else "pearson")
      sum((cc_ref - cc1)^2 * weight)
    })
  }, mc.cores = ncores))
  log(raw[order(unlist(chunks))] * n0 * n0)
}
age_weight <- function(X, age) {
  cc <- apply(X, 2, function(x) cor(age, x, use="pairwise.complete.obs"))
  w <- abs(cc %*% t(cc)); diag(w) <- 0; w <- w/sum(w)
  w
}
run_metric <- function(X, age, use_spearman) {
  ref_idx <- which(age < 75)
  ref <- X[ref_idx,,drop=FALSE]; n0 <- nrow(ref)
  w <- age_weight(X, age)
  mc <- dmc_par(X, ref, w, use_spearman, ncores)
  # 移除常数/近常数列，ridge 正则化协方差避免奇异
  sdx <- apply(ref, 2, sd, na.rm=TRUE)
  keep <- sdx > 1e-8
  Xh <- X[, keep, drop=FALSE]; refh <- ref[, keep, drop=FALSE]
  mu <- colMeans(refh)
  S <- cov(refh) + diag(1e-6, ncol(refh))
  hd <- sqrt(mahalanobis(Xh, mu, S))
  list(mc=mc, hd=hd, nref=n0, p_used=ncol(refh))
}

eval_surv <- function(mc, hd, surv, label) {
  idx <- match(surv$ID, surv$ID)
  d <- data.frame(time=surv$time, status=as.numeric(surv$status), age=surv$age, sex=surv$sex, mc=mc, loghd=log(hd))
  d <- d[d$status %in% c(0,1), ]
  res <- list()
  for (metric in c('mc','loghd')) {
    sd_m <- sd(d[[metric]], na.rm = TRUE)
    f1 <- coxph(as.formula(paste0('Surv(time,status)~',metric)), data=d)
    f2 <- coxph(as.formula(paste0('Surv(time,status)~',metric,'+age+sex')), data=d)
    b1 <- coef(f1)[metric]; se1 <- sqrt(diag(vcov(f1)))[metric]
    b2 <- coef(f2)[metric]; se2 <- sqrt(diag(vcov(f2)))[metric]
    res[[metric]] <- data.frame(
      modality=label, metric=metric, N=f1$n, events=f1$nevent,
      sd_metric=sd_m,
      HR_crude=exp(b1), C_crude=summary(f1)$concordance[1],
      HR_crude_sd=exp(b1*sd_m), lo_crude_sd=exp((b1-1.96*se1)*sd_m), hi_crude_sd=exp((b1+1.96*se1)*sd_m),
      HR_adj=exp(b2), C_adj=summary(f2)$concordance[1],
      HR_adj_sd=exp(b2*sd_m), lo_adj_sd=exp((b2-1.96*se2)*sd_m), hi_adj_sd=exp((b2+1.96*se2)*sd_m))
  }
  do.call(rbind, res)
}

# ---- 数据 ----
pc <- readRDS(file.path(r6,'Rugao.PCs.met.rds'))
mic <- readRDS(file.path(r6,'Rugao.selected.microbe.rds'))
surv_met <- readRDS(file.path(r6,'Rugao.Met.surv.rds'))
surv_mic <- readRDS(file.path(r6,'Rugao.Meta.surv.rds'))
surv_agg <- readRDS(file.path(r6,'Rugao.Omics.surv.rds'))

# metabolomics: 23 PC (Pearson)
X_met <- as.matrix(pc[match(surv_met$ID, pc$ID), -1])
X_met <- apply(X_met, 2, as.numeric)
age_met <- surv_met$age
r_met <- run_metric(X_met, age_met, use_spearman=FALSE)

# microbiome: 27 species (Spearman)
X_mic <- as.matrix(mic[match(surv_mic$ID, mic$ID), -1])
X_mic <- apply(X_mic, 2, as.numeric)
age_mic <- surv_mic$age
r_mic <- run_metric(X_mic, age_mic, use_spearman=TRUE)

# aggregate: clinical(第10列起) + PC + microbe
rugao.phe <- readxl::read_xlsx(file.path(r6,'Analysis.2019.x.xlsx'))
group.phe <- read.csv(file.path(r6,'Group_info.2019.csv'))
rugao.phe <- rugao.phe[,c(colnames(rugao.phe)[1:6], group.phe$label)]
colnames(rugao.phe)[7:ncol(rugao.phe)] <- group.phe$Pheno
rugao.phe1 <- apply(rugao.phe[7:ncol(rugao.phe)], 2, function(x) suppressWarnings(as.numeric(x)))
rugao.phe1 <- apply(rugao.phe1, 2, function(x) Hmisc::impute(as.numeric(x)))
samples <- surv_agg$ID
clin <- rugao.phe1[match(samples, rugao.phe$X2019ID), 10:ncol(rugao.phe1)]
pc_a <- as.matrix(pc[match(samples, pc$ID), -1]); pc_a <- apply(pc_a, 2, as.numeric)
mic_a <- as.matrix(mic[match(samples, mic$ID), -1]); mic_a <- apply(mic_a, 2, as.numeric)
X_agg <- cbind(clin, pc_a, mic_a)
X_agg <- apply(X_agg, 2, as.numeric)
age_agg <- surv_agg$age
r_agg <- run_metric(X_agg, age_agg, use_spearman=TRUE)

# ---- 评估 ----
res <- rbind(
  eval_surv(r_met$mc, r_met$hd, surv_met, 'metabolomics'),
  eval_surv(r_mic$mc, r_mic$hd, surv_mic, 'microbiome'),
  eval_surv(r_agg$mc, r_agg$hd, surv_agg, 'aggregate')
)
res$reference <- 'age<75'
write.csv(res, file.path(out, 'age75_HR_Cindex.csv'), row.names=FALSE)

# manuscript 数字
man <- data.frame(
  modality=c('metabolomics','microbiome','aggregate'),
  HR_ms=c(1.17,1.38,2.49), C_ms=c(0.578,0.610,0.713))
cmp <- merge(res[res$metric=='mc',], man, by='modality')
write.csv(cmp, file.path(out, 'age75_vs_manuscript.csv'), row.names=FALSE)

cat('=== age<75 统一版 HR/C-index ===\n'); print(res, row.names=FALSE)
cat('\n=== 对比 manuscript（mc）===\n'); print(cmp[,c('modality','HR_adj','HR_ms','C_crude','C_ms')], row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo.txt'))
writeLines(c('status=PASS_AGE75_UNIFIED','reference=age<75','legacy_files_modified=false'), file.path(out,'status.txt'))
cat('\nDONE\n')
