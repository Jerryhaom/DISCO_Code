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
disco_init_module("05_organ_network")


# Organ DISCO 扩展协变量 cause-specific 分析（对齐 Supplementary Table 19）
# 最小改动：复用 Result3/Organ/Cal_Organ_spe.R 的 cause-specific 逻辑，
# 仅把 Model2 从 "+age+sex" 改为扩展协变量 "+age+sex+bmi+race+edu+smoke+alcohol"。
# 输出 per-unit HR（exp(coef)）+ 95% CI + P，7 结局 × 18 organ。
# 不改动 Result3 旧文件。

options(stringsAsFactors = FALSE)
suppressMessages(library(survival))
root <- disco_project_root()
out_dir <- disco_results_path("05_organ_network/outputs/organ_expanded_covariates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 数据（同 Cal_Organ_spe.R）----
odis <- readRDS(file.path(root, "Result3/Organ/UKB.odis.rds"))
prot <- readRDS(file.path(root, "Result3/Organ/UKB.prot.imputed.rds"))
eid  <- prot[, 1]  # odis 行序 = prot eid 顺序

# ---- 扩展协变量（按行对齐 odis）----
ext <- readRDS(file.path(root, "UKB_Data/extended_covariates_clinical.rds"))
race <- read.csv(file.path(root, "UKB_Data/ukb_race_21000.csv"), colClasses = "character")
idx <- match(eid, ext$ID)
odis$bmi    <- ext$bmi_raw[idx]
odis$edu    <- factor(ext$edu_high[idx])
odis$smoke  <- factor(ext$smoke_status[idx], levels = c(0, 1, 2))
odis$alcohol<- factor(ext$alcohol_freq[idx], levels = 1:6, ordered = TRUE)
odis$race   <- factor(race$race[match(as.character(eid), race$eid)])
odis$race[odis$race == ""] <- "White"
odis$sex    <- factor(odis$sex)

# ---- 18 organ（主分析锁定）----
organ18 <- c("tongue","salivary.gland","stomach","adipose.tissue","choroid.plexus",
             "retina","lung","pancreas","kidney","heart.muscle","skin","esophagus",
             "skeletal.muscle","intestine","bone.marrow","liver","brain","lymphoid.tissue")
cl <- which(colnames(odis) %in% organ18)
stopifnot(length(cl) == 18)

# ---- cause 顺序（同 Cal_Organ_spe.R）----
cause <- na.omit(unique(unique(odis$cause)))[c(2,1,6,5,3,4)]
# cause = Malignant, Heart, Cerebrovascular, Respiratory, Alzheimer's, Others
cl2 <- which(is.na(odis$cause))  # 存活

cov_str <- "+age+sex+bmi+race+edu+smoke+alcohol"  # 扩展协变量（Model2）

outcome_label <- c("All-cause","Heart disease","Respiratory disease","Cancer",
                   "Cerebrovascular disease","Alzheimer's disease","Others")
# 循环顺序 -> 结局映射（res 行序 = All cause + cause 顺序）
res_order <- c("All-cause","Cancer","Heart disease","Cerebrovascular disease",
               "Respiratory disease","Alzheimer's disease","Others")

res <- list()
for (i in cl) {
  organ <- colnames(odis)[i]
  fl <- paste0("Surv(time,status)~", organ, cov_str)
  hr_vec <- p_vec <- lo_vec <- hi_vec <- n_vec <- ev_vec <- numeric(7)
  names(hr_vec) <- res_order
  # all-cause
  c0 <- coxph(as.formula(paste0("Surv(time,status)~", organ)), data = odis)
  c1 <- coxph(as.formula(fl), data = odis)
  s0 <- summary(c0); s1 <- summary(c1)
  hr_vec["All-cause"] <- s1$coefficients[1,2]; p_vec["All-cause"] <- s1$coefficients[1,5]
  ci <- exp(confint(c1)[1,]); lo_vec["All-cause"] <- ci[1]; hi_vec["All-cause"] <- ci[2]
  n_vec["All-cause"] <- c1$n; ev_vec["All-cause"] <- c1$nevent
  # cause-specific（6 死因，同原始代码：存活 + 该死因）
  for (j in cause) {
    cl1 <- which(odis$cause == j)
    d1 <- odis[c(cl1, cl2), ]
    c <- coxph(as.formula(fl), data = d1)
    s <- summary(c)
    lab <- switch(j, "Heart Diseases"="Heart disease", "Respiratory system"="Respiratory disease",
                  "Malignant Neoplasms"="Cancer", "Cerebrovascular diseases"="Cerebrovascular disease",
                  "Alzheimer's disease"="Alzheimer's disease", "Others"="Others")
    hr_vec[lab] <- s$coefficients[1,2]; p_vec[lab] <- s$coefficients[1,5]
    ci <- exp(confint(c)[1,]); lo_vec[lab] <- ci[1]; hi_vec[lab] <- ci[2]
    n_vec[lab] <- c$n; ev_vec[lab] <- c$nevent
  }
  res[[length(res) + 1]] <- data.frame(
    organ = gsub("[.]", " ", organ),
    outcome = outcome_label,
    events_N = paste0(ev_vec[outcome_label], " / ", n_vec[outcome_label]),
    Model2_HR = round(hr_vec[outcome_label], 2),
    Model2_lo = round(lo_vec[outcome_label], 2),
    Model2_hi = round(hi_vec[outcome_label], 2),
    Model2_P = signif(p_vec[outcome_label], 3),
    stringsAsFactors = FALSE
  )
}
final <- do.call(rbind, res)
write.csv(final, file.path(out_dir, "organ_cause_specific_expanded.csv"), row.names = FALSE)

# 汇总：与 Table 19 旧值（age+sex）对照 brain 等关键 organ
cat("=== brain 扩展协变量 HR（对照 Table 19 旧值 age+sex）===\n")
b <- final[final$organ == "brain", ]
print(b[, c("outcome","Model2_HR","Model2_lo","Model2_hi","Model2_P")], row.names = FALSE)
cat("\nTable 19 旧值（age+sex）：All-cause 1.89, Heart 2.56, Respiratory 2.56, Cancer 1.66, Cerebrovascular 1.85, Alzheimer 1.56\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_cause_specific.txt"))
writeLines(c(
  "status=PASS",
  "model2=age+sex+bmi+race+edu+smoke+alcohol",
  "outcomes=all-cause_6_cause-specific",
  "organs=18",
  "hr_scale=per-unit",
  "reuses=Cal_Organ_spe.R_cause_specific_logic",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status_cause_specific.txt"))
cat("\nDONE\n")
