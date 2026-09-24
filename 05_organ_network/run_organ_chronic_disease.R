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


# Organ DISCO 慢病 incidence 扩展协变量分析（对齐 Fig 4a 的 Chronic disease 面板 / Table 18）
# 8 慢病 × 18 organ，Model2 = 扩展协变量 age+sex+bmi+race+edu+smoke+alcohol，per-unit HR。
# 复用 organ cause-specific 的扩展协变量匹配逻辑，最小改动。

options(stringsAsFactors = FALSE)
suppressMessages(library(survival))
root <- disco_project_root()
out_dir <- disco_results_path("05_organ_network/outputs/organ_expanded_covariates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 数据 + 扩展协变量（同 run_organ_cause_specific.R）----
odis <- readRDS(file.path(root, "Result3/Organ/UKB.odis.rds"))
prot <- readRDS(file.path(root, "Result3/Organ/UKB.prot.imputed.rds"))
eid  <- prot[, 1]
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

# ---- 18 organ ----
organ18 <- c("tongue","salivary.gland","stomach","adipose.tissue","choroid.plexus",
             "retina","lung","pancreas","kidney","heart.muscle","skin","esophagus",
             "skeletal.muscle","intestine","bone.marrow","liver","brain","lymphoid.tissue")
cl <- which(colnames(odis) %in% organ18)

# ---- 8 慢病 incidence（dat_disease.rds）----
dis <- readRDS(file.path(root, "UKB_Data/dat_disease.rds"))
diseases <- c("obesity","hypertension","diabetes","cancer","mi","stroke","copd","dementia")
disease_label <- c("Obesity","Hypertension","Diabetes","Cancer","MI","Stroke","COPD","Dementia")

cov_str <- "+age+sex+bmi+race+edu+smoke+alcohol"

res <- list()
for (i in cl) {
  organ <- colnames(odis)[i]
  for (j in seq_along(diseases)) {
    st <- dis[[paste0("status_", diseases[j])]]
    tm <- dis[[paste0("time_", diseases[j])]]
    d <- data.frame(time = tm[match(eid, dis$eid)], status = st[match(eid, dis$eid)],
                    organ = odis[[organ]], age = odis$age, sex = odis$sex,
                    bmi = odis$bmi, race = odis$race, edu = odis$edu,
                    smoke = odis$smoke, alcohol = odis$alcohol)
    d <- d[!is.na(d$time) & d$time >= 0 & d$status %in% c(0, 1) & !is.na(d$organ), ]
    fl <- paste0("Surv(time,status)~organ", cov_str)
    fit <- coxph(as.formula(fl), data = d)
    s <- summary(fit)
    ci <- exp(confint(fit)[1, ])
    res[[length(res) + 1]] <- data.frame(
      organ = gsub("[.]", " ", organ),
      disease = disease_label[j],
      events_N = paste0(fit$nevent, " / ", fit$n),
      Model2_HR = round(s$coefficients[1, 2], 2),
      Model2_lo = round(ci[1], 2),
      Model2_hi = round(ci[2], 2),
      Model2_P = signif(s$coefficients[1, 5], 3),
      stringsAsFactors = FALSE
    )
  }
}
final <- do.call(rbind, res)
write.csv(final, file.path(out_dir, "organ_chronic_disease_expanded.csv"), row.names = FALSE)

cat("=== brain 慢病扩展协变量 HR ===\n")
print(final[final$organ == "brain", c("disease","Model2_HR","Model2_lo","Model2_hi","Model2_P")], row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_chronic.txt"))
writeLines(c(
  "status=PASS",
  "model2=age+sex+bmi+race+edu+smoke+alcohol",
  "outcomes=8_chronic_diseases",
  "organs=18",
  "hr_scale=per-unit",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status_chronic.txt"))
cat("\nDONE\n")
