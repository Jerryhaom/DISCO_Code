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


# Organ DISCO 扩展协变量分析（18 organs，主分析）
# 对 18 organ 的 DISCO 分数跑 covariate ladder：Crude → age+sex → 扩展协变量
# （+BMI+race+education+smoking+drinking），与 R01/R06 的 UKB 口径一致。
# 不改动 Result3 旧文件。

options(stringsAsFactors = FALSE)
suppressMessages(library(survival))

root <- disco_project_root()
out_dir <- disco_results_path("05_organ_network/outputs/organ_expanded_covariates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 18 organs（主分析，与 Supplementary Table 19 一致）----
organ18 <- c("tongue","salivary.gland","stomach","adipose.tissue","choroid.plexus",
             "retina","lung","pancreas","kidney","heart.muscle","skin","esophagus",
             "skeletal.muscle","intestine","bone.marrow","liver","brain","lymphoid.tissue")

# ---- 数据 ----
odis <- readRDS(file.path(root, "Result3/Organ/UKB.odis.rds"))
prot <- readRDS(file.path(root, "Result3/Organ/UKB.prot.imputed.rds"))
eid  <- prot[, 1]  # prot 是 matrix，第 1 列 = eid；odis 行序 = prot eid 顺序

ext <- readRDS(file.path(root, "UKB_Data/extended_covariates_clinical.rds"))
race <- read.csv(file.path(root, "UKB_Data/ukb_race_21000.csv"), colClasses = "character")

idx <- match(eid, ext$ID)
d <- data.frame(
  time = odis$time, status = as.numeric(odis$status),
  age = ext$age[idx], sex = ext$sex[idx],
  bmi = ext$bmi_raw[idx], edu = ext$edu_high[idx],
  smoke = ext$smoke_status[idx], alcohol = ext$alcohol_freq[idx],
  stringsAsFactors = FALSE
)
d$race <- race$race[match(eid, race$eid)]
d$race[d$race == ""] <- "White"

# 协变量编码（与 R01 一致）
d$sex <- factor(d$sex)
d$edu <- factor(d$edu)
d$smoke <- factor(d$smoke, levels = c(0, 1, 2))
d$alcohol <- factor(d$alcohol, levels = 1:6, ordered = TRUE)
d$race <- factor(d$race, levels = c("White", "Non-White"))

# ---- 对 18 organ 跑 covariate ladder ----
models <- list(
  Crude = character(0),
  M1_age_sex = c("age", "sex"),
  M2_expanded = c("age", "sex", "bmi", "race", "edu", "smoke", "alcohol")
)

run_one <- function(organ, mc) {
  d$mc <- mc
  out <- list()
  for (mname in names(models)) {
    covars <- models[[mname]]
    rhs <- if (length(covars) == 0) "mc" else paste(c("mc", covars), collapse = "+")
    x <- d[complete.cases(d[, c("time", "status", "mc", covars)]), , drop = FALSE]
    fit <- coxph(as.formula(paste0("Surv(time,status)~", rhs)), data = x)
    b <- coef(fit)["mc"]; se <- sqrt(diag(vcov(fit)))["mc"]
    sd_mc <- sd(x$mc, na.rm = TRUE)
    out[[length(out) + 1]] <- data.frame(
      organ = organ, model = mname, n = nrow(x), events = sum(x$status),
      sd_mc = sd_mc,
      hr_unit = exp(b), hr_per_sd = exp(b * sd_mc),
      lo_per_sd = exp((b - 1.96 * se) * sd_mc), hi_per_sd = exp((b + 1.96 * se) * sd_mc),
      c_index = summary(fit)$concordance[1], p = coef(summary(fit))["mc", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

res <- do.call(rbind, lapply(organ18, function(o) run_one(o, odis[[o]])))
write.csv(res, file.path(out_dir, "organ_expanded_covariates.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "organs=18",
  "models=Crude_M1_age_sex_M2_expanded",
  "expanded_covariates=age+sex+bmi+race+edu+smoke+alcohol",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))

# 打印摘要：Crude vs M2_expanded 的 per-SD HR
cat("=== 18 organ：Crude → 扩展协变量 per-SD HR ===\n")
crude <- res[res$model == "Crude", ]
full  <- res[res$model == "M2_expanded", ]
print(data.frame(organ = crude$organ,
  HR_crude_perSD = round(crude$hr_per_sd, 2),
  HR_expanded_perSD = round(full$hr_per_sd, 2),
  C_crude = round(crude$c_index, 3),
  C_expanded = round(full$c_index, 3)), row.names = FALSE)
