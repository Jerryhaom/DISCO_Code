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

# Reproduce PNAS Extended Data Fig. 7 (chronic-disease KM panels), faithful to
# Supple/Result2/UKB.{clinical,prot,clinical.prot}.disease.R (paths adapted).
suppressMessages({library(dplyr); library(survival); library(survminer)})

root <- disco_project_root()
outdir <- disco_results_path("03_multimodal_UKB/outputs/ED7_KM")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

disease <- c("obesity", "hypertension", "diabetes", "cancer", "mi", "stroke", "copd", "dementia")
titles  <- c("Obesity", "Hypertension", "Diabetes", "Cancer",
             "Myocardial infarction", "Stroke", "COPD", "Dementia")

run_modality <- function(disco_path, label) {
  pd <- readRDS(disco_path)
  dd <- readRDS(file.path(root, "UKB_Data/dat_disease.rds"))
  dd <- dd[match(pd$ID, dd$eid), ]
  cat(label, ": unmatched eid =", sum(is.na(dd$eid)), "\n")
  pd <- cbind(pd, dd[-1])

  pd$mcg <- as.numeric(cut(pd$mc, breaks = quantile(pd$mc), include.lowest = TRUE, right = TRUE))
  pd$hdg <- as.numeric(cut(pd$hd, breaks = quantile(pd$hd), include.lowest = TRUE, right = TRUE))
  pd <- pd %>% mutate(Sig = case_when(mcg == 1 & hdg == 1 ~ "LL",
                                      mcg == 4 & hdg == 4 ~ "HH",
                                      mcg == 1 & hdg == 4 ~ "LH",
                                      mcg == 4 & hdg == 1 ~ "HL"))
  cat(label, "Sig table:\n"); print(table(pd$Sig))

  process_disease <- function(time_var, title) {
    pdd <- pd[pd$Sig %in% c("LL", "HH", "LH", "HL"), ]
    pdd <- pdd %>% mutate(Sig_num = case_when(Sig == "LL" ~ 1, Sig == "LH" ~ 2,
                                              Sig == "HL" ~ 3, Sig == "HH" ~ 4))
    pdd1 <- pdd[pdd[[time_var]] >= 0, ]
    pdd1$time1 <- pdd1$time + pdd1$age
    fit1 <- survfit(Surv(age, time1, status) ~ Sig_num, data = pdd1)
    p <- ggsurvplot(fit1, data = pdd1, conf.int = TRUE, censor = FALSE, pval = FALSE,
                    xlab = "Age (Years)", ylab = "Survival probability", legend.title = "",
                    legend.labs = c("DISCO Low Risk, DM Low Risk", "DISCO Low Risk, DM High Risk",
                                    "DISCO High Risk, DM Low Risk", "DISCO High Risk, DM High Risk"),
                    title = title, ylim = c(0, 1), xlim = c(40, max(fit1$time)),
                    pval.coord = c(0.1, 0.5), legend = c(0.4, 0.35), palette = "npg")
    p$plot + theme(legend.background = element_blank())
  }

  pdf(file.path(outdir, paste0(label, ".disease.pdf")), height = 3, width = 4)
  for (j in seq_along(disease)) {
    print(process_disease(paste0("time_", disease[j]), titles[j]))
  }
  dev.off()
  cat("wrote", paste0(label, ".disease.pdf"), "\n")
}

run_modality(file.path(root, "Result1/UKB.Clinical.DISCO.rds"), "UKB.clinical")
run_modality(file.path(root, "Result2/UKB.Prot.DISCO.rds"), "UKB.prot")
run_modality(file.path(root, "Result2/UKB.ClinicalProt.DISCO.rds"), "UKB.clinical.prot")
