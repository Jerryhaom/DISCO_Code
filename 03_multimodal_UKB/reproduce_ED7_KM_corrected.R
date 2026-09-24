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

# Corrected version of PNAS ED7 chronic-disease KM:
# plots DISEASE INCIDENCE (status_<disease>) instead of all-cause mortality (status).
suppressMessages({library(dplyr); library(survival); library(survminer)})

root <- disco_project_root()
outdir <- disco_results_path("03_multimodal_UKB/outputs/ED7_KM")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

disease <- c("obesity", "hypertension", "diabetes", "cancer", "mi", "stroke", "copd", "dementia")
titles  <- c("Obesity", "Hypertension", "Diabetes", "Cancer",
             "Myocardial infarction", "Stroke", "COPD", "Dementia")

run_corrected <- function(disco_path, label) {
  pd <- readRDS(disco_path)
  dd <- readRDS(file.path(root, "UKB_Data/dat_disease.rds"))
  dd <- dd[match(pd$ID, dd$eid), ]
  pd <- cbind(pd, dd[-1])

  pd$mcg <- as.numeric(cut(pd$mc, breaks = quantile(pd$mc), include.lowest = TRUE, right = TRUE))
  pd$hdg <- as.numeric(cut(pd$hd, breaks = quantile(pd$hd), include.lowest = TRUE, right = TRUE))
  pd <- pd %>% mutate(Sig = case_when(mcg == 1 & hdg == 1 ~ "LL",
                                      mcg == 4 & hdg == 4 ~ "HH",
                                      mcg == 1 & hdg == 4 ~ "LH",
                                      mcg == 4 & hdg == 1 ~ "HL"))
  pd <- pd %>% mutate(Sig_num = case_when(Sig == "LL" ~ 1, Sig == "LH" ~ 2,
                                          Sig == "HL" ~ 3, Sig == "HH" ~ 4))
  pd <- pd[pd$Sig %in% c("LL", "HH", "LH", "HL"), ]

  pdf(file.path(outdir, paste0(label, ".disease.CORRECTED.pdf")), height = 3, width = 4)
  for (j in disease) {
    tv <- paste0("time_", j)
    sv <- paste0("status_", j)
    pdd <- pd[pd[[tv]] >= 0, ]
    pdd$start <- pdd$age
    pdd$stop  <- pdd$age + pdd[[tv]]
    pdd$event <- pdd[[sv]]
    fit <- survfit(Surv(start, stop, event) ~ Sig_num, data = pdd)
    p <- ggsurvplot(fit, data = pdd, conf.int = TRUE, censor = FALSE, pval = FALSE,
                    xlab = "Age (Years)", ylab = "Disease-free survival probability",
                    legend.title = "",
                    legend.labs = c("DISCO Low Risk, DM Low Risk", "DISCO Low Risk, DM High Risk",
                                    "DISCO High Risk, DM Low Risk", "DISCO High Risk, DM High Risk"),
                    title = titles[match(j, disease)], ylim = c(0, 1), xlim = c(40, max(fit$time)),
                    legend = c(0.4, 0.35), palette = "npg")
    print(p$plot + theme(legend.background = element_blank()))
  }
  dev.off()
  cat("wrote", paste0(label, ".disease.CORRECTED.pdf"), "\n")
}

run_corrected(file.path(root, "Result1/UKB.Clinical.DISCO.rds"), "UKB.clinical")
run_corrected(file.path(root, "Result2/UKB.Prot.DISCO.rds"), "UKB.prot")
run_corrected(file.path(root, "Result2/UKB.ClinicalProt.DISCO.rds"), "UKB.clinical.prot")
