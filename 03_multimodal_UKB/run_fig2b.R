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


# Regenerate Fig 2b (UKB Clinical+Proteomics KM) with the SAME submitted
# layout/title as before ("UKB, Clinical+Proteomics").
# Replicates the exact plotting block from Result2/Cal_ClinicalProt.R:221-276.
# The ClinicalProt object is already age<=40 reference (unchanged by the
# proteomics-only reference fix), so the original locked object is used.
# No legacy files modified.

suppressMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggsci)
})

root <- disco_project_root()
inp  <- file.path(root, "Result2/UKB.ClinicalProt.DISCO.rds")
out  <- disco_results_path("03_multimodal_UKB/outputs/03_prot_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pd  <- readRDS(inp)
dd3 <- pd

dd4 <- dd3
dd4$Group <- cut(dd4$mc,
                 breaks = quantile(dd4$mc, probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE, right = TRUE)
dd4 <- dd4[as.numeric(dd4$Group) %in% c(1, 4), ]
dd4$Group <- as.numeric(dd4$Group)
dd4$Group[dd4$Group > 1] <- 2

dd5 <- dd3
dd5$Group <- cut(dd5$hd,
                 breaks = quantile(dd5$hd, probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE, right = TRUE)
dd5 <- dd5[as.numeric(dd5$Group) %in% c(1, 4), ]
dd5$Group <- as.numeric(dd5$Group)
dd5$Group[dd5$Group > 1] <- 2
dd5$Group <- dd5$Group + 2

dd4 <- rbind(dd4, dd5)

fit1 <- survfit(Surv(age, age + time, status) ~ Group, data = dd4)

pdf(file.path(out, "Fig2b_UKB_ClinicalProteomics.pdf"), height = 3, width = 4, onefile = FALSE)
p <- ggsurvplot(fit1, data = dd4,
                conf.int = TRUE, censor = FALSE, pval = FALSE,
                xlab = "Age (Years)",
                ylab = "Survival probability",
                legend.title = "",
                legend.labs = c(paste("DISCO", c("Low", "High"), "Risk"),
                                paste("DM/Entropy", c("Low", "High"), "Risk")),
                title = "UKB, Clinical+Proteomics",
                ylim = c(0.3, 1),
                xlim = c(39, max(fit1$time)),
                pval.coord = c(0.1, 0.8),
                legend = c(0.32, 0.35),
                palette = "npg")
print(p$plot + theme(legend.background = element_blank(),
                     legend.text = element_text(size = 12),
                     axis.text = element_text(size = 10),
                     axis.title.y = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     plot.title = element_text(size = 14)))
dev.off()

writeLines(c(
  "status=OK",
  "title=UKB, Clinical+Proteomics",
  "reference=age<=40 (unchanged; ClinicalProt already age<=40)",
  "data=Result2/UKB.ClinicalProt.DISCO.rds (N=47173)",
  "layout=same_as_submitted_Cal_ClinicalProt.R:221-276"
), file.path(out, "Fig2b_status.txt"))
cat("OK: Fig2b_UKB_ClinicalProteomics.pdf (title 'UKB, Clinical+Proteomics')\n")
