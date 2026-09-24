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
disco_init_module("02_primary_clinical")


# Reproduce Fig 1c-d and ED2a (NHANES clinical) from the locked legacy objects.
# Reads Result1/NHANES.Clinical.DISCO.rds and NHANES.Clinical.DISCO2.rds and
# regenerates the same PDFs emitted by Result1/Cal_NHANES.DISCO.R (plotting section).
# No legacy files are modified.

suppressMessages({
  library(ggplot2)
  library(ggpointdensity)
  library(ggthemes)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(ggsci)
})

root <- disco_project_root()
out_dir <- "outputs/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pd  <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO.rds"))
dd2 <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO2.rds"))

# ---- Fig 1c left: DISCO vs Age ------------------------------------------------
pdf(file.path(out_dir, "NHANES.DISCO.Age.pdf"), height = 3, width = 3)
print(
  ggplot(pd, aes(Age, mc)) +
    geom_pointdensity(size = 0.5) +
    scale_color_distiller(palette = "Spectral") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 0.95) +
    scale_y_continuous(limits = c(-5, 10), breaks = c(-5, 0, 5, 10)) +
    labs(title = "NHANES") +
    theme_few() + xlab("Age (Years)") + ylab("DISCO") +
    theme(legend.position = "None")
)
dev.off()

# ---- Fig 1c right: DISCO vs DM (log) ------------------------------------------
pdf(file.path(out_dir, "NHANES.DISCO.HD.pdf"), height = 3, width = 3)
print(
  ggplot(pd, aes(log(hd), mc)) +
    geom_pointdensity(size = 0.5) +
    scale_color_distiller(palette = "Spectral") +
    geom_smooth(method = "lm", color = "black") +
    scale_y_continuous(limits = c(-5, 10), breaks = c(-5, 0, 5, 10)) +
    labs(title = "NHANES") +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 0.95) +
    theme_few() + xlab("DM (Log-transformed)") + ylab("DISCO") +
    theme(legend.position = "None")
)
dev.off()

# ---- ED2a: Entropy distance vs DM (log) ---------------------------------------
pdf(file.path(out_dir, "NHANES.Entropy.HD.pdf"), height = 3, width = 3)
print(
  ggplot(pd, aes(log(hd), en)) +
    geom_pointdensity(size = 0.5) +
    scale_color_distiller(palette = "Spectral") +
    labs(title = "NHANES") +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 0.95) +
    theme_few() + xlab("DM (Log-transformed)") + ylab("Entropy distance") +
    theme(legend.position = "None")
)
dev.off()

# ---- Fig 1d: KM by DISCO/DM-Entropy risk (quartiles) --------------------------
dd3 <- dd2
dd3$age <- pd$Age[match(dd3$id, pd$id)]

dd4 <- dd3
res <- dd4$mc
dd4$Group <- cut(res, breaks = quantile(res, probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE, right = TRUE)
dd4 <- dd4[as.numeric(dd4$Group) %in% c(1, 4), ]
dd4$Group <- as.numeric(dd4$Group)
dd4$Group[dd4$Group > 1] <- 2

dd5 <- dd3
res <- dd5$hd
dd5$Group <- cut(res, breaks = quantile(res, probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE, right = TRUE)
dd5 <- dd5[as.numeric(dd5$Group) %in% c(1, 4), ]
dd5$Group <- as.numeric(dd5$Group)
dd5$Group[dd5$Group > 1] <- 2
dd5$Group <- dd5$Group + 2

dd4 <- rbind(dd4, dd5)
dd4$time1 <- dd4$time / 12 + dd4$age

fit1 <- survfit(Surv(age, time1, status) ~ Group, data = dd4)

pdf(file.path(out_dir, "NHANES.DISCO.all2.pdf"), height = 3, width = 4, onefile = FALSE)
p <- ggsurvplot(fit1, data = dd4,
                conf.int = TRUE, censor = FALSE, pval = FALSE,
                xlab = "Age (Years)", ylab = "Survival probability",
                legend.title = "",
                legend.labs = c(paste("DISCO", c("Low", "High"), "Risk"),
                                paste("DM/Entropy", c("Low", "High"), "Risk")),
                title = "NHANES",
                ylim = c(0, 1), xlim = c(20, max(fit1$time)),
                pval.coord = c(0.1, 0.5), legend = c(0.3, 0.35),
                palette = "npg")
print(p$plot + theme(legend.background = element_blank()))
dev.off()

writeLines(c("status=OK",
             "figures=Fig1c_left,Fig1c_right,ED2a,Fig1d",
             "source_objects=Result1/NHANES.Clinical.DISCO.rds,Result1/NHANES.Clinical.DISCO2.rds"),
           file.path(out_dir, "status.txt"))
cat("OK\n")
