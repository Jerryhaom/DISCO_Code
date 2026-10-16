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


# Reproduce Fig 3a (NHANES clock C-index bars) and Fig 3b (NHANES KM by clock).
# Reads the locked Manuscript/PNAS/cindex_comparison.csv and
# Result1/DNAm/NHANES.DNAm.res.rds; regenerates the same PDFs emitted by
# Result1/DNAm/NHANES.cmp.plot.R and NHANES.surv.R. No legacy files modified.

suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(survival)
  library(survminer)
  library(ggsci)
})

root <- disco_project_root()
out_dir <- "outputs/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Fig 3a: NHANES clock C-index -------------------------------------------
dat <- read.csv(file.path(root, "Manuscript/PNAS/cindex_comparison.csv"),
                stringsAsFactors = FALSE)

label_map <- c(
  mc           = "DISCO",
  hd           = "DM",
  DunedinPoAm  = "DunedinPoAm",
  GrimAgeMort  = "GrimAge",
  PhenoAge     = "PhenoAge",
  HannumAge    = "Hannum",
  HorvathAge   = "Horvath",
  VidalBraloAge = "VidalBralo",
  LinAge       = "Lin",
  ZhangAge     = "Zhang",
  SkinBloodAge = "SkinBlood",
  WeidnerAge   = "Weidner"
)

dat <- dat[dat$metric %in% names(label_map), ]
dat$label <- label_map[dat$metric]

d1 <- dat[order(dat$adj_c), ]
d1$label <- factor(d1$label, levels = d1$label)
ref_disco_adj <- dat$adj_c[dat$metric == "mc"]

p1 <- ggplot(d1, aes(adj_c, label)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.9,
           fill = ifelse(d1$label == "DISCO", "#E64B35", "grey70")) +
  geom_errorbar(aes(xmin = adj_c - adj_se, xmax = adj_c + adj_se), width = 0.4) +
  coord_cartesian(xlim = c(0.72, 0.78)) +
  geom_vline(xintercept = ref_disco_adj, color = "navy", linetype = "dashed") +
  theme_classic() +
  xlab("C-index (age-adjusted)") + ylab(NULL) +
  labs(title = "Original (age-adjusted)") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 11),
        plot.margin = margin(5, 5, 10, 5))

d2 <- dat[order(dat$accel_c), ]
d2$label <- factor(d2$label, levels = d2$label)
ref_disco_accel <- dat$accel_c[dat$metric == "mc"]

p2 <- ggplot(d2, aes(accel_c, label)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.9,
           fill = ifelse(d2$label == "DISCO", "#E64B35", "grey70")) +
  geom_errorbar(aes(xmin = accel_c - accel_se, xmax = accel_c + accel_se), width = 0.4) +
  coord_cartesian(xlim = c(0.50, 0.62)) +
  geom_vline(xintercept = ref_disco_accel, color = "red", linetype = "dashed") +
  theme_classic() +
  xlab("C-index (age-residualized)") + ylab(NULL) +
  labs(title = "Accel") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 11),
        plot.margin = margin(5, 10, 10, 5))

tgrob <- text_grob("NHANES 1999-2002", size = 13, face = "bold")

pdf(file.path(out_dir, "NHANES.cmp.pdf"), height = 5, width = 7.5)
print(ggarrange(tgrob, NULL, NULL,
                p1,   NULL, p2,
                widths = c(1, 0.05, 1),
                heights = c(0.08, 1),
                ncol = 3, nrow = 2))
dev.off()

# ---- Fig 3b: NHANES KM by clock ----------------------------------------------
DNAm  <- readRDS(file.path(root, "Result1/DNAm/NHANES.DNAm.res.rds"))
DNAm1 <- DNAm[!is.na(DNAm$HorvathAge), ]

pd.vd1 <- DNAm1
pd.vd1$group1 <- cut(pd.vd1$GrimAgeMort,
                     breaks = quantile(pd.vd1$GrimAgeMort, probs = seq(0, 1, 0.25)),
                     include.lowest = TRUE, right = TRUE)
pd.vd1$group1 <- as.numeric(pd.vd1$group1)
pd.vd1 <- pd.vd1[pd.vd1$group1 %in% c(1, 4), ]
pd.vd1$group1[pd.vd1$group1 > 1] <- 2
pd.vd1$group1 <- pd.vd1$group1 + 2

pd.vd2 <- DNAm1
pd.vd2$group1 <- cut(pd.vd2$HorvathAge,
                     breaks = quantile(pd.vd2$HorvathAge, probs = seq(0, 1, 0.25)),
                     include.lowest = TRUE, right = TRUE)
pd.vd2$group1 <- as.numeric(pd.vd2$group1)
pd.vd2 <- pd.vd2[pd.vd2$group1 %in% c(1, 4), ]
pd.vd2$group1[pd.vd2$group1 > 1] <- 2
pd.vd2$group1 <- pd.vd2$group1 + 4

pd.vd3 <- DNAm1
pd.vd3$group1 <- cut(pd.vd3$mc,
                     breaks = quantile(pd.vd3$mc, probs = seq(0, 1, 0.25)),
                     include.lowest = TRUE, right = TRUE)
pd.vd3$group1 <- as.numeric(pd.vd3$group1)
pd.vd3 <- pd.vd3[pd.vd3$group1 %in% c(1, 4), ]
pd.vd3$group1[pd.vd3$group1 > 1] <- 2

pd.vd1 <- rbind(pd.vd1, pd.vd2, pd.vd3)
pd.vd1$group1 <- as.factor(pd.vd1$group1)
levels(pd.vd1$group1) <- c(paste("DISCO", c("Low", "High"), "Risk"),
                           paste("GrimAge", c("Low", "High"), "Risk"),
                           paste("HorvathAge", c("Low", "High"), "Risk"))

pd.vd1$age1 <- pd.vd1$age + pd.vd1$time / 12
fit1 <- survfit(Surv(age, age1, status) ~ group1, data = pd.vd1)

p <- ggsurvplot(fit1, data = pd.vd1,
                conf.int = TRUE, conf.int.alpha = 0.1, censor = FALSE, pval = FALSE,
                xlab = "Age (Years)", ylab = "Survival probability",
                legend.title = "",
                legend.labs = c(paste("DISCO", c("Low", "High"), "Risk"),
                                paste("GrimAge", c("Low", "High"), "Risk"),
                                paste("HorvathAge", c("Low", "High"), "Risk")),
                title = "NHANES 1999-2002",
                xlim = c(50, 100), pval.coord = c(0.1, 0.85),
                legend = c(0.27, 0.3), palette = "npg")

pdf(file.path(out_dir, "DMAm.surv.all.pdf"), height = 4, width = 5, onefile = FALSE)
print(p$plot + theme(legend.background = element_blank(),
                     legend.text = element_text(size = 12),
                     axis.title = element_text(size = 14),
                     plot.title = element_text(size = 14)))
dev.off()

writeLines(c("status=OK",
             "figures=Fig3a,Fig3b",
             "source_objects=Manuscript/PNAS/cindex_comparison.csv,Result1/DNAm/NHANES.DNAm.res.rds"),
           file.path(out_dir, "status_fig3ab.txt"))
cat("OK\n")
