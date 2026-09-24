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

# Reproduce PNAS DISCO-aging Fig 2c-f (RuLAS Kaplan-Meier) and ED9f (organ-system C-index).
# Reads LOCKED legacy saved objects from DISCO_DATA_ROOT/Result6.
# Does NOT modify any legacy file.

suppressMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggsci)
  library(ggthemes)
  library(RColorBrewer)
})

LEGACY <- disco_path("Result6")
OUT <- file.path("outputs", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

read_legacy <- function(f) {
  p <- file.path(LEGACY, f)
  if (!file.exists(p)) stop("missing legacy file: ", f)
  readRDS(p)
}

# ---- Common KM theme (from Result6/Plot_surv.R lines 114-120) ----
km_theme <- theme(
  legend.background = element_blank(),
  legend.text       = element_text(size = 18),
  legend.key.height = unit(1.5, "lines"),
  axis.text         = element_text(size = 14),
  axis.title.y      = element_text(size = 20),
  axis.title.x      = element_text(size = 20),
  plot.title        = element_text(size = 22)
)

# ---- Build the 4-curve age-axis KM (Plot_surv.R lines 64-122 / 178-235 / 295-352) ----
# pd must contain columns: age, time, status, mc (DISCO), hd (DM/Entropy)
make_km4 <- function(pd, title, outfile) {
  pd.vd1 <- pd
  pd.vd1$group1 <- cut(pd.vd1$hd,
                       breaks = quantile(pd.vd1$hd, probs = seq(0, 1, 0.25), na.rm = TRUE),
                       include.lowest = TRUE, right = TRUE)
  pd.vd1$group1 <- as.numeric(pd.vd1$group1)
  pd.vd1 <- pd.vd1[pd.vd1$group1 %in% c(1, 4), ]
  pd.vd1$group1[pd.vd1$group1 > 1] <- 2
  pd.vd1$group1 <- pd.vd1$group1 + 2                       # 3 = DM/Entropy Low, 4 = High

  pd.vd2 <- pd
  pd.vd2$group1 <- cut(pd.vd2$mc,
                       breaks = quantile(pd.vd2$mc, probs = seq(0, 1, 0.25), na.rm = TRUE),
                       include.lowest = TRUE, right = TRUE)
  pd.vd2$group1 <- as.numeric(pd.vd2$group1)
  pd.vd2 <- pd.vd2[pd.vd2$group1 %in% c(1, 4), ]
  pd.vd2$group1[pd.vd2$group1 > 1] <- 2                    # 1 = DISCO Low, 2 = High

  pd.vd1 <- rbind(pd.vd1, pd.vd2)

  fit1 <- survfit(Surv(age, age + time / 12, status) ~ group1,
                  data = pd.vd1[pd.vd1$age >= 70, ])

  p <- ggsurvplot(
    fit1,
    data           = pd.vd1,
    conf.int       = TRUE,
    conf.int.alpha = 0.1,
    censor         = TRUE,
    pval           = FALSE,
    xlab           = "Age (years)",
    ylab           = "Survival probability",
    legend.title   = "",
    legend.labs    = c("DISCO Low Risk", "DISCO High Risk",
                       "DM/Entropy Low Risk", "DM/Entropy High Risk"),
    title          = title,
    xlim           = c(70, max(fit1$time)),
    break.x.by     = 5,
    legend         = c(0.28, 0.3),
    palette        = "npg"
  )
  p1 <- p$plot + km_theme

  pdf(outfile, height = 4, width = 6, onefile = FALSE)
  print(p1)
  dev.off()
  invisible(fit1)
}

# ---- Fig 2c: RuLAS, Clinical (Plot_surv.R lines 8-18 + 64-122) ----
pd <- read_legacy("Rugao.Phe.DISCO2.rds")          # has mc + hd
mortality <- read_legacy("Rugao.2019.mortality.rds")
mm <- mortality[match(pd$ID, mortality$x2019ID), ]
mm$time[is.na(mm$status)]   <- 24
mm$status[is.na(mm$status)] <- 0
pd$time   <- mm$time
pd$status <- mm$status
pd <- pd[pd$age >= 60, ]
make_km4(pd, "RuLAS, Clinical", file.path(OUT, "Fig2c_RuLAS_Clinical.pdf"))

# ---- Fig 2d: RuLAS, Metabolomics (Plot_surv.R lines 178-235) ----
pd <- read_legacy("Rugao.Met.surv.rds")            # mc + hd + time + status + age
make_km4(pd, "RuLAS, Metabolomics", file.path(OUT, "Fig2d_RuLAS_Metabolomics.pdf"))

# ---- Fig 2e: RuLAS, Microbiota (Plot_surv.R lines 295-352) ----
pd <- read_legacy("Rugao.Meta.surv.rds")           # mc + hd + time + status + age
make_km4(pd, "RuLAS, Microbiota", file.path(OUT, "Fig2e_RuLAS_Microbiota.pdf"))

# ---- Fig 2f: RuLAS, Aggregate (Plot_subsystems2.R lines 274-367) ----
met   <- read_legacy("Rugao.Met.surv.rds")
meta  <- read_legacy("Rugao.Meta.surv.rds")
disco.sub <- read_legacy("Rugao.subsystems.rds")   # mc per subsystem (ID + All + 8 subsystems)
hd.sub    <- read_legacy("Rugao.subsystems2.rds")  # hd per subsystem (ID, Age, Sex, All + 8 subsystems)

samples <- disco.sub$ID
sub2    <- hd.sub[match(samples, hd.sub$ID), ]     # Age/Sex aligned by ID
met1    <- met[match(samples, met$ID), ]
meta1   <- meta[match(samples, meta$ID), ]
mort    <- mortality[match(samples, mortality$x2019ID), ]
mort$time[is.na(mort$status)]   <- 24
mort$status[is.na(mort$status)] <- 0

# DISCO score = mean of 10 systems (8 clinical subsystems + Met + Meta)  [disco.sub[,3:12]]
disco.sub <- disco.sub[match(samples, disco.sub$ID), ]
disco.sub$Met    <- met1$mc
disco.sub$Meta   <- meta1$mc
disco.sub$Age    <- as.numeric(sub2$Age)
disco.sub$Sex    <- sub2$Sex
disco.sub$status <- mort$status
disco.sub$time   <- mort$time
disco.sub$var    <- apply(disco.sub[, 3:12], 1, function(x) mean(x))
disco.sub$risk1  <- as.numeric(as.factor(cut(disco.sub$var,
                                             quantile(disco.sub$var, right = TRUE,
                                                      include.lowest = TRUE, na.rm = TRUE))))

# DM/Entropy score = mean of 10 systems (8 clinical hd + Met hd + Meta hd) [hd.sub[,5:14]]
hd.sub <- hd.sub[match(samples, hd.sub$ID), ]
hd.sub$Met    <- met1$hd
hd.sub$Meta   <- meta1$hd
hd.sub$Age    <- as.numeric(sub2$Age)
hd.sub$Sex    <- sub2$Sex
hd.sub$status <- mort$status
hd.sub$time   <- mort$time
hd.sub$var    <- apply(hd.sub[, 5:14], 1, function(x) mean(x))
hd.sub$risk1  <- as.numeric(as.factor(cut(hd.sub$var,
                                          quantile(hd.sub$var, right = TRUE,
                                                   include.lowest = TRUE, na.rm = TRUE))))

hd.sub2    <- hd.sub[hd.sub$risk1 %in% c(1, 4), ]
hd.sub2$risk1[hd.sub2$risk1 < 4] <- 3
disco.sub2 <- disco.sub[disco.sub$risk1 %in% c(1, 4), ]
disco.sub2$risk1[disco.sub2$risk1 > 1] <- 2
pd.vd1 <- rbind(disco.sub2, hd.sub2)

fit1 <- survfit(Surv(Age, Age + time / 12, status) ~ risk1,
                data = pd.vd1[pd.vd1$Age >= 70, ])
p <- ggsurvplot(
  fit1,
  data           = pd.vd1,
  conf.int       = TRUE,
  conf.int.alpha = 0.1,
  censor         = TRUE,
  pval           = FALSE,
  xlab           = "Age (years)",
  ylab           = "Survival probability",
  legend.title   = "",
  legend.labs    = c("DISCO Low Risk", "DISCO High Risk",
                     "DM/Entropy Low Risk", "DM/Entropy High Risk"),
  title          = "RuLAS, Aggregate",
  xlim           = c(70, max(fit1$time)),
  break.x.by     = 5,
  legend         = c(0.28, 0.3),
  palette        = "npg"
)
p1 <- p$plot + km_theme
pdf(file.path(OUT, "Fig2f_RuLAS_Aggregate.pdf"), height = 4, width = 6, onefile = FALSE)
print(p1)
dev.off()

# ---- ED9f: organ-system C-index bar (Plot_subsystems.R lines 122-266) ----
sub <- read_legacy("Rugao.subsystems.rds")
samples <- sub$ID
age  <- as.numeric(hd.sub[match(samples, hd.sub$ID), "Age"])
sex  <- hd.sub[match(samples, hd.sub$ID), "Sex"]
mort <- mortality[match(samples, mortality$x2019ID), ]
mort$time[is.na(mort$status)]   <- 24
mort$status[is.na(mort$status)] <- 0

pd <- data.frame(ID = samples, age = age, sex = sex,
                 sub[, 2:ncol(sub)], time = mort$time, status = mort$status)

ci <- matrix(0, ncol = 6, nrow = 13)

# row 1: Age
s <- summary(coxph(Surv(time, status) ~ age, data = pd))
ci[1, 1:2] <- s$concordance

# rows 2-10: All (Clinical) + 8 subsystems
discoll <- data.frame(sub)
k <- 2
for (i in colnames(discoll)[2:ncol(discoll)]) {
  s <- summary(coxph(as.formula(paste("Surv(time,status)~", i)), data = pd))
  ci[k, 1:2] <- s$concordance
  k <- k + 1
}

# row 11: Metabolomics
s <- summary(coxph(Surv(time, status) ~ mc, data = met))
ci[11, 1:2] <- s$concordance

# row 12: Microbiota
s <- summary(coxph(Surv(time, status) ~ mc, data = meta))
ci[12, 1:2] <- s$concordance

# row 13: Aggregate (mean of 10 systems, mc)
disco.sub <- sub
disco.sub$Met    <- met$mc[match(samples, met$ID)]
disco.sub$Meta   <- meta$mc[match(samples, meta$ID)]
disco.sub$Age    <- age
disco.sub$Sex    <- sex
disco.sub$status <- mort$status
disco.sub$time   <- mort$time
disco.sub$var    <- apply(disco.sub[, 3:12], 1, function(x) mean(x))
s <- summary(coxph(Surv(time, status) ~ var, data = disco.sub))
ci[13, 1:2] <- s$concordance

ci <- data.frame(ci)
rownames(ci) <- c("Age", colnames(discoll)[2:ncol(discoll)],
                  "Metabolomics", "Microbiota", "Aggregate")
rownames(ci)[10] <- "Vascular Health"
rownames(ci)[2]  <- "Clinical"

pb <- data.frame(phe = rownames(ci), ci)
pb$phe <- factor(pb$phe, levels = pb$phe[order(pb$X1)])

p <- ggplot(pb, aes(X1, phe, fill = phe)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.9) +
  geom_errorbar(aes(xmin = X1 - X2, xmax = X1 + X2), width = 0.5) +
  scale_y_discrete(labels = gsub("[.]", " ", levels(pb$phe))) +
  coord_cartesian(xlim = c(0.5, 0.78)) +
  scale_fill_manual(values = c(brewer.pal(8, "Set3"), brewer.pal(5, "Set1"))) +
  theme_few() +
  ylab("") + labs(title = "RuLAS") +
  xlab("C-index") +
  theme(legend.position = "None",
        axis.text.y  = element_text(size = 12),
        axis.text.x  = element_text(size = 12),
        axis.title   = element_text(size = 12))

ggsave(file.path(OUT, "ED9f_RuLAS_Cindex.pdf"), p, width = 3, height = 4)

cat("Figures written to", normalizePath(OUT), "\n")
