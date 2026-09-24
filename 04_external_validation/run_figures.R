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
disco_init_module("04_external_validation")


# Reproduce Extended Data figures from LOCKED legacy saved objects.
#   ED3 (CHARLS b / RuLAS c / CLHLS d): DISCO~Age, DISCO~log(DM), KM by risk groups
#   ED6 (HRS): C-index bar charts (Original/Accel) + KM by DISCO/GrimAge/HorvathAge
# Reads only locked legacy objects; writes outputs/figures/*.pdf.
# No legacy file is modified.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(ggpointdensity)
  library(survival)
  library(survminer)
  library(ggsci)
})

root    <- disco_project_root()
out_dir <- disco_results_path("04_external_validation/outputs/figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

missing <- character(0)
record_missing <- function(panel, reason) {
  missing <<- c(missing, sprintf("- %s: %s", panel, reason))
}

save_pdf <- function(plot, path, width, height) {
  pdf(path, width = width, height = height, onefile = FALSE)
  print(plot)
  dev.off()
  cat("WROTE", path, "\n")
}

# ------------------------------------------------------------------------------
# Shared: quartile risk-group builder (Q1 vs Q4), identical to legacy scripts.
# ------------------------------------------------------------------------------
build_km_groups <- function(df) {
  pd.vd1 <- df
  pd.vd1$group1 <- cut(pd.vd1$hd,
                       breaks = quantile(pd.vd1$hd, probs = seq(0, 1, 0.25)),
                       include.lowest = TRUE, right = TRUE)
  pd.vd1$group1 <- as.numeric(pd.vd1$group1)
  pd.vd1 <- pd.vd1[pd.vd1$group1 %in% c(1, 4), ]
  pd.vd1$group1[pd.vd1$group1 > 1] <- 2
  pd.vd1$group1 <- pd.vd1$group1 + 2

  pd.vd2 <- df
  pd.vd2$group1 <- cut(pd.vd2$mc,
                       breaks = quantile(pd.vd2$mc, probs = seq(0, 1, 0.25)),
                       include.lowest = TRUE, right = TRUE)
  pd.vd2$group1 <- as.numeric(pd.vd2$group1)
  pd.vd2 <- pd.vd2[pd.vd2$group1 %in% c(1, 4), ]
  pd.vd2$group1[pd.vd2$group1 > 1] <- 2

  pd.vd1 <- rbind(pd.vd1, pd.vd2)
  pd.vd1$group1 <- as.factor(pd.vd1$group1)
  levels(pd.vd1$group1) <- c(paste("DISCO", c("Low", "High"), "Risk"),
                             paste("MD",    c("Low", "High"), "Risk"))
  pd.vd1
}

# ------------------------------------------------------------------------------
# ED3 helper plots
# ------------------------------------------------------------------------------
disco_age_plot <- function(df, title, x_breaks = NULL, y_breaks = NULL) {
  p <- ggplot(df, aes(age, mc)) +
    geom_pointdensity(size = 1) +
    scale_color_distiller(palette = "Spectral") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "spearman", cor.coef.name = c("rho"), label.y.npc = 0.95) +
    labs(title = title) +
    theme_few() + xlab("Age (Years)") + ylab("DISCO") +
    theme(legend.position = "None")
  if (!is.null(x_breaks)) p <- p + scale_x_continuous(breaks = x_breaks)
  if (!is.null(y_breaks)) p <- p + scale_y_continuous(breaks = y_breaks)
  p
}

disco_hd_plot <- function(df, title, xname = "hd", y_breaks = NULL) {
  p <- ggplot(df, aes(.data[[xname]], mc)) +
    geom_pointdensity(size = 1) +
    scale_color_distiller(palette = "Spectral") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "spearman", cor.coef.name = c("rho"), label.y.npc = 0.95) +
    labs(title = title) +
    theme_few() + xlab("DM (Log-transformed)") + ylab("DISCO") +
    theme(legend.position = "None")
  if (!is.null(y_breaks)) p <- p + scale_y_continuous(breaks = y_breaks)
  p
}

km_age_plot <- function(pd, title, xlim, x_breaks = NULL, legend_pos = c(0.27, 0.25)) {
  fit <- survfit(Surv(age, age + time, status) ~ group1, data = pd)
  p <- ggsurvplot(fit,
                  data = pd,
                  conf.int = TRUE,
                  conf.int.alpha = 0.1,
                  censor = TRUE,
                  pval = FALSE,
                  xlab = "Age (Years)",
                  ylab = "Survival probability",
                  legend.title = "",
                  legend.labs = c(paste("DISCO", c("Low", "High"), "Risk"),
                                  paste("DM/Entropy", c("Low", "High"), "Risk")),
                  title = title,
                  ylim = c(0, 1),
                  xlim = xlim,
                  break.x.by = 10,
                  risk.table = FALSE,
                  legend = legend_pos,
                  palette = "npg")
  p1 <- p$plot + theme(legend.background = element_blank(),
                       legend.text = element_text(size = 14),
                       axis.text = element_text(size = 14),
                       axis.title.y = element_text(size = 18),
                       axis.title.x = element_text(size = 18),
                       plot.title = element_text(size = 20))
  if (!is.null(x_breaks)) p1 <- p1 + scale_x_continuous(breaks = x_breaks)
  p1
}

# ------------------------------------------------------------------------------
# ED3b - CHARLS
# ------------------------------------------------------------------------------
charls_file <- file.path(root, "Result4/CHARLS.DISCO.rds")
if (file.exists(charls_file)) {
  pd.vd <- readRDS(charls_file)
  save_pdf(disco_age_plot(pd.vd, "CHARLS", x_breaks = seq(50, 100, 10)),
           file.path(out_dir, "CHARLS.DISCO.age.pdf"), 3.5, 3)
  save_pdf(disco_hd_plot(transform(pd.vd, hd = log(hd)), "CHARLS"),
           file.path(out_dir, "CHARLS.DISCO.hd.pdf"), 3.5, 3)
  kmd <- build_km_groups(pd.vd)
  save_pdf(km_age_plot(kmd, "CHARLS", xlim = c(45, 105),
                       x_breaks = seq(45, 105, 20), legend_pos = c(0.3, 0.25)),
           file.path(out_dir, "CHARLS.KM.pdf"), 5, 4)
} else {
  record_missing("CHARLS (ED3b)", "CHARLS.DISCO.rds not found")
}

# ------------------------------------------------------------------------------
# ED3c - RuLAS
# ------------------------------------------------------------------------------
rulas_file <- file.path(root, "Result4/RuLAS.DISCO.rds")
if (file.exists(rulas_file)) {
  pd.vd <- readRDS(rulas_file)
  pd.vd <- pd.vd[pd.vd$age > 60, ]
  pd.vd <- pd.vd[pd.vd$status %in% c(0, 1), ]
  save_pdf(disco_age_plot(pd.vd, "RuLAS"),
           file.path(out_dir, "RuLAS.DISCO.age.pdf"), 3.5, 3)
  save_pdf(disco_hd_plot(transform(pd.vd, hd = log(hd)), "RuLAS"),
           file.path(out_dir, "RuLAS.DISCO.hd.pdf"), 3.5, 3)
  kmd <- build_km_groups(pd.vd)
  kmd <- kmd[kmd$age >= 70, ]
  save_pdf(km_age_plot(kmd, "RuLAS", xlim = c(73, 96),
                       x_breaks = seq(75, 95, 5), legend_pos = c(0.27, 0.25)),
           file.path(out_dir, "RuLAS.KM.pdf"), 5, 4)
} else {
  record_missing("RuLAS (ED3c)", "RuLAS.DISCO.rds not found")
}

# ------------------------------------------------------------------------------
# ED3d - CLHLS (status -9 = lost to follow-up, excluded per codebook)
# ------------------------------------------------------------------------------
clhls_file <- file.path(root, "Result4/CLHLS.DISCO.rds")
if (file.exists(clhls_file)) {
  pd.vd <- readRDS(clhls_file)

  # Scatter panels use the full object (all 2546 records), as in the legacy script.
  save_pdf(disco_age_plot(pd.vd, "CLHLS",
                          x_breaks = seq(50, 110, 10), y_breaks = c(-3, 0, 3, 6)),
           file.path(out_dir, "CLHLS.DISCO.age.pdf"), 3.5, 3)

  pd.vd$hd1 <- as.numeric(log(pd.vd$hd))
  pd.vd1 <- pd.vd[pd.vd$mc < 7, ]
  pd.vd1 <- pd.vd1[pd.vd1$hd1 < 0, ]
  save_pdf(disco_hd_plot(pd.vd1, "CLHLS", xname = "hd1", y_breaks = c(-3, 0, 3, 6)),
           file.path(out_dir, "CLHLS.DISCO.hd.pdf"), 3.5, 3)

  # KM panel excludes status -9 (lost to follow-up) and NA, per the codebook.
  pd.vd <- pd.vd[pd.vd$status %in% c(0, 1), ]
  kmd <- build_km_groups(pd.vd)
  save_pdf(km_age_plot(kmd, "CLHLS", xlim = c(45, 115), legend_pos = c(0.27, 0.25)),
           file.path(out_dir, "CLHLS.KM.pdf"), 5, 4)
} else {
  record_missing("CLHLS (ED3d)", "CLHLS.DISCO.rds not found")
}

# ------------------------------------------------------------------------------
# ED6a - HRS C-index comparison (Original + Accel bar charts)
# ------------------------------------------------------------------------------
cindex_se <- function(formula, data) {
  fit <- coxph(as.formula(formula), data = data)
  s   <- summary(fit)
  c(cindex = s$concordance[1], se = s$concordance[2])
}

hrs_cmp_file <- file.path(root, "Manuscript/PNAS/hrs_cindex_comparison.csv")
hrs_data_file <- file.path(root, "HRS/HRS/HRS_data/HRS.data.rds")
if (file.exists(hrs_cmp_file) && file.exists(hrs_data_file)) {
  dat <- read.csv(hrs_cmp_file, stringsAsFactors = FALSE)
  dat$label <- dat$metric
  dat$label[dat$label == "DM (Mahalanobis)"] <- "DM+Age"
  dat$label[dat$label == "DISCO"]            <- "DISCO+Age"

  keep_labels <- c("DISCO+Age", "DM+Age", "DunedinPACE",
                   "GrimAge", "PhenoAge", "Hannum", "Horvath",
                   "VidalBralo", "Lin", "Zhang", "SkinBlood", "Weidner")
  dat <- dat[dat$label %in% keep_labels, ]

  hrs_obj <- readRDS(hrs_data_file)
  tr      <- hrs_obj$Train
  ok      <- !is.na(tr$status) & !is.na(tr$time) & !is.na(tr$age)
  age_ci  <- cindex_se("Surv(time, status) ~ age", tr[ok, ])

  age_row <- data.frame(label = "Age", adj_c = age_ci[1], adj_se = age_ci[2],
                        accel_c = NA, accel_se = NA)

  d <- rbind(
    data.frame(label = dat$label, adj_c = dat$adj_c, adj_se = dat$adj_se,
               accel_c = dat$accel_c, accel_se = dat$accel_se),
    age_row
  )

  highlight <- c("DISCO+Age", "DM+Age")

  d1 <- d[order(d$adj_c), ]
  d1$label <- factor(d1$label, levels = d1$label)
  fill1 <- ifelse(d1$label %in% highlight, "#E64B35",
           ifelse(d1$label == "Age", "#4DBBD5", "grey70"))

  p1 <- ggplot(d1, aes(adj_c, label)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.9, fill = fill1) +
    geom_errorbar(aes(xmin = adj_c - adj_se, xmax = adj_c + adj_se), width = 0.4) +
    coord_cartesian(xlim = c(0.70, 0.80)) +
    theme_classic() + xlab("C-index") + ylab("Metrics") +
    labs(title = "Original") +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(5, 5, 10, 5))
  save_pdf(p1, file.path(out_dir, "HRS.cmp.original.pdf"), 3.5, 5)

  d2 <- d[!is.na(d$accel_c), ]
  d2$label <- gsub("\\+Age$", "", d2$label)
  d2 <- d2[order(d2$accel_c), ]
  d2$label <- factor(d2$label, levels = d2$label)
  fill2 <- ifelse(d2$label %in% c("DISCO", "DM"), "#E64B35", "grey70")

  p2 <- ggplot(d2, aes(accel_c, label)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.9, fill = fill2) +
    geom_errorbar(aes(xmin = accel_c - accel_se, xmax = accel_c + accel_se), width = 0.4) +
    coord_cartesian(xlim = c(0.50, 0.65)) +
    theme_classic() + xlab("C-index") + ylab("Metrics") +
    labs(title = "Accel") +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          plot.margin = margin(5, 10, 10, 5))
  save_pdf(p2, file.path(out_dir, "HRS.cmp.accel.pdf"), 3.5, 5)
} else {
  if (!file.exists(hrs_cmp_file)) record_missing("HRS cmp (ED6a)", "hrs_cindex_comparison.csv not found")
  if (!file.exists(hrs_data_file)) record_missing("HRS cmp Age row (ED6a)", "HRS.data.rds not found")
}

# ------------------------------------------------------------------------------
# ED6b - HRS KM by DISCO / GrimAge / HorvathAge
# ------------------------------------------------------------------------------
hrs_res_file <- file.path(root, "Result1/DNAm/HRS.DNAm.res.rds")
if (file.exists(hrs_res_file)) {
  DNAm  <- readRDS(hrs_res_file)
  DNAm1 <- DNAm[!is.na(DNAm$DNAmcliage), ]

  pd.vd1 <- DNAm1
  pd.vd1$group1 <- cut(pd.vd1$DNAMGRIMAGE,
                       breaks = quantile(pd.vd1$DNAMGRIMAGE, probs = seq(0, 1, 0.25)),
                       include.lowest = TRUE, right = TRUE)
  pd.vd1$group1 <- as.numeric(pd.vd1$group1)
  pd.vd1 <- pd.vd1[pd.vd1$group1 %in% c(1, 4), ]
  pd.vd1$group1[pd.vd1$group1 > 1] <- 2
  pd.vd1$group1 <- pd.vd1$group1 + 2

  pd.vd2 <- DNAm1
  pd.vd2$group1 <- cut(pd.vd2$HORVATH_DNAMAGE,
                       breaks = quantile(pd.vd2$HORVATH_DNAMAGE, probs = seq(0, 1, 0.25)),
                       include.lowest = TRUE, right = TRUE)
  pd.vd2$group1 <- as.numeric(pd.vd2$group1)
  pd.vd2 <- pd.vd2[pd.vd2$group1 %in% c(1, 4), ]
  pd.vd2$group1[pd.vd2$group1 > 1] <- 2
  pd.vd2$group1 <- pd.vd2$group1 + 4

  pd.vd3 <- DNAm1
  pd.vd3$group1 <- cut(pd.vd3$DNAmcliage,
                       breaks = quantile(pd.vd3$DNAmcliage, probs = seq(0, 1, 0.25)),
                       include.lowest = TRUE, right = TRUE)
  pd.vd3$group1 <- as.numeric(pd.vd3$group1)
  pd.vd3 <- pd.vd3[pd.vd3$group1 %in% c(1, 4), ]
  pd.vd3$group1[pd.vd3$group1 > 1] <- 2

  pd.vd1 <- rbind(pd.vd1, pd.vd2, pd.vd3)
  pd.vd1$group1 <- as.factor(pd.vd1$group1)
  levels(pd.vd1$group1) <- c(paste("DISCO", c("Low", "High"), "Risk"),
                             paste("GrimAge", c("Low", "High"), "Risk"),
                             paste("HorvathAge", c("Low", "High"), "Risk"))

  pd.vd1$age1 <- pd.vd1$age + pd.vd1$time

  fit1 <- survfit(Surv(age, age1, status) ~ group1, data = pd.vd1)

  p <- ggsurvplot(fit1,
                  data = pd.vd1,
                  conf.int = TRUE,
                  conf.int.alpha = 0.1,
                  censor = FALSE,
                  pval = FALSE,
                  xlab = "Age (Years)",
                  ylab = "Survival probability",
                  legend.title = "",
                  legend.labs = c(paste("DISCO", c("Low", "High"), "Risk"),
                                  paste("GrimAge", c("Low", "High"), "Risk"),
                                  paste("HorvathAge", c("Low", "High"), "Risk")),
                  title = "HRS",
                  xlim = c(50, 100),
                  pval.coord = c(0.1, 0.85),
                  legend = "right",
                  palette = "npg")

  p1 <- p$plot + theme(legend.background = element_blank(),
                       legend.text = element_text(size = 10),
                       axis.title = element_text(size = 14),
                       plot.title = element_text(size = 14))
  save_pdf(p1, file.path(out_dir, "HRS.KM.pdf"), 7, 4)
} else {
  record_missing("HRS KM (ED6b)", "HRS.DNAm.res.rds not found")
}

# ------------------------------------------------------------------------------
# Missing-panel report
# ------------------------------------------------------------------------------
if (length(missing) > 0) {
  writeLines(c("# Missing panels", "", missing),
             file.path(out_dir, "MISSING.md"))
  cat("\nMISSING PANELS:\n", paste(missing, collapse = "\n"), "\n", sep = "")
} else {
  cat("\nALL PANELS WRITTEN SUCCESSFULLY\n")
}

cat("\nOutputs in:", out_dir, "\n")
