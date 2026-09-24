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
disco_init_module("06_frailty_function")


# =============================================================================
# Reproduce Extended Data Fig 8 (CHARLS frailty) from LOCKED legacy objects.
#
#   (a) forest plot + table: frailty incidence OR (95% CI, P) by DISCO quartile
#       for Wave1-3 and Wave3-4 cohorts
#   (b) Kaplan-Meier: "CHARLS, DISCO + Frailty" (Frail/Prefrail/Robust x High/Low)
#   (c) violin: DISCO by Robust/Prefrail/Frail (Wave 1 to 3, Wave 3 to 4)
#   (d) ROC: Age+DISCO / Age+FI / Age / DISCO+FI / DISCO / FI (Wave 1-3, Wave 3-4)
#
# Source plotting logic copied from:
#   Result5/Asso_DISCO_FI.R   (forest + KM + violin)
#   Result5/Plot_ROC.R        (ROC; curves taken precomputed from Result5/ROC.rds)
#
# Legacy files are read-only; nothing under Result5/ is modified.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(ggthemes)
  library(ggsci)
  library(forestploter)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(pROC)
})

# ---- paths -----------------------------------------------------------------
root    <- disco_project_root()
r5      <- file.path(root, "Result5")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

missing <- character()          # records of missing objects/columns

note_missing <- function(msg) {
  missing <<- c(missing, msg)
  cat("[MISSING] ", msg, "\n", sep = "")
}

# ---- load locked objects ---------------------------------------------------
load_required <- function(path, label) {
  if (file.exists(path)) return(readRDS(path))
  note_missing(paste(label, "->", path))
  NULL
}

pd2011 <- load_required(file.path(r5, "CHARLS.2011.disco.rds"), "CHARLS.2011.disco.rds")
pd2015 <- load_required(file.path(r5, "CHARLS.2015.disco.rds"), "CHARLS.2015.disco.rds")
roc_obj <- load_required(file.path(r5, "ROC.rds"), "ROC.rds")

fi2011 <- load_required(file.path(r5, "Frailty", "CHARLS.FI.2011.rds"), "CHARLS.FI.2011.rds")
fi2013 <- load_required(file.path(r5, "Frailty", "CHARLS.FI.2013.rds"), "CHARLS.FI.2013.rds")
fi2015 <- load_required(file.path(r5, "Frailty", "CHARLS.FI.2015.rds"), "CHARLS.FI.2015.rds")
fi2018 <- load_required(file.path(r5, "Frailty", "CHARLS.FI.2018.rds"), "CHARLS.FI.2018.rds")

# ---- build frailty-incidence cohorts (exact data flow from Asso_DISCO_FI.R) ----
if (!is.null(pd2011) && !is.null(pd2015) &&
    !is.null(fi2011) && !is.null(fi2013) && !is.null(fi2015) && !is.null(fi2018)) {

  # Wave1-3 incidence cohort
  fi1  <- fi2011[match(pd2011$id, fi2011$ID), ]
  fi21 <- fi2013[match(fi1$ID1, fi2013$ID), ]
  fi31 <- fi2015[match(fi1$ID1, fi2015$ID), ]
  pd1  <- data.frame(y1 = as.numeric(fi1$scoreg),
                     y2 = as.numeric(fi21$scoreg),
                     y3 = as.numeric(fi31$scoreg))
  pd1$y <- rep(0, nrow(fi1))
  pd1$y[pd1$y3 == 3] <- 1
  pd1$y[pd1$y2 == 3] <- 1
  pd1$y[is.na(pd1$y3) & is.na(pd1$y2)] <- NA
  pd1$y[is.na(pd1$y1)] <- NA
  pd1$y[pd1$y1 == 3] <- NA
  pd1 <- data.frame(pd1, pd2011)
  pd1$edu <- as.numeric(pd1$edu)
  pd1$gg  <- as.numeric(cut(pd1$mc, breaks = quantile(pd1$mc),
                            right = TRUE, include.lowest = TRUE))

  # Wave3-4 incidence cohort (exclude wave1 re-interviewed individuals)
  pd2015b <- pd2015[!pd2015$id %in% pd2011$id1, ]
  fi3  <- fi2015[match(pd2015b$id, fi2015$ID), ]
  fi41 <- fi2018[match(fi3$ID, fi2018$ID), ]
  pd2  <- data.frame(y1 = as.numeric(fi3$scoreg),
                     y2 = as.numeric(fi41$scoreg))
  pd2$y <- rep(0, nrow(fi3))
  pd2$y[pd2$y2 == 3] <- 1
  pd2$y[is.na(pd2$y2)] <- NA
  pd2$y[is.na(pd2$y1)] <- NA
  pd2$y[pd2$y1 == 3] <- NA
  pd2 <- data.frame(pd2, pd2015b)
  pd2$edu <- as.numeric(pd2$edu)
  pd2$gg  <- as.numeric(cut(pd2$mc, breaks = quantile(pd2$mc),
                            right = TRUE, include.lowest = TRUE))
}

# ---- panel (a) forest ------------------------------------------------------
compute_forest_block <- function(pd) {
  o1 <- matrix(0, ncol = 4, nrow = 5)
  no <- matrix(0, ncol = 1, nrow = 5)
  for (i in 0:4) {
    cc <- c()
    j  <- i
    if (i == 0) {                                   # "All" (continuous DISCO)
      tb <- as.numeric(table(pd$y))
      no[1, ] <- paste(tb[2], sum(tb), sep = "/")
      l  <- glm(data = pd, y ~ mc + age + sex + smoke + drink + edu,
                family = "binomial")
      s  <- summary(l)
      or <- round(c(exp(s$coefficients[2, 1]),
                    exp(confint(l)[2, 1]), exp(confint(l)[2, 2])), 2)
      cc <- c(cc, or, s$coefficients[2, 4])
      o1[i + 1, ] <- cc
      next
    }
    pd11      <- pd[pd$gg %in% c(1, j), ]
    pd11$gg1  <- pd11$gg
    if (i == 1) {                                   # Q1 (reference)
      o1[i + 1, ] <- rep("", ncol(o1))
      next
    }
    if (i > 1) pd11$gg1 <- as.numeric(as.factor(pd11$gg))
    tb <- table(data.frame(pd11$gg, pd11$y))
    no[2, ]     <- paste(tb[1, 2], sum(tb[1, ]), sep = "/")
    no[j + 1, ] <- paste(tb[2, 2], sum(tb[2, ]), sep = "/")
    l  <- glm(data = pd11, y ~ gg1 + age + sex + smoke + drink + edu,
              family = "binomial")
    s  <- summary(l)
    or <- round(c(exp(s$coefficients[2, 1]),
                  exp(confint(l)[2, 1]), exp(confint(l)[2, 2])), 2)
    cc <- c(cc, or, s$coefficients[2, 4])
    o1[i + 1, ] <- cc
  }
  cbind(no, o1)
}

forest_done <- FALSE
if (exists("pd1") && exists("pd2")) {
  res1 <- compute_forest_block(pd1)
  res2 <- compute_forest_block(pd2)

  # overall trend P values (continuous quartile) used for the group labels
  p1 <- summary(glm(y ~ gg + age + sex + smoke + drink + edu,
                    data = pd1, family = "binomial"))$coefficients[2, 4]
  p2 <- summary(glm(y ~ gg + age + sex + smoke + drink + edu,
                    data = pd2, family = "binomial"))$coefficients[2, 4]
  fmt_p <- function(p) ifelse(is.na(p) | p < 0.001, "<0.001",
                              sprintf("%.3f", round(p, 3)))

  tm <- forest_theme(
    base_size = 12,
    core = list(bg_params = list(fill = c("white"))),
    ci_pch = 15,
    ci_fill = "#374E55FF",
    ci_col = "black",
    ci_lty = 1,
    ci_lwd = 1.5,
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "grey20",
    footnote_cex = 0.6,
    footnote_fontface = "italic",
    footnote_col = "red4"
  )

  res <- rbind(c(rep("", 5)), res1, c(rep("", 5)), res2)
  x  <- res[, 1]
  x[x == ""] <- "./."
  x1 <- unlist(strsplit(x, "/"))
  x1 <- as.numeric(x1)
  x2 <- sprintf("%.2f", x1[seq(1, length(x1), 2)] * 100 / x1[seq(2, length(x1), 2)])
  x2 <- paste(x2, "%", sep = "")
  x2 <- paste(x1[seq(1, length(x1), 2)], " (", x2, ")", sep = "")
  x2[x == "./."] <- ""
  res[, 1] <- x2

  res <- data.frame(
    gg = c(c("Wave1-3", "All", "Q1", paste("Q", seq(2, 4), sep = "")),
           c("Wave3-4", "All", "Q1", paste("Q", seq(2, 4), sep = ""))),
    res)
  colnames(res) <- c("Category", "Events (%)", "OR", "low", "up", "p")
  res$p <- round(as.numeric(res$p), 3)
  res$p[res$p < 0.001] <- "<0.001"
  res$lb <- paste(sprintf("%.2f", as.numeric(res$OR)), " (",
                  sprintf("%.2f", as.numeric(res$low)), "-",
                  sprintf("%.2f", as.numeric(res$up)), ")", sep = "")
  res$lb[c(1, 7)] <- ""
  res$lb[c(3, 9)] <- "1 (Ref)"

  res3 <- data.frame(res[, c(1, 2, 7, 6)],
                     P = rep(paste(rep(" ", 40), collapse = ""), nrow(res)))
  res3[res3 == "NA"] <- ""
  res3[is.na(res3)] <- ""
  res3 <- as.matrix(res3)
  colnames(res3)[2:4] <- c("Events (%)", "OR (95% CI)", "P value")
  colnames(res3)[5]   <- ""
  res4 <- res3[c(1:6, 7, 7:12, 7), ]
  res4[c(7, 14), 4] <- c(fmt_p(p1), fmt_p(p2))
  res  <- res[c(1:6, 7, 7:12, 7), ]

  p <- forest(res4,
              est = as.numeric(res$OR),
              lower = as.numeric(res$low),
              upper = as.numeric(res$up),
              sizes = 0.4,
              ci_column = c(5),
              ref_line = 1,
              x_trans = "log",
              xlim = c(0.3, 4),
              ticks_digits = 1,
              ticks_at = c(0.5, 1, 2, 3),
              title = "CHARLS, Frailty incidence",
              theme = tm)

  p <- add_border(p, row = 1:14, col = c(1:4), where = "top",
                  gp = gpar(lwd = 0.5))
  p <- add_border(p, part = "header", col = c(1:4),
                  where = "bottom", gp = gpar(lwd = 1.5))
  p <- edit_plot(p, row = c(1, 8), col = 1, gp = gpar(fontface = "bold"))
  p <- edit_plot(p, row = c(1:30)[-c(1, 8)], col = 1,
                 hjust = 0, x = unit(0.2, "npc"))

  pdf(file.path(out_dir, "frailty_forest.pdf"), height = 5, width = 7)
  print(p)
  dev.off()
  forest_done <- TRUE
  cat("[OK] frailty_forest.pdf (a)\n")
} else {
  note_missing("panel (a) forest: missing pd1/pd2 incidence cohorts")
}

# ---- panel (b) Kaplan-Meier --------------------------------------------------
km_done <- FALSE
if (!is.null(pd2011) && all(c("mc", "fl", "age", "time", "status") %in% names(pd2011))) {
  pdv <- pd2011
  pdv$g <- as.numeric(cut(pdv$mc, breaks = quantile(pdv$mc)))
  pdv   <- pdv[!is.na(pdv$fl), ]
  pdv1  <- pdv[pdv$g %in% c(1, 4), ]
  pdv1$flg <- paste(pdv1$fl, pdv1$g)
  pdv1$flg <- factor(pdv1$flg,
                     levels = names(table(pdv1$flg))[c(2, 1, 4, 3, 6, 5)])

  fit1 <- survfit(Surv(age, age + time, status) ~ flg, data = pdv1)

  p_km <- ggsurvplot(fit1,
                     data = pdv1,
                     conf.int = TRUE,
                     conf.int.alpha = 0.1,
                     censor = TRUE,
                     pval = FALSE,
                     xlab = "Age (years)",
                     ylab = "Survival probability",
                     legend.title = "",
                     legend.labs = c(paste("Frail, DISCO", c("High", "Low"), "Risk"),
                                     paste("Prefrail, DISCO", c("High", "Low"), "Risk"),
                                     paste("Robust, DISCO", c("High", "Low"), "Risk")),
                     title = "CHARLS, DISCO + Frailty",
                     xlim = c(45, max(fit1$time)),
                     break.x.by = 5,
                     pval.coord = c(16, 0.82),
                     legend = c(0.3, 0.26),
                     palette = "npg")

  p_km <- p_km$plot + theme(legend.background = element_blank(),
                            legend.text = element_text(size = 12),
                            legend.key.height = unit(0.8, "lines"),
                            axis.text = element_text(size = 12),
                            axis.title = element_text(size = 16),
                            plot.title = element_text(size = 16))

  ggsave(file.path(out_dir, "frailty_KM.pdf"), p_km, height = 4, width = 5)
  km_done <- TRUE
  cat("[OK] frailty_KM.pdf (b)\n")
} else {
  note_missing("panel (b) KM: pd2011 missing or lacks mc/fl/age/time/status")
}

# ---- panel (c) violin --------------------------------------------------------
violin_plot <- function(pd, title) {
  ggplot(pd[!is.na(pd$fl), ], aes(fl, mc)) +
    geom_violin(aes(fill = fl)) +
    geom_boxplot(width = 0.2, fill = "white") +
    labs(title = title) +
    stat_compare_means(comparisons = list(c("Robust", "Prefrail"),
                                          c("Prefrail", "Frail"),
                                          c("Robust", "Frail"))) +
    theme_few() + xlab("") + ylab("DISCO") +
    scale_fill_manual(values = rev(pal_npg("nrc")(3))) +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          legend.position = "None") +
    scale_y_continuous(expand = c(0.1, 0.1))
}

violin_done <- FALSE
if (!is.null(pd2011) && !is.null(pd2015) &&
    all(c("fl", "mc") %in% names(pd2011)) && all(c("fl", "mc") %in% names(pd2015))) {
  v1 <- violin_plot(pd2011, "CHARLS, Wave 1 to 3")
  v2 <- violin_plot(pd2015, "CHARLS, Wave 3 to 4")
  p_v <- ggarrange(v1, v2, ncol = 2, nrow = 1)
  ggsave(file.path(out_dir, "frailty_violin.pdf"), p_v, height = 3.5, width = 8)
  violin_done <- TRUE
  cat("[OK] frailty_violin.pdf (c)\n")
} else {
  note_missing("panel (c) violin: pd2011/pd2015 missing or lacks fl/mc")
}

# ---- panel (d) ROC -----------------------------------------------------------
roc_done <- FALSE
if (!is.null(roc_obj) && is.list(roc_obj) && length(roc_obj) == 2) {
  d1 <- roc_obj[[1]]
  d2 <- roc_obj[[2]]
  if (all(c("x", "y", "gg") %in% names(d1)) && all(c("x", "y", "gg") %in% names(d2))) {
    d1$tt <- "CHARLS, Wave 1 to 3"
    d2$tt <- "CHARLS, Wave 3 to 4"
    dd <- rbind(d1, d2)

    ls   <- unlist(strsplit(as.character(dd$gg), ":"))
    dd$gg1 <- ls[seq(1, length(ls), 2)]
    dd$gg1 <- factor(dd$gg1,
                     levels = c("Age+DISCO", "Age+Frailty index", "Age",
                                "DISCO+Frailty index", "DISCO", "Frailty index"))
    dd$yy <- 0.45 - as.numeric(dd$gg1) * 0.09

    p_roc <- ggplot(dd, aes(x = x, y = y, color = gg1)) +
      geom_path(size = 0.8) +
      scale_color_npg() +
      facet_wrap(tt ~ ., ncol = 2) +
      geom_text(aes(label = gg, y = yy + 0.1), x = 0.3, size = 3,
                check_overlap = TRUE, hjust = 0) +
      ylab("Sensitivity") + xlab("1-Specificity") +
      theme_few() +
      theme(legend.position = "None") +
      theme(strip.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12))

    ggsave(file.path(out_dir, "frailty_ROC.pdf"), p_roc, height = 4.5, width = 9)
    roc_done <- TRUE
    cat("[OK] frailty_ROC.pdf (d)\n")
  } else {
    note_missing("panel (d) ROC: ROC.rds data.frames lack x/y/gg columns")
  }
} else {
  note_missing("panel (d) ROC: ROC.rds not a list of 2 (precomputed curves unavailable)")
}

# ---- write MISSING.md ---------------------------------------------------------
writeLines(
  if (length(missing) == 0) {
    c("# Missing objects / columns",
      "",
      "None. All required legacy objects and columns were present.",
      "")
  } else {
    c("# Missing objects / columns", "", missing, "")
  },
  file.path(out_dir, "MISSING.md")
)

cat("\n===== summary =====\n",
    "forest:  ", if (forest_done)  "OK" else "FAIL", "\n",
    "KM:      ", if (km_done)      "OK" else "FAIL", "\n",
    "violin:  ", if (violin_done)  "OK" else "FAIL", "\n",
    "ROC:     ", if (roc_done)     "OK" else "FAIL", "\n",
    "missing: ", length(missing), " item(s)\n", sep = "")
