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

# Reproduce publication figures from LOCKED legacy saved objects.
#   Fig 3c  : 3-modality (ClinicalProt / Prot / Clinical) UKB chronic-disease
#             hazard-ratio forest plot  (forestploter)  -> Forest.BioAge.hd.pdf
#   ED7     : 8 panels cumulative-incidence KM (DISCO Low/High + DM Low/High)
#             for UKB chronic diseases (ClinicalProt).
#
# NOTE: Fig 2a-b (proteomics KM) already produced in outputs/03_prot_analysis/
#       and is intentionally skipped here.
#
# All source objects are READ-ONLY (never modified).

options(stringsAsFactors = FALSE)
suppressMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggsci)
  library(forestploter)
  library(grid)
  library(dplyr)
  library(ggpubr)
})

root <- disco_project_root()
wkd  <- disco_results_path("03_multimodal_UKB")
out  <- file.path(wkd, "outputs", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

res3  <- file.path(root, "Result3", "Result3")
supple2 <- file.path(root, "Supple", "Result2")

###############################################################################
# Fig 3c — 3-modality forest plot (submitted format, Forest.BioAge.hd.pdf)
###############################################################################

forest11 <- readRDS(file.path(res3, "forest.rds"))   # 7-col display table
incres   <- readRDS(file.path(res3, "inc.res.rds"))  # + HR / low / up

# forest_theme copied verbatim from Result3/Result3/forest.R
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

p3c <- forest(forest11,
              est   = as.numeric(incres$HR),
              lower = as.numeric(incres$low),
              upper = as.numeric(incres$up),
              sizes = 0.4,
              ci_column = 7,
              ref_line = 1,
              x_trans = "log",
              xlim = c(0.7, 2.3),
              ticks_digits = 1,
              ticks_at = c(0.9, 1.3, 1.7, 2.3),
              theme = tm)

nrow3 <- nrow(forest11)  # 24 rows (6 diseases x 4 rows each)
p3c <- add_border(p3c, row = 1:nrow3, col = 1:6, where = "top",
                  gp = gpar(lwd = 0.5))
p3c <- add_border(p3c, part = "header", col = 1:6,
                  where = "bottom", gp = gpar(lwd = 1.5))
p3c <- edit_plot(p3c, row = seq(1, nrow3, 4), col = 1,
                 gp = gpar(fontface = "bold"))
p3c <- edit_plot(p3c, row = c(1:nrow3)[-seq(1, nrow3, 4)], col = 1,
                 hjust = 0, x = unit(0.2, "npc"))

pdf(file.path(out, "Fig3c_forest_3modality.pdf"), height = 8, width = 10)
print(p3c)
dev.off()

###############################################################################
# ED7 — 8-panel cumulative-incidence KM (ClinicalProt, DISCO vs DM)
###############################################################################

pd  <- readRDS(file.path(root, "Result2", "UKB.ClinicalProt.DISCO.rds"))
dis <- readRDS(file.path(root, "UKB_Data", "dat_disease.rds"))

dis2 <- dis[match(pd$ID, dis$eid), ]
pd <- cbind(pd, dis2[, setdiff(colnames(dis2), "eid")])

pd$mcg <- as.numeric(cut(pd$mc, breaks = quantile(pd$mc, na.rm = TRUE),
                         include.lowest = TRUE, right = TRUE))
pd$hdg <- as.numeric(cut(pd$hd, breaks = quantile(pd$hd, na.rm = TRUE),
                         include.lowest = TRUE, right = TRUE))

nn     <- c("obesity", "hypertension", "diabetes", "cancer",
            "mi", "stroke", "copd", "dementia")
titles <- c("Obesity", "Hypertension", "Diabetes", "Cancer",
            "Myocardial infarction", "Stroke", "COPD", "Dementia")

npg_cols <- c("DISCO Low Risk"  = "#E64B35",
              "DISCO High Risk" = "#4DBBD5",
              "DM Low Risk"     = "#00A087",
              "DM High Risk"    = "#3C5488")

process_disease <- function(status_var, time_var, title) {
  tm <- pd[[time_var]]; st <- pd[[status_var]]
  keep <- !is.na(tm) & tm >= 0 & !is.na(st) & st %in% c(0, 1) & !is.na(pd$mc)
  p <- pd[keep, ]
  p$time   <- p[[time_var]]
  p$status <- p[[status_var]]
  p$time1  <- p$time + p$age
  # delayed entry requires stop > start; nudge ties at baseline
  p$time1  <- pmax(p$time1, p$age + 1e-4)

  # unadjusted C-index (matches Table 17 / ED7 annotation)
  c_dis <- summary(coxph(Surv(time, status) ~ mc,     data = p))$concordance[1]
  c_dm  <- summary(coxph(Surv(time, status) ~ log(hd), data = p))$concordance[1]

  d1 <- p[p$mcg %in% c(1, 4), ]
  d1$grp <- ifelse(d1$mcg == 1, "DISCO Low Risk", "DISCO High Risk")
  d2 <- p[p$hdg %in% c(1, 4), ]
  d2$grp <- ifelse(d2$hdg == 1, "DM Low Risk", "DM High Risk")

  fit1 <- survfit(Surv(age, time1, status) ~ grp, data = d1)
  fit2 <- survfit(Surv(age, time1, status) ~ grp, data = d2)

  get_surv <- function(fit) {
    s <- summary(fit)
    data.frame(time = s$time, surv = s$surv,
               lower = s$lower, upper = s$upper,
               grp = sub("^grp=", "", s$strata), stringsAsFactors = FALSE)
  }

  dd <- rbind(get_surv(fit1), get_surv(fit2))
  dd$grp <- factor(dd$grp, levels = names(npg_cols))

  # dynamic y-axis limit so no curve / confidence band is clipped
  yl  <- max(1 - dd$lower, na.rm = TRUE) * 1.12
  ann <- sprintf("DISCO C-index: %.3f; DM C-index: %.3f", c_dis, c_dm)

  gg <- ggplot(dd, aes(x = time, y = 1 - surv, color = grp, fill = grp)) +
    geom_ribbon(aes(ymax = 1 - lower, ymin = 1 - upper), alpha = 0.3,
                linewidth = 0, color = "white") +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = npg_cols) +
    scale_fill_manual(values = npg_cols) +
    theme_classic() +
    labs(x = "Age (Years)", y = "Cumulative incidence", title = title) +
    ylim(0, yl) +
    annotate("text", x = 42, y = yl * 0.90, label = ann,
             hjust = 0, size = 3.0, color = "black") +
    theme(axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 12, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.position = c(0.30, 0.72),
          legend.background = element_rect(fill = "transparent"))
  gg
}

pl <- list()
for (i in seq_along(nn)) {
  pl[[i]] <- process_disease(paste0("status_", nn[i]), paste0("time_", nn[i]),
                             titles[i])
}

combined <- ggarrange(
  plotlist = pl, nrow = 4, ncol = 2,
  common.legend = TRUE, legend = "right",
  align = "hv"
)

pdf(file.path(out, "ED7_disease_KM.pdf"), height = 12, width = 10)
print(combined)
dev.off()

cat("Fig3c ->", file.path(out, "Fig3c_forest_3modality.pdf"), "\n")
cat("ED7   ->", file.path(out, "ED7_disease_KM.pdf"), "\n")
writeLines(capture.output(sessionInfo()),
           file.path(out, "sessionInfo_figures.txt"))
cat("DONE\n")
