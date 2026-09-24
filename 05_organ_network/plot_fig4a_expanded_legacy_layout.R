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

# Fig 4a: expanded-covariate HR numbers plotted in the ORIGINAL submission layout.
# Layout copied verbatim from Result3/Organ/Plot_Fig4C.R Organ.outcomes.pdf block
# (facet_wrap two horizontal panels, strip on top, log HR, Reds distiller, theme_few).
# Only the numbers are replaced with expanded-covariate Model2 HR.

suppressPackageStartupMessages({
  library(ggplot2); library(ggthemes); library(reshape2)
})

root <- disco_project_root()
src <- disco_results_path("05_organ_network/outputs/organ_expanded_covariates")

# --- legacy factor levels (verbatim from Fig4a_legacy_plot_data.rds) ---
legacy_var1 <- c("brain","esophagus","intestine","lymphoid tissue","choroid plexus",
                 "liver","adipose tissue","heart muscle","lung","bone marrow",
                 "salivary gland","kidney","pancreas","skin","stomach",
                 "skeletal muscle","tongue","retina")
legacy_var2 <- c("Dementia","COPD","Stroke","MI","Diabetes","Hypertension","Obesity",
                 "Others","Alzheimer's disease","Cerebrovascular disease","Cancer",
                 "Respiratory disease","Heart disease","All-cause")

stars <- function(p) ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))

# --- expanded-covariate data ---
m <- read.csv(file.path(src, "organ_cause_specific_expanded.csv"), stringsAsFactors = FALSE)
m$tt <- "Mortality"; m$Var2 <- m$outcome; m$Var1 <- m$organ
m$value <- m$Model2_HR; m$p <- stars(m$Model2_P)
m <- m[, c("Var1","Var2","value","p","tt")]

d <- read.csv(file.path(src, "organ_chronic_disease_expanded.csv"), stringsAsFactors = FALSE)
d$tt <- "Chronic disease incidence"; d$Var2 <- d$disease; d$Var1 <- d$organ
d$value <- d$Model2_HR; d$p <- stars(d$Model2_P)
d <- d[, c("Var1","Var2","value","p","tt")]

rr <- rbind(m, d)
rr$Var1 <- factor(rr$Var1, levels = legacy_var1)
rr$Var2 <- factor(rr$Var2, levels = legacy_var2)
rr$tt <- factor(rr$tt, levels = c("Mortality","Chronic disease incidence"))

# --- original submission layout, transposed: x=organ, y=outcome, nrow=2 ---
pdf(file.path(src, "Fig4a_expanded_legacy_layout.pdf"), height = 6, width = 12)
print(
  ggplot(rr, aes(Var1, Var2)) +
    geom_tile(aes(fill = value), color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.5, hjust = 0.7) +
    geom_text(aes(label = p), hjust = -0.8, vjust = 0.7) +
    facet_wrap(tt ~ ., nrow = 2, scales = "free_y") +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_few() + labs(fill = "HR") +
    xlab("Organ and tissues") + ylab("Outcomes") +
    theme(panel.spacing.y = unit(0.6, "lines"),
          axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12, color = "black"),
          strip.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 14))
)
dev.off()
cat("saved Fig4a_expanded_legacy_layout.pdf\n")
