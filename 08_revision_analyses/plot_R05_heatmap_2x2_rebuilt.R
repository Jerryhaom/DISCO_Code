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
disco_init_module("08_revision_analyses")

# R05 heatmap (3/3): 2x2 group HR by tertile, events<5 greyed, diseases on top.

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

wkd <- disco_results_path("08_revision_analyses")
src <- file.path(wkd, "outputs/R05_genetic_interaction_rebuilt")
out <- file.path(src, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(src, "R05_2x2_group_hr.csv"), stringsAsFactors = FALSE)

lbl <- c(af = "AF", cad = "CAD", hypertension = "Hypertension",
         t2d = "Type 2 diabetes", osteoporosis = "Osteoporosis",
         vte = "VTE", prostate_ca = "Prostate cancer", amd = "AMD")

d$label <- factor(lbl[d$disease], levels = unname(lbl))
d$disco <- factor(d$disco, levels = c("Low", "High"))
d$prs   <- factor(d$prs,   levels = c("Low", "High"))
d$version <- ifelse(d$version == "clinical", "Clinical DISCO", "Proteomics DISCO")

d$fillv <- ifelse(d$events < 5, NA, log(d$hr))
d$star <- ifelse(is.na(d$p), "",
                 ifelse(d$p < 0.001, "***",
                 ifelse(d$p < 0.01, "**",
                 ifelse(d$p < 0.05, "*", ""))))
d$lab <- ifelse(d$events < 5, "n<5",
                ifelse(d$disco == "Low" & d$prs == "Low", "ref",
                       paste0(sprintf("%.1f", d$hr), d$star)))

gg <- ggplot(d, aes(x = disco, y = prs, fill = fillv)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 2.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey85",
                       limits = log(c(1/8, 8)), oob = scales::squish,
                       breaks = log(c(0.5, 1, 2, 4, 8)),
                       labels = c("0.5", "1", "2", "4", "8"),
                       name = "HR") +
  facet_grid(version ~ label) +
  labs(x = "DISCO (median split)", y = "PRS (median split)",
       title = "Incident-disease HR by DISCO and genetic-risk strata (reference = Low/Low; * P<0.05, ** P<0.01, *** P<0.001)") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        strip.text.x = element_text(size = 7, face = "bold"),
        strip.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(out, "R05_2x2_heatmap.pdf"), gg, width = 11, height = 4.5)
ggsave(file.path(out, "R05_2x2_heatmap.png"), gg, width = 11, height = 4.5, dpi = 300)
cat("saved 2x2 heatmap (rebuilt)\n")
