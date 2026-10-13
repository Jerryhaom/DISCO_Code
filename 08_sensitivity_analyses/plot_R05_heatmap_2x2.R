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
disco_init_module("08_sensitivity_analyses")

# R05 visualization (3/3): 2x2 group-effect heatmap (DISCO median x PRS median).

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

wkd <- disco_results_path("08_sensitivity_analyses")
out <- file.path(wkd, "outputs/R05_genetic_interaction/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(wkd, "outputs/R05_genetic_interaction/R05_2x2_group_hr.csv"),
              stringsAsFactors = FALSE)

lbl <- c(
  af = "Atrial fibrillation", cad = "Coronary artery disease",
  hypertension = "Hypertension", t2d = "Type 2 diabetes",
  osteoporosis = "Osteoporosis", vte = "Venous thromboembolism",
  prostate_ca = "Prostate cancer", amd = "Age-related macular degeneration")

d$label <- factor(lbl[d$disease], levels = unname(lbl))
d$disco <- factor(d$disco, levels = c("Low", "High"))
d$prs   <- factor(d$prs, levels = c("Low", "High"))
d$version <- ifelse(d$version == "clinical", "Clinical DISCO", "Proteomics DISCO")

gg <- ggplot(d, aes(x = disco, y = prs, fill = log2(hr))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(is.na(hr), "ref", sprintf("%.1f", hr))),
            colour = "black", size = 2.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "HR (log2)") +
  facet_grid(label ~ version) +
  labs(x = "DISCO", y = "PRS",
       title = "Incident-disease HR by DISCO and genetic-risk strata (reference = Low/Low)") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(size = 7, face = "bold"),
        strip.text.x = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(out, "R05_2x2_heatmap.pdf"), gg, width = 7, height = 9)
ggsave(file.path(out, "R05_2x2_heatmap.png"), gg, width = 7, height = 9, dpi = 300)
cat("saved 2x2 heatmap\n")
