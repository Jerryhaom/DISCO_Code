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
disco_init_module("01_method_definition")

# R27 combined into a single panel (no facet): model = colour, weight = linetype.
suppressMessages({library(ggplot2); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
res <- readRDS(file.path(out, "R27_versions_full.rds"))
res$model <- factor(res$model, levels = c("crude", "age+sex", "accel"))
refsizes <- sort(unique(res$refsize))

p <- ggplot(res, aes(refsize, ci, colour = model, linetype = weight)) +
  geom_point(size = 1, position = position_dodge(0.12)) +
  geom_smooth(aes(group = interaction(model, weight)), colour = "black", se = FALSE, linewidth = 0.5) +
  scale_x_log10(breaks = refsizes) +
  scale_colour_manual(values = c("crude" = "#4DBBD5FF", "age+sex" = "#E64B35FF", "accel" = "#00A087FF")) +
  scale_linetype_manual(values = c("age-weighted" = "solid", "uniform" = "dashed")) +
  theme_bw() +
  labs(x = "Reference size", y = "C-index for mortality", colour = NULL, linetype = NULL,
       title = "UKB") +
  theme(legend.position = "right",
        axis.title = element_text(size = 12))

ggsave(file.path(out, "R27_models_vs_refsize_combined.pdf"), p, width = 6.5, height = 4.5)
cat("wrote R27_models_vs_refsize_combined.pdf\n")
