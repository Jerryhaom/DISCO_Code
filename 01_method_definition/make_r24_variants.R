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

# R2.4 nonlinear detection — clearer variants.
# Option A: deviation-from-healthy grouped bars (0 = healthy baseline).
# Option B: deviation heatmap (scenario x metric).

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(ggthemes)})

root <- disco_project_root()
csv  <- disco_results_path("08_sensitivity_analyses/outputs/R03_correlation_robustness/R24_monotone_nonlinear.csv")
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

r24 <- read.csv(csv)
r24$scenario <- factor(r24$scenario,
  levels = c("healthy","linear","quadratic","cubic","sigmoid","funnel"),
  labels = c("Healthy (baseline)","Linear (weakened)","Quadratic (U-shape)",
             "Cubic (monotone)","Sigmoid (monotone)","Funnel (heterosced.)"))

# deviation from healthy for each metric
h <- r24[r24$scenario == "Healthy (baseline)", ]
dev <- data.frame(
  scenario = r24$scenario,
  `DISCO (Pearson)` = r24$disco_pearson - h$disco_pearson,
  `DISCO (rank)`    = r24$disco_rank    - h$disco_rank,
  `DM`              = r24$dm            - h$dm,
  check.names = FALSE
)
devl <- dev %>% pivot_longer(-scenario, names_to = "metric", values_to = "dev")
devl$metric <- factor(devl$metric, levels = c("DISCO (Pearson)","DISCO (rank)","DM"))

# ---------- Option A: grouped bar, deviation from healthy ----------
pA <- ggplot(devl, aes(scenario, dev, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_hline(yintercept = 0, linewidth = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("DISCO (Pearson)" = "#E64B35FF",
                               "DISCO (rank)" = "#4DBBD5FF",
                               "DM" = "#00A087FF")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, colour = "black"),
        legend.position = "top", legend.title = element_blank()) +
  labs(x = NULL, y = "Deviation from healthy baseline",
       title = "R2.4  nonlinear-dependence detection",
       subtitle = "Positive = detected as dysregulated; negative = looks 'healthier' than baseline (missed)")
ggsave(file.path(out, "R24_optA_deviation_bars.pdf"), pA, width = 8, height = 5)

# ---------- Option B: heatmap ----------
pB <- ggplot(devl, aes(metric, scenario, fill = dev)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", dev)), colour = "black", size = 3.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Deviation\nfrom healthy") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "R2.4  nonlinear-dependence detection (heatmap)",
       subtitle = "Red = detected; blue = missed (looks healthier than baseline)")
ggsave(file.path(out, "R24_optB_heatmap.pdf"), pB, width = 7, height = 5)

cat("wrote optA + optB\n")
