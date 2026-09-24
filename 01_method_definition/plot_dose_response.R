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

# Final R2.4 dose-response figure: nonlinearity strength -> correlation & DISCO detection.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(ggthemes); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
res  <- readRDS(file.path(out, "dose_response_results.rds"))

res$type_lab <- factor(res$type, levels = c("monotone","nonmonotone"),
                       labels = c("Monotone nonlinear (e.g. cubic)",
                                  "Non-monotone nonlinear (e.g. U-shape)"))

# correlation vs alpha
corl <- res %>% select(type_lab, alpha, `Pearson r` = pearson_r, `Spearman rho` = spearman_r) %>%
  pivot_longer(c(`Pearson r`, `Spearman rho`), names_to = "metric", values_to = "r")
corl$metric <- factor(corl$metric, levels = c("Pearson r", "Spearman rho"))

p_cor <- ggplot(corl, aes(alpha, r, colour = metric, shape = metric)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ type_lab) +
  scale_colour_manual(values = c("Pearson r" = "#E64B35FF", "Spearman rho" = "#4DBBD5FF")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(limits = c(-0.1, 1.05)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "Nonlinearity strength", y = "Correlation")

# DISCO detection (deviation from linear alpha=0)
d0 <- res %>% group_by(type) %>% summarise(base_p = disco_pearson[alpha == 0], base_r = disco_rank[alpha == 0])
res2 <- merge(res, d0, by = "type")
res2$dev_pearson <- res2$disco_pearson - res2$base_p
res2$dev_rank    <- res2$disco_rank    - res2$base_r
det <- res2 %>% select(type_lab, alpha, `DISCO (Pearson)` = dev_pearson, `DISCO (rank)` = dev_rank) %>%
  pivot_longer(c(`DISCO (Pearson)`, `DISCO (rank)`), names_to = "metric", values_to = "dev")
det$metric <- factor(det$metric, levels = c("DISCO (Pearson)", "DISCO (rank)"))

p_det <- ggplot(det, aes(alpha, dev, colour = metric)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ type_lab) +
  scale_colour_manual(values = c("DISCO (Pearson)" = "#E64B35FF", "DISCO (rank)" = "#4DBBD5FF")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "Nonlinearity strength", y = "DISCO detection (deviation from linear)")

p <- (p_cor / p_det) + plot_annotation(
  tag_levels = "a",
  title = "R2.4: DISCO(Pearson) misses monotone nonlinear dependence",
  subtitle = "Monotone nonlinear: Pearson r stays high -> DISCO does not flag (score falls below linear).\nNon-monotone: Pearson r collapses -> DISCO flags, but only because correlation vanished.")
ggsave(file.path(out, "R24_dose_response.pdf"), p, width = 8, height = 6.5)
cat("wrote R24_dose_response.pdf\n")
