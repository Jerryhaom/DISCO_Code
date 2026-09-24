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

# Method sensitivity figures: R2.4 (nonlinear), R2.8 (age changepoint), R4.3 (high-dim).
# Reads existing R03/R24/R28/R43 CSVs + NHANES locked object. No legacy files modified.

suppressMessages({
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  library(tidyr)
  library(ggpointdensity)
})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
csv  <- disco_results_path("08_revision_analyses/outputs/R03_correlation_robustness")

# =====================================================================
# R2.4  nonlinear-dependence detection
# =====================================================================
r24 <- read.csv(file.path(csv, "R24_monotone_nonlinear.csv"))
r24$scenario <- factor(r24$scenario,
  levels = c("healthy","linear","quadratic","cubic","sigmoid","funnel"),
  labels = c("Healthy","Linear","Quadratic","Cubic (monotone)","Sigmoid (monotone)","Funnel"))
r24l <- r24 %>%
  select(scenario, DISCO_Pearson = disco_pearson, DISCO_rank = disco_rank, DM = dm) %>%
  pivot_longer(-scenario, names_to = "metric", values_to = "score")

p24 <- ggplot(r24l, aes(scenario, score, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  scale_fill_manual(values = c("DISCO_Pearson" = "#E64B35FF",
                               "DISCO_rank" = "#4DBBD5FF",
                               "DM" = "#00A087FF"),
                    labels = c("DISCO (Pearson)", "DISCO (rank)", "DM")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
        legend.position = "top", legend.title = element_blank()) +
  labs(x = "Simulated dependence", y = "Median DISCO / DM score",
       title = "Detection of nonlinear dependence") +
  annotate("text", x = 4, y = max(r24l$score), hjust = 0.5, size = 3,
           label = "rank-DISCO detects monotone nonlinear (Cubic/Sigmoid);\nPearson-DISCO misses it")
ggsave(file.path(out, "R24_nonlinear_detection.pdf"), p24, width = 7.5, height = 4.5)

# =====================================================================
# R2.8  age changepoint (NHANES DISCO ~ age hinge)
# =====================================================================
nh <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO.rds"))
d  <- nh[nh$Age >= 20 & nh$Age <= 85 & is.finite(nh$mc), c("Age", "mc")]
age <- d$Age; mc <- d$mc

grid <- seq(25, 65, 0.5)
sse  <- sapply(grid, function(cp) {
  fit <- lm(mc ~ age + pmax(age - cp, 0))
  sum(residuals(fit)^2)
})
cp  <- grid[which.min(sse)]
fit <- lm(mc ~ age + pmax(age - cp, 0))
s   <- summary(fit)
b1  <- coef(fit)[2]                      # slope before
b2  <- coef(fit)[2] + coef(fit)[3]       # slope after
# hinge significance (difference in slope) via contrast
p_hinge <- summary(lm(mc ~ age + pmax(age - cp, 0)))$coefficients[3, 4]

# binned means for clean rendering
db <- d %>% mutate(age_bin = round(age)) %>%
  group_by(age_bin) %>% summarise(mc_mean = mean(mc), n = n())
fitline <- data.frame(age = seq(20, 85, 0.5))
fitline$fit <- predict(fit, newdata = data.frame(age = fitline$age))

p28 <- ggplot() +
  geom_pointdensity(data = d[sample(nrow(d), 5000), ], aes(Age, mc), size = 0.4) +
  scale_color_distiller(palette = "Spectral") +
  geom_line(data = fitline, aes(age, fit), colour = "black", linewidth = 1) +
  geom_vline(xintercept = cp, linetype = "dashed", colour = "navy") +
  theme_classic(base_size = 12) +
  labs(x = "Age (Years)", y = "DISCO",
       title = "Age acceleration of DISCO (NHANES)") +
  theme(legend.position = "none") +
  annotate("text", x = 62, y = -2.5, size = 3.5, hjust = 0,
           label = sprintf("changepoint = %.0f yr\nslope %.3f -> %.3f /yr\nhinge P = %.1e",
                           cp, b1, b2, p_hinge))
ggsave(file.path(out, "R28_age_changepoint.pdf"), p28, width = 5.5, height = 4.5)

# =====================================================================
# R4.3  high-dimensional estimation error + shrinkage
# =====================================================================
r43  <- read.csv(file.path(csv, "R43_highdim.csv"))
r43s <- read.csv(file.path(csv, "R43_shrinkage.csv"))

# left: p/n_ref -> cor estimation error + AUC
r43l <- r43 %>% select(p_n_ref, auc, cor_est_err) %>%
  pivot_longer(-p_n_ref, names_to = "metric", values_to = "value")
p43a <- ggplot(r43l, aes(p_n_ref, value, colour = metric)) +
  geom_point(size = 2.5) + geom_line() +
  scale_colour_manual(values = c("auc" = "#E64B35FF", "cor_est_err" = "#4DBBD5FF"),
                      labels = c("AUC (detection)", "Correlation-matrix error")) +
  theme_classic(base_size = 12) +
  labs(x = "p / n_ref", y = "Value",
       title = "High-dimensional estimation error") +
  theme(legend.position = "top", legend.title = element_blank())

# right: shrinkage recovery
r43s$p <- factor(r43s$p, levels = c(200, 500), labels = c("p = 200", "p = 500"))
p43b <- ggplot(r43s, aes(lambda, auc, colour = p)) +
  geom_point(size = 2.5) + geom_line() +
  theme_classic(base_size = 12) +
  labs(x = expression(lambda ~ "(ridge shrinkage)"), y = "AUC",
       title = "Ridge shrinkage recovers detection") +
  theme(legend.position = "top", legend.title = element_blank()) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50")

library(patchwork)
p43 <- p43a + p43b + plot_annotation(tag_levels = "a")
ggsave(file.path(out, "R43_highdim.pdf"), p43, width = 9, height = 4)

cat("wrote R24, R28, R43 figures\n")
cat(sprintf("R2.8: cp=%.1f slope %.4f->%.4f hinge_p=%.2e\n", cp, b1, b2, p_hinge))
