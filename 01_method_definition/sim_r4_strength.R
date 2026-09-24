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

# Weighted degree (strength = sum |cor|), no threshold. Power-law check + strength vs perturbation.
suppressMessages({library(ggplot2); library(dplyr); library(parallel)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; p <- ncol(X)

young <- which(age <= 40 & !is.na(age))
old   <- which(age >= 65 & !is.na(age))
Xref <- X[young, , drop = FALSE]; n <- nrow(Xref)
mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)

# ---- weighted degree (strength) = sum |cor| ----
strength <- colSums(abs(Rref)) - 1   # minus self (diag=1)
cat("strength range:", round(range(strength), 1), " mean:", round(mean(strength), 1), "\n")

# survival function (complementary CDF) on log-log
s_sorted <- sort(strength, decreasing = TRUE)
ccdf <- data.frame(s = s_sorted, P = seq_along(s_sorted) / length(s_sorted))  # P(strength > s)
# fit power-law on the middle range (exclude extremes)
fit <- lm(log(P) ~ log(s), data = ccdf[ccdf$P > 0.01 & ccdf$P < 0.99, ])
gamma <- -coef(fit)[2]   # P(S>s) ~ s^(-(gamma-1)), so slope = -(gamma-1)
cat("CCDF log-log slope =", round(coef(fit)[2], 3), " -> gamma =", round(gamma, 3), "\n")

# ---- strength vs node perturbation (C) ----
c1 <- (n-1)/n; c2 <- 1/(n+1)
set.seed(42); sample_idx <- sample(old, 300)
perturb_one <- function(i) {
  dd <- X[i, ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Delta <- abs(Rref - Rnew); diag(Delta) <- 0; Delta
}
Dlist <- mclapply(sample_idx, perturb_one, mc.cores = detectCores())
Delta_mean <- Reduce(`+`, Dlist) / length(Dlist)
node_pert <- colMeans(Delta_mean)
cor_dp <- cor(strength, node_pert, method = "spearman")
cat("Spearman(strength, node perturbation) =", round(cor_dp, 4), "\n")

saveRDS(list(strength = strength, ccdf = ccdf, gamma = gamma, cor_dp = cor_dp, node_pert = node_pert),
        file.path(out, "R4_strength.rds"))

# figures
p1 <- ggplot(ccdf, aes(s, P)) +
  geom_point(colour = "#4DBBD5FF") +
  scale_x_log10() + scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(x = "Strength (sum |cor|)", y = "P(strength > s)",
       title = sprintf("Weighted degree distribution (slope %.2f)", coef(fit)[2])) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.5)
ggsave(file.path(out, "R4_strength_distribution.pdf"), p1, width = 6, height = 4.5)

p2 <- ggplot(data.frame(s = strength, p = node_pert), aes(s, p)) +
  geom_point(colour = "#E64B35FF") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30") +
  theme_classic(base_size = 14) +
  labs(x = "Strength (sum |cor|)", y = "Mean perturbation",
       title = sprintf("Strength vs perturbation (rho=%.2f)", cor_dp))
ggsave(file.path(out, "R4_strength_vs_perturbation.pdf"), p2, width = 6, height = 4.5)
cat("DONE\n")
