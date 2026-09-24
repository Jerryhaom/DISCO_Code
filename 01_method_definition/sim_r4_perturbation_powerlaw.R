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

# Perturbation CC(N+1)-CC(N): is its distribution also power-law (structured)?
# Parallel mclapply over a sample of individuals. Also edge-strength vs perturbation (A)
# and degree vs perturbation (C).
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
c1 <- (n-1)/n; c2 <- 1/(n+1)

set.seed(42); sample_idx <- sample(old, 300)

perturb_one <- function(i) {
  dd <- X[i, ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Delta <- abs(Rref - Rnew); diag(Delta) <- 0
  Delta
}
ncores <- detectCores()
Dlist <- mclapply(sample_idx, perturb_one, mc.cores = ncores)

# aggregate |Delta| distribution across all edges and individuals
Delta_all <- unlist(lapply(Dlist, function(M) M[lower.tri(M)]))
Delta_all <- Delta_all[Delta_all > 0]

# reference edge strength (|r_ref| off-diagonal)
off_r <- abs(Rref[lower.tri(Rref)])

# ---- power-law check of perturbation distribution (log-log) ----
h <- hist(Delta_all, breaks = 100, plot = FALSE)
keep <- h$counts > 0 & h$mids > 0
fit <- lm(log(h$counts[keep]) ~ log(h$mids[keep]))
gamma_p <- -coef(fit)[2]
cat("perturbation |Delta| log-log slope =", round(gamma_p, 3), "\n")

# ---- A: edge strength vs perturbation ----
# for each edge, average perturbation across individuals
Delta_mean <- Reduce(`+`, Dlist) / length(Dlist)
edge <- data.frame(r_ref = off_r, delta = Delta_mean[lower.tri(Rref)])
edge <- edge %>% mutate(bin = cut(r_ref, breaks = 10, include.lowest = TRUE))
edgeA <- edge %>% group_by(bin) %>% summarise(r = mean(r_ref), delta = mean(delta))
cor_rd <- cor(edge$r_ref, edge$delta, method = "spearman")
cat("Spearman(edge strength, perturbation) =", round(cor_rd, 4), "\n")

# ---- C: degree vs perturbation ----
deg <- colSums(abs(Rref)) - 1
node_pert <- colMeans(Reduce(`+`, Dlist) / length(Dlist))
cor_dp <- cor(deg, node_pert, method = "spearman")
cat("Spearman(degree, node perturbation) =", round(cor_dp, 4), "\n")

saveRDS(list(Delta_all = Delta_all, gamma_p = gamma_p, edgeA = edgeA,
             cor_rd = cor_rd, cor_dp = cor_dp),
        file.path(out, "R4_perturbation_powerlaw.rds"))

# figures
p1 <- ggplot(data.frame(d = Delta_all), aes(d)) +
  geom_histogram(bins = 100, fill = "#E64B35FF", colour = NA) +
  scale_x_log10() + scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(x = "|CC(N+1) - CC(N)| per edge", y = "Count (log)",
       title = "Dysregulation perturbation distribution")
ggsave(file.path(out, "R4_perturbation_distribution.pdf"), p1, width = 6, height = 4.5)

p2 <- ggplot(edgeA, aes(r, delta)) +
  geom_point(colour = "#4DBBD5FF") + geom_smooth(method = "lm", se = FALSE, colour = "grey30") +
  theme_classic(base_size = 14) +
  labs(x = "Edge strength |r_ref|", y = "Mean perturbation",
       title = sprintf("Perturbation vs co-expression strength (rho=%.2f)", cor_rd))
ggsave(file.path(out, "R4_edge_strength.pdf"), p2, width = 6, height = 4.5)
cat("DONE\n")
