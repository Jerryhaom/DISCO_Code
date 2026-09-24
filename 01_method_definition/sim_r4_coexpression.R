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

# R4 co-expression: show DISCO's perturbation is concentrated on co-expression network hub edges.
# Build reference co-expression network (correlation matrix). For each individual, compute
# edge-wise perturbation (r_ref - r_new)^2, and test whether high-degree (hub) edges are
# perturbed more. Uses proteomics block only (367 proteins) for a clean network.
suppressMessages({library(dplyr); library(igraph)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age

ref_age <- which(age <= 40 & !is.na(age))
Xref <- X[ref_age, , drop = FALSE]; n <- nrow(Xref)
mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)
c1 <- (n-1)/n; c2 <- 1/(n+1)

# reference co-expression network: degree = sum |r| per node
deg <- colSums(abs(Rref)) - 1   # node degree (sum of |cor|, minus self)
hub_rank <- rank(deg)           # higher = more connected

# sample of individuals (old + young) to compute edge perturbation
ok <- which(!is.na(age)); 
set.seed(42); sample_idx <- sample(ok, 2000)

# for each individual, per-edge perturbation |r_ref - r_new| for each edge; aggregate by hub degree
edge_pert <- matrix(NA, nrow = length(sample_idx), ncol = ncol(X))  # per-node aggregate
for (i in seq_along(sample_idx)) {
  dd <- X[sample_idx[i], ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Diff <- abs(Rref - Rnew)           # per-edge |change|
  diag(Diff) <- 0
  edge_pert[i, ] <- colSums(Diff)     # total perturbation incident to each node
}

# average perturbation per node (across individuals)
node_pert <- colMeans(edge_pert)
res <- data.frame(
  protein = colnames(X),
  degree = deg,
  hub_rank = hub_rank,
  node_perturbation = node_pert
)
# correlation between node degree and node perturbation (are hubs perturbed more?)
cor_deg_pert <- cor(res$degree, res$node_perturbation, method = "spearman")
cat("Spearman(degree, perturbation) =", round(cor_deg_pert, 4), "\n")

# binned view: bin proteins by degree decile, average perturbation
res <- res %>% arrange(degree) %>% mutate(bin = cut(seq_len(n()), breaks = 10, labels = FALSE))
binres <- res %>% group_by(bin) %>% summarise(deg = mean(degree), pert = mean(node_perturbation))
cat("\n=== degree decile vs perturbation ===\n"); print(binres)

saveRDS(list(res = res, cor_deg_pert = cor_deg_pert, binres = binres),
        file.path(out, "R4_coexpression.rds"))
cat("DONE\n")
