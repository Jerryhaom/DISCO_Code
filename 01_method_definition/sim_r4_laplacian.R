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

# Intuitive: (1) Laplacian embedding of protein co-expression network -> clusters?
#            (2) cluster individual perturbations CC(N+1)-CC(N) -> dysregulation types?
suppressMessages({library(igraph); library(ggplot2); library(dplyr); library(parallel)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names
age <- d$age; p <- ncol(X)

young <- which(age <= 40 & !is.na(age))
old   <- which(age >= 65 & !is.na(age))
Xref <- X[young, , drop = FALSE]; n <- nrow(Xref)
mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)

# ================= Part 1: Laplacian embedding of protein network =================
A <- abs(Rref); diag(A) <- 0          # weighted adjacency (co-expression strength)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
comm <- cluster_louvain(g, weights = E(g)$weight)
cat("Louvain communities:", max(membership(comm)), "\n")

# normalized graph Laplacian eigenmap (Fiedler + next -> 2D)
D <- rowSums(A); Dinv <- 1/sqrt(D)
Ln <- diag(p) - (Dinv %*% t(Dinv)) * A
ev <- eigen(Ln, symmetric = TRUE)
coord <- ev$vectors[, c(p, p-1)]      # smallest two non-trivial eigenvectors
colnames(coord) <- c("Fiedler", "EV2")

df1 <- data.frame(coord, protein = colnames(X), module = factor(membership(comm)))
p1 <- ggplot(df1, aes(Fiedler, EV2, colour = module)) +
  geom_point(size = 1.8) +
  theme_classic(base_size = 14) +
  labs(colour = "module", title = "Protein co-expression network (Laplacian embedding)") +
  theme(legend.position = "right")
ggsave(file.path(out, "R4_laplacian_network.pdf"), p1, width = 7, height = 5.5)

# ================= Part 2: cluster individual perturbations =================
c1 <- (n-1)/n; c2 <- 1/(n+1)
set.seed(42); sample_idx <- sample(old, 400)
perturb_vec <- function(i) {
  dd <- X[i, ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Delta <- Rref - Rnew               # perturbation matrix (signed)
  Delta[lower.tri(Delta)]            # vectorize
}
Pmat <- do.call(rbind, mclapply(sample_idx, perturb_vec, mc.cores = detectCores()))
# PCA on the perturbation vectors
pca <- prcomp(Pmat, center = TRUE, scale. = FALSE)
km <- kmeans(pca$x[, 1:3], centers = 3, nstart = 25)
cat("perturbation clusters (k=3) sizes:", paste(table(km$cluster), collapse=", "), "\n")

df2 <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], cluster = factor(km$cluster))
p2 <- ggplot(df2, aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 1.8) +
  theme_classic(base_size = 14) +
  labs(colour = "type", title = "Individual perturbation CC(N+1)-CC(N) (clustered)") +
  theme(legend.position = "right")
ggsave(file.path(out, "R4_perturbation_clusters.pdf"), p2, width = 7, height = 5.5)

saveRDS(list(coord = df1, comm = comm, pca = pca$x, km = km),
        file.path(out, "R4_laplacian_cluster.rds"))
cat("DONE\n")
