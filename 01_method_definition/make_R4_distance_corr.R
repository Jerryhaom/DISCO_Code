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

# 1) correlation heatmap of DISCO scores under 7 covariance distances
# 2) standalone distance-metric benchmark (NHANES, larger fonts)
# 3) top protein loadings of modes 1/2/3 (how the three aging axes are annotated)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
genes <- expl$names
X <- as.matrix(prot[, -1]); colnames(X) <- genes; X <- scale(X)
age <- d$age; disco <- d$mc
young <- which(age <= 40); old <- which(age >= 65)

# ---------------- distance metrics (p=100, 500 samples) ----------------
set.seed(1)
test <- sample(setdiff(seq_len(nrow(X)), young), 500)
p <- 100; sel <- sample(ncol(X), p)
Xr <- X[young, sel, drop = FALSE]
n <- nrow(Xr); mu <- colMeans(Xr); S <- cov(Xr); R1 <- cor(Xr)
c1 <- (n-1)/n; c2 <- 1/(n+1)

logdet <- function(A) sum(log(eigen(A, symmetric = TRUE)$values))
matlog <- function(A){e<-eigen(A, symmetric = TRUE); e$vectors %*% diag(log(pmax(e$values,1e-12))) %*% t(e$vectors)}
matsqrt <- function(A){e<-eigen(A, symmetric = TRUE); e$vectors %*% diag(sqrt(pmax(e$values,1e-12))) %*% t(e$vectors)}
minvsqrt <- function(A){e<-eigen(A, symmetric = TRUE); e$vectors %*% diag(1/sqrt(pmax(e$values,1e-12))) %*% t(e$vectors)}

dist_euc  <- function(R2) sqrt(sum((R2 - R1)^2))
dist_chol <- function(R2) sqrt(sum((chol(R2) - chol(R1))^2))
dist_root <- function(R2) sqrt(sum((matsqrt(R2) - matsqrt(R1))^2))
dist_ler  <- function(R2) sqrt(sum((matlog(R2) - matlog(R1))^2))
dist_airm <- function(R2){ M <- minvsqrt(R1) %*% R2 %*% minvsqrt(R1); sqrt(sum(log(pmax(eigen(M, symmetric=TRUE)$values,1e-12))^2)) }
dist_jbld <- function(R2) logdet((R1+R2)/2) - 0.5*(logdet(R1) + logdet(R2))
dist_kl   <- function(R2){ p0 <- nrow(R1)
  kl12 <- 0.5*(sum(diag(solve(R2) %*% R1)) - p0 + logdet(R2) - logdet(R1))
  kl21 <- 0.5*(sum(diag(solve(R1) %*% R2)) - p0 + logdet(R1) - logdet(R2))
  0.5*(kl12 + kl21) }

dists <- c("Euclidean", "Cholesky", "Root-Euclidean", "Log-Euclidean", "Affine-invariant", "Jensen-Bregman", "KL")
fun <- list(dist_euc, dist_chol, dist_root, dist_ler, dist_airm, dist_jbld, dist_kl)
M <- matrix(NA, length(test), length(dists)); colnames(M) <- dists
for (i in seq_along(test)) {
  dd <- X[test[i], sel] - mu
  Sn <- c1*S + c2*(dd %*% t(dd)); s <- sqrt(diag(Sn)); R2 <- Sn/(s %*% t(s))
  for (j in seq_along(fun)) M[i, j] <- fun[[j]](R2)
}
C <- cor(M, method = "spearman")
saveRDS(C, file.path(out, "R4_distance_corr.rds"))

# ---- heatmap ----
ord <- dists
hm <- data.frame(x = rep(ord, each = length(ord)), y = rep(ord, length(ord)),
                 r = as.numeric(C))
hm$y <- factor(hm$y, levels = rev(ord)); hm$x <- factor(hm$x, levels = ord)
ph <- ggplot(hm, aes(x, y, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.5, colour = "black") +
  scale_fill_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF",
                       midpoint = median(hm$r), limits = c(min(hm$r), 1)) +
  theme_classic(base_size = 15) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13), legend.position = "none") +
  ggtitle("Spearman correlation of DISCO scores across distance metrics")
ggsave(file.path(out, "R4_distance_corr_heatmap.pdf"), ph, width = 8.5, height = 7)
cat("wrote R4_distance_corr_heatmap.pdf\n")
print(round(C, 3))

# ---- standalone benchmark (NHANES) ----
bench <- readRDS(file.path(root, "Result1/Cov.dis.rds"))
names_map <- c(AIRM = "Affine-invariant Riemannian", Bhattacharyya = "Bhattacharyya",
               Cholesky = "Cholesky", Euclidean = "Euclidean (unweighted)",
               Hellinger = "Hellinger", JBLD = "Jensen-Bregman log-det",
               KLDM = "Kullback-Leibler", LERM = "Log-Euclidean Riemannian",
               "Procrustes.SS" = "Procrustes (size-shape)", "Procrustes.Full" = "Procrustes (scale)")
bench$name <- names_map[bench$d]
bench$cls  <- ifelse(bench$d %in% c("LERM", "AIRM"), "Riemannian",
                     ifelse(bench$d == "Euclidean", "DISCO (Euclidean)", "Other"))
bench$name <- factor(bench$name, levels = bench$name[order(bench$ci)])
pb <- ggplot(bench, aes(ci, name, colour = cls)) +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = ci - se, xmax = ci + se), height = 0.35, linewidth = 0.6) +
  scale_colour_manual(values = c("Riemannian" = "#4DBBD5FF", "DISCO (Euclidean)" = "#E64B35FF", "Other" = "grey55"),
                      name = NULL) +
  coord_cartesian(xlim = c(0.735, 0.78)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 17),
        plot.title = element_text(size = 19, face = "bold"), legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", colour = "grey80")) +
  labs(x = "C-index (NHANES mortality)", y = NULL, title = "Distance metric robustness")
ggsave(file.path(out, "R4_distance_benchmark.pdf"), pb, width = 8.5, height = 6)
cat("wrote R4_distance_benchmark.pdf\n")

# ---- top protein loadings of modes 1/2/3 ----
Rref <- cor(X[young, , drop = FALSE])
E <- eigen(Rref, symmetric = TRUE); U <- E$vectors
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
delta <- delta/sqrt(sum(delta^2))
for (k in 1:3) {
  l <- U[, k]
  cat(sprintf("\nmode %d (eigenvalue %.1f, |aging load| %.3f):\n", k, E$values[k], abs(sum(l*delta))))
  top <- head(order(-abs(l)), 12)
  cat("  up:  ", paste(sprintf("%s(%.2f)", genes[top[l[top]>0]], l[top[l[top]>0]]), collapse=", "), "\n")
  cat("  down:", paste(sprintf("%s(%.2f)", genes[top[l[top]<0]], l[top[l[top]<0]]), collapse=", "), "\n")
}
