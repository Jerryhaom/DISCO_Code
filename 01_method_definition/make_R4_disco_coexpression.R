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

# DISCO <-> co-expression geometry (2x2, R43 style, larger fonts).
# a: DISCO = co-expression perturbation (definitional)
# b: Riemannian (log-Euclidean) ~ Euclidean (linearization valid)
# c: distance-metric robustness benchmark (NHANES)
# d: DISCO aligns with the aging structure in co-expression space
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; disco <- d$mc
young <- which(age <= 40); old <- which(age >= 65)

# ---------- geometry data (a & b): p=100, 800 test samples ----------
set.seed(1)
test <- sample(setdiff(seq_len(nrow(X)), young), 800)
p <- 100; sel <- sample(ncol(X), p)
Xr <- X[young, sel, drop = FALSE]
n  <- nrow(Xr); mu <- colMeans(Xr); S <- cov(Xr); Rref <- cor(Xr)
c1 <- (n - 1)/n; c2 <- 1/(n + 1)
matlog <- function(A){ e <- eigen(A, symmetric = TRUE)
  e$vectors %*% diag(log(pmax(e$values, 1e-12))) %*% t(e$vectors) }
Lref <- matlog(Rref)
eucl <- ler <- numeric(length(test))
for (i in seq_along(test)) {
  dd <- X[test[i], sel] - mu
  Sn <- c1*S + c2*(dd %*% t(dd))
  s  <- sqrt(diag(Sn)); Rn <- Sn/(s %*% t(s))
  eucl[i] <- sqrt(sum((Rn - Rref)^2))
  ler[i]  <- sqrt(sum((matlog(Rn) - Lref)^2))
}
disco_score <- log(n^2 * eucl^2)

# ---------- d: aging loading vs DISCO correlation ----------
Rref_full <- cor(X[young, , drop = FALSE])
E <- eigen(Rref_full, symmetric = TRUE); U <- E$vectors
PC <- X %*% U
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
delta <- delta/sqrt(sum(delta^2))
ck <- as.numeric(t(U) %*% delta)
rk <- sapply(seq_len(ncol(X)), function(k) cor(PC[, k], disco, method = "spearman"))

# ---------- theme (larger fonts) ----------
thm <- theme_classic(base_size = 17) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15, colour = "black"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15), legend.title = element_blank())

# ---------- a: DISCO = perturbation ----------
pa <- ggplot(data.frame(x = eucl, y = disco_score), aes(x, y)) +
  geom_point(colour = "#E64B35FF", size = 1.1, alpha = 0.5) +
  scale_x_log10() + thm +
  labs(x = "Co-expression perturbation ||R(N+1)-R(N)|| (log)", y = "DISCO score",
       title = "a. DISCO is the co-expression perturbation") +
  annotate("text", x = max(eucl)*0.5, y = min(disco_score)*0.9, size = 5, hjust = 0,
           label = "DISCO = log(n²·||ΔR||²)", colour = "grey20")

# ---------- b: Riemannian ~ Euclidean ----------
rho_geom <- cor(eucl, ler)
pb <- ggplot(data.frame(x = eucl, y = ler), aes(x, y)) +
  geom_point(colour = "#4DBBD5FF", size = 1.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
  thm + labs(x = "Euclidean perturbation", y = "Riemannian (log-Euclidean) distance",
             title = "b. Euclidean ≈ Riemannian") +
  annotate("text", x = max(eucl)*0.5, y = min(ler), size = 5, hjust = 0,
           label = sprintf("r = %.2f", rho_geom), colour = "grey20")

# ---------- c: distance benchmark (NHANES) ----------
bench <- readRDS(file.path(root, "Result1/Cov.dis.rds"))
names_map <- c(AIRM = "Affine-invariant Riemannian", Bhattacharyya = "Bhattacharyya",
               Cholesky = "Cholesky", Euclidean = "Euclidean (unweighted)",
               Hellinger = "Hellinger", JBLD = "Jensen-Bregman log-det",
               KLDM = "Kullback-Leibler", LERM = "Log-Euclidean Riemannian",
               "Procrustes.SS" = "Procrustes (size-shape)", "Procrustes.Full" = "Procrustes (scale)")
bench$name <- names_map[bench$d]
bench$cls  <- ifelse(bench$d %in% c("LERM", "AIRM"), "Riemannian",
                     ifelse(bench$d == "Euclidean", "Euclidean", "Other"))
bench$name <- factor(bench$name, levels = bench$name[order(bench$ci)])
pc <- ggplot(bench, aes(ci, name, colour = cls)) +
  geom_point(size = 2.2) +
  geom_errorbarh(aes(xmin = ci - se, xmax = ci + se), height = 0.3) +
  scale_colour_manual(values = c(Riemannian = "#4DBBD5FF", Euclidean = "#E64B35FF", Other = "grey55")) +
  coord_cartesian(xlim = c(0.73, 0.78)) +
  thm + labs(x = "C-index", y = NULL, title = "c. Distance metric (NHANES)")

# ---------- d: DISCO aligns with aging ----------
top3 <- order(-ck^2)[1:3]
dfd <- data.frame(ck = ck, rk = rk, k = seq_along(ck), top = seq_along(ck) %in% top3)
rho_d <- cor(ck, rk, method = "spearman")
pd <- ggplot(dfd, aes(ck, rk)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = top), size = 1.4) +
  scale_colour_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "#4DBBD5FF"), guide = "none") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
  ggrepel::geom_text_repel(data = dfd[dfd$top, ], aes(label = k), size = 5,
                           nudge_x = 0.02, nudge_y = 0.03) +
  thm + labs(x = "Aging loading (signed)", y = "DISCO correlation (Spearman)",
             title = "d. DISCO aligns with aging structure") +
  annotate("text", x = min(ck), y = max(rk), size = 5, hjust = 0, vjust = 1,
           label = sprintf("Spearman = %.2f", rho_d))

fig <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20, face = "bold"))
ggsave(file.path(out, "R4_disco_coexpression.pdf"), fig, width = 11, height = 9)
cat("wrote R4_disco_coexpression.pdf | geom r =", round(rho_geom, 3), "| disco rho =", round(rho_d, 3), "\n")
