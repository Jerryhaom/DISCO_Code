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

# Top-10 eigenvalues of proteomics correlation matrix + correlation of each eigenmode (PC) with DISCO.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; disco <- d$mc

young <- which(age <= 40)
Rref <- cor(X[young, , drop = FALSE])
eig <- eigen(Rref, symmetric = TRUE)
lam <- eig$values; U <- eig$vectors

# top 10 eigenvalues
top10 <- data.frame(k = 1:10, lambda = lam[1:10], var_explained = lam[1:10]/sum(lam))

# PC scores (project full data onto eigenvectors) and correlation with DISCO
PC <- X %*% U
r_disco <- sapply(1:10, function(k) cor(PC[, k], disco, method = "spearman"))
r_age   <- sapply(1:10, function(k) cor(PC[, k], age, method = "spearman"))
pctab <- data.frame(k = 1:10, r_disco = r_disco, r_age = r_age)

cat("=== top 10 eigenvalues ===\n"); print(top10, row.names = FALSE)
cat("\n=== eigenmode (PC) vs DISCO / age ===\n"); print(pctab, row.names = FALSE)

saveRDS(list(top10 = top10, pctab = pctab), file.path(out, "R4_eigen_disco.rds"))

p1 <- ggplot(top10, aes(factor(k), lambda)) +
  geom_col(fill = "#4DBBD5FF", colour = "black", width = 0.6) +
  geom_text(aes(label = round(lambda)), vjust = -0.3, size = 3.5) +
  theme_classic(base_size = 13) +
  labs(x = "Eigenmode rank", y = "Eigenvalue", title = "Top-10 eigenvalues (global factor)")

p2 <- ggplot(pctab, aes(factor(k), r_disco, fill = r_disco > 0)) +
  geom_col(colour = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r_disco)), vjust = ifelse(pctab$r_disco>0,-0.3,1.3), size = 3.5) +
  scale_fill_manual(values = c("TRUE"="#E64B35FF","FALSE"="#4DBBD5FF"), guide="none") +
  geom_hline(yintercept = 0, colour = "grey40") +
  theme_classic(base_size = 13) +
  labs(x = "Eigenmode rank", y = "Correlation with DISCO", title = "Eigenmode vs DISCO")

p <- p1 | p2
ggsave(file.path(out, "R4_eigen_disco.pdf"), p, width = 10, height = 4.5)
cat("wrote R4_eigen_disco.pdf\n")
