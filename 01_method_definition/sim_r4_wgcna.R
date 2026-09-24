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

# Part 3 (simplified WGCNA): co-expression modules -> module eigengene vs DISCO.
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco_score <- d$mc   # DISCO (already computed)

# ---- simplified WGCNA: soft-threshold adjacency + TOM + hierarchical clustering ----
beta <- 6
R <- cor(X)
A <- abs(R)^beta                      # soft-thresholded adjacency
# topological overlap matrix (TOM)
k <- colSums(A)
TOM <- matrix(0, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X)))
for (i in 1:ncol(X)) {
  for (j in i:ncol(X)) {
    if (i == j) next
    num <- sum(A[i, ] * A[, j]) + A[i, j]
    den <- min(k[i], k[j]) + 1 - A[i, j]
    tom <- num / den
    TOM[i, j] <- TOM[j, i] <- tom
  }
}
hc <- hclust(as.dist(1 - TOM), method = "average")
K <- 12
module <- cutree(hc, k = K)
cat("module sizes:", paste(table(module), collapse=", "), "\n")

# ---- module eigengene (first PC) and its correlation with DISCO / age ----
mod_eig <- lapply(sort(unique(module)), function(k) {
  M <- X[, module == k, drop = FALSE]
  prcomp(M, center = TRUE, scale. = TRUE)$x[, 1]   # eigengene
})
names(mod_eig) <- paste0("ME", sort(unique(module)))

res <- data.frame(
  module = paste0("ME", sort(unique(module))),
  n = as.numeric(table(module)),
  r_disco = sapply(mod_eig, function(me) cor(me, disco_score, method = "spearman")),
  r_age   = sapply(mod_eig, function(me) cor(me, age, method = "spearman"))
)
res <- res[order(-abs(res$r_disco)), ]
cat("\n=== module eigengene vs DISCO / age ===\n"); print(res, row.names = FALSE)

saveRDS(list(module = module, eig = mod_eig, res = res, X = X),
        file.path(out, "R4_wgcna.rds"))

p <- ggplot(res, aes(r_age, r_disco, size = n)) +
  geom_point(colour = "#E64B35FF") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_text(aes(label = module), vjust = -1, size = 3.5) +
  theme_classic(base_size = 14) +
  labs(x = "Correlation with age (Spearman)", y = "Correlation with DISCO (Spearman)",
       size = "module size", title = "Co-expression modules vs DISCO / age")
ggsave(file.path(out, "R4_wgcna_eigengene.pdf"), p, width = 7, height = 5.5)
cat("wrote R4_wgcna_eigengene.pdf\n")
