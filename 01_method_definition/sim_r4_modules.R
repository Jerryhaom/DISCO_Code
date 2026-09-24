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

# Co-expression modules + biological-process annotation + aging deviation.
suppressMessages({library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names   # gene symbols
age <- d$age; p <- ncol(X)

young <- which(age <= 40 & !is.na(age))
old   <- which(age >= 65 & !is.na(age))
Rref <- cor(X[young, , drop = FALSE])

# ---- co-expression modules (hierarchical clustering on 1-|r|) ----
hc <- hclust(as.dist(1 - abs(Rref)), method = "complete")
K <- 12
module <- cutree(hc, k = K)
cat("module sizes:\n"); print(table(module))

# ---- aging deviation (mean old - mean young) ----
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
names(delta) <- colnames(X)

# ---- gene -> ENTREZ mapping ----
cls <- readRDS(file.path(root, "Result3/Organ/UKB.Prot.class.rds"))
g2e <- setNames(cls$ENTREZID, cls$Gene)
entrez <- unname(g2e[colnames(X)])

# ---- GO enrichment per module ----
reactome <- readRDS(file.path(root, "UKB_Data/reactome_pathways.rds"))

mod_go <- list(); mod_react <- list()
for (k in sort(unique(module))) {
  genes <- colnames(X)[module == k]
  e <- entrez[module == k]; e <- e[!is.na(e)]
  # GO BP
  ego <- tryCatch(enrichGO(gene = e, universe = entrez[!is.na(entrez)],
                           OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                           qvalueCutoff = 0.2, readable = TRUE),
                  error = function(e) NULL)
  topgo <- if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) head(as.data.frame(ego)$Description, 3) else character(0)
  # Reactome (hypergeometric)
  rp <- sapply(reactome, function(pw) {
    ov <- sum(genes %in% pw$proteins)
    if (ov == 0) return(1)
    phyper(ov - 1, pw$n, p - pw$n, length(genes), lower.tail = FALSE)
  })
  rp <- p.adjust(rp, "BH")
  topreact <- names(reactome)[order(rp)][1:3]
  topreact <- sapply(topreact, function(nm) reactome[[nm]]$name)
  mod_go[[as.character(k)]] <- topgo
  mod_react[[as.character(k)]] <- topreact
  cat(sprintf("module %d (n=%d) aging|delta|mean=%.3f | GO: %s\n",
              k, sum(module==k), mean(abs(delta[module==k])),
              paste(topgo, collapse="; ")))
}

# ---- aging-hit proteins (top |delta|) ----
top_n <- 25
ord <- order(abs(delta), decreasing = TRUE)
top <- data.frame(
  protein = colnames(X)[ord][1:top_n],
  delta = delta[ord][1:top_n],
  module = module[ord][1:top_n],
  tissue = cls$specific.tissue[match(colnames(X)[ord][1:top_n], cls$Gene)]
)
cat("\n=== top aging-hit proteins ===\n"); print(top, row.names = FALSE)

saveRDS(list(Rref = Rref, module = module, delta = delta, mod_go = mod_go,
             mod_react = mod_react, top = top, hc = hc),
        file.path(out, "R4_modules_annotation.rds"))
cat("DONE\n")
