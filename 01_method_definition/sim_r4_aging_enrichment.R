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

# Enrich the aging-dysregulation signature (top |delta| proteins) as a group.
suppressMessages({library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
r <- readRDS(file.path(out, "R4_modules_annotation.rds"))
cls <- readRDS(file.path(root, "Result3/Organ/UKB.Prot.class.rds"))
reactome <- readRDS(file.path(root, "UKB_Data/reactome_pathways.rds"))

delta <- r$delta; genes <- names(delta); p <- length(genes)
top <- genes[order(abs(delta), decreasing = TRUE)][1:40]   # aging signature

# GO enrichment (ENTREZ)
g2e <- setNames(cls$ENTREZID, cls$Gene)
univ_e <- unname(g2e[genes]); univ_e <- univ_e[!is.na(univ_e)]
top_e <- unname(g2e[top]); top_e <- top_e[!is.na(top_e)]
ego <- enrichGO(gene = top_e, universe = univ_e, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5,
                minGSSize = 3, readable = TRUE)
cat("=== GO BP enrichment (aging signature) ===\n")
print(head(as.data.frame(ego)[, c("Description","pvalue","p.adjust","Count")], 15), row.names = FALSE)

# Reactome enrichment
pv <- sapply(reactome, function(pw) {
  ov <- sum(top %in% pw$proteins); m <- sum(pw$proteins %in% genes)
  if (ov == 0 || m == 0) 1 else phyper(ov - 1, m, p - m, length(top), lower.tail = FALSE)
})
pv <- p.adjust(pv, "BH")
ord <- order(pv)
cat("\n=== Reactome enrichment (aging signature) ===\n")
for (i in ord[1:12]) cat(sprintf("%s (ov=%d, P=%.2e)\n", reactome[[i]]$name,
                                 sum(top %in% reactome[[i]]$proteins), pv[i]))

# tissue annotation of aging signature
cat("\n=== aging-signature proteins: tissue ===\n")
df <- data.frame(protein = top, delta = round(delta[top], 3),
                 tissue = cls$specific.tissue[match(top, cls$Gene)])
print(df, row.names = FALSE)

saveRDS(list(ego = as.data.frame(ego), reactome = data.frame(pathway=names(reactome)[ord[1:12]], pv=pv[ord[1:12]]), top = df),
        file.path(out, "R4_aging_signature_enrichment.rds"))
cat("DONE\n")
