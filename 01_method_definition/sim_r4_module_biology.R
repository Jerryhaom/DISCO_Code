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

# Biological annotation of the top DISCO-driving modules (ME8, ME7, ME5, ME3, ME4, ME1).
suppressMessages({library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
cls <- readRDS(file.path(root, "Result3/Organ/UKB.Prot.class.rds"))
reactome <- readRDS(file.path(root, "UKB_Data/reactome_pathways.rds"))

module <- r$module
genes <- expl$names
delta <- colMeans(scale(readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))[, -1])[d$age >= 65, , drop = FALSE]) -
        colMeans(scale(readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))[, -1])[d$age <= 40, , drop = FALSE])
names(delta) <- genes

g2e <- setNames(cls$ENTREZID, cls$Gene)
univ_e <- unname(g2e[genes]); univ_e <- univ_e[!is.na(univ_e)]

targets <- c("ME8","ME7","ME5","ME3","ME4","ME1")
for (tm in targets) {
  m <- as.integer(sub("ME","",tm))
  g <- genes[module == m]
  cat("\n========== ", tm, " (n=", length(g), ") ==========\n", sep="")
  # top proteins by |delta| with tissue
  dd <- data.frame(protein = g, delta = delta[g], tissue = cls$specific.tissue[match(g, cls$Gene)])
  dd <- dd[order(-abs(dd$delta)), ]
  cat("top proteins:", paste(head(dd$protein, 15), collapse=", "), "\n")
  cat("tissues:", paste(names(sort(table(dd$tissue), decreasing=TRUE))[1:5], collapse=", "), "\n")
  # GO enrichment
  e <- unname(g2e[g]); e <- e[!is.na(e)]
  ego <- tryCatch(enrichGO(gene=e, universe=univ_e, OrgDb=org.Hs.eg.db, ont="BP",
                           pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.5,
                           minGSSize=3, readable=TRUE), error=function(e) NULL)
  if (!is.null(ego) && nrow(as.data.frame(ego))>0) {
    cat("GO:", paste(head(as.data.frame(ego)$Description, 3), collapse=" | "), "\n")
  } else cat("GO: (none significant)\n")
  # Reactome enrichment (correct hypergeometric)
  pv <- sapply(reactome, function(pw) {
    ov <- sum(g %in% pw$proteins); mm <- sum(pw$proteins %in% genes)
    if (ov==0 || mm==0) 1 else phyper(ov-1, mm, length(genes)-mm, length(g), lower.tail=FALSE)
  })
  pv <- p.adjust(pv, "BH"); ord <- order(pv)
  tops <- head(ord[order(pv)[1:3]], 3)
  cat("Reactome:", paste(sapply(ord[1:3], function(i) sprintf("%s(P=%.1e)", reactome[[i]]$name, pv[i])), collapse=" | "), "\n")
}
