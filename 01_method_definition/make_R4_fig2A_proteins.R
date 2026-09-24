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

# New Fig2A: dendrogram + protein heatmap combined (pheatmap), important proteins annotated.
suppressMessages({library(pheatmap); library(dplyr); library(dendextend)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
module <- r$module

X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age

# aging deviation for importance ranking
delta <- colMeans(X[age >= 65, , drop = FALSE]) - colMeans(X[age <= 40, , drop = FALSE])
important <- names(sort(abs(delta), decreasing = TRUE))[1:25]

# protein labels: only show important ones
lab <- colnames(X); lab[!lab %in% important] <- ""

# module annotation (with biology labels)
bio_map <- c("1"="", "2"="", "3"="Senescence", "4"="Bone/fibrolysis", "5"="Cardiac",
             "6"="", "7"="Immune", "8"="Inflammation", "9"="", "10"="", "11"="", "12"="")
mod_ann <- data.frame(Module = factor(module), row.names = colnames(X))

# reorder proteins by module for cleaner heatmap, keep dendrogram
mat <- cor(X)
# color palette: use module colors for the annotation
ann_colors <- list(Module = setNames(rep(c("#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF","#8491B4FF","grey80"), each=2),
                                     levels(mod_ann$Module)))

pdf(file.path(out, "R4_fig2A_proteins.pdf"), width = 10, height = 9)
pheatmap(mat,
         labels_row = lab, labels_col = lab,
         annotation_row = mod_ann, annotation_col = mod_ann,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 5, fontsize_col = 5,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         border_color = NA,
         treeheight_row = 60, treeheight_col = 60,
         main = "Co-expression modules (important aging proteins annotated)")
dev.off()
cat("wrote R4_fig2A_proteins.pdf\n")
cat("important proteins:", paste(important, collapse=", "), "\n")
