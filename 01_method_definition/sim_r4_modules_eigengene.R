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

# Balanced co-expression modules (hierarchical clustering on |cor|, K=12) -> eigengene vs DISCO.
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco <- d$mc

# balanced modules: hierarchical clustering on 1-|cor|, complete linkage, K=12
hc <- hclust(as.dist(1 - abs(cor(X))), method = "complete")
module <- cutree(hc, k = 12)
cat("module sizes:", paste(table(module), collapse=", "), "\n")

# eigengene per module
mods <- sort(unique(module))
eig <- lapply(mods, function(m) prcomp(X[, module == m], center=TRUE, scale.=TRUE)$x[,1])
names(eig) <- paste0("ME", mods)

res <- data.frame(
  module = paste0("ME", mods),
  n = as.numeric(table(module)),
  r_disco = sapply(eig, function(e) cor(e, disco, method="spearman")),
  r_age   = sapply(eig, function(e) cor(e, age, method="spearman"))
)
res <- res[order(-abs(res$r_disco)), ]
cat("\n=== module eigengene vs DISCO / age ===\n"); print(res, row.names = FALSE)

saveRDS(list(module = module, eig = eig, res = res), file.path(out, "R4_modules_eigengene.rds"))

# figure 1: scatter r_age vs r_disco (sized by n)
p1 <- ggplot(res, aes(r_age, r_disco, size = n)) +
  geom_point(colour = "#E64B35FF") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_text(aes(label = module), vjust = -1, size = 4) +
  theme_classic(base_size = 14) +
  labs(x = "Eigengene correlation with age", y = "Eigengene correlation with DISCO",
       size = "module size", title = "Co-expression modules vs DISCO / age")
ggsave(file.path(out, "R4_modules_eigengene_scatter.pdf"), p1, width = 7, height = 5.5)

# figure 2: sorted bar of |r_disco|
res2 <- res[order(-abs(res$r_disco)), ]
res2$module <- factor(res2$module, levels = res2$module)
p2 <- ggplot(res2, aes(module, r_disco, fill = r_disco > 0)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r_disco)), vjust = ifelse(res2$r_disco>0, -0.3, 1.3), size = 3.5) +
  scale_fill_manual(values = c("TRUE"="#E64B35FF","FALSE"="#4DBBD5FF"), guide="none") +
  theme_classic(base_size = 14) +
  labs(x = "Module", y = "Correlation with DISCO (Spearman)",
       title = "Which co-expression modules drive DISCO")
ggsave(file.path(out, "R4_modules_eigengene_bar.pdf"), p2, width = 7, height = 4.5)
cat("DONE\n")
