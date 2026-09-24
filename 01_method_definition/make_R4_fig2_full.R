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

# Fig2 update: module biology bar + module eigengene correlation heatmap + variance contribution.
# Plus supplementary table of proteins per module (CSV).
suppressMessages({library(ggplot2); library(dplyr); library(reshape2); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
cls <- readRDS(file.path(root, "Result3/Organ/UKB.Prot.class.rds"))

module <- r$module; res <- r$res; eig <- r$eig
genes <- expl$names
X <- as.matrix(prot[, -1]); colnames(X) <- genes; X <- scale(X)
age <- d$age

bio_map <- c("ME8"="Inflammation / TNF-kB","ME7"="Immune / cytokines","ME5"="Cardiac / natriuretic",
             "ME3"="Senescence / neuro","ME4"="Bone / fibrinolysis")
bio <- function(m) ifelse(m %in% names(bio_map), unname(bio_map[m]), "Mixed")

# order modules by |r_disco|
ord <- res$module[order(-abs(res$r_disco))]

# ---- module eigengene correlation heatmap ----
mods <- names(eig)
E <- do.call(cbind, eig); colnames(E) <- mods
Me <- cor(E)[ord, ord]
Me[upper.tri(Me)] <- NA
hm <- reshape2::melt(Me, na.rm = TRUE)
colnames(hm) <- c("m1","m2","cor")
hm$m1 <- factor(hm$m1, levels = ord); hm$m2 <- factor(hm$m2, levels = rev(ord))
lab <- sapply(ord, bio)
ph <- ggplot(hm, aes(m1, m2, fill = cor)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "cor") +
  scale_x_discrete(labels = lab, position = "top") +
  scale_y_discrete(labels = rev(lab)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
        axis.text.y = element_text(size = 9), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "Module eigengene correlations")

# ---- module variance contribution ----
varfr <- sapply(mods, function(m) {
  M <- X[, module == as.integer(sub("ME","",m)), drop = FALSE]
  eigen(cor(M), symmetric = TRUE, only.values = TRUE)$values[1] / ncol(M)
})
names(varfr) <- mods
vd <- data.frame(module = ord, vf = varfr[ord]); vd$module <- factor(vd$module, levels = ord)
vd$bio <- sapply(ord, bio)
pv <- ggplot(vd, aes(module, vf, fill = bio)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", vf)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Inflammation / TNF-kB"="#E64B35FF","Immune / cytokines"="#4DBBD5FF",
                               "Cardiac / natriuretic"="#00A087FF","Senescence / neuro"="#F39B7FFF",
                               "Bone / fibrinolysis"="#8491B4FF","Mixed"="grey80"), guide="none") +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Variance explained by eigengene", title = "Module variance contribution")

# ---- module biology bar ----
res$bio <- sapply(res$module, bio)
res <- res[order(-abs(res$r_disco)), ]; res$module <- factor(res$module, levels = res$module)
pb <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = bio), angle = 90, hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c("Inflammation / TNF-kB"="#E64B35FF","Immune / cytokines"="#4DBBD5FF",
                               "Cardiac / natriuretic"="#00A087FF","Senescence / neuro"="#F39B7FFF",
                               "Bone / fibrinolysis"="#8491B4FF","Mixed"="grey80"), guide="none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_classic(base_size = 13) +
  labs(x = "Co-expression module", y = "|Correlation with DISCO|", title = "DISCO driven by biological modules")

fig2 <- (pb | ph) / (pv | plot_spacer())
ggsave(file.path(out, "R4_fig2_modules.pdf"), fig2, width = 12, height = 8)
cat("wrote R4_fig2_modules.pdf\n")

# ---- supplementary table: proteins per module ----
delta <- colMeans(X[d$age >= 65, , drop = FALSE]) - colMeans(X[d$age <= 40, , drop = FALSE])
tab <- data.frame()
for (m in sort(unique(module))) {
  g <- genes[module == m]
  tab <- rbind(tab, data.frame(
    module = paste0("ME", m),
    r_disco = res$r_disco[match(paste0("ME", m), res$module)],
    protein = g,
    tissue = cls$specific.tissue[match(g, cls$Gene)],
    aging_delta = round(delta[g], 3)
  ))
}
tab <- tab[order(-abs(tab$r_disco), -abs(tab$aging_delta)), ]
write.csv(tab, file.path(out, "R4_module_proteins.csv"), row.names = FALSE)
cat("wrote R4_module_proteins.csv (", nrow(tab), "proteins)\n")
