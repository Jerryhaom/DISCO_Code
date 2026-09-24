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

# Complete R4_fig2: (A) module derivation dendrogram, (B) module relationships heatmap,
# (C) module-DISCO bar, (D) module biology (protein-level heatmap).
suppressMessages({library(ggplot2); library(dplyr); library(reshape2); library(patchwork); library(dendextend)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
module <- r$module; res <- r$res; eig <- r$eig
genes <- expl$names
X <- as.matrix(prot[, -1]); colnames(X) <- genes; X <- scale(X)

bio_map <- c("ME8"="Inflammation / TNF-kB","ME7"="Immune / cytokines","ME5"="Cardiac / natriuretic",
             "ME3"="Senescence / neuro","ME4"="Bone / fibrinolysis")
bio <- function(m) ifelse(m %in% names(bio_map), unname(bio_map[m]), "Mixed")
colpal <- c("Inflammation / TNF-kB"="#E64B35FF","Immune / cytokines"="#4DBBD5FF",
            "Cardiac / natriuretic"="#00A087FF","Senescence / neuro"="#F39B7FFF",
            "Bone / fibrinolysis"="#8491B4FF","Mixed"="grey80")

ord <- res$module[order(-abs(res$r_disco))]

# ---- A: dendrogram (how modules derived) ----
hc <- hclust(as.dist(1 - abs(cor(X))), method = "complete")
dend <- as.dendrogram(hc) %>% color_branches(k = 12) %>% set("labels_cex", 0.05)
ggd <- as.ggdend(dend)
pA <- ggplot(ggd$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               colour = ifelse(is.na(ggd$segments$col), "grey60", ggd$segments$col), linewidth = 0.3) +
  coord_flip() + scale_y_reverse(expand = c(0.1, 0)) +
  theme_void() +
  labs(title = "A. Module derivation (hierarchical clustering)")

# ---- B: module relationships (eigengene heatmap) ----
mods <- names(eig)
E <- do.call(cbind, eig); colnames(E) <- mods
Me <- cor(E)[ord, ord]
Me[upper.tri(Me)] <- NA
hm <- reshape2::melt(Me, na.rm = TRUE); colnames(hm) <- c("m1","m2","cor")
hm$m1 <- factor(hm$m1, levels = ord); hm$m2 <- factor(hm$m2, levels = rev(ord))
lab <- sapply(ord, bio)
pB <- ggplot(hm, aes(m1, m2, fill = cor)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "cor") +
  scale_x_discrete(labels = lab, position = "top") +
  scale_y_discrete(labels = rev(lab)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_text(size = 8), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "B. Module relationships")

# ---- C: module-DISCO bar ----
res$bio <- sapply(res$module, bio)
res <- res[order(-abs(res$r_disco)), ]; res$module <- factor(res$module, levels = res$module)
pC <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = bio), angle = 90, hjust = -0.05, size = 3) +
  scale_fill_manual(values = colpal, guide = "none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_classic(base_size = 12) +
  labs(x = "Module", y = "|Correlation with DISCO|", title = "C. Module vs DISCO")

# ---- D: protein-level heatmap (module biology) ----
prot_order <- order(module)
Rsub <- cor(X)[prot_order, prot_order]
boundaries <- cumsum(as.numeric(table(module)))
Rsub_disp <- Rsub
Rsub_disp[upper.tri(Rsub_disp)] <- NA
hm2 <- reshape2::melt(Rsub_disp, na.rm = TRUE); colnames(hm2) <- c("i","j","cor")
pD <- ggplot(hm2, aes(i, j, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "cor") +
  geom_vline(xintercept = boundaries + 0.5, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = boundaries + 0.5, colour = "black", linewidth = 0.2) +
  coord_fixed() +
  theme_void(base_size = 10) +
  theme(legend.position = "right") +
  labs(title = "D. Module biology (proteins grouped by module)")
ggsave(file.path(out, "R4_fig2_modules.pdf"), (pA | pB) / (pC | pD), width = 14, height = 10)
cat("wrote R4_fig2_modules.pdf\n")
