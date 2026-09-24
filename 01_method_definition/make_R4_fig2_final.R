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

# Final Figure 2 (4 panels): UMAP clustering, module relationships, module-DISCO, important proteins.
suppressMessages({library(ggplot2); library(dplyr); library(reshape2); library(patchwork); library(uwot)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
module <- r$module; res <- r$res; eig <- r$eig
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco <- d$mc

bio_map <- c("ME8"="Inflammation","ME7"="Immune","ME5"="Cardiac","ME3"="Senescence","ME4"="Bone")
bio <- function(m) ifelse(paste0("ME",m) %in% names(bio_map), unname(bio_map[paste0("ME",m)]), "Other")
colpal <- c("Inflammation"="#E64B35FF","Immune"="#4DBBD5FF","Cardiac"="#00A087FF",
            "Senescence"="#F39B7FFF","Bone"="#8491B4FF","Other"="grey80")

tt <- theme(plot.title = element_text(size = 14, face = "bold"))
ax <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11, colour = "black"))

# ---- A: UMAP after removing global PC ----
PC1 <- prcomp(X, center=TRUE, scale.=TRUE)$x[,1]
Xr <- residuals(lm(X ~ PC1))
Rr <- cor(Xr)
set.seed(1); um <- uwot::umap(Rr, n_neighbors = 15, n_components = 2, n_epochs = 500)
dfA <- data.frame(x = um[,1], y = um[,2], module = factor(module))
dfA$bio <- sapply(dfA$module, bio)
cent <- dfA %>% group_by(module) %>% summarise(x = mean(x), y = mean(y), bio = first(bio))
pA <- ggplot(dfA, aes(x, y, colour = bio)) +
  geom_point(size = 1.3, alpha = 0.6) +
  geom_text(data = cent, aes(x, y, label = bio), size = 4, fontface = "bold", colour = "black", vjust = -0.7) +
  scale_colour_manual(values = colpal, guide = "none") +
  theme_classic() + tt + ax +
  labs(x = "UMAP 1", y = "UMAP 2", title = "A. Co-expression modules (UMAP)")

# ---- B: module relationships (eigengene heatmap) ----
ord <- res$module[order(-abs(res$r_disco))]
mods <- names(eig); E <- do.call(cbind, eig); colnames(E) <- mods
Me <- cor(E)[ord, ord]; Me[upper.tri(Me)] <- NA
hm <- reshape2::melt(Me, na.rm = TRUE); colnames(hm) <- c("m1","m2","cor")
hm$m1 <- factor(hm$m1, levels = ord); hm$m2 <- factor(hm$m2, levels = rev(ord))
lab <- sapply(ord, function(m) bio(as.integer(sub("ME","",m))))
pB <- ggplot(hm, aes(m1, m2, fill = cor)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "cor") +
  scale_x_discrete(labels = lab, position = "top") + scale_y_discrete(labels = rev(lab)) +
  theme_minimal(base_size = 10) + tt +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8), axis.text.y = element_text(size = 8),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "B. Module relationships")

# ---- C: module vs DISCO ----
res$bio <- sapply(res$module, function(m) bio(as.integer(sub("ME","",m))))
res <- res[order(-abs(res$r_disco)), ]; res$module <- factor(res$module, levels = res$module)
pC <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", abs(r_disco))), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = colpal, guide = "none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_classic() + tt + ax +
  labs(x = "Module", y = "|Correlation with DISCO|", title = "C. Which modules drive DISCO")

# ---- D: important proteins ----
delta <- colMeans(X[age >= 65, , drop = FALSE]) - colMeans(X[age <= 40, , drop = FALSE])
top <- names(sort(abs(delta), decreasing = TRUE))[1:25]
dfD <- data.frame(protein = factor(top, levels = rev(top)), delta = delta[top],
                  bio = sapply(module[top], bio))
pD <- ggplot(dfD, aes(protein, delta, fill = bio)) +
  geom_col(width = 0.7, colour = "black") + coord_flip() +
  scale_fill_manual(values = colpal, name = "Module") +
  theme_classic() + tt + ax +
  labs(x = NULL, y = "Aging deviation", title = "D. Important proteins")

fig2 <- (pA | pB) / (pC | pD)
ggsave(file.path(out, "R4_fig2_modules.pdf"), fig2, width = 13, height = 9.5)
cat("wrote R4_fig2_modules.pdf\n")
