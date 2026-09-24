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

# Clean Figure 2 (3 panels):
#  2A: dimension-reduction scatter of proteins (MDS), clusters colored by module
#  2B: module vs DISCO (bar)
#  2C: important proteins per module (top by aging deviation, colored by module)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
module <- r$module; res <- r$res

X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco <- d$mc

bio_map <- c("ME8"="Inflammation","ME7"="Immune","ME5"="Cardiac","ME3"="Senescence","ME4"="Bone")
bio <- function(m) ifelse(paste0("ME",m) %in% names(bio_map), unname(bio_map[paste0("ME",m)]), "Other")
colpal <- c("Inflammation"="#E64B35FF","Immune"="#4DBBD5FF","Cardiac"="#00A087FF",
            "Senescence"="#F39B7FFF","Bone"="#8491B4FF","Other"="grey80")

tt <- theme(plot.title = element_text(size = 14, face = "bold"))
ax <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11, colour = "black"))

# ---- 2A: MDS scatter (clusters) ----
R <- cor(X)
D <- 1 - abs(R); diag(D) <- 0
mds <- cmdscale(D, k = 2)
dfA <- data.frame(x = mds[,1], y = mds[,2], module = factor(module), protein = colnames(X))
dfA$bio <- sapply(dfA$module, bio)
# module centroids for labels
cent <- dfA %>% group_by(module) %>% summarise(x = mean(x), y = mean(y), bio = first(bio))
pA <- ggplot(dfA, aes(x, y, colour = bio)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text(data = cent, aes(x, y, label = bio), size = 4, fontface = "bold", colour = "black",
            vjust = -0.7) +
  scale_colour_manual(values = colpal, guide = "none") +
  theme_classic() + tt + ax +
  labs(x = "MDS dimension 1", y = "MDS dimension 2", title = "A. Co-expression modules")

# ---- 2B: module vs DISCO ----
res$bio <- sapply(res$module, function(m) bio(as.integer(sub("ME","",m))))
res <- res[order(-abs(res$r_disco)), ]; res$module <- factor(res$module, levels = res$module)
pB <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", abs(r_disco))), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = colpal, guide = "none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_classic() + tt + ax +
  labs(x = "Module", y = "|Correlation with DISCO|", title = "B. Which modules drive DISCO")

# ---- 2C: important proteins per module ----
delta <- colMeans(X[age >= 65, , drop = FALSE]) - colMeans(X[age <= 40, , drop = FALSE])
top <- names(sort(abs(delta), decreasing = TRUE))[1:25]
dfC <- data.frame(protein = factor(top, levels = rev(top)), delta = delta[top],
                  module = factor(module[top]), bio = sapply(module[top], bio))
pC <- ggplot(dfC, aes(protein, delta, fill = bio)) +
  geom_col(width = 0.7, colour = "black") +
  coord_flip() +
  scale_fill_manual(values = colpal, name = "Module") +
  theme_classic() + tt + ax +
  labs(x = NULL, y = "Aging deviation (old - young)", title = "C. Important proteins")

fig2 <- (pA | pB) / (pC + plot_spacer())
ggsave(file.path(out, "R4_fig2_modules.pdf"), fig2, width = 11, height = 9)
cat("wrote R4_fig2_modules.pdf\n")
