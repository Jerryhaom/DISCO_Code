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

# Unified, optimized R4 figures (consistent theme). 3 core figures:
#  Fig1: structure is in the difference (modular reference + structured dysregulation)
#  Fig2: DISCO driven by specific biological modules
#  Fig3: block-wise (cross-modality carries independent signal)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco <- d$mc

thm <- theme_classic(base_size = 13) +
  theme(axis.title = element_text(size = 13), plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(colour = "black"))

# ================= Fig 1: structure is in the difference =================
young <- which(age <= 40); old <- which(age >= 65)
Rref <- cor(X[young, , drop = FALSE])
ev <- eigen(Rref, symmetric = TRUE)$values
# (a) reference modularity
df_ev <- data.frame(k = 1:20, lambda = ev[1:20])
pa <- ggplot(df_ev, aes(k, lambda)) + geom_col(fill = "#4DBBD5FF", colour = "black", width = 0.6) +
  geom_text(aes(label = round(lambda)), vjust = -0.3, size = 3) +
  thm + labs(x = "Co-expression mode", y = "Eigenvalue", title = "Reference is modular (low-rank)") +
  annotate("text", x = 14, y = max(ev[1:20])*0.9, size = 3.5, hjust = 0,
           label = sprintf("effective rank: %d / %d", sum(ev > 1), length(ev)))
# (b) dysregulation structured
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
delta <- delta / sqrt(sum(delta^2))
U <- eigen(Rref, symmetric = TRUE)$vectors
c2 <- (t(U) %*% delta)^2; c2 <- c2 / sum(c2)
cum_real <- cumsum(c2)
set.seed(1); cum_null <- matrix(NA, 200, length(ev))
for (rep in 1:200) { dn <- sample(delta); c2n <- (t(U) %*% dn)^2; c2n <- c2n/sum(c2n); cum_null[rep,] <- cumsum(c2n) }
df_c <- data.frame(k = 1:length(ev), real = cum_real, null = colMeans(cum_null))
pb <- ggplot(df_c, aes(k)) +
  geom_line(aes(y = real, colour = "Aging deviation")) +
  geom_line(aes(y = null, colour = "Random (shuffled)")) +
  scale_colour_manual(values = c("Aging deviation"="#E64B35FF", "Random (shuffled)"="#4DBBD5FF")) +
  coord_cartesian(xlim = c(0, 60)) +
  thm + labs(x = "Top-K co-expression modes", y = "Cumulative deviation explained", colour = NULL,
             title = "Dysregulation is structured, not random") +
  theme(legend.position = c(0.75, 0.3), legend.background = element_blank()) +
  annotate("text", x = 30, y = 0.4, size = 3.5,
           label = sprintf("~%.0f modes vs %.0f random", 1/sum(c2^2), length(ev)))
fig1 <- pa | pb
ggsave(file.path(out, "R4_fig1_structure.pdf"), fig1, width = 9.5, height = 4.3)

# ================= Fig 2: module biology =================
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
res <- r$res
bio <- c("ME8"="Inflammation / TNF-kB", "ME7"="Immune / cytokines", "ME5"="Cardiac / natriuretic",
         "ME3"="Senescence / neuro", "ME4"="Bone / fibrinolysis", "ME1"="Mixed")
res$bio <- bio[res$module]; res$bio[is.na(res$bio)] <- "Mixed"
res <- res[order(-abs(res$r_disco)), ]
res$module <- factor(res$module, levels = res$module)
colpal <- c("Inflammation / TNF-kB"="#E64B35FF", "Immune / cytokines"="#4DBBD5FF",
            "Cardiac / natriuretic"="#00A087FF", "Senescence / neuro"="#F39B7FFF",
            "Bone / fibrinolysis"="#8491B4FF", "Mixed"="grey80")
fig2 <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = bio), angle = 90, hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = colpal, guide = "none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  thm + labs(x = "Co-expression module", y = "|Correlation with DISCO|",
             title = "DISCO is driven by specific biological modules")
ggsave(file.path(out, "R4_fig2_modules.pdf"), fig2, width = 8.5, height = 4.8)

# ================= Fig 3: block-wise =================
cindex <- data.frame(block = c("within-clinical","within-proteomics","cross","full"),
                     ci = c(0.6317, 0.6872, 0.6701, 0.6882))
cindex$block <- factor(cindex$block, levels = c("within-clinical","within-proteomics","cross","full"))
fig3 <- ggplot(cindex, aes(block, ci, fill = block)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = sprintf("%.3f", ci)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("within-clinical"="#E64B35FF","within-proteomics"="#4DBBD5FF",
                               "cross"="#00A087FF","full"="grey40"), guide="none") +
  coord_cartesian(ylim = c(0.60, 0.71)) +
  thm + labs(x = NULL, y = "Mortality C-index", title = "Cross-modality block carries signal") +
  annotate("text", x = 4.5, y = 0.635, size = 3.5, hjust = 0,
           label = "cross independent: LRT p=0.015")
ggsave(file.path(out, "R4_fig3_blockwise.pdf"), fig3, width = 6.5, height = 4.8)
cat("DONE\n")
