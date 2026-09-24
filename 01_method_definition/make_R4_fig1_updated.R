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

# Figure 1 (5 panels, larger fonts): A = co-expression continuous spectrum,
# B = aging deviation concentration (cumulative), C = ~3 modes carry the signal,
# D = DISCO aligned with the aging structure, E = top-20 protein loadings of modes 1/2/3.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
genes <- expl$names
X <- as.matrix(prot[, -1]); colnames(X) <- genes; X <- scale(X)
age <- d$age; sex <- d$sex; disco <- d$mc
young <- which(age <= 40); old <- which(age >= 65)

tt <- theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 10)))
ax <- theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12, colour = "black"),
            plot.margin = margin(14, 6, 16, 6))

# ---- shared decomposition ----
Rref  <- cor(X[young, , drop = FALSE])
E     <- eigen(Rref, symmetric = TRUE)
U     <- E$vectors; ev <- E$values
PC    <- X %*% U
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
delta <- delta / sqrt(sum(delta^2))
ck    <- as.numeric(t(U) %*% delta)
c2    <- ck^2
rk    <- sapply(seq_len(ncol(X)), function(k) cor(PC[, k], disco, method = "spearman"))

# ---- A: eigenvalue spectrum ----
pA <- ggplot(data.frame(k = seq_along(ev), lambda = ev), aes(k, lambda)) +
  geom_line(colour = "#4DBBD5FF") + geom_point(size = 0.6, colour = "#4DBBD5FF") +
  scale_y_log10() + theme_classic() + tt + ax +
  labs(x = "Co-expression mode", y = "Eigenvalue (log)") +
  annotate("text", x = 140, y = max(ev), size = 5, hjust = 0,
           label = sprintf("effective rank %d/%d", sum(ev > 1), length(ev)))

# ---- B: cumulative concentration (raw / age+sex-adjusted / random) ----
Xadj <- apply(X, 2, function(y) residuals(lm(y ~ age + sex))); colnames(Xadj) <- colnames(X)
conc <- function(M) { Uu <- eigen(cor(M), symmetric = TRUE)$vectors
  cc <- (t(Uu) %*% delta)^2; cc <- cc / sum(cc)
  list(cum = cumsum(cc), eff = 1/sum(cc^2)) }
raw <- conc(X); adj <- conc(Xadj)
set.seed(1); cum_null <- matrix(NA, 50, ncol(X))
Uadj <- eigen(cor(Xadj), symmetric = TRUE)$vectors
for (rep in 1:50) { dn <- sample(delta); cc <- (t(Uadj) %*% dn)^2; cc <- cc/sum(cc); cum_null[rep, ] <- cumsum(cc) }
dfB <- data.frame(k = seq_len(ncol(X)), raw = raw$cum, adj = adj$cum, null = colMeans(cum_null))
pB <- ggplot(dfB, aes(k)) +
  geom_line(aes(y = raw, colour = "Raw")) +
  geom_line(aes(y = adj, colour = "Age+sex-adjusted")) +
  geom_line(aes(y = null, colour = "Random"), linetype = "dashed") +
  scale_colour_manual(values = c("Raw" = "#E64B35FF", "Age+sex-adjusted" = "#00A087FF", "Random" = "grey50")) +
  coord_cartesian(xlim = c(0, 60)) +
  theme_classic() + tt + ax +
  theme(legend.position = c(0.8, 0.5), legend.background = element_blank(), legend.text = element_text(size = 13)) +
  labs(x = "Top-K co-expression modes", y = "Cumulative aging deviation", colour = NULL)

# ---- C: ~3 modes (bar) ----
nshow <- 15
idx   <- order(-c2)[1:nshow]
dfC   <- data.frame(k = idx, frac = 100 * c2[idx], top = rank(-c2[idx]) <= 3)
pC <- ggplot(dfC, aes(factor(k, levels = k), frac, fill = top)) +
  geom_col() +
  geom_hline(yintercept = 100/367, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "grey70"), guide = "none") +
  theme_classic() + tt + ax +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)) +
  labs(x = "Mode (by aging load)", y = "Aging deviation (%)") +
  annotate("text", x = 6, y = max(dfC$frac) * 0.95, size = 5, hjust = 0,
           label = sprintf("top-3 = %.0f%%\nrandom ~%.1f%%/mode", 100*sum(c2[1:3]), 100/367))

# ---- D: DISCO alignment ----
top3 <- order(-c2)[1:3]
dfD  <- data.frame(ck = ck, rk = rk, k = seq_along(ck), top = seq_along(ck) %in% top3)
rho  <- cor(ck, rk, method = "spearman")
pD <- ggplot(dfD, aes(ck, rk)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = top), size = 1.5) +
  scale_colour_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "#4DBBD5FF"), guide = "none") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
  ggrepel::geom_text_repel(data = dfD[dfD$top, ], aes(label = k), size = 5,
                           nudge_x = 0.02, nudge_y = 0.03) +
  theme_classic() + tt + ax +
  labs(x = "Aging loading (signed)", y = "DISCO correlation") +
  annotate("text", x = min(ck), y = max(rk), size = 5, hjust = 0, vjust = 1,
           label = sprintf("Spearman = %.2f", rho))

# ---- E: top-20 protein loadings heatmap ----
cat_map <- c(
  IL6="Inflammation", GDF15="Inflammation", ADM="Inflammation", KYNU="Inflammation",
  HGF="Inflammation", TNFRSF10B="Inflammation", WARS="Inflammation", LGMN="Inflammation",
  TNFRSF12A="Inflammation", TFF3="Inflammation", ANGPTL4="Inflammation", CTF1="Inflammation",
  TNFRSF11A="Inflammation",
  BCAN="Neuronal", NCAM1="Neuronal", DNER="Neuronal", CA14="Neuronal", SPOCK1="Neuronal",
  BOC="Neuronal", CNTN1="Neuronal", NPTXR="Neuronal", NELL2="Neuronal", SEZ6L="Neuronal",
  SLITRK6="Neuronal", RGMB="Neuronal", SCG3="Neuronal",
  CD55="Immune/complement", CD59="Immune/complement", TNFRSF9="Immune/complement",
  FLI1="Immune/complement", CX3CL1="Immune/complement", ALCAM="Immune/complement",
  NFKBIE="Immune/complement",
  CCDC80="Stromal/ECM", SMOC1="Stromal/ECM", PRSS8="Stromal/ECM", S100A4="Stromal/ECM",
  LTBP2="Stromal/ECM", ITGB6="Stromal/ECM",
  OGFR="Other", MAD1L1="Other", BAG6="Other", SORT1="Other", DNAJB1="Other",
  EPS8L2="Other", MANSC1="Other", ABHD14B="Other", CERT="Other", FOLR1="Other",
  ATOX1="Other", SULT1A1="Other", AKT3="Other", FUS="Other", FABP4="Other",
  PRDX1="Other", DCTN2="Other")

blocks <- list()
for (k in 1:3) {
  o <- order(-abs(U[, k]))[1:10]
  blocks[[k]] <- data.frame(gene = genes[o], loading = U[o, k], mode = k, stringsAsFactors = FALSE)
}
dfE <- do.call(rbind, blocks)
dfE$cat <- ifelse(dfE$gene %in% names(cat_map), cat_map[dfE$gene], "Other")
dfE$mode_f <- factor(c("M1", "M2", "M3")[dfE$mode], levels = c("M1", "M2", "M3"))
ord <- dfE %>% arrange(mode, desc(loading)) %>% pull(gene)
dfE$gene <- factor(dfE$gene, levels = unique(ord))

cat_cols <- c(Inflammation="#E64B35FF", Neuronal="#4DBBD5FF",
              "Immune/complement"="#F39B7FFF", "Stromal/ECM"="#00A087FF", Other="grey65")
gene_col <- cat_cols[as.character(dfE$cat[match(levels(dfE$gene), dfE$gene)])]
gene_lab <- paste0("<span style='color:", gene_col, "'>", levels(dfE$gene), "</span>")

p_heat <- ggplot(dfE, aes(gene, mode_f)) +
  geom_tile(aes(fill = loading), colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF",
                       midpoint = 0, limits = c(-0.15, 0.15), name = "Loading") +
  geom_point(aes(colour = cat), size = 0, stroke = 0) +
  scale_colour_manual(values = cat_cols, name = "Biology") +
  scale_x_discrete(labels = gene_lab) +
  scale_y_discrete(limits = rev(levels(dfE$mode_f))) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = ggtext::element_markdown(size = 11, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.spacing.x = unit(0.8, "cm"),
        plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
        plot.margin = margin(2, 2, 2, 2)) +
  guides(fill = guide_colorbar(barwidth = 14, barheight = 0.8),
         colour = guide_legend(override.aes = list(size = 3), nrow = 1))

pE <- p_heat

fig <- (pA | pB) / (pC | pD) / pE + plot_layout(heights = c(1.2, 1.2, 0.9))
fig <- fig + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold", vjust = 0.2))
ggsave(file.path(out, "R4_fig1_structure.pdf"), fig, width = 12, height = 12)
cat("wrote R4_fig1_structure.pdf | signed Spearman(ck, rk) =", round(rho, 3), "\n")
