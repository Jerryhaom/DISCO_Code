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

# Heatmap of top-20 protein loadings for modes 1/2/3, annotated by biology.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
genes <- expl$names
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; young <- which(age <= 40)
R <- cor(X[young, , drop = FALSE]); E <- eigen(R, symmetric = TRUE); U <- E$vectors

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
  o <- order(-abs(U[, k]))[1:20]
  blocks[[k]] <- data.frame(gene = genes[o], loading = U[o, k], mode = k,
                            stringsAsFactors = FALSE)
}
df <- do.call(rbind, blocks)
df$cat <- ifelse(df$gene %in% names(cat_map), cat_map[df$gene], "Other")
df$block_lab <- c("M1 (global)", "M2 (stromal/immune)", "M3 (inflammation vs neuronal)")[df$mode]

# order genes: by block, then by loading within block
ord <- df %>% arrange(mode, desc(loading)) %>% pull(gene)
df$gene <- factor(df$gene, levels = rev(unique(ord)))
df$mode_f <- factor(df$block_lab, levels = c("M1 (global)", "M2 (stromal/immune)", "M3 (inflammation vs neuronal)"))

cat_cols <- c(Inflammation="#E64B35FF", Neuronal="#4DBBD5FF",
              "Immune/complement"="#F39B7FFF", "Stromal/ECM"="#00A087FF", Other="grey65")

# loading heatmap (wide: protein x mode)
w <- df %>% select(gene, mode_f, loading) %>%
  tidyr::pivot_wider(names_from = mode_f, values_from = loading)
w <- w[order(match(w$gene, levels(df$gene))), ]

p_heat <- ggplot(df, aes(mode_f, gene, fill = loading)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF",
                       midpoint = 0, limits = c(-0.15, 0.15)) +
  scale_x_discrete(position = "top") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 0, face = "bold"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(2, 2, 2, 2)) +
  labs(fill = "Loading")

# category strip
df_cat <- df %>% distinct(gene, cat)
df_cat$gene <- factor(df_cat$gene, levels = levels(df$gene))
p_cat <- ggplot(df_cat, aes("", gene, fill = cat)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = cat_cols, name = "Biology") +
  scale_x_discrete(position = "top") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), axis.text.x = element_text(size = 11),
        legend.position = "none", plot.margin = margin(2, 2, 2, 2))

fig <- p_heat + p_cat + plot_layout(widths = c(1, 0.09))
ggsave(file.path(out, "R4_mode_protein_heatmap.pdf"), fig, width = 8.5, height = 10)

# supplementary table
write.csv(df[, c("mode", "block_lab", "gene", "loading", "cat")],
          file.path(out, "R4_mode_top20_proteins.csv"), row.names = FALSE)
cat("wrote R4_mode_protein_heatmap.pdf and R4_mode_top20_proteins.csv\n")
cat("unique proteins:", length(unique(df$gene)), "\n")
print(table(df$cat))
