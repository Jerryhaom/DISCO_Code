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

# Visualize: (1) module -> biological process (Reactome + GO, lenient);
#            (2) aging deviation per protein (which proteins are hit).
suppressMessages({library(dplyr); library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
r <- readRDS(file.path(out, "R4_modules_annotation.rds"))
cls <- readRDS(file.path(root, "Result3/Organ/UKB.Prot.class.rds"))
reactome <- readRDS(file.path(root, "UKB_Data/reactome_pathways.rds"))

module <- r$module; delta <- r$delta; Rref <- r$Rref
p <- ncol(Rref)
genes <- names(delta)

# ---- Reactome enrichment per module (hypergeometric, correct universe) ----
universe <- genes
mod_react <- sapply(sort(unique(module)), function(k) {
  g <- genes[module == k]
  pv <- sapply(reactome, function(pw) {
    ov <- sum(g %in% pw$proteins)
    m  <- sum(pw$proteins %in% universe)   # pathway genes present in the 367
    if (ov == 0 || m == 0) 1 else phyper(ov - 1, m, length(universe) - m, length(g), lower.tail = FALSE)
  })
  pv <- p.adjust(pv, "BH")
  ord <- order(pv)
  top <- ord[1:2]
  paste(sapply(top, function(i) sprintf("%s (P=%.1e)", reactome[[i]]$name, pv[i])), collapse=" | ")
})
names(mod_react) <- paste0("M", sort(unique(module)))

# ---- GO BP per module (lenient) ----
g2e <- setNames(cls$ENTREZID, cls$Gene)
entrez <- unname(g2e[genes])
mod_go <- sapply(sort(unique(module)), function(k) {
  e <- entrez[module == k]; e <- e[!is.na(e)]
  ego <- tryCatch(enrichGO(gene = e, universe = entrez[!is.na(entrez)],
                           OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                           pvalueCutoff = 0.05, qvalueCutoff = 0.5, minGSSize = 3, readable = TRUE),
                  error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return("")
  paste(head(as.data.frame(ego)$Description, 2), collapse = " | ")
})
names(mod_go) <- paste0("M", sort(unique(module)))

# module aging loading
mod_load <- sapply(sort(unique(module)), function(k) mean(abs(delta[module == k])))
names(mod_load) <- paste0("M", sort(unique(module)))

modtab <- data.frame(
  module = names(mod_load), n = as.numeric(table(module)),
  aging_loading = round(mod_load, 3),
  Reactome = unname(mod_react), GO = unname(mod_go)
)
cat("=== module -> biological process ===\n"); print(modtab, row.names = FALSE)

saveRDS(list(modtab = modtab, module = module, delta = delta),
        file.path(out, "R4_modules_table.rds"))

# ---- figure: aging deviation waterfall (top proteins labeled) ----
ord <- order(delta, decreasing = TRUE)
df <- data.frame(rank = seq_along(ord), delta = delta[ord],
                 protein = genes[ord], module = factor(module[ord]))
top_idx <- c(1:15, (nrow(df)-14):nrow(df))   # top up + top down
df$lab <- ifelse(df$rank %in% top_idx, df$protein, "")

modcols <- scales::hue_pal()(length(unique(module)))
p1 <- ggplot(df, aes(rank, delta, fill = module)) +
  geom_col(width = 1) +
  scale_fill_manual(values = modcols, guide = "none") +
  geom_text(data = df[df$lab != "", ], aes(label = lab), angle = 90, hjust = -0.1, size = 2.6) +
  theme_classic(base_size = 14) +
  labs(x = "Proteins (ranked by aging deviation)", y = "Aging deviation (old - young)",
       title = "Which proteins aging dysregulation hits") +
  coord_cartesian(ylim = c(-1.5, 1.5))
ggsave(file.path(out, "R4_aging_deviation_waterfall.pdf"), p1, width = 11, height = 5)
cat("wrote figures\n")
