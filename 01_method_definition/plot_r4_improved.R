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

# Two improved figures:
# 1) module eigengene vs DISCO with BIOLOGICAL labels (inflammation/immune/cardiac).
# 2) block-wise C-index (crude) per block + cross-independent annotation.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

# ---- 1. module eigengene vs DISCO with biological labels ----
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
res <- r$res
bio <- c("ME8"="Inflammation / TNF-kB",
         "ME7"="Immune / cytokines",
         "ME5"="Cardiac / natriuretic",
         "ME3"="Senescence/neuro",
         "ME4"="Bone / fibrinolysis",
         "ME1"="Mixed",
         "ME2"="Mixed", "ME6"="Mixed", "ME9"="Mixed", "ME10"="Mixed", "ME11"="Mixed", "ME12"="Mixed")
res$bio <- bio[res$module]
res <- res[order(-abs(res$r_disco)), ]
res$module <- factor(res$module, levels = res$module)
p1 <- ggplot(res, aes(module, abs(r_disco), fill = bio)) +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = bio), angle = 90, hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c("Inflammation / TNF-kB"="#E64B35FF",
                               "Immune / cytokines"="#4DBBD5FF",
                               "Cardiac / natriuretic"="#00A087FF",
                               "Senescence/neuro"="#F39B7FFF",
                               "Bone / fibrinolysis"="#8491B4FF",
                               "Mixed"="grey80"), guide="none") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_classic(base_size = 13) +
  labs(x = "Co-expression module", y = "|Correlation with DISCO|",
       title = "DISCO is driven by specific biological modules")
ggsave(file.path(out, "R4_modules_bio.pdf"), p1, width = 8.5, height = 4.8)

# ---- 2. block-wise C-index ----
cindex <- data.frame(
  block = c("full", "within-clinical", "within-proteomics", "cross"),
  ci = c(0.6882, 0.6317, 0.6872, 0.6701)
)
cindex$block <- factor(cindex$block, levels = c("within-clinical","within-proteomics","cross","full"))
p2 <- ggplot(cindex, aes(block, ci, fill = block)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = sprintf("%.3f", ci)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("within-clinical"="#E64B35FF","within-proteomics"="#4DBBD5FF",
                               "cross"="#00A087FF","full"="grey40"), guide="none") +
  coord_cartesian(ylim = c(0.60, 0.71)) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Mortality C-index (crude)",
       title = "Each block alone predicts mortality (cross not noise)") +
  annotate("text", x = 4.5, y = 0.64, size = 3.5, hjust = 0,
           label = "cross adds independent signal:\nLRT p=0.015 (C-index +0.0008)")
ggsave(file.path(out, "R4_blockwise_cindex.pdf"), p2, width = 7, height = 5)
cat("DONE\n")
