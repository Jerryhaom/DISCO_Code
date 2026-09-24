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

# Symmetric distance-metric correlation heatmap (Reds palette), placed in ED4.
suppressMessages({library(ggplot2)})

root <- disco_project_root()
C <- readRDS(disco_results_path("01_method_definition/outputs/sensitivity/R4_distance_corr.rds"))
C <- (C + t(C)) / 2   # 显式对称化
out <- disco_results_path("01_method_definition/reproduction/ED4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ord <- rownames(C)
hm <- expand.grid(x = ord, y = ord)
hm$r <- as.numeric(C)
hm$x <- factor(hm$x, levels = rev(ord))
hm$y <- factor(hm$y, levels = rev(ord))

ph <- ggplot(hm, aes(x, y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5, alpha = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4, colour = "black") +
  scale_fill_distiller(palette = "Blues", direction = 1,
                       limits = c(0.8, 1), name = "Spearman\ncorrelation") +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        legend.position = "right") +
  ggtitle("Distance-metric correlation")

ggsave(file.path(out, "R4_distance_corr_heatmap.pdf"), ph, width = 8, height = 6.5)
cat("wrote", file.path(out, "R4_distance_corr_heatmap.pdf"), "\n")
