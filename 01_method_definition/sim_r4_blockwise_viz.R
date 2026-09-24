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

# Part 4: block-wise visualization of ClinicalProt correlation matrix.
# within-clinical / within-proteomics / cross (clinical x proteomics).
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
age <- d$age

m <- match(prot$eid, clin$sampleID); keep <- which(!is.na(m))
Xc <- scale(as.matrix(clin[m[keep], -1])); Xp <- scale(as.matrix(prot[keep, -1]))
X <- cbind(Xc, Xp); nc <- ncol(Xc); np <- ncol(Xp)
young <- which(age[keep] <= 40)
R <- cor(X[young, , drop = FALSE])

# block summaries (mean |cor| off-diagonal)
ic <- 1:nc; ip <- (nc+1):(nc+np)
within_c <- mean(abs(R[ic, ic][lower.tri(R[ic, ic])]))
within_p <- mean(abs(R[ip, ip][lower.tri(R[ip, ip])]))
cross    <- mean(abs(R[ic, ip]))
cat(sprintf("mean|cor| within-clinical=%.3f  within-proteomics=%.3f  cross=%.3f\n", within_c, within_p, cross))

# block heatmap (downsample proteomics for readability)
set.seed(1); psub <- sort(sample(ip, 150))
sel <- c(ic, psub)
Rsub <- R[sel, sel]
block <- data.frame(
  idx = seq_along(sel),
  type = factor(ifelse(sel <= nc, "clinical", "proteomics"), levels = c("clinical","proteomics"))
)
# reorder: clinical first, then proteomics
ord <- order(block$type)
Rsub <- Rsub[ord, ord]; block <- block[ord, ]
block$type <- factor(block$type, levels = c("clinical","proteomics"))

# melt for heatmap
Rsub[upper.tri(Rsub)] <- NA
hm <- reshape2::melt(Rsub, na.rm = TRUE)
p1 <- ggplot(hm, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "cor") +
  theme_void(base_size = 12) +
  theme(legend.position = "right") +
  labs(title = "ClinicalProt correlation (block structure)")
ggsave(file.path(out, "R4_blockwise_heatmap.pdf"), p1, width = 6, height = 5.5)

# block summary bar
s <- data.frame(
  block = c("within-clinical","within-proteomics","cross"),
  mean_abs_cor = c(within_c, within_p, cross)
)
p2 <- ggplot(s, aes(block, mean_abs_cor, fill = block)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", mean_abs_cor)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("within-clinical"="#E64B35FF","within-proteomics"="#4DBBD5FF","cross"="#00A087FF")) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Mean |correlation|", title = "Correlation strength by block") +
  theme(legend.position = "none")
ggsave(file.path(out, "R4_blockwise_summary.pdf"), p2, width = 5.5, height = 4.5)

saveRDS(list(within_c = within_c, within_p = within_p, cross = cross), file.path(out, "R4_blockwise.rds"))
cat("DONE\n")
