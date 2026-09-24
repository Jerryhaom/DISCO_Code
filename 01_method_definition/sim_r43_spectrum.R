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

# Eigenvalue spectrum / effective rank of the real UKB proteomics correlation matrix.
# Supports: high correlation (low effective rank) -> estimation error stays at noise floor.
suppressMessages({library(ggplot2); library(dplyr); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; n <- nrow(X); p <- ncol(X)
R <- cor(X)
saveRDS(R, file.path(out, "R_true_full.rds"))   # reuse for Experiment 2

ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
ev <- sort(ev, decreasing = TRUE)

# effective rank: eigenvalues > 1 (variance above a single standardized variable)
eff_rank_gt1 <- sum(ev > 1)
# 95% variance PCs
cumvar <- cumsum(ev) / sum(ev)
pcs95 <- which(cumvar >= 0.95)[1]

# Marchenko-Pastur bulk edge for p/n ratio (upper edge of noise band)
q <- p / n
mp_edge <- (1 + sqrt(q))^2
n_above_mp <- sum(ev > mp_edge)

cat(sprintf("n=%d p=%d\n", n, p))
cat(sprintf("eigenvalues > 1: %d / %d = %.1f%%\n", eff_rank_gt1, p, 100*eff_rank_gt1/p))
cat(sprintf("PCs for 95%% variance: %d / %d = %.1f%%\n", pcs95, p, 100*pcs95/p))
cat(sprintf("Marchenko-Pastur bulk edge = %.3f; eigenvalues above it: %d\n", mp_edge, n_above_mp))

# scree plot
df <- data.frame(idx = seq_along(ev), ev = ev)
pfig <- ggplot(df, aes(idx, ev)) +
  geom_point(size = 0.4, colour = "#4DBBD5FF") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#E64B35FF") +
  geom_hline(yintercept = mp_edge, linetype = "dotted", colour = "grey40") +
  annotate("text", x = p*0.7, y = 1.1, colour = "#E64B35FF", size = 3.2,
           label = "eigenvalue = 1") +
  annotate("text", x = p*0.7, y = mp_edge + 0.3, colour = "grey30", size = 3.2,
           label = sprintf("MP bulk edge (%.2f)", mp_edge)) +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  labs(x = "Eigenvalue rank", y = "Eigenvalue (log scale)",
       title = "Real UKB proteomics: correlation-matrix eigenvalue spectrum",
       subtitle = sprintf("Effective rank: %d eigenvalues > 1 (%.0f%%); %d PCs explain 95%% variance",
                          eff_rank_gt1, 100*eff_rank_gt1/p, pcs95))
ggsave(file.path(out, "R43_eigenvalue_spectrum.pdf"), pfig, width = 7.5, height = 5)
cat("wrote R43_eigenvalue_spectrum.pdf\n")
