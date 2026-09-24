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

# Does DISCO capture the main variance? Scatter: eigenvalue (x) vs cor(PC, DISCO) (y), all 367 modes.
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; disco <- d$mc

young <- which(age <= 40)
Rref <- cor(X[young, , drop = FALSE])
eig <- eigen(Rref, symmetric = TRUE)
lam <- eig$values; U <- eig$vectors

PC <- X %*% U
r_disco <- sapply(1:ncol(X), function(k) cor(PC[, k], disco, method = "spearman"))
r_age   <- sapply(1:ncol(X), function(k) cor(PC[, k], age, method = "spearman"))

df <- data.frame(k = 1:ncol(X), lambda = lam, r_disco = r_disco, r_age = r_age)
df$ab <- abs(df$r_disco)

# is |r_disco| concentrated on large-lambda modes?
cat("Spearman(lambda, |r_disco|) =", round(cor(df$lambda, df$ab, method="spearman"), 3), "\n")
cat("top-1 mode r_disco =", round(r_disco[1], 3), " | mean |r_disco| =", round(mean(df$ab), 3), "\n")

p <- ggplot(df, aes(lambda, r_disco)) +
  geom_point(colour = "#4DBBD5FF", alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(x = "Eigenvalue (log scale)", y = "Correlation with DISCO",
       title = "Does DISCO capture the main variance?")
ggsave(file.path(out, "R4_eigen_disco_scatter.pdf"), p, width = 7, height = 5)

# second: |r_disco| vs lambda, with loess
p2 <- ggplot(df, aes(lambda, abs(r_disco))) +
  geom_point(colour = "#E64B35FF", alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, colour = "grey30") +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(x = "Eigenvalue (log scale)", y = "|Correlation with DISCO|",
       title = "|DISCO correlation| vs eigenvalue")
ggsave(file.path(out, "R4_eigen_disco_scatter_abs.pdf"), p2, width = 7, height = 5)

cat("DONE\n")
