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

# Age+sex-adjusted co-expression structure: does the "aging deviation concentrated in few modes" hold?
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; sex <- d$sex
young <- which(age <= 40); old <- which(age >= 65)

# ---- age+sex residualize (each protein regressed on age + sex) ----
Xadj <- apply(X, 2, function(y) residuals(lm(y ~ age + sex)))
colnames(Xadj) <- colnames(X)

# ---- 1. co-expression structure: raw vs adjusted ----
eff_rank <- function(M) sum(eigen(cor(M), symmetric=TRUE, only.values=TRUE)$values > 1)
cat(sprintf("effective rank: raw = %d, age+sex-adjusted = %d (of %d)\n",
            eff_rank(X), eff_rank(Xadj), ncol(X)))

# ---- 2. aging deviation concentration: raw vs adjusted modes ----
aging_dev <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
aging_dev <- aging_dev / sqrt(sum(aging_dev^2))

conc <- function(M) {
  U <- eigen(cor(M), symmetric = TRUE)$vectors
  c2 <- (t(U) %*% aging_dev)^2; c2 <- c2 / sum(c2)
  list(eff_modes = 1 / sum(c2^2), cum = cumsum(c2), c2 = c2)
}
raw <- conc(X)
adj <- conc(Xadj)
cat(sprintf("aging-deviation effective modes: raw = %.1f, age+sex-adjusted = %.1f\n",
            raw$eff_modes, adj$eff_modes))

# null (random) reference
set.seed(1)
null_eff <- replicate(50, { dn <- sample(aging_dev); U <- eigen(cor(Xadj), symmetric=TRUE)$vectors;
  c2 <- (t(U) %*% dn)^2; c2 <- c2/sum(c2); 1/sum(c2^2) })
cat(sprintf("null (random) effective modes = %.0f\n", mean(null_eff)))

# ---- figure: cumulative concentration raw vs adjusted vs null ----
set.seed(1); cum_null <- matrix(NA, 50, ncol(X))
for (rep in 1:50) { dn <- sample(aging_dev); U <- eigen(cor(Xadj), symmetric=TRUE)$vectors;
  c2 <- (t(U)%*%dn)^2; c2 <- c2/sum(c2); cum_null[rep,] <- cumsum(c2) }
df <- data.frame(k = 1:ncol(X), raw = raw$cum, adj = adj$cum, null = colMeans(cum_null))
p <- ggplot(df, aes(k)) +
  geom_line(aes(y = raw, colour = "Raw modes")) +
  geom_line(aes(y = adj, colour = "Age+sex-adjusted modes")) +
  geom_line(aes(y = null, colour = "Random"), linetype = "dashed") +
  scale_colour_manual(values = c("Raw modes"="#E64B35FF","Age+sex-adjusted modes"="#00A087FF","Random"="grey50")) +
  coord_cartesian(xlim = c(0, 80)) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.75, 0.3), legend.background = element_blank()) +
  labs(x = "Top-K co-expression modes", y = "Cumulative aging deviation", colour = NULL,
       title = "Aging deviation concentration (raw vs age+sex-adjusted)")
ggsave(file.path(out, "R4_fig1_age_sex_adjusted.pdf"), p, width = 7, height = 5)

saveRDS(list(raw = raw, adj = adj, null_eff = null_eff), file.path(out, "R4_age_sex_adjusted.rds"))
cat("wrote R4_fig1_age_sex_adjusted.pdf\n")
