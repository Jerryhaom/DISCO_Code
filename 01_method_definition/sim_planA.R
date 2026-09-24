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

# Plan A: nonlinear-dependence simulation using REAL UKB protein pair (FLI1-AKT3, r=0.81).
# Same scenario framework as the synthetic R24, but X uses the real FLI1 marginal.
# For each scenario compute per-individual DISCO(Pearson), DISCO(rank), DM; report median.

suppressMessages({library(MASS); library(ggplot2); library(dplyr); library(tidyr); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
exp  <- readRDS(file.path(out, "proteomics_explore.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X    <- as.matrix(prot[, -1]); colnames(X) <- exp$names
age  <- d$age; ref_idx <- exp$ref_idx

# ---- pick pair: FLI1 (X) and AKT3 (Y) ----
P1 <- "FLI1"; P2 <- "AKT3"
x_all <- X[, P1]; y_all <- X[, P2]
xr <- scale(x_all[ref_idx]); yr <- scale(y_all[ref_idx])
ref <- cbind(xr, yr)                 # reference N=565
cat(sprintf("reference cor(%s,%s) = %.3f\n", P1, P2, cor(xr, yr)))

nref <- nrow(ref); set.seed(20260918); N <- 1000

# scenario generators (X = real FLI1 marginal, standardized)
make_scen <- function() {
  xb <- scale(x_all)                  # full FLI1 (real marginal)
  xs <- sample(xb, N, replace = TRUE)
  list(
    healthy   = cbind(xs, (0.81 * xs + sqrt(1 - 0.81^2) * rnorm(N))),
    linear    = cbind(xs, (0.30 * xs + sqrt(1 - 0.30^2) * rnorm(N))),
    cubic     = { y <- xs + 0.6 * xs^3; cbind(xs, scale(y)) },
    sigmoid   = { y <- tanh(1.5 * xs) + 0.3 * rnorm(N); cbind(xs, scale(y)) },
    quadratic = { y <- xs^2 - 1; cbind(xs, scale(y)) },
    funnel    = { y <- xs * rnorm(N); cbind(xs, scale(y)) }
  )
}
scen <- make_scen()

# ---- DISCO per-individual score (uniform weight, 2x2 corr) ----
disco_one <- function(x, Mref, Rcor) {
  ca <- cor(rbind(Mref, x))
  ds <- sum((Rcor - ca)^2)
  log(ds * nrow(Mref)^2)
}
disco_group <- function(M) {
  Rc <- cor(M)
  vapply(seq_len(nrow(M)), function(i) disco_one(M[i, , drop = FALSE], M, Rc), numeric(1))
}
# NOTE: DISCO is computed against the REFERENCE correlation, not self-reference.
disco_against_ref <- function(M) {
  Rc <- cor(ref)
  vapply(seq_len(nrow(M)), function(i) {
    ca <- cor(rbind(ref, M[i, , drop = FALSE]))
    log(sum((Rc - ca)^2) * nref^2)
  }, numeric(1))
}
dm_group <- function(M) {
  mu <- colMeans(ref); Sig <- cov(ref)
  vapply(seq_len(nrow(M)), function(i) as.numeric(mahalanobis(M[i, , drop = FALSE], mu, Sig)), numeric(1))
}

# rank-DISCO: rank transform (globally) then same computation
Rall <- apply(rbind(ref, do.call(rbind, scen)), 2, rank)
Rref <- Rall[seq_len(nref), , drop = FALSE]
rank_disco_against_ref <- function(M, start) {
  Rc <- cor(Rref)
  vapply(seq_len(nrow(M)), function(i) {
    ca <- cor(rbind(Rref, Rall[start + i, , drop = FALSE]))
    log(sum((Rc - ca)^2) * nref^2)
  }, numeric(1))
}

# ---- compute ----
res <- data.frame(scenario = names(scen))
for (nm in names(scen)) {
  M <- scen[[nm]]
  res$disco_pearson[res$scenario == nm] <- median(disco_against_ref(M))
  res$dm[res$scenario == nm] <- median(dm_group(M))
}
# rank disco per scenario (need global rank offset)
offset <- nref
for (nm in names(scen)) {
  M <- scen[[nm]]
  res$disco_rank[res$scenario == nm] <- median(rank_disco_against_ref(M, offset))
  offset <- offset + N
}
res$pearson_r <- vapply(scen, function(M) cor(M)[1,2], numeric(1))
res$spearman_r <- vapply(scen, function(M) cor(M, method="spearman")[1,2], numeric(1))

print(res)
saveRDS(res, file.path(out, "planA_results.rds"))

# ---- figure: deviation from healthy ----
h <- res[res$scenario == "healthy", ]
dev <- data.frame(
  scenario = factor(res$scenario, levels = c("healthy","linear","cubic","sigmoid","quadratic","funnel"),
                    labels = c("Healthy","Linear (weakened)","Cubic (monotone)","Sigmoid (monotone)",
                               "Quadratic (U-shape)","Funnel (heterosced.)")),
  `DISCO (Pearson)` = res$disco_pearson - h$disco_pearson,
  `DISCO (rank)`    = res$disco_rank - h$disco_rank,
  `DM`              = res$dm - h$dm, check.names = FALSE)
devl <- dev %>% pivot_longer(-scenario, names_to = "metric", values_to = "dev")
devl$metric <- factor(devl$metric, levels = c("DISCO (Pearson)","DISCO (rank)","DM"))

pA <- ggplot(devl, aes(scenario, dev, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey30") +
  scale_fill_manual(values = c("DISCO (Pearson)" = "#E64B35FF",
                               "DISCO (rank)" = "#4DBBD5FF", "DM" = "#00A087FF")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, colour = "black"),
        legend.position = "top", legend.title = element_blank()) +
  labs(x = NULL, y = "Deviation from healthy (median score)",
       title = paste0("R2.4 nonlinear-dependence detection — real UKB pair (", P1, "-", P2, ")"),
       subtitle = "Positive = detected as dysregulated; negative = missed")
ggsave(file.path(out, "R24_planA_protein_pair.pdf"), pA, width = 8, height = 5)
cat("wrote planA figure\n")
