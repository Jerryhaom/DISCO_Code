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

# R4 block-wise: decompose DISCO perturbation into within-clinical, within-proteomics, cross blocks.
# ClinicalProt (58 clinical + 367 proteomics). Reference age<=40. Mortality C-index per block.
suppressMessages({library(survival); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))

m <- match(prot$eid, clin$sampleID)
keep <- which(!is.na(m))
Xc <- as.matrix(clin[m[keep], -1]); Xp <- as.matrix(prot[keep, -1])
Xc <- scale(Xc); Xp <- scale(Xp)
X <- cbind(Xc, Xp)
nc <- ncol(Xc); np <- ncol(Xp)
age <- d$age[keep]; time <- d$time[keep]; status <- d$status[keep]
cat("ClinicalProt N =", length(keep), " clinical =", nc, " proteomics =", np, "\n")

ref_age <- which(age <= 40 & !is.na(age))
cat("reference (age<=40) N =", length(ref_age), "\n")

# block masks
ic <- 1:nc; ip <- (nc+1):(nc+np)

# reference stats
Xref <- X[ref_age, , drop = FALSE]; n <- nrow(Xref)
mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)
c1 <- (n-1)/n; c2 <- 1/(n+1)

# analysis subset
ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 5000)

res <- matrix(NA, length(analysis), 4)  # S_full, S_withinC, S_withinP, S_cross
for (i in seq_along(analysis)) {
  dd <- X[analysis[i], ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Diff <- (Rref - Rnew)^2
  S_withinC <- sum(Diff[ic, ic])
  S_withinP <- sum(Diff[ip, ip])
  S_cross   <- 2 * sum(Diff[ic, ip])
  res[i, ] <- c(S_withinC + S_withinP + S_cross, S_withinC, S_withinP, S_cross)
}
colnames(res) <- c("full", "within_clinical", "within_proteomics", "cross")

# per-pair average (normalize by number of off-diagonal pairs in each block)
n_pairs <- c(
  within_clinical = nc*(nc-1)/2,
  within_proteomics = np*(np-1)/2,
  cross = nc*np
)
perp <- data.frame(
  block = names(n_pairs),
  n_pairs = as.numeric(n_pairs),
  mean_per_pair = c(mean(res[,2])/n_pairs[1], mean(res[,3])/n_pairs[2], mean(res[,4])/n_pairs[3])
)
cat("\n=== per-pair average contribution ===\n"); print(perp)

# sub-DISCO = log(n^2 * S_block); mortality C-index + age correlation
sub <- log(n^2 * res)  # log per block
dd0 <- data.frame(tm = time[analysis], st = status[analysis], age = age[analysis])
cindex <- sapply(seq_len(ncol(sub)), function(k) {
  f <- coxph(Surv(tm, st) ~ sub[,k], data = dd0)
  summary(f)$concordance[1]
})
names(cindex) <- colnames(sub)
cat("\n=== mortality C-index per block (sub-DISCO) ===\n"); print(round(cindex, 4))
agecor <- sapply(seq_len(ncol(sub)), function(k) cor(sub[,k], dd0$age, method = "spearman"))
names(agecor) <- colnames(sub)
cat("\n=== Spearman cor(sub-DISCO, age) per block ===\n"); print(round(agecor, 4))

saveRDS(list(res = res, perp = perp, cindex = cindex, agecor = agecor, n = n, nc = nc, np = np),
        file.path(out, "R4_blockwise.rds"))
cat("DONE\n")
