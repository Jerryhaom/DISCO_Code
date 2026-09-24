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

# R2.1/R2.2 correlation-metric robustness figures.
# Recomputes DISCO (Pearson / rank-Spearman / winsorized) on the 50k UKB clinical
# subset, SAVES raw scores, and plots: (a) score scatter, (b) C-index, (c) HR per SD.
# Parallelized with mclapply. No legacy files modified.

suppressMessages({
  library(survival)
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  library(tidyr)
  library(ggpointdensity)
})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scores_rds <- file.path(out, "R21_scores.rds")

if (file.exists(scores_rds)) {
  scores <- readRDS(scores_rds)
} else {
  clin  <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
  disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
  stopifnot(identical(clin$sampleID, disco$ID))
  X <- as.matrix(clin[, -1]); age <- disco$age; tm <- disco$time; st <- disco$status
  ref_idx <- which(age <= 40)
  ref <- X[ref_idx, , drop = FALSE]; nref <- nrow(ref)
  cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
  weight <- abs(cc %*% t(cc)); diag(weight) <- 0; weight <- weight / sum(weight)

  set.seed(20260914)
  idx <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)
  Xsub <- X[idx, , drop = FALSE]

  ref_cor <- cor(ref, method = "pearson")
  ncores <- parallel::detectCores()
  disco_spec <- function(M) {
    M <- as.matrix(M)
    parallel::mclapply(seq_len(nrow(M)), function(i) {
      cc1 <- cor(rbind(ref, M[i, , drop = FALSE]), method = "pearson")
      sum((ref_cor - cc1)^2 * weight, na.rm = TRUE)
    }, mc.cores = ncores)
  }

  cat("pearson...\n");  s_p <- log(unlist(disco_spec(Xsub)) * nref^2)
  cat("spearman...\n")
  Xr <- apply(X, 2, rank); Rsub <- Xr[idx, , drop = FALSE]; Rref <- Xr[ref_idx, , drop = FALSE]
  ref_cor_r <- cor(Rref, method = "pearson")
  s_s <- log(unlist(parallel::mclapply(seq_len(nrow(Rsub)), function(i) {
      cc1 <- cor(rbind(Rref, Rsub[i, , drop = FALSE]), method = "pearson")
      sum((ref_cor_r - cc1)^2 * weight, na.rm = TRUE)
    }, mc.cores = ncores)) * nref^2)
  cat("winsorized...\n")
  q <- apply(ref, 2, quantile, probs = c(0.01, 0.99), na.rm = TRUE)
  Xw <- Xsub
  for (j in seq_len(ncol(Xw))) Xw[, j] <- pmax(q[1, j], pmin(q[2, j], Xw[, j]))
  Wref <- ref; for (j in seq_len(ncol(Wref))) Wref[, j] <- pmax(q[1, j], pmin(q[2, j], Wref[, j]))
  ref_cor_w <- cor(Wref, method = "pearson")
  s_w <- log(unlist(parallel::mclapply(seq_len(nrow(Xw)), function(i) {
      cc1 <- cor(rbind(Wref, Xw[i, , drop = FALSE]), method = "pearson")
      sum((ref_cor_w - cc1)^2 * weight, na.rm = TRUE)
    }, mc.cores = ncores)) * nref^2)

  scores <- data.frame(pearson = s_p, spearman = s_s, winsorized = s_w)
  saveRDS(scores, scores_rds)
}

# ---- performance (C-index / HR per SD) ---------------------------------
perf <- lapply(scores, function(s) {
  ok <- is.finite(s)
  d <- data.frame(tm = disco$time[idx][ok], st = disco$status[idx][ok],
                  age = disco$age[idx][ok], sex = disco$sex[idx][ok], s = s[ok])
  f <- coxph(Surv(tm, st) ~ s + age + sex, data = d)
  sm <- summary(f)
  data.frame(C_index = unname(sm$concordance[1]), se = unname(sm$concordance[2]),
             HR = exp(coef(f)["s"] * sd(s[ok])))
})
perf <- do.call(rbind, perf); perf$spec <- rownames(perf)
perf$spec <- factor(perf$spec, levels = c("pearson","spearman","winsorized"),
                    labels = c("Pearson","Spearman (rank)","Winsorized 1-99%"))

# ---- figures -----------------------------------------------------------
# (a) scatter: Pearson vs Spearman, Pearson vs winsorized
sc1 <- ggplot(scores, aes(pearson, spearman)) +
  geom_pointdensity(size = 0.4) + scale_color_distiller(palette = "Spectral") +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  labs(x = "DISCO (Pearson)", y = "DISCO (Spearman)") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.5,
           label = sprintf("r = %.2f", cor(scores$pearson, scores$spearman)))
sc2 <- ggplot(scores, aes(pearson, winsorized)) +
  geom_pointdensity(size = 0.4) + scale_color_distiller(palette = "Spectral") +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  labs(x = "DISCO (Pearson)", y = "DISCO (Winsorized)") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.5,
           label = sprintf("r = %.2f", cor(scores$pearson, scores$winsorized)))

# (b) C-index bar
p_ci <- ggplot(perf, aes(spec, C_index)) +
  geom_col(fill = "#4DBBD5FF", colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = C_index - se, ymax = C_index + se), width = 0.25) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "C-index (age+sex)", title = "Mortality discrimination") +
  coord_cartesian(ylim = c(0.70, 0.76))

# (c) HR per SD bar
p_hr <- ggplot(perf, aes(spec, HR)) +
  geom_col(fill = "#E64B35FF", colour = "black", width = 0.6) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "HR per SD", title = "Mortality hazard") +
  coord_cartesian(ylim = c(1.4, 1.8))

library(patchwork)
p <- (sc1 | sc2) / (p_ci | p_hr) + plot_annotation(tag_levels = "a",
  title = "DISCO robustness to correlation metric (R2.1/R2.2)")
ggsave(file.path(out, "R21_correlation_robustness.pdf"), p, width = 8, height = 6.5)
cat("wrote R21_correlation_robustness.pdf\n"); print(perf)
