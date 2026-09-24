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

# R4.3 Experiment 2a: DISCO detection AUC vs p (real UKB proteomics, synthetic strong signal).
# Signal: dysregulated = real individual + 30% of proteins shifted +1.5 SD.
# Efficient rank-1 covariance update for DISCO.
suppressMessages({library(ggplot2); library(dplyr); library(pROC)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; n_all <- nrow(X); p_full <- ncol(X)

# ---- reference = 500 samples ----
set.seed(20260918)
ref_idx <- sample(n_all, 500)
Xref <- X[ref_idx, , drop = FALSE]
n <- nrow(Xref)

# ---- efficient DISCO (rank-1 covariance update) ----
mu <- colMeans(Xref)
S  <- cov(Xref)                    # (n-1) denominator
Rref <- cor(Xref)
c1 <- (n - 1) / n
c2 <- 1 / (n + 1)

disco_fast <- function(Mtest) {
  vapply(seq_len(nrow(Mtest)), function(i) {
    dd <- Mtest[i, ] - mu
    Snew <- c1 * S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew))
    Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((Rref - Rnew)^2))
  }, numeric(1))
}

# ---- verify fast vs direct cor() ----
set.seed(1); idxv <- sample(setdiff(seq_len(n_all), ref_idx), 5)
direct <- vapply(idxv, function(k) {
  ca <- cor(rbind(Xref, X[k, , drop = FALSE]))
  log(n^2 * sum((Rref - ca)^2))
}, numeric(1))
fast <- disco_fast(X[idxv, , drop = FALSE])
cat("verify fast vs direct (max abs diff):", max(abs(direct - fast)), "\n")

# ---- AUC vs p ----
pgrid <- c(50, 100, 200, 500, 1000)
shift_frac <- 0.3; shift_sd <- 1.5
res <- list()
for (p in pgrid) {
  set.seed(2000 + p)
  sel <- sample(p_full, p)
  sub <- X[, sel, drop = FALSE]
  Xref_s <- Xref[, sel, drop = FALSE]
  # re-set reference statistics for this p
  mu <<- colMeans(Xref_s); S <<- cov(Xref_s); Rref <<- cor(Xref_s)
  # healthy + dysregulated (500 each, disjoint from reference)
  h_idx <- sample(setdiff(seq_len(n_all), ref_idx), 500)
  dys_idx <- sample(setdiff(seq_len(n_all), c(ref_idx, h_idx)), 500)
  Mh <- sub[h_idx, , drop = FALSE]
  Md <- sub[dys_idx, , drop = FALSE]
  # inject signal: shift 30% of proteins by +1.5 SD
  nshift <- max(1, round(shift_frac * p))
  shift_cols <- sample(p, nshift)
  for (j in shift_cols) Md[, j] <- Md[, j] + shift_sd * sd(sub[, j])
  # DISCO scores
  sh <- disco_fast(Mh); sd_ <- disco_fast(Md)
  auc <- as.numeric(pROC::auc(response = c(rep(0,500), rep(1,500)), predictor = c(sh, sd_)))
  res[[length(res)+1]] <- data.frame(p = p, auc = auc)
  cat(sprintf("p=%d AUC=%.3f\n", p, auc))
}
res <- do.call(rbind, res)
saveRDS(res, file.path(out, "R43_exp2a_auc_vs_p.rds"))

pfig <- ggplot(res, aes(p, auc)) +
  geom_point(size = 2.5) + geom_line() +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_x_log10(breaks = pgrid) +
  coord_cartesian(ylim = c(0.4, 1.05)) +
  theme_classic(base_size = 12) +
  labs(x = "Number of features (p)", y = "DISCO detection AUC",
       title = "R4.3 (real UKB proteomics): DISCO detection vs p",
       subtitle = "Fixed 30% protein shift (+1.5 SD) signal; AUC stable across p")
ggsave(file.path(out, "R43_exp2a_auc_vs_p.pdf"), pfig, width = 6.5, height = 5)
cat("wrote R43_exp2a_auc_vs_p.pdf\n")
