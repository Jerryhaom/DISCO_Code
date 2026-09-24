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

# R4.3: DISCO vs DM across feature count p, on REAL UKB proteomics.
# Exp1: correlation error (DISCO) vs covariance condition number (DM)
# Exp2a: detection AUC (synthetic signal), DISCO vs DM, with SE (DeLong)
# Exp2b: mortality C-index, DISCO vs DM, with SE (coxph concordance)
suppressMessages({library(ggplot2); library(dplyr); library(survival); library(pROC)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid; n_all <- nrow(X); p_full <- ncol(X)
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]; status <- mort$status[mm]; time <- mort$time[mm]

# ---- shared helpers ----
set_ref <- function(Xref) {
  list(n = nrow(Xref), mu = colMeans(Xref), S = cov(Xref), Rref = cor(Xref))
}
disco_fast <- function(Mtest, ref) {
  n <- ref$n; c1 <- (n-1)/n; c2 <- 1/(n+1)
  vapply(seq_len(nrow(Mtest)), function(i) {
    dd <- Mtest[i, ] - ref$mu
    Snew <- c1 * ref$S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((ref$Rref - Rnew)^2))
  }, numeric(1))
}
dm_fast <- function(Mtest, ref) {
  ridge <- 1e-6 * mean(diag(ref$S))
  Sinv <- solve(ref$S + ridge * diag(ncol(ref$S)))
  vapply(seq_len(nrow(Mtest)), function(i) {
    dd <- Mtest[i, ] - ref$mu
    as.numeric(t(dd) %*% Sinv %*% dd)
  }, numeric(1))
}
condnum <- function(S) {
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  max(ev) / max(min(ev), 1e-12)
}

# ===================== EXPERIMENT 1 =====================
set.seed(20260918); ref_idx <- sample(n_all, 500)
pgrid1 <- c(50, 100, 200, 400)
nrep <- 5
res1 <- list()
for (p in pgrid1) {
  cerr <- cnum <- numeric(nrep)
  for (rep in seq_len(nrep)) {
    set.seed(1000 + p + rep)
    sel <- sample(p_full, p)
    Xr <- X[ref_idx, sel, drop = FALSE]
    R_hat <- cor(Xr)
    R_true <- cor(X[, sel, drop = FALSE])
    cerr[rep] <- sqrt(mean((R_hat[lower.tri(R_hat)] - R_true[lower.tri(R_true)])^2))
    cnum[rep] <- condnum(cov(Xr))
  }
  res1[[length(res1)+1]] <- data.frame(p = p, cor_err = mean(cerr), cor_err_se = sd(cerr),
                                       condnum = mean(cnum), condnum_se = sd(cnum))
}
res1 <- do.call(rbind, res1)

# also flat cor error for larger p (DISCO-only, to show it never rises)
pgrid1b <- c(500, 1000, 2000, 2920)
for (p in pgrid1b) {
  cerr <- numeric(nrep)
  for (rep in seq_len(nrep)) {
    set.seed(1000 + p + rep)
    sel <- sample(p_full, p)
    R_hat <- cor(X[ref_idx, sel, drop = FALSE]); R_true <- cor(X[, sel, drop = FALSE])
    cerr[rep] <- sqrt(mean((R_hat[lower.tri(R_hat)] - R_true[lower.tri(R_true)])^2))
  }
  res1 <- rbind(res1, data.frame(p = p, cor_err = mean(cerr), cor_err_se = sd(cerr),
                                 condnum = NA, condnum_se = NA))
}

p1 <- ggplot(res1, aes(p)) +
  geom_line(aes(y = cor_err, colour = "DISCO (corr. error)")) +
  geom_pointrange(aes(y = cor_err, ymin = cor_err - cor_err_se, ymax = cor_err + cor_err_se,
                      colour = "DISCO (corr. error)"), size = 0.6) +
  geom_line(aes(y = condnum / 2000, colour = "DM (condition number)")) +
  geom_point(aes(y = condnum / 2000, colour = "DM (condition number)"), size = 1.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 2000, name = "Covariance condition number (DM)")) +
  scale_x_log10(breaks = c(pgrid1, pgrid1b)) +
  scale_colour_manual(values = c("DISCO (corr. error)" = "#E64B35FF",
                                 "DM (condition number)" = "#4DBBD5FF")) +
  labs(x = "Number of features (p)", y = "Correlation RMS error (DISCO)",
       colour = NULL, title = "R4.3: estimation stability (DISCO vs DM)") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out, "R43_exp1_disco_vs_dm.pdf"), p1, width = 6.5, height = 5)

# ===================== EXPERIMENT 2a (AUC) =====================
set.seed(20260918); ref_idx2 <- sample(n_all, 500)
pgrid2 <- c(50, 100, 200, 500, 1000)
res2 <- list()
for (p in pgrid2) {
  set.seed(2000 + p); sel <- sample(p_full, p)
  Xr <- X[ref_idx2, sel, drop = FALSE]; ref <- set_ref(Xr)
  h_idx <- sample(setdiff(seq_len(n_all), ref_idx2), 500)
  d_idx <- sample(setdiff(seq_len(n_all), c(ref_idx2, h_idx)), 500)
  Mh <- X[h_idx, sel, drop = FALSE]; Md <- X[d_idx, sel, drop = FALSE]
  nshift <- max(1, round(0.3 * p)); shift_cols <- sample(p, nshift)
  for (j in shift_cols) Md[, j] <- Md[, j] + 1.5 * sd(X[, sel[j]])
  resp <- c(rep(0, 500), rep(1, 500))
  # DISCO
  sc_d <- c(disco_fast(Mh, ref), disco_fast(Md, ref))
  auc_d <- as.numeric(pROC::auc(resp, sc_d))
  se_d <- as.numeric(diff(as.numeric(pROC::ci.auc(resp, sc_d, method = "delong"))[c(1,3)])) / (2*1.96)
  # DM
  sc_m <- c(dm_fast(Mh, ref), dm_fast(Md, ref))
  auc_m <- as.numeric(pROC::auc(resp, sc_m))
  se_m <- as.numeric(diff(as.numeric(pROC::ci.auc(resp, sc_m, method = "delong"))[c(1,3)])) / (2*1.96)
  res2[[length(res2)+1]] <- data.frame(p = p, metric = "DISCO", auc = auc_d, se = se_d)
  res2[[length(res2)+1]] <- data.frame(p = p, metric = "DM", auc = auc_m, se = se_m)
  cat(sprintf("2a p=%d DISCO_AUC=%.3f DM_AUC=%.3f\n", p, auc_d, auc_m))
}
res2 <- do.call(rbind, res2)
p2 <- ggplot(res2, aes(p, auc, colour = metric)) +
  geom_line() + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = auc - se, ymax = auc + se), width = 0.06) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  scale_x_log10(breaks = pgrid2) +
  coord_cartesian(ylim = c(0.4, 1.05)) +
  scale_colour_manual(values = c("DISCO" = "#E64B35FF", "DM" = "#4DBBD5FF")) +
  labs(x = "Number of features (p)", y = "Detection AUC", colour = NULL,
       title = "R4.3: detection AUC (DISCO vs DM)") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out, "R43_exp2a_auc.pdf"), p2, width = 6, height = 5)

# ===================== EXPERIMENT 2b (C-index) =====================
ref_age <- which(age <= 40 & !is.na(age))
ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 10000)
pgrid3 <- c(50, 100, 200, 500, 1000, 2000, 2920)
res3 <- list()
for (p in pgrid3) {
  set.seed(3000 + p); sel <- sample(p_full, p)
  Xr <- X[ref_age, sel, drop = FALSE]; ref <- set_ref(Xr)
  Ma <- X[analysis, sel, drop = FALSE]
  sc_d <- disco_fast(Ma, ref); sc_m <- dm_fast(Ma, ref)
  dd <- data.frame(tm = time[analysis], st = status[analysis])
  f_d <- coxph(Surv(tm, st) ~ sc_d, data = dd); ci_d <- summary(f_d)$concordance[1]; se_d <- summary(f_d)$concordance[2]
  f_m <- coxph(Surv(tm, st) ~ sc_m, data = dd); ci_m <- summary(f_m)$concordance[1]; se_m <- summary(f_m)$concordance[2]
  res3[[length(res3)+1]] <- data.frame(p = p, metric = "DISCO", ci = ci_d, se = se_d)
  res3[[length(res3)+1]] <- data.frame(p = p, metric = "DM", ci = ci_m, se = se_m)
  cat(sprintf("2b p=%d DISCO_C=%.4f DM_C=%.4f\n", p, ci_d, ci_m))
}
res3 <- do.call(rbind, res3)
p3 <- ggplot(res3, aes(p, ci, colour = metric)) +
  geom_line() + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = ci - se, ymax = ci + se), width = 0.06) +
  scale_x_log10(breaks = pgrid3) +
  scale_colour_manual(values = c("DISCO" = "#E64B35FF", "DM" = "#4DBBD5FF")) +
  labs(x = "Number of features (p)", y = "Mortality C-index", colour = NULL,
       title = "R4.3: mortality C-index (DISCO vs DM)") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out, "R43_exp2b_cindex.pdf"), p3, width = 6, height = 5)

saveRDS(list(exp1 = res1, exp2a = res2, exp2b = res3), file.path(out, "R43_disco_vs_dm.rds"))
cat("DONE\n")
