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

# R4.3: DISCO-only, real UKB proteomics, p = {50,200,500,1000,2000,2920}, equal x-spacing.
# Exp1: correlation RMS error. Exp2a: detection AUC. Exp2b: mortality C-index.
suppressMessages({library(ggplot2); library(dplyr); library(survival); library(pROC)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid; n_all <- nrow(X); p_full <- ncol(X)
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]; status <- mort$status[mm]; time <- mort$time[mm]

PGRID <- c(50, 200, 500, 1000, 2000, 2920)
PGRIDF <- factor(PGRID, levels = PGRID)

disco_vec <- function(Ma, ref) {
  n <- ref$n; c1 <- (n-1)/n; c2 <- 1/(n+1)
  vapply(seq_len(nrow(Ma)), function(i) {
    dd <- Ma[i, ] - ref$mu
    Snew <- c1 * ref$S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((ref$Rref - Rnew)^2))
  }, numeric(1))
}

# ---------- Exp1: correlation RMS error ----------
set.seed(20260918); ref_idx <- sample(n_all, 500)
R_true <- cor(X)
res1 <- list()
for (p in PGRID) {
  cerr <- numeric(5)
  for (rep in 1:5) {
    set.seed(1000 + p + rep); sel <- sample(p_full, p)
    R_hat <- cor(X[ref_idx, sel, drop = FALSE]); R_t <- R_true[sel, sel, drop = FALSE]
    cerr[rep] <- sqrt(mean((R_hat[lower.tri(R_hat)] - R_t[lower.tri(R_t)])^2))
  }
  res1[[length(res1)+1]] <- data.frame(p = p, val = mean(cerr), se = sd(cerr))
}
res1 <- do.call(rbind, res1); res1$p <- factor(res1$p, levels = PGRID)
p1 <- ggplot(res1, aes(p, val)) +
  geom_pointrange(aes(ymin = val - se, ymax = val + se), size = 0.7) + geom_line(group = 1) +
  theme_classic(base_size = 13) +
  labs(x = "Number of features (p)", y = "Correlation RMS error",
       title = "Estimation error of correlation matrix")
ggsave(file.path(out, "R43_exp1_disco.pdf"), p1, width = 6, height = 4.5)

# ---------- Exp2a: detection AUC ----------
set.seed(20260918); ref_idx2 <- sample(n_all, 500)
res2 <- list()
for (p in PGRID) {
  set.seed(2000 + p); sel <- sample(p_full, p)
  Xr <- X[ref_idx2, sel, drop = FALSE]
  ref <- list(n = nrow(Xr), mu = colMeans(Xr), S = cov(Xr), Rref = cor(Xr))
  hi <- sample(setdiff(seq_len(n_all), ref_idx2), 500)
  di <- sample(setdiff(seq_len(n_all), c(ref_idx2, hi)), 500)
  Mh <- X[hi, sel, drop = FALSE]; Md <- X[di, sel, drop = FALSE]
  ns <- max(1, round(0.3 * p)); sc <- sample(p, ns)
  for (j in sc) Md[, j] <- Md[, j] + 1.5 * sd(X[, sel[j]])
  resp <- c(rep(0, 500), rep(1, 500))
  sc_ <- c(disco_vec(Mh, ref), disco_vec(Md, ref))
  auc <- as.numeric(pROC::auc(resp, sc_))
  se <- diff(as.numeric(pROC::ci.auc(resp, sc_, method = "delong"))[c(1,3)]) / (2*1.96)
  res2[[length(res2)+1]] <- data.frame(p = p, val = auc, se = se)
  cat(sprintf("2a p=%d AUC=%.3f\n", p, auc))
}
res2 <- do.call(rbind, res2); res2$p <- factor(res2$p, levels = PGRID)
p2 <- ggplot(res2, aes(p, val)) +
  geom_pointrange(aes(ymin = val - se, ymax = val + se), size = 0.7) + geom_line(group = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  coord_cartesian(ylim = c(0.4, 1.05)) +
  theme_classic(base_size = 13) +
  labs(x = "Number of features (p)", y = "Detection AUC", title = "DISCO detection vs p")
ggsave(file.path(out, "R43_exp2a_disco.pdf"), p2, width = 6, height = 4.5)

# ---------- Exp2b: mortality C-index ----------
ref_age <- which(age <= 40 & !is.na(age))
ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 3000)
res3 <- list()
for (p in PGRID) {
  set.seed(3000 + p); sel <- sample(p_full, p)
  Xr <- X[ref_age, sel, drop = FALSE]
  ref <- list(n = nrow(Xr), mu = colMeans(Xr), S = cov(Xr), Rref = cor(Xr))
  Ma <- X[analysis, sel, drop = FALSE]
  sc_ <- disco_vec(Ma, ref)
  dd <- data.frame(tm = time[analysis], st = status[analysis])
  f <- coxph(Surv(tm, st) ~ sc_, data = dd)
  res3[[length(res3)+1]] <- data.frame(p = p, val = summary(f)$concordance[1], se = summary(f)$concordance[2])
  cat(sprintf("2b p=%d C=%.4f\n", p, summary(f)$concordance[1]))
}
res3 <- do.call(rbind, res3); res3$p <- factor(res3$p, levels = PGRID)
p3 <- ggplot(res3, aes(p, val)) +
  geom_pointrange(aes(ymin = val - se, ymax = val + se), size = 0.7) + geom_line(group = 1) +
  theme_classic(base_size = 13) +
  labs(x = "Number of features (p)", y = "Mortality C-index", title = "DISCO mortality C-index vs p")
ggsave(file.path(out, "R43_exp2b_disco.pdf"), p3, width = 6, height = 4.5)

saveRDS(list(exp1 = res1, exp2a = res2, exp2b = res3), file.path(out, "R43_disco_only.rds"))
cat("DONE\n")
