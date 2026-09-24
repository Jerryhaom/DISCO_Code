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

# Revised 4-panel R4.3 figure:
#  a: correlation error (DISCO)
#  b: score stability (two refs)
#  c: 5-year mortality AUC (DISCO vs DM)
#  d: mortality C-index (DISCO vs DM)
suppressMessages({library(ggplot2); library(dplyr); library(survival); library(pROC); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid; n_all <- nrow(X); p_full <- ncol(X)
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]; status <- mort$status[mm]; time <- mort$time[mm]

disco_vec <- function(Ma, ref) {
  n <- ref$n; c1 <- (n-1)/n; c2 <- 1/(n+1)
  vapply(seq_len(nrow(Ma)), function(i) {
    dd <- Ma[i, ] - ref$mu
    Snew <- c1 * ref$S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((ref$Rref - Rnew)^2))
  }, numeric(1))
}
dm_vec <- function(Ma, ref) {
  ridge <- 1e-6 * mean(diag(ref$S))
  Sinv <- solve(ref$S + ridge * diag(ncol(ref$S)))
  D <- sweep(Ma, 2, ref$mu, "-"); rowSums(D * (D %*% Sinv))
}
set_ref <- function(Xr) list(n = nrow(Xr), mu = colMeans(Xr), S = cov(Xr), Rref = cor(Xr))

ref_age <- which(age <= 40 & !is.na(age))
ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 10000)

y5 <- status[analysis] == 1 & time[analysis] <= 5     # 5-year death
lab <- as.integer(y5)

PGRID <- c(50, 200, 500, 1000, 2000, 2920)
resc <- resd <- list()
for (p in PGRID) {
  set.seed(3000 + p); sel <- sample(p_full, p)
  ref <- set_ref(X[ref_age, sel, drop = FALSE])
  Ma <- X[analysis, sel, drop = FALSE]
  sc_d <- disco_vec(Ma, ref); sc_m <- dm_vec(Ma, ref)
  # panel c: 5-year mortality AUC
  auc_d <- as.numeric(pROC::auc(lab, sc_d)); se_d <- diff(as.numeric(pROC::ci.auc(lab, sc_d, method="delong"))[c(1,3)])/(2*1.96)
  auc_m <- as.numeric(pROC::auc(lab, sc_m)); se_m <- diff(as.numeric(pROC::ci.auc(lab, sc_m, method="delong"))[c(1,3)])/(2*1.96)
  resc[[length(resc)+1]] <- data.frame(p = p, metric = "DISCO", val = auc_d, se = se_d)
  resc[[length(resc)+1]] <- data.frame(p = p, metric = "DM", val = auc_m, se = se_m)
  # panel d: mortality C-index
  dd <- data.frame(tm = time[analysis], st = status[analysis])
  fd <- coxph(Surv(tm, st) ~ sc_d, data = dd); fm <- coxph(Surv(tm, st) ~ sc_m, data = dd)
  resd[[length(resd)+1]] <- data.frame(p = p, metric = "DISCO", val = summary(fd)$concordance[1], se = summary(fd)$concordance[2])
  resd[[length(resd)+1]] <- data.frame(p = p, metric = "DM", val = summary(fm)$concordance[1], se = summary(fm)$concordance[2])
  cat(sprintf("p=%d AUC DISCO=%.3f DM=%.3f | C DISCO=%.3f DM=%.3f\n", p, auc_d, auc_m, summary(fd)$concordance[1], summary(fm)$concordance[1]))
}
resc <- do.call(rbind, resc); resc$p <- factor(resc$p, levels = PGRID)
resd <- do.call(rbind, resd); resd$p <- factor(resd$p, levels = PGRID)
saveRDS(list(auc = resc, cindex = resd), file.path(out, "R43_mortality_disco_dm.rds"))

# load a (correlation error) and b (stability)
ra <- readRDS(file.path(out, "R43_disco_only.rds"))
sb <- readRDS(file.path(out, "R43_score_stability.rds"))

thm_a <- theme_classic(base_size = 15) + theme(legend.position = "none",
  axis.title = element_text(size = 16), plot.title = element_text(size = 18))
thm_l <- theme_classic(base_size = 15) +
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.title = element_blank(), legend.key.size = unit(0.7, "cm"),
        legend.text = element_text(size = 14), axis.title = element_text(size = 16),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
        legend.margin = margin(3, 3, 3, 3), plot.title = element_text(size = 18))

pa <- ggplot(ra$exp1, aes(p, val)) +
  geom_pointrange(aes(ymin = val - se, ymax = val + se), size = 0.6) + geom_line(group = 1) +
  thm_a + labs(x = "Number of features (p)", y = "Correlation RMS error", title = "Estimation error")

pb <- ggplot(sb, aes(p, rho, colour = ref, group = ref)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se), width = 0.08) +
  coord_cartesian(ylim = c(0.5, 1.02)) +
  scale_colour_manual(values = c("DISCO(50)" = "#E64B35FF", "DISCO(full)" = "#4DBBD5FF")) +
  thm_l + labs(x = "Number of features (p)", y = "Spearman correlation", title = "Score stability")

pc <- ggplot(resc, aes(p, val, colour = metric, group = metric)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = val - se, ymax = val + se), width = 0.08) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  scale_colour_manual(values = c("DISCO" = "#E64B35FF", "DM" = "#4DBBD5FF")) +
  thm_l + labs(x = "Number of features (p)", y = "5-year mortality AUC", title = "5-year mortality")

pd <- ggplot(resd, aes(p, val, colour = metric, group = metric)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = val - se, ymax = val + se), width = 0.08) +
  coord_cartesian(ylim = c(0.5, 0.63)) +
  scale_colour_manual(values = c("DISCO" = "#E64B35FF", "DM" = "#4DBBD5FF")) +
  thm_l + labs(x = "Number of features (p)", y = "Mortality C-index", title = "Mortality prediction")

p <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a") +
  plot_layout(widths = c(1, 1), heights = c(1, 1))
ggsave(file.path(out, "R43_sensitivity_combined.pdf"), p, width = 9.5, height = 8)
cat("wrote R43_sensitivity_combined.pdf\n")
