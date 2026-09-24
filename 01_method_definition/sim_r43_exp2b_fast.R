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

# R4.3 Exp2b (optimized): DISCO vs DM mortality C-index vs p.
suppressMessages({library(ggplot2); library(dplyr); library(survival)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid; n_all <- nrow(X); p_full <- ncol(X)
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]; status <- mort$status[mm]; time <- mort$time[mm]

ref_age <- which(age <= 40 & !is.na(age))
ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 3000)

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
  D <- sweep(Ma, 2, ref$mu, "-")
  rowSums(D * (D %*% Sinv))
}

pgrid <- c(50, 100, 200, 500, 1000, 2000, 2920)
res <- list()
for (p in pgrid) {
  set.seed(3000 + p); sel <- sample(p_full, p)
  Xr <- X[ref_age, sel, drop = FALSE]
  ref <- list(n = nrow(Xr), mu = colMeans(Xr), S = cov(Xr), Rref = cor(Xr))
  Ma <- X[analysis, sel, drop = FALSE]
  sc_d <- disco_vec(Ma, ref); sc_m <- dm_vec(Ma, ref)
  dd <- data.frame(tm = time[analysis], st = status[analysis])
  f_d <- coxph(Surv(tm, st) ~ sc_d, data = dd); ci_d <- summary(f_d)$concordance[1]; se_d <- summary(f_d)$concordance[2]
  f_m <- coxph(Surv(tm, st) ~ sc_m, data = dd); ci_m <- summary(f_m)$concordance[1]; se_m <- summary(f_m)$concordance[2]
  res[[length(res)+1]] <- data.frame(p = p, metric = "DISCO", ci = ci_d, se = se_d)
  res[[length(res)+1]] <- data.frame(p = p, metric = "DM", ci = ci_m, se = se_m)
  cat(sprintf("2b p=%d DISCO_C=%.4f DM_C=%.4f\n", p, ci_d, ci_m))
}
res <- do.call(rbind, res)
saveRDS(res, file.path(out, "R43_exp2b_disco_dm.rds"))

p3 <- ggplot(res, aes(p, ci, colour = metric)) +
  geom_line() + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = ci - se, ymax = ci + se), width = 0.06) +
  scale_x_log10(breaks = pgrid) +
  scale_colour_manual(values = c("DISCO" = "#E64B35FF", "DM" = "#4DBBD5FF")) +
  labs(x = "Number of features (p)", y = "Mortality C-index", colour = NULL,
       title = "R4.3: mortality C-index (DISCO vs DM)") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out, "R43_exp2b_cindex.pdf"), p3, width = 6, height = 5)
cat("DONE\n")
