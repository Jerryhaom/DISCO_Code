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

# R2.7 (final): reference size -> crude (unadjusted) mortality C-index.
suppressMessages({library(survival); library(ggplot2)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
X <- as.matrix(clin[, -1]); X <- scale(X)
age <- disco$age; time <- disco$time; status <- disco$status

cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
weight <- abs(cc %*% t(cc)); diag(weight) <- 0; weight <- weight / sum(weight)

ref_idx <- which(age <= 40)
set.seed(20260914)
analysis <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)

refsizes <- c(50, 100, 200, 500, 1000, 2000, 5282)
res <- data.frame(refsize = refsizes, c_index = NA)
for (k in seq_along(refsizes)) {
  nref <- refsizes[k]
  set.seed(123); ref_sub <- sample(ref_idx, nref)
  Xr <- X[ref_sub, , drop = FALSE]
  mu <- colMeans(Xr); S <- cov(Xr); Rref <- cor(Xr)
  c1 <- (nref - 1) / nref; c2 <- 1 / (nref + 1)
  sc <- vapply(analysis, function(i) {
    dd <- X[i, ] - mu
    Snew <- c1 * S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(nref^2 * sum(weight * (Rref - Rnew)^2))
  }, numeric(1))
  fit <- coxph(Surv(time[analysis], status[analysis]) ~ sc)
  res$c_index[k] <- summary(fit)$concordance[1]
  cat(sprintf("refsize=%d crude C-index=%.4f\n", nref, res$c_index[k]))
}
saveRDS(res, file.path(out, "R27_crude_cindex.rds"))

p <- ggplot(res, aes(refsize, c_index)) +
  geom_pointrange(aes(ymin = c_index - 0.003, ymax = c_index + 0.003), size = 0.7) +
  geom_line() +
  scale_x_log10(breaks = res$refsize) +
  coord_cartesian(ylim = c(0.64, 0.70)) +
  theme_classic(base_size = 14) +
  labs(x = "Reference size", y = "Crude mortality C-index",
       title = "DISCO C-index vs reference size")
ggsave(file.path(out, "R27_refsize_cindex.pdf"), p, width = 6.5, height = 4.5)
cat("wrote R27_refsize_cindex.pdf\n")
