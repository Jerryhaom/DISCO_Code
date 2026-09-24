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

# R4.3 Experiment 1: empirical correlation-matrix estimation error vs p, on REAL UKB proteomics.
# "True" correlation = full-sample (N=53014) matrix; "empirical" = n_ref=500 subsample.
# Frobenius/RMS off-diagonal error, normalized per element, vs feature count p.
suppressMessages({library(ggplot2); library(dplyr); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; n_all <- nrow(X); p_full <- ncol(X)
cat("full data:", n_all, "x", p_full, "\n")

# ---- true correlation (full sample) ----
R_true <- cor(X)                       # 2920 x 2920

# ---- fixed 500-sample reference ----
set.seed(20260918)
ref_idx <- sample(n_all, 500)

pgrid <- c(10, 50, 100, 200, 500, 1000, 2000, 2900)
nrep <- 5

offdiag_rms <- function(Rhat, Rtrue) {
  sqrt(mean((Rhat[lower.tri(Rhat)] - Rtrue[lower.tri(Rtrue)])^2))
}

res <- list()
for (p in pgrid) {
  errs <- numeric(nrep)
  for (rep in seq_len(nrep)) {
    set.seed(1000 + p + rep)
    sel <- sample(p_full, min(p, p_full))
    R_hat <- cor(X[ref_idx, sel, drop = FALSE])
    R_t   <- R_true[sel, sel, drop = FALSE]
    errs[rep] <- offdiag_rms(R_hat, R_t)
  }
  res[[length(res) + 1]] <- data.frame(p = p, err_mean = mean(errs), err_sd = sd(errs))
}
res <- do.call(rbind, res)
print(res)
saveRDS(res, file.path(out, "R43_real_ukb_estimation_error.rds"))

# noise floor for independent features: SE of sample correlation ~ 1/sqrt(n-1)
noise_floor <- 1 / sqrt(499)

pfig <- ggplot(res, aes(p, err_mean)) +
  geom_point(size = 2.5) + geom_line() +
  geom_errorbar(aes(ymin = err_mean - err_sd, ymax = err_mean + err_sd), width = 0.08) +
  geom_hline(yintercept = noise_floor, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 20, y = noise_floor + 0.01, hjust = 0, size = 3,
           label = sprintf("noise floor for independent features (%.3f)", noise_floor)) +
  scale_x_log10(breaks = pgrid) +
  theme_classic(base_size = 12) +
  labs(x = "Number of features (p)", y = "RMS off-diagonal correlation error",
       title = "R4.3 (real UKB proteomics): estimation error of correlation matrix",
       subtitle = "Empirical (n=500) vs full-sample (n=53014) correlation, across p")
ggsave(file.path(out, "R43_real_ukb_estimation_error.pdf"), pfig, width = 7, height = 5)
cat("wrote R43_real_ukb_estimation_error.pdf\n")
