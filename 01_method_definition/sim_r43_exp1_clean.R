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

# R4.3 Exp1 (clean): DISCO correlation error (flat) vs DM covariance condition number (exploding).
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds")); X <- d$X; n_all <- nrow(X); p_full <- ncol(X)
set.seed(20260918); ref_idx <- sample(n_all, 500)

R_true <- cor(X)   # full-sample "true" correlation (compute once)

nrep <- 5
pgrid_all <- c(50, 100, 200, 400, 500, 1000, 2000, 2920)
res <- list()
for (p in pgrid_all) {
  cerr <- numeric(nrep); cnum <- numeric(nrep)
  for (rep in seq_len(nrep)) {
    set.seed(1000 + p + rep); sel <- sample(p_full, p)
    Xr <- X[ref_idx, sel, drop = FALSE]
    R_hat <- cor(Xr); R_t <- R_true[sel, sel, drop = FALSE]
    cerr[rep] <- sqrt(mean((R_hat[lower.tri(R_hat)] - R_t[lower.tri(R_t)])^2))
    if (p < 500) { ev <- eigen(cov(Xr), symmetric = TRUE, only.values = TRUE)$values
                   cnum[rep] <- max(ev) / max(min(ev), 1e-12) } else cnum[rep] <- NA
  }
  res[[length(res)+1]] <- data.frame(p = p, cor_err = mean(cerr), cor_err_se = sd(cerr),
                                     condnum = mean(cnum, na.rm = TRUE))
}
res <- do.call(rbind, res)
print(res)
saveRDS(res, file.path(out, "R43_exp1_disco_dm.rds"))

p1 <- ggplot(res, aes(p)) +
  geom_pointrange(aes(y = cor_err, ymin = cor_err - cor_err_se, ymax = cor_err + cor_err_se),
                  colour = "#E64B35FF", size = 0.6) +
  geom_line(aes(y = cor_err), colour = "#E64B35FF", linewidth = 0.8) +
  geom_point(aes(y = log10(condnum)), colour = "#4DBBD5FF", size = 1.8) +
  geom_line(aes(y = log10(condnum)), colour = "#4DBBD5FF", linewidth = 0.8) +
  scale_x_log10(breaks = pgrid_all) +
  scale_y_continuous(name = "Correlation RMS error (DISCO)",
                     limits = c(0, 5),
                     sec.axis = sec_axis(~ 10^., name = "Covariance condition number (DM)",
                                         breaks = c(10, 1e2, 1e3, 1e4, 1e5))) +
  theme_classic(base_size = 12) +
  labs(x = "Number of features (p)", title = "R4.3: estimation stability")
ggsave(file.path(out, "R43_exp1_disco_vs_dm.pdf"), p1, width = 6.5, height = 5)
cat("wrote R43_exp1_disco_vs_dm.pdf\n")
