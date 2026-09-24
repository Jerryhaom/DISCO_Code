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

# R2.7 + weight sensitivity: reference size -> C-index, across 4 specs:
#   weight (age-weighted / uniform) x model (crude / age+sex).
suppressMessages({library(survival); library(ggplot2); library(dplyr); library(tidyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
X <- as.matrix(clin[, -1]); X <- scale(X)
age <- disco$age; sex <- disco$sex; time <- disco$time; status <- disco$status

# weights
cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
w_age <- abs(cc %*% t(cc)); diag(w_age) <- 0; w_age <- w_age / sum(w_age)
w_unif <- matrix(1, ncol(X), ncol(X)); diag(w_unif) <- 0; w_unif <- w_unif / sum(w_unif)

ref_idx <- which(age <= 40)
set.seed(20260914)
analysis <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)

refsizes <- c(50, 100, 200, 500, 1000, 2000, 5282)
res <- list()
for (wname in c("age-weighted", "uniform")) {
  w <- if (wname == "age-weighted") w_age else w_unif
  for (nref in refsizes) {
    set.seed(123); ref_sub <- sample(ref_idx, nref)
    Xr <- X[ref_sub, , drop = FALSE]
    mu <- colMeans(Xr); S <- cov(Xr); Rref <- cor(Xr)
    c1 <- (nref - 1) / nref; c2 <- 1 / (nref + 1)
    sc <- vapply(analysis, function(i) {
      dd <- X[i, ] - mu
      Snew <- c1 * S + c2 * (dd %*% t(dd))
      snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
      log(nref^2 * sum(w * (Rref - Rnew)^2))
    }, numeric(1))
    f_crude <- coxph(Surv(time[analysis], status[analysis]) ~ sc)
    f_adj   <- coxph(Surv(time[analysis], status[analysis]) ~ sc + age[analysis] + sex[analysis])
    res[[length(res)+1]] <- data.frame(weight = wname, refsize = nref, model = "crude", c_index = summary(f_crude)$concordance[1])
    res[[length(res)+1]] <- data.frame(weight = wname, refsize = nref, model = "age+sex", c_index = summary(f_adj)$concordance[1])
  }
}
res <- do.call(rbind, res)
res$spec <- paste(res$weight, res$model, sep=", ")
print(res %>% pivot_wider(names_from = spec, values_from = c_index))
saveRDS(res, file.path(out, "R27_versions.rds"))

p <- ggplot(res, aes(refsize, c_index, colour = weight, linetype = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_x_log10(breaks = refsizes) +
  scale_colour_manual(values = c("age-weighted" = "#E64B35FF", "uniform" = "#4DBBD5FF")) +
  theme_classic(base_size = 14) +
  labs(x = "Reference size", y = "Mortality C-index", colour = NULL, linetype = NULL,
       title = "DISCO C-index vs reference size (weight x model)")
ggsave(file.path(out, "R27_refsize_versions.pdf"), p, width = 7, height = 5)
cat("wrote R27_refsize_versions.pdf\n")
