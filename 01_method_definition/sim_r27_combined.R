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

# Combined: crude / age+sex / accel C-index vs reference size, age-weighted vs uniform.
suppressMessages({library(survival); library(ggplot2); library(dplyr); library(tidyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
X <- as.matrix(clin[, -1]); X <- scale(X)
age <- disco$age; sex <- disco$sex; time <- disco$time; status <- disco$status

cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
w_age <- abs(cc %*% t(cc)); diag(w_age) <- 0; w_age <- w_age / sum(w_age)
w_unif <- matrix(1, ncol(X), ncol(X)); diag(w_unif) <- 0; w_unif <- w_unif / sum(w_unif)

ref_idx <- which(age <= 40)
set.seed(20260914); analysis <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)
refsizes <- c(50, 100, 200, 500, 1000, 2000, 5282)

res <- list()
for (wname in c("age-weighted", "uniform")) {
  w <- if (wname == "age-weighted") w_age else w_unif
  for (nref in refsizes) {
    set.seed(123); ref_sub <- sample(ref_idx, nref)
    Xr <- X[ref_sub, , drop = FALSE]
    mu <- colMeans(Xr); S <- cov(Xr); Rref <- cor(Xr)
    c1 <- (nref-1)/nref; c2 <- 1/(nref+1)
    sc <- vapply(analysis, function(i){ dd<-X[i,]-mu; Sn<-c1*S+c2*(dd%*%t(dd)); sn<-sqrt(diag(Sn)); Rn<-Sn/(sn%*%t(sn)); log(nref^2*sum(w*(Rref-Rn)^2)) }, numeric(1))
    accel <- residuals(lm(sc ~ age[analysis]))
    f_crude <- coxph(Surv(time[analysis],status[analysis]) ~ sc)
    f_adj   <- coxph(Surv(time[analysis],status[analysis]) ~ sc + age[analysis] + sex[analysis])
    f_accel <- coxph(Surv(time[analysis],status[analysis]) ~ accel)
    res[[length(res)+1]] <- data.frame(weight=wname, refsize=nref, model="crude", ci=summary(f_crude)$concordance[1])
    res[[length(res)+1]] <- data.frame(weight=wname, refsize=nref, model="age+sex", ci=summary(f_adj)$concordance[1])
    res[[length(res)+1]] <- data.frame(weight=wname, refsize=nref, model="accel", ci=summary(f_accel)$concordance[1])
  }
}
res <- do.call(rbind, res)
res$model <- factor(res$model, levels = c("crude","age+sex","accel"))
saveRDS(res, file.path(out, "R27_versions_full.rds"))

p <- ggplot(res, aes(refsize, ci, colour = weight)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~ model, nrow = 1) +
  scale_x_log10(breaks = refsizes) +
  scale_colour_manual(values = c("age-weighted" = "#E64B35FF", "uniform" = "#4DBBD5FF")) +
  theme_classic(base_size = 13) +
  labs(x = "Reference size", y = "Mortality C-index", colour = NULL,
       title = "DISCO C-index vs reference size (model x weight)") +
  theme(legend.position = "top")
ggsave(file.path(out, "R27_models_vs_refsize.pdf"), p, width = 10.5, height = 4.2)
cat("wrote R27_models_vs_refsize.pdf\n")
print(res %>% pivot_wider(id_cols = c(model, refsize), names_from = weight, values_from = ci))
