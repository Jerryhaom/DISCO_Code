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

# Accel (age-residualized) DISCO C-index: age-weighted vs uniform.
suppressMessages({library(survival)})

root <- disco_project_root()
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
X <- as.matrix(clin[, -1]); X <- scale(X)
age <- disco$age; time <- disco$time; status <- disco$status

cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
w_age <- abs(cc %*% t(cc)); diag(w_age) <- 0; w_age <- w_age / sum(w_age)
w_unif <- matrix(1, ncol(X), ncol(X)); diag(w_unif) <- 0; w_unif <- w_unif / sum(w_unif)

ref_idx <- which(age <= 40); nref <- length(ref_idx)
Xr <- X[ref_idx, ]; mu <- colMeans(Xr); S <- cov(Xr); Rref <- cor(Xr)
c1 <- (nref-1)/nref; c2 <- 1/(nref+1)
set.seed(20260914); analysis <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)

for (wname in c("age-weighted","uniform")) {
  w <- if (wname=="age-weighted") w_age else w_unif
  sc <- vapply(analysis, function(i){ dd<-X[i,]-mu; Sn<-c1*S+c2*(dd%*%t(dd)); sn<-sqrt(diag(Sn)); Rn<-Sn/(sn%*%t(sn)); log(nref^2*sum(w*(Rref-Rn)^2)) }, numeric(1))
  accel <- residuals(lm(sc ~ age[analysis]))
  f <- coxph(Surv(time[analysis], status[analysis]) ~ accel)
  cat(sprintf("%s Accel C-index = %.4f\n", wname, summary(f)$concordance[1]))
}
