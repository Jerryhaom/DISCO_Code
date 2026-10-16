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
disco_init_module("08_sensitivity_analyses")


# R2.4: nonlinear-dependence simulation.
# Does DISCO, which uses Pearson (linear) correlation, still flag "dysregulated"
# individuals whose biomarker dependence is nonlinear? Compare DISCO-Pearson,
# DISCO-rank (Spearman-type) and Mahalanobis distance (DM) across scenarios whose
# marginal mean/variance are held equal so DM cannot win by a trivial mean shift.

options(stringsAsFactors = FALSE)
out_dir <- "outputs/R03_correlation_robustness"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260914)
n_ref <- 1000
n_test <- 1000

# --- reference & healthy: bivariate normal, rho = 0.5 -------------------
rho <- 0.5
S <- matrix(c(1, rho, rho, 1), 2, 2)
ref <- MASS::mvrnorm(n_ref, c(0, 0), S)
healthy <- MASS::mvrnorm(n_test, c(0, 0), S)

# --- dysregulated scenarios (same marginal mean 0, var 1) ---------------
# linear: weakened correlation
linear <- MASS::mvrnorm(n_test, c(0, 0), matrix(c(1, 0.3, 0.3, 1), 2, 2))
# quadratic (nonlinear monotone): X2 = a*(X1^2 - 1) + noise, Var(X2)=1
a <- 0.6
X1_q <- rnorm(n_test)
quadratic <- cbind(X1_q, a * (X1_q^2 - 1) + rnorm(n_test, 0, sqrt(1 - 2*a^2)))
# funnel / heteroscedastic: X2 = X1 * e, Var(X2)=1
X1_f <- rnorm(n_test)
funnel <- cbind(X1_f, X1_f * rnorm(n_test))

# --- metrics ------------------------------------------------------------
# uniform weights (no age in simulation): mean squared off-diagonal correlation change
disco_pearson <- function(x) {
  rc <- cor(ref)
  ca <- cor(rbind(ref, x))
  log(n_ref^2 * mean((rc - ca)^2))
}
# rank-based DISCO: Pearson on global ranks (monotone/Spearman-type)
rank_all <- function(M) apply(M, 2, rank)
R_all <- rank_all(rbind(ref, healthy, linear, quadratic, funnel))
R_ref <- R_all[1:n_ref, ]
offset <- n_ref
rank_disco <- function(x, Rmat, k) {
  rc <- cor(R_ref)
  ca <- cor(rbind(R_ref, Rmat[k, , drop = FALSE]))
  log(n_ref^2 * mean((rc - ca)^2))
}
mahalanobis <- function(x) {
  mu <- colMeans(ref); Sig <- cov(ref)
  as.numeric(stats::mahalanobis(matrix(x, nrow = 1), mu, Sig))
}

metrics <- function(M, label, Rmat = NULL) {
  if (is.null(Rmat)) {
    d <- sapply(seq_len(nrow(M)), function(i) disco_pearson(M[i, , drop = FALSE]))
  } else {
    d <- sapply(seq_len(nrow(M)), function(i) rank_disco(M[i,,drop=FALSE], Rmat, i))
  }
  dm <- sapply(seq_len(nrow(M)), function(i) mahalanobis(M[i, , drop = FALSE]))
  data.frame(scenario = label,
             disco_pearson = if (is.null(Rmat)) d else disco_pearson_vec(M),
             disco_rank = d,
             dm = dm)
}
# pearson vector (compute once per scenario)
disco_pearson_vec <- function(M) sapply(seq_len(nrow(M)), function(i) disco_pearson(M[i,,drop=FALSE]))

# healthy (need rank offset)
get_rank_disco <- function(M, Rmat, start) {
  sapply(seq_len(nrow(M)), function(i) {
    rc <- cor(R_ref); ca <- cor(rbind(R_ref, Rmat[start + i, , drop=FALSE]))
    log(n_ref^2 * mean((rc - ca)^2))
  })
}

res <- data.frame(
  scenario = c("healthy","linear","quadratic","funnel"),
  disco_pearson = NA_real_, disco_rank = NA_real_, dm = NA_real_
)
res$disco_pearson[1] <- median(disco_pearson_vec(healthy))
res$disco_pearson[2] <- median(disco_pearson_vec(linear))
res$disco_pearson[3] <- median(disco_pearson_vec(quadratic))
res$disco_pearson[4] <- median(disco_pearson_vec(funnel))
res$disco_rank[1] <- median(get_rank_disco(healthy, R_all, n_ref))
res$disco_rank[2] <- median(get_rank_disco(linear, R_all, n_ref + n_test))
res$disco_rank[3] <- median(get_rank_disco(quadratic, R_all, n_ref + 2*n_test))
res$disco_rank[4] <- median(get_rank_disco(funnel, R_all, n_ref + 3*n_test))
res$dm[1] <- median(sapply(seq_len(nrow(healthy)), function(i) mahalanobis(healthy[i,,drop=FALSE])))
res$dm[2] <- median(sapply(seq_len(nrow(linear)), function(i) mahalanobis(linear[i,,drop=FALSE])))
res$dm[3] <- median(sapply(seq_len(nrow(quadratic)), function(i) mahalanobis(quadratic[i,,drop=FALSE])))
res$dm[4] <- median(sapply(seq_len(nrow(funnel)), function(i) mahalanobis(funnel[i,,drop=FALSE])))

write.csv(res, file.path(out_dir, "R24_nonlinear_simulation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_R24.txt"))
cat("R2.4 done\n")
res_num <- res; res_num$disco_pearson <- round(res_num$disco_pearson, 3)
res_num$disco_rank <- round(res_num$disco_rank, 3); res_num$dm <- round(res_num$dm, 3)
print(res_num)
