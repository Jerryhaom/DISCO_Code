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

# Dose-response: gradually increase nonlinearity strength (alpha 0->1), track DISCO.
# Two nonlinearity types: monotone (cubic) and non-monotone (quadratic U).
# Real UKB pair FLI1-AKT3 (reference r=0.81).
suppressMessages({library(MASS); library(ggplot2); library(dplyr); library(tidyr); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
exp  <- readRDS(file.path(out, "proteomics_explore.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X    <- as.matrix(prot[, -1]); colnames(X) <- exp$names
age  <- d$age; ref_idx <- exp$ref_idx

P1 <- "FLI1"; P2 <- "AKT3"
x_all <- scale(X[, P1])
ref <- cbind(scale(X[ref_idx, P1]), scale(X[ref_idx, P2]))
nref <- nrow(ref); Rc <- cor(ref)
set.seed(20260918); N <- 2000

disco_p <- function(M) {  # against reference, Pearson
  vapply(seq_len(nrow(M)), function(i) {
    ca <- cor(rbind(ref, M[i, , drop = FALSE]))
    log(sum((Rc - ca)^2) * nref^2)
  }, numeric(1))
}
dm_g <- function(M) {
  mu <- colMeans(ref); Sig <- cov(ref)
  vapply(seq_len(nrow(M)), function(i) as.numeric(mahalanobis(M[i,,drop=FALSE], mu, Sig)), numeric(1))
}

# global rank transform over ref + all test points
alphas <- seq(0, 1, 0.1)
gen <- function(a, type) {
  xs <- sample(x_all, N, replace = TRUE)
  if (type == "monotone") {
    y <- (1 - a) * 0.81 * xs + a * (xs + 0.6 * xs^3) + 0.3 * rnorm(N)
  } else {
    y <- (1 - a) * 0.81 * xs + a * (xs^2 - 1) + 0.3 * rnorm(N)
  }
  cbind(xs, scale(y))
}
allM <- lapply(c("monotone","nonmonotone"), function(ty) lapply(alphas, function(a) gen(a, ty)))
all_stack <- do.call(rbind, unlist(allM, recursive = FALSE))
Rall <- apply(rbind(ref, all_stack), 2, rank)
Rref <- Rall[seq_len(nref), , drop = FALSE]
Rc_r <- cor(Rref)
disco_r <- function(M, start) {
  vapply(seq_len(nrow(M)), function(i) {
    ca <- cor(rbind(Rref, Rall[start + i, , drop = FALSE]))
    log(sum((Rc_r - ca)^2) * nref^2)
  }, numeric(1))
}

res <- list()
offset <- nref
for (ty in c("monotone","nonmonotone")) {
  for (k in seq_along(alphas)) {
    M <- allM[[if (ty == "monotone") 1 else 2]][[k]]
    res[[length(res)+1]] <- data.frame(
      type = ty, alpha = alphas[k],
      disco_pearson = median(disco_p(M)),
      disco_rank    = median(disco_r(M, offset)),
      dm            = median(dm_g(M)),
      pearson_r     = cor(M)[1,2], spearman_r = cor(M, method="spearman")[1,2])
    offset <- offset + N
  }
}
res <- do.call(rbind, res)
print(res)
saveRDS(res, file.path(out, "dose_response_results.rds"))
cat("saved dose_response_results.rds\n")
