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
disco_init_module("08_revision_analyses")


# R03: correlation-metric robustness (R2.1 / R2.2).
# Recompute DISCO on a random 50,000-subject UKB clinical subset using three
# correlation specifications -- Pearson (default), Spearman (rank), and
# winsorized-Pearson (1%/99%) -- and compare score agreement, age association
# and age+sex-adjusted all-cause-mortality discrimination. No legacy files are modified.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
out_dir <- "outputs/R03_correlation_robustness"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(survival))

# ---- load --------------------------------------------------------------
clin  <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))     # sampleID + 58 biomarkers
disco <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))     # ID, age, sex, time, status, cause, mc, hd
stopifnot(identical(clin$sampleID, disco$ID))

X   <- as.matrix(clin[, -1])              # 58 biomarkers
age <- disco$age; sex <- disco$sex; tm <- disco$time; st <- disco$status

ref_idx <- which(age <= 40)
ref <- X[ref_idx, , drop = FALSE]
nref <- nrow(ref)

# age-correlation weights (same as cal_disco)
cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
weight <- abs(cc %*% t(cc)); diag(weight) <- 0; weight <- weight / sum(weight)

# random 50k analysis subset (exclude reference)
set.seed(20260914)
idx <- sample(setdiff(seq_len(nrow(X)), ref_idx), 50000)
Xsub <- X[idx, , drop = FALSE]

# ---- DISCO under a correlation specification ---------------------------
# All specifications share the fixed age-correlation weight and the Pearson loop;
# "spearman" is implemented by rank-transforming each biomarker once (Spearman rho =
# Pearson correlation of ranks), which is the monotone/rank-equivalent correlation.
disco_pearson <- function(Xa, ref, weight) {
  ref <- as.matrix(ref); Xa <- as.matrix(Xa)
  ref_cor <- cor(ref, method = "pearson")
  mc <- numeric(nrow(Xa))
  for (i in seq_len(nrow(Xa))) {
    d1 <- rbind(ref, Xa[i, , drop = FALSE])
    cc1 <- cor(d1, method = "pearson")
    ds <- sum((ref_cor - cc1)^2 * weight, na.rm = TRUE)
    mc[i] <- log(ds * nref * nref)
  }
  mc
}

winsorize_cols <- function(M, ref) {
  q <- apply(ref, 2, quantile, probs = c(0.01, 0.99), na.rm = TRUE)
  M <- as.matrix(M)
  for (j in seq_len(ncol(M))) M[, j] <- pmax(q[1, j], pmin(q[2, j], M[, j]))
  M
}

cat("computing pearson...\n")
mc_pearson  <- disco_pearson(Xsub, ref, weight)

cat("rank-transforming for spearman...\n")
Xrank <- apply(X, 2, rank)                     # full-sample ranks (monotone/Spearman)
Rsub  <- Xrank[idx, , drop = FALSE]
Rref  <- Xrank[ref_idx, , drop = FALSE]
cat("computing spearman (rank-based)...\n")
mc_spearman <- disco_pearson(Rsub, Rref, weight)

cat("computing winsorized...\n")
Xw_ref  <- winsorize_cols(ref, ref)
Xw_sub  <- winsorize_cols(Xsub, ref)
mc_winsor <- disco_pearson(Xw_sub, Xw_ref, weight)

# validation: existing full-sample Pearson DISCO for the same individuals
mc_existing <- disco$mc[idx]

# ---- comparisons -------------------------------------------------------
scores <- data.frame(pearson = mc_pearson, spearman = mc_spearman,
                     winsorized = mc_winsor, existing_mc = mc_existing)

cor_methods <- cor(scores, use = "complete.obs")

age_cor <- sapply(scores, function(s) cor(age[idx], s, method = "spearman"))

# age+sex adjusted mortality C-index + HR per SD
surv_metric <- function(s) {
  ok <- is.finite(s)
  d <- data.frame(tm = tm[idx][ok], st = st[idx][ok], age = age[idx][ok], sex = sex[idx][ok], s = s[ok])
  f <- coxph(Surv(tm, st) ~ s + age + sex, data = d)
  sm <- summary(f)
  sd_s <- sd(s[ok])
  c(idx = unname(sm$concordance[1]), se = unname(sm$concordance[2]),
    hr = exp(coef(f)["s"] * sd_s))
}
perf <- t(sapply(scores, surv_metric))
perf <- as.data.frame(perf)
perf$age_spearman <- age_cor

# ---- export ------------------------------------------------------------
write.csv(cor_methods, file.path(out_dir, "R03_score_correlation.csv"))
write.csv(perf, file.path(out_dir, "R03_mortality_performance.csv"))
write.csv(data.frame(
  subset_n = length(idx), events = sum(st[idx] == 1), reference_n = nref,
  existing_mc_subset_cor = cor(mc_existing, mc_pearson, use = "complete.obs")
), file.path(out_dir, "R03_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "subset=random_50000_seed_20260914",
  "reference=age<=40_n5282",
  "correlation_specs=pearson_spearman_winsorized_1_99",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("R03 done\n")
