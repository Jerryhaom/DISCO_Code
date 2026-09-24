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
disco_init_module("03_multimodal_UKB")


# Regenerate UKB proteomics DISCO (mc) and DM (hd) from scratch.
# reference unified to age<=40 (consistent with clinical/clinical+prot/metabolomics).
#   - DISCO = log(sum(w * (cor_ref - cor_new)^2) * n_ref^2)
#   - DM = sqrt(Mahalanobis squared distance)  (standard Mahalanobis distance)
# No legacy files are modified.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
od <- disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration")
r2 <- file.path(root, "Result2")
dir.create(od, recursive = TRUE, showWarnings = FALSE)

prot <- readRDS(file.path(od, "prot_367_selected.rds"))
meta <- readRDS(file.path(od, "meta.rds"))
X <- as.matrix(prot[, -1])
age <- meta$age
n_total <- nrow(X)

ref_idx <- which(age <= 40)
ref <- X[ref_idx, , drop = FALSE]
n0 <- nrow(ref)
cat("reference age<=40 N =", n0, "\n")

mu <- colMeans(ref)
S <- cov(ref)
ref_cor <- cor(ref)
cc <- apply(X, 2, function(x) cor(age, x, use = "complete.obs"))
weight <- abs(cc %*% t(cc)); diag(weight) <- 0; weight <- weight / sum(weight)

# --- DISCO (exact incremental update) ---
disco_incremental <- function(Xa, mu, S, n0, ref_cor, weight) {
  n <- nrow(Xa)
  out <- numeric(n)
  for (k in seq_len(n)) {
    d <- Xa[k, ] - mu
    sum_new <- (n0 - 1) * S + (n0 / (n0 + 1)) * (d %*% t(d))
    S_new <- sum_new / n0
    sd_new <- sqrt(diag(S_new))
    cor_new <- S_new / (sd_new %*% t(sd_new))
    ds <- sum((ref_cor - cor_new)^2 * weight)
    out[k] <- log(ds * n0 * n0)
  }
  out
}

cat("computing DISCO for all", n_total, "individuals...\n")
mc <- disco_incremental(X, mu, S, n0, ref_cor, weight)

# --- DM = standard Mahalanobis distance (no scaling) ---
m2 <- mahalanobis(X, mu, S)
hd <- sqrt(m2)

new_obj <- data.frame(
  ID = meta$ID, age = meta$age, sex = meta$sex,
  time = meta$time, status = meta$status, cause = meta$cause,
  mc = mc, hd = hd, stringsAsFactors = FALSE
)
saveRDS(new_obj, file.path(od, "UKB.Prot.DISCO.regenerated.rds"))
saveRDS(mc, file.path(od, "DISCO.UKB.prot.regenerated.rds"))

# --- compare with legacy ---
old <- readRDS(file.path(r2, "UKB.Prot.DISCO.rds"))
stopifnot(identical(old$ID, new_obj$ID))
cmp <- data.frame(
  metric = c("DISCO (mc)", "DM (hd)"),
  pearson_cor = c(cor(mc, old$mc), cor(hd, old$hd)),
  max_abs_diff = c(max(abs(mc - old$mc)), max(abs(hd - old$hd))),
  stringsAsFactors = FALSE
)
write.csv(cmp, file.path(od, "regeneration_comparison.csv"), row.names = FALSE)

summary_tbl <- data.frame(
  item = c("analysis_N", "events", "reference_age_cutoff", "reference_N",
           "n_features", "DM_scale_constant"),
  value = c(n_total, sum(meta$status == 1), "<=40", n0, ncol(X), "1.0 (standard)")
)
write.csv(summary_tbl, file.path(od, "regeneration_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo_regenerate.txt"))
writeLines(c(
  "status=PASS_UNIFIED_REFERENCE_AGE40",
  "disco_reference=age<=40",
  "dm_reference=age<=40",
  "reference_unified_with_clinical_and_other_modalities=true",
  "legacy_files_modified=false"
), file.path(od, "regeneration_status.txt"))

cat("\n=== comparison vs legacy ===\n")
print(cmp)
cat("\n=== summary ===\n")
print(summary_tbl)
cat("\nDONE\n")
