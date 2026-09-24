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


# R02: organ-DISCO versus raw protein/panel predictors (R1.3 / R4.1).
# For each of the 18 locked organ panels, compare the organ-DISCO score against
# (a) the single best raw protein, (b) the panel mean and (c) the panel first
# principal component, all in age+sex-adjusted Cox models for all-cause mortality.
# No legacy files are modified.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
organ_dir <- file.path(root, "Result3/Organ")
out_dir <- "outputs/R02_raw_biomarker_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(survival))

# ---- load --------------------------------------------------------------
odis  <- readRDS(file.path(organ_dir, "UKB.odis.rds"))          # 36 organ DISCO + outcomes
prot  <- readRDS(file.path(organ_dir, "UKB.prot.imputed.rds"))  # eid + imputed raw proteins
cls   <- readRDS(file.path(organ_dir, "UKB.Prot.class.rds"))    # protein -> specific.tissue
sel18 <- readRDS(file.path(organ_dir, "Tissue.selected.rds"))   # 18 locked organs (spaces)

stopifnot(nrow(odis) == nrow(prot))
odis$eid <- prot[, "eid"]                                       # row-aligned (Cal_Organ_spe.R)
stopifnot(identical(odis$eid, prot[, "eid"]))

organ_dot <- gsub(" ", ".", sel18)                             # odis uses dot names
stopifnot(all(organ_dot %in% colnames(odis)))

age <- odis$age; sex <- odis$sex; tm <- odis$time; st <- odis$status
P <- prot[, -1]                                                 # 1553 proteins, no missing
pnames <- colnames(P)

# ---- helpers -----------------------------------------------------------
# age+sex-adjusted Cox for one continuous predictor x; returns c(index, se, hr_per_sd)
fit_one_pred <- function(x, label) {
  ok <- is.finite(x) & !is.na(x)
  if (sum(ok) < 100 || length(unique(x[ok])) < 2) return(c(NA, NA, NA))
  dd <- data.frame(tm = tm[ok], st = st[ok], age = age[ok], sex = sex[ok], x = x[ok])
  f <- tryCatch(coxph(Surv(tm, st) ~ x + age + sex, data = dd), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA, NA))
  s <- summary(f)
  sd_x <- sd(x[ok], na.rm = TRUE)
  c(idx = unname(s$concordance[1]), se = unname(s$concordance[2]),
    hr = exp(coef(f)["x"] * sd_x))
}

# panel first PC score (prcomp on the scaled panel); returns scores or NULL
panel_pc1 <- function(M) {
  if (ncol(M) < 2) return(M[, 1])
  fit <- tryCatch(prcomp(M, center = TRUE, scale. = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  fit$x[, 1]
}

# ---- loop over organs --------------------------------------------------
res <- list()
for (i in seq_along(sel18)) {
  organ_space <- sel18[i]
  organ_col   <- organ_dot[i]

  panel <- unique(cls$Gene[cls$specific.tissue == organ_space])
  panel <- panel[panel %in% pnames]
  if (length(panel) == 0) next
  M <- P[, panel, drop = FALSE]

  # organ-DISCO
  d_disco <- odis[[organ_col]]
  disco <- fit_one_pred(d_disco, "disco")

  # panel mean (of z-scored proteins)
  Z <- scale(M)
  mean_score <- rowMeans(Z)
  mean_p <- fit_one_pred(mean_score, "mean")

  # panel PC1
  pc1 <- panel_pc1(M)
  pc1_p <- if (is.null(pc1)) c(NA, NA, NA) else fit_one_pred(pc1, "pc1")

  # best single protein (max age+sex C-index; selection is optimistic)
  prot_c <- sapply(panel, function(p) fit_one_pred(M[, p], p)["idx"])
  best_c <- max(prot_c, na.rm = TRUE)
  best_p <- panel[which.max(prot_c)]
  best_hr <- fit_one_pred(M[, best_p], best_p)["hr"]

  res[[length(res) + 1L]] <- data.frame(
    organ = organ_space,
    n_proteins = length(panel),
    disco_c = disco["idx"], disco_se = disco["se"], disco_hr_per_sd = disco["hr"],
    best_protein_c = best_c, best_protein = best_p, best_protein_hr_per_sd = best_hr,
    panel_mean_c = mean_p["idx"], panel_mean_hr_per_sd = mean_p["hr"],
    panel_pc1_c = pc1_p["idx"], panel_pc1_hr_per_sd = pc1_p["hr"],
    disco_minus_best_protein = disco["idx"] - best_c,
    disco_minus_panel_mean = disco["idx"] - mean_p["idx"],
    disco_minus_panel_pc1 = disco["idx"] - pc1_p["idx"],
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, res)

# ---- summary -----------------------------------------------------------
summ <- data.frame(
  n_organs = nrow(res),
  disco_beats_best_protein = sum(res$disco_minus_best_protein > 0, na.rm = TRUE),
  disco_beats_panel_mean   = sum(res$disco_minus_panel_mean > 0, na.rm = TRUE),
  disco_beats_panel_pc1    = sum(res$disco_minus_panel_pc1 > 0, na.rm = TRUE),
  median_disco_minus_best  = median(res$disco_minus_best_protein, na.rm = TRUE),
  median_disco_minus_mean  = median(res$disco_minus_panel_mean, na.rm = TRUE),
  median_disco_minus_pc1   = median(res$disco_minus_panel_pc1, na.rm = TRUE),
  median_panel_n           = median(res$n_proteins),
  stringsAsFactors = FALSE
)

# ---- export ------------------------------------------------------------
write.csv(res, file.path(out_dir, "R02_organ_disco_vs_raw_protein.csv"), row.names = FALSE)
write.csv(summ, file.path(out_dir, "R02_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "panel_source=HPA_tissue_specific_proteins",
  "raw_protein_matrix=UKB.prot.imputed.rds_1553_proteins_no_missing",
  "organ_disco=UKB.odis.rds_age_sex_adjusted",
  "best_protein_is_selection_optimistic=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("R02 done:", nrow(res), "organs\n")
