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


# Prepare the UKB proteomics feature matrix for DISCO/DM regeneration.
# Reproduces the QC + impute + LASSO 367-protein selection from Cal_Prot.R.
# No legacy files are modified; outputs go to outputs/03_prot_regeneration/.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
out_dir <- disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(survival); library(Hmisc); library(readr); library(glmnet)})

info <- read_csv(file.path(root, "UKB_Data/UKB.Mortality.csv"), show_col_types = FALSE)
load(file.path(root, "UKB_Data/20240321_prot_instance0.RData"))
prot <- prot.dat
cat("raw prot dim:", dim(prot), "\n")

ns <- apply(prot, 1, function(x) length(x[is.na(x)]))
nc <- apply(prot, 2, function(x) length(x[is.na(x)]))
prot1 <- prot[ns < 0.5 * ncol(prot), nc < 0.1 * nrow(prot)]
cat("after QC dim:", dim(prot1), "\n")

prot1 <- apply(prot1, 2, Hmisc::impute)
prot1 <- data.frame(prot1)
info1 <- info[match(prot1$eid, info$id), ]
info1 <- info1[which(info1$time > 0), ]
prot1 <- prot1[match(info1$id, prot1$eid), ]
cat("after mortality align dim:", dim(prot1), "\n")

# LASSO 367 proteins (lambda.min)
fit_cv <- readRDS(file.path(root, "Result2/Lasso.model.prot.rds"))
coef_min <- coef(fit_cv, s = "lambda.min")
bb <- rownames(coef_min)[coef_min[, 1] != 0]
cat("selected proteins (lambda.min):", length(bb), "\n")

prot2 <- data.frame(prot1[, c("eid", bb)])
stopifnot(nrow(prot2) == nrow(prot1))

meta <- data.frame(
  ID = info1$id, age = info1$Age, sex = info1$Gender,
  time = info1$time, status = info1$status, cause = info1$cause
)
stopifnot(nrow(meta) == nrow(prot2))

saveRDS(prot2, file.path(out_dir, "prot_367_selected.rds"))
saveRDS(meta, file.path(out_dir, "meta.rds"))
saveRDS(bb, file.path(out_dir, "selected_proteins.rds"))

cat("reference age<=40 N:", sum(meta$age <= 40), "\n")
cat("analysis N:", nrow(meta), "\n")
cat("events:", sum(meta$status == 1), "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_prepare.txt"))
cat("PREPARE DONE\n")
