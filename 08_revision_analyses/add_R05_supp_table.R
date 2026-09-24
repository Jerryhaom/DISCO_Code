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

# Add Supplementary Table 17 (PRS + genetic interaction) to Nature Aging Supplementary tables.xlsx.
suppressPackageStartupMessages({library(openxlsx)})

root <- disco_project_root()
src <- disco_results_path("08_revision_analyses/outputs/R05_genetic_interaction_rebuilt")
template <- file.path(root, "Manuscript/NA_submission/Supplementary tables.xlsx")
tgt <- disco_results_path("08_revision_analyses/outputs/R05_genetic_interaction_rebuilt/Supplementary tables.xlsx")
dir.create(dirname(tgt), recursive=TRUE, showWarnings=FALSE)

d <- read.csv(file.path(src, "R05_genetic_interaction.csv"), stringsAsFactors = FALSE)

lbl <- c(
  amd = "Age-related macular degeneration", ad = "Alzheimer's disease",
  asthma = "Asthma", af = "Atrial fibrillation", bipolar = "Bipolar disorder",
  bowel_ca = "Bowel cancer", breast_ca = "Breast cancer", cvd = "Cardiovascular disease",
  coeliac = "Coeliac disease", cad = "Coronary artery disease", crohn = "Crohn's disease",
  ovarian_ca = "Ovarian cancer", hypertension = "Hypertension",
  ischaemic_stroke = "Ischaemic stroke", melanoma = "Melanoma", ms = "Multiple sclerosis",
  osteoporosis = "Osteoporosis", parkinson = "Parkinson's disease",
  glaucoma = "Glaucoma", prostate_ca = "Prostate cancer", psoriasis = "Psoriasis",
  ra = "Rheumatoid arthritis", schizophrenia = "Schizophrenia", sle = "Systemic lupus",
  t1d = "Type 1 diabetes", t2d = "Type 2 diabetes",
  ulcerative_colitis = "Ulcerative colitis", vte = "Venous thromboembolism")

sys <- c(
  af = "Cardiovascular", cad = "Cardiovascular", cvd = "Cardiovascular",
  hypertension = "Cardiovascular", ischaemic_stroke = "Cardiovascular", vte = "Cardiovascular",
  t1d = "Metabolic", t2d = "Metabolic",
  ad = "Neurological", bipolar = "Neurological", ms = "Neurological",
  parkinson = "Neurological", schizophrenia = "Neurological",
  bowel_ca = "Cancer", breast_ca = "Cancer", melanoma = "Cancer",
  ovarian_ca = "Cancer", prostate_ca = "Cancer",
  coeliac = "Autoimmune", crohn = "Autoimmune", psoriasis = "Autoimmune",
  ra = "Autoimmune", sle = "Autoimmune", ulcerative_colitis = "Autoimmune",
  asthma = "Respiratory", amd = "Ophthalmic", glaucoma = "Ophthalmic",
  osteoporosis = "Musculoskeletal")

sys_order <- c("Cardiovascular", "Metabolic", "Neurological", "Cancer",
               "Autoimmune", "Respiratory", "Ophthalmic", "Musculoskeletal")

star <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
       ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p == 0) return("<1e-300")
  if (p < 0.001) return(sprintf("%.1e", p))
  sprintf("%.3f", p)
}
ci <- function(lo, hi) sprintf("%.2f\u2013%.2f", lo, hi)

d$System <- factor(sys[d$disease], levels = sys_order)
d$Disease <- lbl[d$disease]
d$Cohort <- ifelse(d$version == "clinical", "Clinical DISCO", "Proteomics DISCO")
d <- d[order(d$System, d$disease, d$version), ]

tbl <- data.frame(
  Disease = d$Disease,
  Cohort = d$Cohort,
  N = d$n,
  Events = d$events,
  `DISCO HR (95% CI)` = paste0(sprintf("%.2f", d$disco_hr), star(d$disco_p),
                                " (", ci(d$disco_lo, d$disco_hi), ")"),
  `DISCO P` = sapply(d$disco_p, fmt_p),
  `PRS HR (95% CI)` = paste0(sprintf("%.2f", d$prs_hr), star(d$prs_p),
                              " (", ci(d$prs_lo, d$prs_hi), ")"),
  `PRS P` = sapply(d$prs_p, fmt_p),
  `Interaction HR (95% CI)` = paste0(sprintf("%.3f", d$interaction_hr), star(d$interaction_p),
                                      " (", ci(d$interaction_lo, d$interaction_hi), ")"),
  `Interaction P` = sapply(d$interaction_p, fmt_p),
  check.names = FALSE
)

title <- "Supplementary Table 17. Associations of DISCO and disease-specific polygenic risk scores (PRS) with incident chronic diseases and their interaction in UKB. DISCO and PRS were standardized per SD; Cox models adjusted for age and sex (sex-specific cancers: age only). * P<0.05, ** P<0.01, *** P<0.001."

wb <- loadWorkbook(template)
if ("Supplementary Table 17" %in% names(wb)) removeWorksheet(wb, "Supplementary Table 17")
addWorksheet(wb, "Supplementary Table 17")
writeData(wb, "Supplementary Table 17", title, startRow = 1, startCol = 1)
mergeCells(wb, "Supplementary Table 17", cols = 1:ncol(tbl), rows = 1)
writeData(wb, "Supplementary Table 17", tbl, startRow = 3, colNames = TRUE)
hs <- createStyle(textDecoration = "bold")
addStyle(wb, "Supplementary Table 17", hs, rows = 1, cols = 1, gridExpand = FALSE)
addStyle(wb, "Supplementary Table 17", hs, rows = 3, cols = 1:ncol(tbl), gridExpand = TRUE)
setColWidths(wb, "Supplementary Table 17", cols = 1:ncol(tbl),
             widths = c(28, 16, 10, 8, 20, 10, 20, 10, 22, 12))

backup <- paste0(tgt, ".bak_R05")
if (file.exists(tgt) && !file.exists(backup)) file.copy(tgt, backup)
saveWorkbook(wb, tgt, overwrite = TRUE)

cat("Added Supplementary Table 17 (", nrow(tbl), "rows) to:\n", tgt, "\n", sep = "")
cat("Original template preserved:", template, "\n")
