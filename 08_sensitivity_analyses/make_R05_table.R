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

# R05 interaction results -> formatted supplementary table (CSV + XLSX).
suppressPackageStartupMessages({library(openxlsx)})
source(disco_code_path('approved/prs_joint.R'))

wkd <- disco_results_path("08_sensitivity_analyses")
src <- file.path(wkd, "outputs/R05_genetic_interaction_rebuilt")
out <- src
dir.create(out, recursive = TRUE, showWarnings = FALSE)

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

star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fmt2 <- function(x) sprintf("%.2f", x)
ci <- function(lo, hi) sprintf("%.2f\u2013%.2f", lo, hi)

d$Disease <- lbl[d$disease]
d$System <- factor(sys[d$disease], levels = sys_order)
d$Version <- ifelse(d$version == "clinical", "Clinical DISCO", "Proteomics DISCO")
d$multiplicative <- d$disco_hr * d$prs_hr
# Legacy CSV main effects come from f0, interaction from f1: mixed multiplication
# is not an estimand. Only an explicitly supplied same-f1 Cox object is accepted.
d$joint <- d$joint_lo <- d$joint_hi <- NA_real_
d$joint_note <- 'Unavailable: f0 main effects and f1 interaction cannot be mixed; same-f1 beta/vcov required'
joint_path <- Sys.getenv('DISCO_R05_JOINT_MODELS_RDS')
if(nzchar(joint_path)){
  fits <- readRDS(joint_path)
  if(!is.list(fits) || is.null(names(fits)))stop('Named same-f1 Cox model list required')
  for(i in seq_len(nrow(d))){
    key <- paste(d$disease[i],d$version[i],sep='__')
    if(!key%in%names(fits))next
    if(!inherits(fits[[key]],'coxph') || fits[[key]]$n!=d$n[i] || fits[[key]]$nevent!=d$events[i])
      stop('Supplied f1 does not match the table evaluation N/events; choose and verify the branch explicitly')
    h <- disco_joint_from_f1(fits[[key]])
    d$joint[i]<-h$HR;d$joint_lo[i]<-h$lower95;d$joint_hi[i]<-h$upper95;d$joint_note[i]<-h$note
  }
}

d <- d[order(d$System, d$disease, d$version), ]

tbl <- data.frame(
  System = as.character(d$System),
  Disease = d$Disease,
  Version = d$Version,
  N = d$n,
  Events = d$events,
  `DISCO HR (per SD)` = round(d$disco_hr, 2),
  `DISCO 95% CI` = ci(d$disco_lo, d$disco_hi),
  `DISCO P` = signif(d$disco_p, 2),
  `PRS HR (per SD)` = round(d$prs_hr, 2),
  `PRS 95% CI` = ci(d$prs_lo, d$prs_hi),
  `PRS P` = signif(d$prs_p, 2),
  `Interaction HR` = round(d$interaction_hr, 3),
  `Interaction 95% CI` = ci(d$interaction_lo, d$interaction_hi),
  `Interaction P` = signif(d$interaction_p, 2),
  `Main-effect contrast from f0 (not interaction model)` = round(d$multiplicative, 2),
  `Joint HR from one f1` = round(d$joint, 2),
  `Joint 95% CI from full f1 covariance` = ifelse(is.na(d$joint),NA_character_,ci(d$joint_lo,d$joint_hi)),
  `Joint source note` = d$joint_note,
  `C-index` = round(d$c_index, 3),
  check.names = FALSE
)

write.csv(tbl, file.path(out, "R05_genetic_interaction_table.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "R05_interaction")
writeData(wb, "R05_interaction", tbl)
# bold header
hs <- createStyle(textDecoration = "bold")
addStyle(wb, "R05_interaction", hs, rows = 1, cols = 1:ncol(tbl), gridExpand = TRUE)
freezePane(wb, "R05_interaction", firstRow = TRUE)
setColWidths(wb, "R05_interaction", cols = 1:ncol(tbl),
             widths = c(14, 26, 16, 9, 8, 15, 14, 9, 14, 14, 9, 12, 15, 11, 20, 16, 20, 48, 8))
saveWorkbook(wb, file.path(out, "R05_genetic_interaction_table.xlsx"), overwrite = TRUE)

cat("saved table: ", nrow(tbl), " rows (clinical",
    sum(d$version == "clinical"), "+ proteomics", sum(d$version == "proteomics"), ")\n")
