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
disco_init_module("09_figures_tables")

# Evaluate ORIGINAL source expressions verbatim. No legacy output is overwritten.
# Statistical reproduction and visual approval are distinct: see the audit report.
options(warn = 1)
root <- disco_project_root()
out <- Sys.getenv("DISCO_REPLAY_OUT",disco_results_path("09_figures_tables/outputs/original_block_replay"))
allowed_output <- paste0(disco_results_path("09_figures_tables/outputs"),"/")
stopifnot(startsWith(out,allowed_output),!grepl("(^|/)\\.\\.(/|$)",out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out,mustWork=TRUE)
stopifnot(startsWith(out,allowed_output))
options(device=function(...) grDevices::pdf(file=NULL, ...))

suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(ggthemes); library(ggsci)
  library(survival); library(survminer); library(reshape2); library(dplyr)
})
requested <- commandArgs(trailingOnly=TRUE)
strict_capture <- identical(Sys.getenv("DISCO_REPLAY_STRICT"),"1")
if (strict_capture) {
  stopifnot(length(requested)==1)
  source(disco_code_path("09_figures_tables", "strict", "capture_runtime.R"))
}
provenance <- if (length(requested) && file.exists(file.path(out,"source_blocks.rds"))) readRDS(file.path(out,"source_blocks.rds")) else list()
new_context <- function(tag, wd) {
  e <- new.env(parent = globalenv())
  e$.tag <- tag
  e$pdf <- function(file, ...) grDevices::pdf(file.path(out, paste0(tag, "__", basename(file))), ...)
  # Host lacks libXrender needed by cairo.so. ONLY the vector-PDF backend changes;
  # width/height, ggplot expressions, data and theme remain untouched.
  e$cairo_pdf <- function(filename, ...) grDevices::pdf(file.path(out, paste0(tag, "__", basename(filename))), ...)
  e$read.csv <- function(file, ...) {
    f <- path.expand(file)
    old <- file.path(root,"Manuscript/PNAS")
    if (!file.exists(f) && dirname(f)==old && basename(f) %in% c("cindex_comparison.csv","hrs_cindex_comparison.csv")) {
      f <- file.path(old,"_archive",basename(f))
      cat("INPUT_PATH_ADAPTATION",file,"=>",f,"\n")
    }
    utils::read.csv(f,...)
  }
  e$saveRDS <- e$write.csv <- e$write.table <- function(...) stop("Legacy write blocked")
  setwd(file.path(root, wd))
  e
}
replay_block <- function(id, e) {
  full <- disco_code_path("09_figures_tables", "replay_blocks", paste0(id, ".R"))
  text <- readLines(full, warn = FALSE)
  provenance[[length(provenance) + 1L]] <<- list(tag=e$.tag, block=id, source=full,
    md5=unname(tools::md5sum(full)), original_text=text)
  cat("BUNDLED_SOURCE", id, "\n")
  if (strict_capture) strict_record_input(full, "source_block")
  for (expr in parse(text=text)) {
    result <- withVisible(eval(expr, envir=e))
    if (result$visible && inherits(result$value, c("ggplot", "ggarrange", "patchwork", "ggsurvplot"))) {
      if (dev.cur() > 1) print(result$value)
    }
  }
}
status <- list()
if (length(requested) && file.exists(file.path(out,"execution_status.txt"))) {
  old_status <- read.delim(file.path(out,"execution_status.txt"), header=FALSE)
  status <- as.list(setNames(old_status[[2]],old_status[[1]]))
}
run <- function(tag, f) {
  if (length(requested) && !tag %in% requested) return(invisible(NULL))
  if (strict_capture) .strict_active <<- TRUE
  cat("START", tag, "\n")
  status[[tag]] <<- tryCatch({f(); "EXECUTED_NOT_YET_VISUALLY_APPROVED"},
    error=function(e) paste("FAILED", conditionMessage(e)))
  while (dev.cur() > 1) dev.off()
  cat("END", tag, status[[tag]], "\n")
  writeLines(unlist(Map(function(k,v) paste(k,v,sep="\t"),names(status),status)),file.path(out,"execution_status.txt"))
  saveRDS(provenance,file.path(out,"source_blocks.rds"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  if (strict_capture) {
    strict_finish(tag,status[[tag]])
    .strict_active <<- FALSE
  }
}
run("Fig3a", function() {
  e <- new_context("Fig3a", "Result1/DNAm")
  replay_block("Result1_DNAm_NHANES_cmp_v2_R_1_102", e)
})
run("ED6a", function() {
  e <- new_context("ED6a", "Result1/DNAm")
  replay_block("Result1_DNAm_HRS_cmp_v2_R_1_92", e)
})
run("Fig2a", function() {
  e <- new_context("Fig2a", "Result2")
  replay_block("Result2_Cal_Prot_R_183_239", e)
})
run("Fig2b", function() {
  e <- new_context("Fig2b", "Result2")
  replay_block("Result2_Cal_ClinicalProt_R_220_276", e)
})
run("ED4a", function() {
  e <- new_context("ED4a", "Result1")
  replay_block("Result1_DNAm_Cmp_with_other_R_371_372", e)
  replay_block("Result1_DNAm_Cmp_with_other_R_395_412", e)
})
run("Fig3b_age", function() {
  if (!requireNamespace("ggpp", quietly=TRUE)) stop("Original geom_text_npc needs ggpp; no replacement used")
  e <- new_context("Fig3b_age", "Result1/DNAm")
  replay_block("Result1_DNAm_Cmp_with_other_R_129_130", e)
  replay_block("Result1_DNAm_Cmp_with_other_R_240_306", e)
})
run("ED1", function() {
  e <- new_context("ED1", "Simulation")
  replay_block("Simulation_Sim_Distances_R_3_101", e)
  # Original global initialization (10-dimensional mean) required by plot_KL3.
  e$mean_vector <- rep(0, 10)
  replay_block("Simulation_Sim_Distances_R_681_798", e)
})
run("Fig4a_legacy", function() {
  e <- new_context("Fig4a_legacy", "Result3/Organ")
  e$or2 <- readRDS("Tissue.order.rds")
  e$sel.t <- e$or2
  e$odis <- readRDS("UKB.odis.rds")
  replay_block("Result3_Organ_Plot_Fig4E_R_127_128", e)
  replay_block("Result3_Organ_Plot_Fig4E_R_265_314", e)
  replay_block("Result3_Organ_Plot_Fig4E_R_323_361", e)
  saveRDS(e$rr,file.path(out,"Fig4a_legacy_plot_data.rds"))
})
run("Fig4a_expanded", function() {
  e <- new_context("Fig4a_expanded", "Result3/Organ")
  legacy <- readRDS(file.path(out,"Fig4a_legacy_plot_data.rds"))
  e$rr <- legacy
  input <- disco_results_path("05_organ_network/outputs/organ_expanded_covariates")
  m <- read.csv(file.path(input,"organ_cause_specific_expanded.csv"))
  d <- read.csv(file.path(input,"organ_chronic_disease_expanded.csv"))
  names(d)[names(d)=="disease"] <- "outcome"
  m$tt <- "Mortality"; d$tt <- "Chronic disease incidence"
  cols <- c("organ","outcome","tt","Model2_HR","Model2_P")
  updated <- rbind(m[,cols],d[,cols])
  key <- function(a,b,c) paste(a,b,c,sep="|")
  idx <- match(key(e$rr$Var1,e$rr$Var2,e$rr$tt),key(updated$organ,updated$outcome,updated$tt))
  stopifnot(nrow(e$rr)==270, !anyNA(idx), !anyDuplicated(idx))
  # ONLY replace HR and significance in the original plotting data.
  e$rr$value <- updated$Model2_HR[idx]
  p <- updated$Model2_P[idx]
  e$rr$p <- ifelse(p<0.01,"**",ifelse(p<0.05," *",""))
  stopifnot(all(is.finite(e$rr$value)),all(is.finite(p)))
  replay_block("Result3_Organ_Plot_Fig4E_R_341_361", e)
  saveRDS(e$rr,file.path(out,"Fig4a_expanded_plot_data.rds"))
})
run("Fig4e_legacy", function() {
  suppressPackageStartupMessages({library(igraph);library(ggnetwork)})
  e <- new_context("Fig4e_legacy", "Result3/Organ")
  e$col <- readRDS("Tissue.col.rds")
  # Skip only library(Cairo): no Cairo-package symbol is used by this block.
  replay_block("Result3_Organ_Cal_Organ_spe_bk_R_539_541", e)
  replay_block("Result3_Organ_Cal_Organ_spe_bk_R_543_555", e)
  saveRDS(e$g2,file.path(out,"Fig4e_legacy_plot_data.rds"))
})
run("ED3a", function() {
  e <- new_context("ED3a", "Result1")
  replay_block("Result1_Cal_UKB_DISCO_R_140_140", e)
  replay_block("Result1_Cal_UKB_DISCO_R_162_246", e)
})
run("Fig4c_CANDIDATE", function() {
  suppressPackageStartupMessages({library(readr);library(glmnet);library(igraph);library(ggrepel)})
  e <- new_context("Fig4c_CANDIDATE", "Result2")
  replay_block("Result2_Proteins_Plot_Fig5AB_R_54_97", e)
  # Validate original positional join BEFORE interpreting a reconstructed plot.
  stopifnot(identical(as.character(e$common),as.character(e$pd$ID)))
  stopifnot(identical(as.character(e$common),as.character(e$prot.data$eid)))
  replay_block("Result2_Proteins_Plot_Fig5AB_R_100_143", e)
  setwd(file.path(root,"Result2/Proteins"))
  replay_block("Result2_Proteins_Plot_Fig5AB_R_224_224", e)
  replay_block("Result2_Proteins_Plot_Fig5AB_R_301_303", e)
  replay_block("Result2_Proteins_Plot_Fig5AB_R_306_335", e)
  saveRDS(e$pd1,file.path(out,"Fig4c_CANDIDATE_plot_data.rds"))
  capture.output({
    cat("Status: candidate only; original network-generation provenance remains to verify\n")
    cat("Rows:",nrow(e$pd1),"\n")
    cat("Missing degree:",sum(is.na(e$pd1$dg)),"\n")
    print(cor.test(abs(e$pd1$cc),e$pd1$dg,method="spearman"))
    print(e$pd1[e$pd1$lb1!="",c("cc","dg","lb1")])
  },file=file.path(out,"Fig4c_CANDIDATE_checks.txt"))
})
run("Fig4c_chain", function() {
  suppressPackageStartupMessages({library(readr);library(glmnet);library(igraph)})
  e <- new_context("Fig4c_chain", "Result2")
  replay_block("Result2_Proteins_Plot_Fig5AB_R_54_97", e)
  stopifnot(identical(as.character(e$common),as.character(e$pd$ID)))
  replay_block("Result2_Proteins_Plot_Fig5AB_R_100_143", e)
  setwd(file.path(root,"Result2/Proteins"))
  replay_block("Result2_Proteins_Plot_Fig5AB_R_199_214", e)
  rebuilt <- e$ig
  saved <- readRDS("Prot.net.rds")
  edges <- function(g) {
    x <- igraph::as_edgelist(g,names=TRUE)
    sort(apply(x,1,paste,collapse=" -> "))
  }
  pd <- readRDS(file.path(out,"Fig4c_CANDIDATE_plot_data.rds"))
  missing <- rownames(pd)[is.na(pd$dg)]
  absent <- setdiff(rownames(pd),V(saved)$name)
  capture.output({
    cat("Rebuilt: absolute Pearson correlation > .2 and < 1, as source lines199–214\n")
    cat("Rebuilt vertices/edges/directed:",vcount(rebuilt),ecount(rebuilt),is_directed(rebuilt),"\n")
    cat("Saved vertices/edges/directed:",vcount(saved),ecount(saved),is_directed(saved),"\n")
    cat("Identical edge multiset:",identical(edges(rebuilt),edges(saved)),"\n")
    cat("Only rebuilt:",length(setdiff(edges(rebuilt),edges(saved))),"\n")
    cat("Only saved:",length(setdiff(edges(saved),edges(rebuilt))),"\n")
    cat("Missing degree exactly absent vertices:",setequal(missing,absent),"\n")
    cat("Absent plot rows:",paste(missing,collapse=", "),"\n")
    cat("Absent rows present in rebuilt:",paste(intersect(missing,V(rebuilt)$name),collapse=", "),"\n")
  },file=file.path(out,"Fig4c_network_chain_audit.txt"))
})
run("ED4b", function() {
  e <- new_context("ED4b", "Result1/Robust")
  replay_block("Result1_Robust_Robust_DISCO_R_294_395", e)
})
run("ED4c", function() {
  e <- new_context("ED4c", "Result1/Robust")
  replay_block("Result1_Robust_Robust_DISCO_R_225_260", e)
})
run("ED4d", function() {
  e <- new_context("ED4d", "Result1")
  replay_block("Result1_Robust_Test_Distance_R_130_131", e)
  # Lines132–134 estimate an unused age-only comparator; dd1 is not saved here.
  # Its only plotting use (geom_vline) is commented out in the original source.
  replay_block("Result1_Robust_Test_Distance_R_135_177", e)
})
run("ED4e", function() {
  e <- new_context("ED4e", "Result1/Num")
  replay_block("Result1_Num_Test2_R_166_207", e)
})
run("ED4f", function() {
  e <- new_context("ED4f", "Result1/Leaveone")
  replay_block("Result1_Leaveone_Plot_leave_R_2_25", e)
})
run("ED2b", function() {
  suppressPackageStartupMessages(library(ggpointdensity))
  e <- new_context("ED2b", "Result1")
  # The source saves the entropy-bearing frame as DISCO2 at line137.
  # DISCO.rds has a DIFFERENT hd vector and must not substitute for it.
  e$pd <- readRDS("UKB.Clinical.DISCO2.rds")
  stopifnot(isTRUE(all.equal(e$pd$en,readRDS("En.UKB2.rds"))),
            isTRUE(all.equal(e$pd$hd,readRDS("HD.UKB2.rds"))))
  replay_block("Result1_Cal_UKB_DISCO_R_109_123", e)
})
run("ED2c", function() {
  suppressPackageStartupMessages(library(ggpointdensity))
  e <- new_context("ED2c", "Result2")
  e$en <- readRDS("En.UKB2.rds")
  replay_block("Result2_Cal_ClinicalProt_R_132_133", e)
  stopifnot(length(e$en)==nrow(e$pd),cor(e$pd$hd,e$en,method="spearman")>.999)
  replay_block("Result2_Cal_ClinicalProt_R_135_149", e)
})

# Original entropy mathematics, serial execution only. The source requests
# matrixcalc but calls no symbol from it; all matrix operations are base R.
entropy_context <- function(tag) {
  e <- new_context(tag,"Result4")
  e$makeCluster <- function(...) NULL
  e$registerDoParallel <- e$clusterExport <- e$stopCluster <- function(...) invisible(NULL)
  e$foreach <- function(..., .packages=NULL) {
    args <- match.call(expand.dots=FALSE)$...
    eval(as.call(c(list(quote(foreach::foreach)),args)),envir=parent.frame())
  }
  e$`%dopar%` <- foreach::`%do%`
  e$detectCores <- function(...) 2L
  replay_block("Result4_Validate_Clhls_R_15_79", e)
  cat("BACKEND_ADAPTATION: original entropy function evaluated serially; unused matrixcalc request omitted\n")
  e
}
entropy_checks <- function(e,reference,ids) {
  # Use the historical saved hd, never the locally installed BioAge stub.
  stopifnot(length(e$hd)==nrow(e$biomarker1),length(ids)==length(e$hd),
            all(is.finite(e$hd)),all(e$hd>0),!anyDuplicated(ids))
  x <- as.matrix(e$biomarker1[,e$bb])
  r <- as.matrix(reference)
  sq <- stats::mahalanobis(x,colMeans(r),cov(r))
  rho <- cor(e$hd,sq,method="spearman")
  stopifnot(rho>0.9999)
  ratio <- e$hd/sqrt(sq)
  capture.output({
    cat("Saved historical hd, not reconstructed/scaled by adapter\n")
    cat("Sample count:",nrow(x),"Reference count:",nrow(r),"Features:",ncol(x),"\n")
    cat("hd vs Mahalanobis-squared Spearman:",rho,"\n")
    cat("Observed saved hd / sqrt(Mahalanobis-squared) range:",range(ratio),"\n")
    cat("log(hd) range:",range(log(e$hd)),"\n")
    cat("Entropy range:",range(e$pn$en),"\n")
    print(cor.test(e$hd,e$pn$en,method="spearman"))
  },file=file.path(out,paste0(e$.tag,"_input_checks.txt")))
  saveRDS(list(id=ids,data=e$pn,reference_n=nrow(r),features=e$bb),
          file.path(out,paste0(e$.tag,"_plot_data.rds")))
}
run("ED2d",function() {
  suppressPackageStartupMessages(library(ggpointdensity))
  e <- entropy_context("ED2d")
  replay_block("Result4_Validate_Charls_R_38_52", e)
  replay_block("Result4_Validate_Charls_R_63_65", e)
  saved <- readRDS("CHARLS.DISCO.rds")
  stopifnot(identical(as.character(e$biomarker1$sampleID),as.character(saved$id)),
            isTRUE(all.equal(e$info$age,saved$age)))
  e$hd <- saved$hd
  replay_block("Result4_Validate_Charls_R_72_74", e)
  entropy_checks(e,e$biomarker1[which(e$info$age<=50),e$bb],saved$id)
  replay_block("Result4_Validate_Charls_R_76_90", e)
})
run("ED2e",function() {
  suppressPackageStartupMessages(library(ggpointdensity))
  e <- entropy_context("ED2e")
  replay_block("Result4_Validate_RuLAS_R_20_35", e)
  replay_block("Result4_Validate_RuLAS_R_62_68", e)
  saved <- readRDS("RuLAS.DISCO.rds")
  stopifnot(identical(as.character(e$biomarker1$sampleID),as.character(saved$ID)),
            isTRUE(all.equal(e$biomarker$age2016,saved$age)))
  e$hd <- saved$hd
  replay_block("Result4_Validate_RuLAS_R_75_78", e)
  entropy_checks(e,e$biomarker1[which(e$phe$age2016<=75),e$bb],saved$ID)
  replay_block("Result4_Validate_RuLAS_R_80_94", e)
})
run("ED2f",function() {
  suppressPackageStartupMessages(library(ggpointdensity))
  e <- entropy_context("ED2f")
  replay_block("Result4_Validate_Clhls_R_87_93", e)
  replay_block("Result4_Validate_Clhls_R_104_106", e)
  saved <- readRDS("CLHLS.DISCO.rds")
  info <- readRDS("CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds")
  idx <- match(e$biomarker$id,info$id)
  stopifnot(!anyNA(idx),length(idx)==nrow(saved),
            identical(as.character(unclass(saved$age)),as.character(unclass(info$trueage[idx]))),
            identical(as.character(unclass(saved$sex)),as.character(unclass(info$a1[idx]))))
  e$hd <- saved$hd
  replay_block("Result4_Validate_Clhls_R_113_115", e)
  entropy_checks(e,e$biomarker1[which(e$biomarker$trueage<=70),e$bb],e$biomarker$id)
  replay_block("Result4_Validate_Clhls_R_117_131", e)
})
run("Fig4d",function() {
  suppressPackageStartupMessages(library(readr))
  e <- new_context("Fig4d","Result3/Organ")
  e$prot.cla <- readRDS("UKB.Prot.class.rds")
  replay_block("Result3_Organ_Plot_Fig5C_R_230_232", e)
  e$prot1 <- readRDS("UKB.prot.imputed.rds")
  mortality <- readr::read_csv("../../UKB_Data/UKB.Mortality.csv",show_col_types=FALSE)
  e$info1 <- mortality[match(e$prot1[,1],mortality$id),]
  old <- readRDS("UKB.odis.rds")
  stopifnot(isTRUE(all.equal(as.numeric(e$info1$Age),old$age)),
            identical(as.character(e$info1$Gender),as.character(old$sex)))
  e$or2 <- readRDS("Tissue.order.rds")
  replay_block("Result3_Organ_Plot_Fig5C_R_253_292", e)
  replay_block("Result3_Organ_Plot_Fig5C_R_305_306", e)
  saved <- readRDS(disco_results_path("05_organ_network/outputs/figures/Fig4d_pcorr_matrix.rds"))
  stopifnot(isTRUE(all.equal(e$pcorr1,saved,tolerance=1e-12)))
  saveRDS(e$pcorr1,file.path(out,"Fig4d_original_pcorr_matrix.rds"))
  replay_block("Result3_Organ_Plot_Fig5C_R_313_351", e)
})
run("ED11a_overlap",function() {
  suppressPackageStartupMessages(library(readxl))
  e <- new_context("ED11a_overlap","Result3/Organ")
  # Explicit cross-script context: the submission is explained by the overlap-
  # ADJUSTED matrix from Plot_Fig5C, not the noverlap matrix loaded elsewhere
  # in Plot_Organ_Net. Eight displayed r and P values independently agree.
  e$pcorr1 <- readRDS(file.path(out,"Fig4d_original_pcorr_matrix.rds"))
  replay_block("Result3_Organ_Plot_Organ_Net_R_304_329", e)
  replay_block("Result3_Organ_Plot_Organ_Net_R_356_381", e)
  saveRDS(e$dci1,file.path(out,"ED11a_overlap_plot_data.rds"))
  checks <- do.call(rbind,lapply(split(e$dci1,e$dci1$di),function(x) {
    ct <- cor.test(x$ci,x$de)
    data.frame(disease=as.character(x$di[1]),r=unname(ct$estimate),p=ct$p.value)
  }))
  capture.output(print(checks,digits=10),file=file.path(out,"ED11a_overlap_checks.txt"))
})
run("ED9a_before",function() {
  suppressPackageStartupMessages(library(readxl))
  e <- new_context("ED9a_before","Result3/Organ")
  replay_block("Result3_Organ_Plot_Fig4A_R_6_22", e)
})
run("ED9a_selected",function() {
  suppressPackageStartupMessages(library(readxl))
  e <- new_context("ED9a_selected","Result3/Organ")
  # The same prot.cla symbol denotes TWO different objects in the legacy files.
  # The unexpanded Excel loaded by Plot_Fig4A contains one row per protein.
  # Expanded UKB.Prot.class.rds is appropriate for overlaps, NOT this count.
  replay_block("Result3_Organ_Plot_Fig4A_R_6_12", e)
  stopifnot(!anyDuplicated(e$prot.cla$Gene))
  e$sel <- readRDS("Lasso.select.prot.rds")
  e$odis <- readRDS("UKB.odis.rds")
  e$cc <- cor(e$odis[,1:36])
  e$cl <- readRDS("Selected.tissue.rds")
  e$sel.t1 <- colnames(e$cc)[e$cl]
  e$overlapped <- matrix(0,nrow=length(e$cl),ncol=length(e$cl),
                        dimnames=list(NULL,gsub("[.]"," ",colnames(e$cc)[e$cl])))
  replay_block("Result3_Organ_Plot_overlapped_R_186_219", e)
  counts <- reshape2::dcast(e$count,tissue~variable,value.var="value")
  counts$selected_n <- vapply(e$sel[e$cl],length,integer(1))[match(as.character(counts$tissue),colnames(e$overlapped))]
  stopifnot(identical(as.integer(counts$Enhanced+counts$Enriched),as.integer(counts$selected_n)))
  saveRDS(e$count,file.path(out,"ED9a_selected_plot_data.rds"))
  capture.output({
    cat("Input: ukb_simple_specific_1553.xlsx, 1553 unique Gene rows\n")
    cat("Original Plot_overlapped count/plot expressions unchanged; explicit prot.cla binding\n")
    print(counts)
  },file=file.path(out,"ED9a_selected_checks.txt"))
})
run("ED9b",function() {
  e <- new_context("ED9b","Result3/Organ")
  e$sel <- readRDS("Lasso.select.prot.rds")
  replay_block("Result3_Organ_Plot_overlapped_R_9_36", e)
  replay_block("Result3_Organ_Plot_overlapped_R_67_91", e)
  replay_block("Result3_Organ_Plot_overlapped_R_95_98", e)
  replay_block("Result3_Organ_Plot_overlapped_R_110_173", e)
  saveRDS(list(before=e$overlapped3,after=e$overlapped11),file.path(out,"ED9b_overlap_counts.rds"))
})
run("ED8c",function() {
  e <- new_context("ED8c","Result5")
  e$pd.vd <- readRDS("CHARLS.2011.disco.rds")
  replay_block("Result5_Asso_DISCO_FI_R_161_175", e)
  e$pd.vd <- readRDS("CHARLS.2015.disco.rds")
  replay_block("Result5_Asso_DISCO_FI_R_348_362", e)
})
run("ED10c_ED11cde",function() {
  e <- new_context("ED10c_ED11cde","Result3/Sensitivity")
  e$col <- readRDS("../Organ/Tissue.col.rds")
  replay_block("Result3_Sensitivity_plot_sensitivity_R_4_25", e)
  # Skip 200-replicate DISCO recomputation: the source explicitly reloads the
  # retained correlation object at line95, before all three following plots.
  replay_block("Result3_Sensitivity_plot_sensitivity_R_95_190", e)
})
if (strict_capture && (is.null(status[[requested]]) || startsWith(status[[requested]],"FAILED"))) {
  quit(save="no",status=1)
}
