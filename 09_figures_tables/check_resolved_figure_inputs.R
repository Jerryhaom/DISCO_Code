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

# Regression checks for the content-based object disambiguations of 2026-09-17.
# Only reads scientific inputs; report is written inside Results_final.
root <- disco_project_root()
out <- disco_results_path("09_figures_tables/outputs/original_block_replay")
checks <- character()
ok <- function(label, condition) {
  if (!isTRUE(condition)) stop(label)
  checks <<- c(checks,paste("PASS",label,sep="\t"))
}
p <- readRDS(file.path(out,"source_blocks.rds"))
ok("original source files unchanged since captured blocks",all(vapply(p,function(x)
  identical(unname(tools::md5sum(x$source)),x$md5),logical(1))))

m <- readRDS(file.path(out,"Fig4d_original_pcorr_matrix.rds"))
old <- readRDS(disco_results_path("05_organ_network/outputs/figures/Fig4d_pcorr_matrix.rds"))
ok("Fig4d original matrix vs independent module replay",isTRUE(all.equal(m,old,tolerance=1e-12)))
ok("Fig4d is symmetric, finite, 18 by 18",identical(dim(m),c(18L,18L)) &&
     all(is.finite(m)) && max(abs(m-t(m)))<1e-14 && all(diag(m)==1))

d <- readRDS(file.path(out,"ED11a_overlap_plot_data.rds"))
target_r <- c(.58,.69,.65,.59,.67,.69,.69,.69)
target_p <- c(.012,.0016,.0033,.0094,.0023,.0015,.0016,.0017)
parts <- split(d,d$di)
tests <- lapply(parts,function(x) cor.test(x$ci,x$de))
ok("ED11a all eight printed Pearson r values",identical(
  unname(vapply(tests,function(x) round(unname(x$estimate),2),numeric(1))),target_r))
ok("ED11a all eight printed P values",identical(
  unname(vapply(tests,function(x) signif(x$p.value,2),numeric(1))),target_p))
ok("ED11a degree sourced from rebuilt Fig4d matrix",max(abs(d$de-(rowSums(m)-1)[as.character(d$tissue)]))<1e-12)

c <- readRDS(file.path(out,"ED9a_selected_plot_data.rds"))
b <- readRDS(file.path(out,"ED9b_overlap_counts.rds"))
counts <- tapply(c$value,c$tissue,sum)
ok("ED9a selected counts equal ED9b after-selection diagonal for all 18 organs",
   all(counts==diag(b$after)[match(names(counts),rownames(b$after))]))
ok("ED9a highlighted counts 158/143/128",identical(as.numeric(counts[c("lymphoid tissue","brain","liver")]),c(158,143,128)))
for (tag in c("ED2d","ED2e","ED2f")) {
  x <- readRDS(file.path(out,paste0(tag,"_plot_data.rds")))
  ok(paste(tag,"unique aligned IDs and finite hd/en"),!anyDuplicated(x$id) &&
       length(x$id)==nrow(x$data) && all(is.finite(as.matrix(x$data))) && all(x$data$hd>0))
  ok(paste(tag,"DM/entropy rank agreement"),cor(x$data$hd,x$data$en,method="spearman")>.999999)
}
writeLines(checks,file.path(out,"resolved_input_regression_checks.txt"))
cat(paste(checks,collapse="\n"),"\n")
