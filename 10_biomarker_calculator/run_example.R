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
disco_init_module("10_biomarker_calculator")


# Fictional inputs only: no participant records or controlled data needed.
if (!requireNamespace("DISCO4HD", quietly=TRUE) || packageVersion("DISCO4HD") < "0.3.0")
  stop("Install DISCO4HD >= 0.3.0 from Jerryhaom/DISCO4HD")
x <- data.frame(example=c("A","B"),CRP=c(1,8),GLU=c(5,7),CHOL=c(4.5,6),
                HDL=c(1.3,.8),CREA=c(80,130),UA=c(300,480),WBC=c(6,12),
                UREA=c(5,9),RBC=c(4.8,3.4),ALB=c(42,30))
refs <- DISCO4HD::disco10_references()
scores <- do.call(rbind,lapply(refs$cohort,function(cohort) DISCO4HD::cal_disco10(x,reference=cohort)))
relative <- DISCO4HD::cal_disco10(x,reference="UKB",comparison="CHARLS")
stopifnot(all(is.finite(scores$DISCO)), all(is.finite(scores$DM)),
          abs(scores$DISCO[1] - (-2.118126394406057)) < 1e-9,
          abs(scores$DM[1] - 2.152268597777491) < 1e-9,
          all(relative$comparison_marker_count==8))
dir.create("outputs",showWarnings=FALSE)
write.csv(refs,"outputs/references.csv",row.names=FALSE)
write.csv(scores,"outputs/example_scores.csv",row.names=FALSE)
write.csv(relative,"outputs/example_comparison.csv",row.names=FALSE)
writeLines(c("PASS: frozen references and web fixture",capture.output(sessionInfo())),"outputs/validation.txt")
print(scores[c("example","reference_cohort","marker_count","DISCO","DM")],row.names=FALSE)
cat("PASS: 10-biomarker reference example\n")
