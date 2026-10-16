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


# Small, fully synthetic end-to-end package example.
if (!requireNamespace("DISCO4HD", quietly=TRUE) || packageVersion("DISCO4HD") < "0.3.1")
  stop("Install DISCO4HD >= 0.3.1 first")
set.seed(123)
age <- seq(20,79)
d <- data.frame(age, a=age/10+rnorm(60), b=age/30+rnorm(60), c=rnorm(60), e=age/15+rnorm(60))
v <- names(d)[-1]; idx <- which(age <= 34)
a <- DISCO4HD::cal_disco(d,v,d[idx,])
b <- DISCO4HD::cal_disco(d,v,d[idx,],cpp=TRUE)
z <- DISCO4HD::disco_modality_decompose(d[v],list(first=1:2,second=3:4),idx,age=age)
stopifnot(max(abs(a$DISCO-b$DISCO)) < 1e-9,
          max(abs(z$sub_disco[,"full"]-a$DISCO[z$analysis_idx])) < 1e-9,
          max(abs(z$decomposition[,"full"]-rowSums(z$decomposition[,-1]))) < 1e-9)
dir.create("outputs/synthetic_demo",recursive=TRUE,showWarnings=FALSE)
write.csv(a,"outputs/synthetic_demo/synthetic_scores.csv",row.names=FALSE)
writeLines(c("PASS: synthetic backend and decomposition checks", capture.output(sessionInfo())),
           "outputs/synthetic_demo/validation.txt")
cat("PASS: synthetic demo\n")
