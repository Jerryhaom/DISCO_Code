#!/usr/bin/env Rscript
# Correct only demonstrable PCA inverse and duplicated-row bugs in original blocks.
options(warn=1, device=function(...) grDevices::pdf(file=NULL,...))
args <- commandArgs(TRUE)
stopifnot(length(args)==1L, args[1] %in% c("pass_1","pass_2"))
base <- Sys.getenv('DISCO_RUN_OUTPUT')
out <- file.path(base,args[1]); dir.create(out, recursive=TRUE, showWarnings=FALSE)
src <- file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/geometry_original_blocks.R')
src_md5 <- unname(tools::md5sum(src)); txt <- readLines(src,warn=FALSE)
blocks <- list(); values <- list(); pca_checks <- list()
e <- new.env(parent=globalenv())
e$pdf <- function(file,...) grDevices::pdf(file.path(out,basename(file)),...)
e$saveRDS <- e$write.csv <- e$write.table <- function(...) stop("Original writes forbidden")
e$setwd <- function(...) stop("Original setwd forbidden")
e$.inverse_pca <- function(pcs,pca) {
  z <- as.matrix(pcs) %*% t(pca$rotation)
  x <- sweep(sweep(z,2,pca$scale,"*"),2,pca$center,"+")
  back <- predict(pca,newdata=x)
  legacy <- pca$center + z
  check <- data.frame(rho=e$rho, n=nrow(x), p=ncol(x),
    max_roundtrip_error=max(abs(back-as.matrix(pcs))),
    max_legacy_feature_error=max(abs(x-legacy)))
  stopifnot(check$max_roundtrip_error < 1e-12)
  pca_checks[[length(pca_checks)+1L]] <<- check
  values[[paste0("PCA_",e$rho)]] <<- list(pca=pca, pcs=pcs, corrected=x, legacy=legacy)
  x
}
block <- function(first,last,transform=identity) {
  original <- paste(txt[first:last],collapse="\n")
  revised <- transform(original)
  blocks[[length(blocks)+1L]] <<- list(first=first,last=last,source=src,original=original,revised=revised)
  for(expr in parse(text=revised)) {
    v <- withVisible(eval(expr,e))
    if(v$visible && inherits(v$value,c("ggplot","ggarrange","patchwork")) && dev.cur()>1) print(v$value)
  }
}
block(3,101)
suppressPackageStartupMessages(library(ggpubr))
block(131,133) # original standardize function
block(241,434,function(x) {
  old <- "pca$center + as.matrix(pc_full_grid) %*% t(pca$rotation)"
  stopifnot(length(gregexpr(old,x,fixed=TRUE)[[1]])==2L)
  x <- gsub(old,".inverse_pca(pc_full_grid,pca)",x,fixed=TRUE)
  gsub("scale_x_continuous(expand = c(0,0))",
    "scale_x_continuous(expand = c(0,0), breaks=function(x) pretty(x,n=3))",x,fixed=TRUE)
})
layers <- function(p) lapply(p$layers,function(x)x$data)
values$S2 <- list(low=layers(e$p_contour_low),high=layers(e$p_contour_high))
e$.capture <- function(name,rho,x) {
  stopifnot(nrow(x)==1000L, !anyNA(x), all(is.finite(as.matrix(x))))
  values[[paste0(name,"_",rho)]] <<- x
}
fix_rows <- function(x) {
  x <- gsub("data <- mvrnorm(n_samples, mean_vector, cov_matrix)",
     "mean_vector <- rep(0,n_dim)\n  data <- mvrnorm(n_samples, mean_vector, cov_matrix)",x,fixed=TRUE)
  x <- gsub("    PC1 = pc_grid$PC1,\n    PC2 = pc_grid$PC2,\n", "",x,fixed=TRUE)
  x
}
block(442,503,function(x) {
  x <- fix_rows(x)
  # S4 only: original decimal KL labels overlap in narrow right-hand panels.
  # Change tick formatting/count, not data, limits, font, layout, or dimensions.
  stopifnot(length(gregexpr("theme_classic()",x,fixed=TRUE)[[1]])==3L)
  gsub("theme_classic()",
    "theme_classic() + scale_x_continuous(breaks=scales::breaks_pretty(n=2), labels=scales::label_scientific(digits=2))",
    x,fixed=TRUE)
})
trace("plot_KL",where=e,exit=quote(.capture("S4",rho,grid_distances)),print=FALSE)
block(505,512); untrace("plot_KL",where=e)
block(514,575,fix_rows)
trace("plot_KLL",where=e,exit=quote(.capture("S3",rho,grid_distances)),print=FALSE)
block(576,583); untrace("plot_KLL",where=e)
while(dev.cur()>1) dev.off()
old <- readRDS(file.path(Sys.getenv('DISCO_CODE_ROOT'),'resources/geometry_prior_objects.rds'))
checks <- do.call(rbind,lapply(grep("^S[34]_",names(values),value=TRUE),function(nm) {
  new <- values[[nm]]; orig <- old[[nm]]
  cols <- colnames(new)
  err <- max(abs(as.matrix(new)-as.matrix(orig[1:1000,cols])))
  cr <- cor(new,method="spearman"); oldcr <- cor(orig[,cols],method="spearman")
  stopifnot(err==0, max(abs(cr-oldcr))<1e-12)
  data.frame(panel=nm,n_original=nrow(orig),n_corrected=nrow(new),
    max_distance_change=err,max_spearman_change=max(abs(cr-oldcr)),
    disco_kl=cr["DISCO","KL"],dm_kl=cr["Mahalanobis","KL"],entropy_kl=cr["Entropy","KL"])
}))
saveRDS(values,file.path(out,"numeric_objects.rds")); saveRDS(blocks,file.path(out,"source_blocks.rds"))
write.csv(checks,file.path(out,"S3_S4_dedup_checks.csv"),row.names=FALSE)
write.csv(do.call(rbind,pca_checks),file.path(out,"S2_inverse_checks.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
stopifnot(src_md5==unname(tools::md5sum(src)))
jsonlite::write_json(list(source=src,source_md5=src_md5,original_unchanged=TRUE,
  changes=c("S2 scale and columnwise centering restored in PCA inverse", "S2 fewer x ticks to avoid collisions at facet boundaries, same dimensions and font sizes", "S3/S4 1000 genuine rows, remove unrelated 10000-row PC grid and define mean locally", "S4 shorter scientific KL ticks with fewer breaks; data, fonts, dimensions and layout unchanged"),
  unchanged=c("original RNG seed and draw order", "original DISCO and distance definitions", "original ggplot designs/fonts/dimensions"),
  scope="New corrected simulation, not recovery of historical final source"),file.path(out,"evidence.json"),pretty=TRUE,auto_unbox=TRUE)
cat("COMPLETE",args[1],"\n")
