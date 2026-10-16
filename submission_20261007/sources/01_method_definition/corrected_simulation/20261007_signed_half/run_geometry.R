#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Minimal corrections to the actual original plotting blocks. Do not choose abs
# entropy or change the scientific model. The source and old outputs stay intact.
options(warn=1,device=function(...) grDevices::pdf(file=NULL,...))
args <- commandArgs(TRUE)
pass <- if(length(args)) args[1] else 'pass_1'
stopifnot(pass %in% c('pass_1','pass_2'))
base <- .disco_path("results","01_method_definition/corrected_simulation/20261007_signed_half")
out <- file.path(base,pass)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- .disco_path("results","09_figures_tables/checklist_actions_20261006/outputs/S1_S5_redraw_05/Sim_Distances_entropy05.R")
source_md5 <- unname(tools::md5sum(src))
txt <- readLines(src,warn=FALSE)
stopifnot(length(txt)==798L)
blocks <- list(); values <- list(); refs <- list(); pca_checks <- list()
e <- new.env(parent=globalenv())
e$pdf <- function(file,...) grDevices::pdf(file.path(out,basename(file)),...)
e$setwd <- function(...) stop('Original setwd forbidden')
e$saveRDS <- e$write.csv <- e$write.table <- function(...) stop('Original writes forbidden')
e$.capture_grid <- function(x) {
  pc <- 'PC1' %in% colnames(x)
  nm <- if(pc) paste0('S2_',e$rho) else paste0('S1_',round(e$cov_matrix[1,2],1))
  stopifnot(nrow(x)==10000L,all(is.finite(as.matrix(x))))
  values[[nm]] <<- x
  refs[[nm]] <<- e$reference_data
}
e$.inverse_pca <- function(pcs,pca) {
  z <- as.matrix(pcs) %*% t(pca$rotation)
  x <- sweep(sweep(z,2,pca$scale,'*'),2,pca$center,'+')
  check <- data.frame(rho=e$rho,n=nrow(x),p=ncol(x),
    max_roundtrip_error=max(abs(predict(pca,newdata=x)-as.matrix(pcs))))
  stopifnot(check$max_roundtrip_error<1e-12)
  pca_checks[[length(pca_checks)+1L]] <<- check
  values[[paste0('PCA_',e$rho)]] <<- list(pca=pca,pcs=pcs,corrected=x)
  x
}
e$.capture <- function(name,rho,x,ref) {
  stopifnot(nrow(x)==1000L,all(is.finite(as.matrix(x))))
  nm <- paste0(name,'_',rho)
  values[[nm]] <<- x; refs[[nm]] <<- ref
}
e$.capture_noise <- function(x,ref) {
  stopifnot(nrow(x)==11000L,all(is.finite(as.matrix(x))),
            all(table(x$noise)==1000L),length(unique(x$subject))==1000L)
  values$S5 <<- x; refs$S5 <<- ref
}
block <- function(first,last,transform=identity) {
  original <- paste(txt[first:last],collapse='\n')
  revised <- transform(original)
  blocks[[length(blocks)+1L]] <<- list(first=first,last=last,original=original,revised=revised)
  for(expr in parse(text=revised)) {
    v <- withVisible(eval(expr,e))
    if(v$visible && inherits(v$value,c('ggplot','ggarrange','patchwork')) && dev.cur()>1) print(v$value)
  }
}
block(3,101)
stopifnot(grepl('0.5 * (log(new_det) - log(ref_det))',paste(txt,collapse='\n'),fixed=TRUE))
suppressPackageStartupMessages(library(ggpubr))
block(105,237,function(x) {
  gsub('grid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)',
    '.capture_grid(grid_distances)\ngrid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)',x,fixed=TRUE)
})
block(241,434,function(x) {
  old <- 'pca$center + as.matrix(pc_full_grid) %*% t(pca$rotation)'
  stopifnot(length(gregexpr(old,x,fixed=TRUE)[[1]])==2L)
  x <- gsub(old,'.inverse_pca(pc_full_grid,pca)',x,fixed=TRUE)
  x <- gsub('grid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)',
    '.capture_grid(grid_distances)\ngrid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)',x,fixed=TRUE)
  gsub('scale_x_continuous(expand = c(0,0))',
    'scale_x_continuous(expand = c(0,0), breaks=function(x) pretty(x,n=3))',x,fixed=TRUE)
})
fix_rows <- function(x) {
  x <- gsub('data <- mvrnorm(n_samples, mean_vector, cov_matrix)',
    'mean_vector <- rep(0,n_dim)\n  data <- mvrnorm(n_samples, mean_vector, cov_matrix)',x,fixed=TRUE)
  gsub('    PC1 = pc_grid$PC1,\n    PC2 = pc_grid$PC2,\n','',x,fixed=TRUE)
}
block(442,503,function(x) {
  x <- fix_rows(x)
  x <- gsub('  p1=ggplot(grid_distances,aes(KL,DISCO))',
    '  .capture("S4",rho,grid_distances,reference_data)\n  p1=ggplot(grid_distances,aes(KL,DISCO))',x,fixed=TRUE)
  gsub('theme_classic()',
    'theme_classic() + scale_x_continuous(breaks=scales::breaks_pretty(n=2), labels=scales::label_scientific(digits=2))',x,fixed=TRUE)
})
block(505,512)
block(514,575,function(x) {
  x <- fix_rows(x)
  gsub('  p1=ggplot(grid_distances,aes(Mahalanobis,DISCO))',
    '  .capture("S3",rho,grid_distances,reference_data)\n  p1=ggplot(grid_distances,aes(Mahalanobis,DISCO))',x,fixed=TRUE)
})
block(576,583)
block(584,680,function(x) {
  x <- fix_rows(x)
  x <- gsub('    grid_distances = rbind(grid_distances,distances)',
    '    distances$subject <- seq_len(nrow(reference_data))\n    grid_distances = rbind(grid_distances,distances)',x,fixed=TRUE)
  x <- gsub('  p1=ggplot(grid_distances,aes(as.factor(noise),DISCO))',
    '  .capture_noise(grid_distances,reference_data)\n  p1=ggplot(grid_distances,aes(as.factor(noise),DISCO))',x,fixed=TRUE)
  x <- gsub('labs(fill="Var")','labs(fill="SD")',x,fixed=TRUE)
  # Keep the descriptive rank correlation, not the invalid independent-row P
  # annotation for 11 repeated perturbations of the same 1000 reference people.
  for(v in c('DISCO','Mahalanobis','Entropy','KL')) {
    x <- gsub(paste0('stat_cor(aes(noise,',v,')'),
      paste0('stat_cor(aes(noise,',v,', label=after_stat(r.label))'),x,fixed=TRUE)
  }
  x
})
while(dev.cur()>1) dev.off()
# Verify against the previously independently validated 9/26 numerical branch.
old <- readRDS(.disco_path("results","01_method_definition/corrected_simulation/20260926/pass_2/numeric_objects.rds"))
checks <- do.call(rbind,lapply(grep('^S[34]_',names(values),value=TRUE),function(nm) {
  new <- values[[nm]]; expected <- old[[nm]]
  expected$Entropy <- .5*expected$Entropy
  distance_error <- max(abs(as.matrix(new)-as.matrix(expected[,colnames(new)])))
  rho_error <- max(abs(cor(new,method='spearman')-cor(expected[,colnames(new)],method='spearman')))
  stopifnot(distance_error<1e-13,rho_error<1e-12,nrow(new)==1000L)
  data.frame(panel=nm,n=nrow(new),max_distance_error=distance_error,
    max_spearman_error=rho_error,negative_entropy=sum(new$Entropy<0),
    disco_entropy=cor(new$DISCO,new$Entropy,method='spearman'),
    disco_kl=cor(new$DISCO,new$KL,method='spearman'))
}))
for(nm in grep('^PCA_',names(values),value=TRUE)) {
  stopifnot(max(abs(values[[nm]]$corrected-old[[nm]]$corrected))<1e-12)
}
saveRDS(values,file.path(out,'numeric_objects.rds'))
saveRDS(refs,file.path(out,'reference_data.rds'))
saveRDS(blocks,file.path(out,'source_blocks.rds'))
writeLines(unlist(lapply(blocks,function(b)c(paste('# Source lines',b$first,b$last),b$revised))),file.path(out,'executed_source.R'))
write.csv(checks,file.path(out,'S3_S4_verification.csv'),row.names=FALSE)
write.csv(do.call(rbind,pca_checks),file.path(out,'S2_inverse_verification.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
stopifnot(source_md5==unname(tools::md5sum(src)))
jsonlite::write_json(list(source=src,source_md5=source_md5,original_preserved=TRUE,
  entropy_definition='signed 0.5*(logdet(new)-logdet(reference)); not abs distance',
  kl_direction='KL(new||reference), unchanged',
  changes=c('correct inverse PCA scale and columnwise centering','1000 genuine S3/S4 rows without recycled PC coordinates',
    'local 10-dimensional zero mean','fewer non-overlapping S2/S4 axis ticks','S5 SD label',
    'S5 rho only, no independent-row P annotation; subject ID stored'),
  pending=c('author confirmation of signed vs abs entropy','theoretical conditions and KL formula in manuscript',
    'general singular covariance/logDISCO domain handling','independent simulation validation; same-sample reference design retained')),
  file.path(out,'evidence.json'),pretty=TRUE,auto_unbox=TRUE)
cat('COMPLETE',pass,'; technical checks passed; theoretical definition still requires author confirmation.\n')
