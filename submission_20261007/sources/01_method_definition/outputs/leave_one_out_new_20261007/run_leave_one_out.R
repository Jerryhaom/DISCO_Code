#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Transparent new ED4i pipeline. The historical upper loop/seed is not recovered.
# Full mode requires an explicit endpoint CSV and reference sample size.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(Rcpp);library(survival);library(data.table);library(digest);library(jsonlite)})
args <- commandArgs(TRUE)
opt <- function(name,default=NULL) {
  i <- match(paste0('--',name),args)
  if(is.na(i)) return(default)
  if(i==length(args)) stop('Missing value for --',name)
  args[i+1L]
}
base <- opt('code-dir',.disco_path("results","01_method_definition/outputs/leave_one_out_new_20261007"))
mode <- opt('mode','smoke'); stopifnot(mode %in% c('smoke','full'))
out <- opt('out',file.path(base,mode))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
clin_path <- opt('clinical',.disco_path("legacy","Result1/Clean.UKB.Clinical.rds"))
meta_path <- opt('metadata',.disco_path("legacy","Result1/UKB.Clinical.DISCO.rds"))
ep_path <- opt('endpoint')
nref <- as.integer(opt('n-ref',if(mode=='smoke') '500' else NA_character_))
neval <- as.integer(opt('n-eval',if(mode=='smoke') '500' else '50000'))
nrep <- as.integer(opt('reps',if(mode=='smoke') '1' else '50'))
private_sampling <- identical(tolower(opt('save-private-sampling','false')),'true')
if(mode=='full' && (is.null(ep_path)||is.na(nref))) stop('Full run requires --endpoint and --n-ref explicitly')
stopifnot(nref>=3L,neval>=20L,nrep>=1L)
sourceCpp(file.path(base,'pearson_leave_one.cpp'),cacheDir=file.path(base,'compile_cache'),rebuild=FALSE)
md5_before <- tools::md5sum(c(clin_path,meta_path,if(!is.null(ep_path)) ep_path))
clin <- readRDS(clin_path); meta <- readRDS(meta_path)
stopifnot(all(c('sampleID') %in% names(clin)),all(c('ID','age','sex','time','status') %in% names(meta)),
  !anyDuplicated(clin$sampleID),!anyDuplicated(meta$ID),setequal(clin$sampleID,meta$ID))
# Stable ID matching, never positional attachment of outcomes.
ix <- match(clin$sampleID,meta$ID);stopifnot(!anyNA(ix))
meta <- meta[ix,];stopifnot(identical(as.character(clin$sampleID),as.character(meta$ID)))
X <- as.matrix(clin[,setdiff(names(clin),'sampleID'),drop=FALSE]);rm(clin)
stopifnot(ncol(X)==58L,all(is.finite(X)))
markers <- colnames(X)
endpoint_identity <- 'SMOKE ONLY: saved historical metadata endpoints'
eligible <- complete.cases(meta[,c('age','sex','time','status')]) & meta$time>0 & meta$status %in% 0:1
if(!is.null(ep_path)) {
  ep <- fread(ep_path); id <- opt('endpoint-id','eid');tm <- opt('endpoint-time','time');ev <- opt('endpoint-event','status')
  stopifnot(all(c(id,tm,ev) %in% names(ep)),!anyDuplicated(ep[[id]]))
  j <- match(as.character(meta$ID),as.character(ep[[id]]));hit <- !is.na(j)
  meta$time <- ep[[tm]][j];meta$status <- ep[[ev]][j]
  eligible <- hit & complete.cases(meta[,c('age','sex','time','status')]) &
    is.finite(meta$time) & meta$time>0 & meta$status %in% 0:1
  if('reason' %in% names(ep)) eligible <- eligible & (is.na(ep$reason[j])|ep$reason[j]=='')
  endpoint_identity <- paste('CONDITIONAL current endpoint CSV, time unit:',opt('endpoint-time-unit','months'))
  rm(ep)
}
# Preserve the existing clinical sensitivity full-panel age weighting. These
# are fixed weights, not independent training/test weights; omission only
# subsets and re-normalises them. That interpretation is documented.
cc <- apply(X,2,function(x)cor(x,meta$age,use='complete.obs'))
w <- abs(cc %o% cc);diag(w)<-0;w <- w/sum(w)
pool <- which(meta$age<=40 & is.finite(meta$age))
candidates <- setdiff(which(eligible),pool)
stopifnot(nref<=length(pool),neval<=length(candidates))
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(20260914);evaluation <- sample(candidates,neval)
if(private_sampling) {
  saveRDS(meta$ID[evaluation],file.path(out,'PRIVATE_evaluation_IDs.rds'))
}
Xa <- X[evaluation,,drop=FALSE]
# Direct, separately evaluated literal cor(rbind(...)) checks.
literal <- function(ref,q,w) {
  rb <- cor(ref);rn <- cor(rbind(ref,q))
  log(nrow(ref)^2*sum((rn-rb)^2*w)/sum(w))
}
set.seed(123);validation_ref <- X[sample(pool,nref),,drop=FALSE]
vf <- pearson_all_omissions(validation_ref,Xa[seq_len(min(5L,nrow(Xa))),,drop=FALSE],w)
validate_cols <- c(1L,match('CRP',markers)+1L,match('bmi',markers)+1L)
errors <- vapply(validate_cols,function(k) {
  keep <- if(k==1L) seq_len(ncol(X)) else setdiff(seq_len(ncol(X)),k-1L)
  v <- vapply(seq_len(nrow(vf)),function(i)literal(validation_ref[,keep,drop=FALSE],Xa[i,keep],w[keep,keep,drop=FALSE]),numeric(1))
  max(abs(v-vf[,k]))
},numeric(1))
stopifnot(max(errors)<1e-7)
rows <- list();timing <- list();score_digests <- list()
groups <- c('All',markers)
selected <- if(mode=='smoke') validate_cols else seq_along(groups)
for(rep in seq_len(nrep)) {
  set.seed(123L+(rep-1L));ri <- sample(pool,nref)
  if(private_sampling) saveRDS(meta$ID[ri],file.path(out,sprintf('PRIVATE_reference_IDs_rep%02d.rds',rep)))
  started <- proc.time()[['elapsed']]
  scores <- pearson_all_omissions(X[ri,,drop=FALSE],Xa,w)
  stopifnot(all(is.finite(scores)))
  score_digests[[as.character(rep)]] <- digest(scores,algo='sha256')
  for(k in selected) {
    dat <- data.frame(tm=meta$time[evaluation],ev=meta$status[evaluation],score=scores[,k])
    warns <- character()
    fit <- withCallingHandlers(coxph(Surv(tm,ev)~score,data=dat),warning=function(e){warns<<-c(warns,conditionMessage(e));invokeRestart('muffleWarning')})
    s <- summary(fit)
    rows[[length(rows)+1L]] <- data.frame(rep=rep,phe=groups[k],N=fit$n,deaths=fit$nevent,
      C_index=unname(s$concordance[1]),C_SE=unname(s$concordance[2]),
      warning_count=length(warns),warning=paste(unique(warns),collapse=' | '))
  }
  timing[[rep]] <- data.frame(rep=rep,n_eval=neval,n_ref=nref,models_run=length(selected),elapsed_sec=proc.time()[['elapsed']]-started)
  cat('Completed rep',rep,'/',nrep,';',length(selected),'Cox fits\n')
}
res <- do.call(rbind,rows);tim <- do.call(rbind,timing)
write.csv(res,file.path(out,'leave_one_out_Cindex.csv'),row.names=FALSE)
if(mode=='full') {
  suppressPackageStartupMessages(library(ggplot2))
  # Retain the actual historical Plot_leave.R plotting design; only map the
  # new explicitly keyed C-index column and exclude the non-biomarker ID.
  med <- aggregate(C_index~phe,res,median)
  order_phe <- unique(c('All',med$phe[order(med$C_index)]))
  plot_res <- res;plot_res$phe <- factor(plot_res$phe,levels=order_phe)
  p <- ggplot(plot_res,aes(phe,C_index))+geom_boxplot()+
    geom_hline(yintercept=mean(plot_res$C_index[plot_res$phe!='All']),
      linetype='dashed',color='red',linewidth=1)+
    ylab('DISCO C-index')+xlab('Biomarkers')+theme_classic()+
    labs(title='Leave-one-out analysis of DISCO')+
    theme(axis.text.x=element_text(size=8,angle=60,hjust=1,color='black'))
  ggsave(file.path(out,'ED4i_leave_one_out_CONDITIONAL_NEW.pdf'),p,height=2.5,width=8)
}
write.csv(tim,file.path(out,'timing.csv'),row.names=FALSE)
write.csv(data.frame(marker=markers,age_correlation=cc),file.path(out,'marker_weights.csv'),row.names=FALSE)
write.csv(data.frame(group=groups[validate_cols],literal_max_error=errors),file.path(out,'literal_validation.csv'),row.names=FALSE)
stopifnot(identical(md5_before,tools::md5sum(names(md5_before))))
write_json(list(mode=mode,status=if(mode=='full') 'CONDITIONAL_NEW_FULL_COMPLETED' else 'SMOKE_ONLY_NOT_FULL',
  endpoint=endpoint_identity,input_N=nrow(X),biomarkers=58,reference_pool_N=length(pool),reference_sample_N=nref,
  evaluation_N=neval,events=sum(meta$status[evaluation]),replicates=nrep,
  all_group_plus_omission_models_run=nrow(res),full_target_models=50L*59L,
  reference_rule='age <= 40; sampling without replacement; seed=123+(rep-1)',evaluation_seed=20260914,
  historical_ED4i_seed_recovered=FALSE,input_md5=as.list(md5_before),inputs_preserved=TRUE,
  private_sampling_lists_saved=private_sampling,
  maximum_literal_score_error=max(errors),selection_sha256=digest(meta$ID[evaluation],algo='sha256'),
  scores_sha256=score_digests,
  interpretation='Crude apparent C-index on a fixed evaluation sample. Fixed full-panel age-correlation weights subset and renormalised for each omission. No independent test-set or historical numeric identity claim.'),
  file.path(out,'verification.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('COMPLETE:',mode,'; individual IDs, score rows, and Cox model objects were not exported.\n')
