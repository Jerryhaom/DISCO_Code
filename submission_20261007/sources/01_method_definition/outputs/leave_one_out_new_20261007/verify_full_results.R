# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
base <- .disco_path("results","01_method_definition/outputs/leave_one_out_new_20261007")
out <- file.path(base,'full_new_endpoint_ref500')
x <- read.csv(file.path(out,'leave_one_out_Cindex.csv'))
v <- jsonlite::read_json(file.path(out,'verification.json'),simplifyVector=TRUE)
stopifnot(nrow(x)==2950L,length(unique(x$phe))==59L,!anyDuplicated(paste(x$rep,x$phe)),
  setequal(x$rep,1:50),all(x$N==50000L),all(x$deaths==4452L),all(x$warning_count==0),
  all(is.finite(x$C_index)),all(x$C_index>=0 & x$C_index<=1),all(is.finite(x$C_SE)),
  v$private_sampling_lists_saved,v$inputs_preserved,v$maximum_literal_score_error<1e-10,
  identical(tools::md5sum(file.path(base,'run_leave_one_out.R'))[[1]],
    tools::md5sum(file.path(out,'executed_run_leave_one_out.R'))[[1]]))
all_c <- x$C_index[x$phe=='All'][match(1:50,x$rep[x$phe=='All'])]
groups <- unique(x$phe)
tab <- do.call(rbind,lapply(groups,function(g) {
  y <- x[x$phe==g,];stopifnot(setequal(y$rep,1:50),nrow(y)==50)
  y <- y[match(1:50,y$rep),];delta <- y$C_index-all_c
  data.frame(omitted_marker=g,n_rep=50L,Cmean=mean(y$C_index),Cmedian=median(y$C_index),
    C_reference_SD=sd(y$C_index),Cmin=min(y$C_index),Cmax=max(y$C_index),
    paired_delta_mean=mean(delta),paired_delta_reference_SD=sd(delta))
}))
write.csv(tab,file.path(out,'summary_by_omitted_marker.csv'),row.names=FALSE)
tim <- read.csv(file.path(out,'timing.csv'))
writeLines(c('PASS: 2950 unique (rep,marker) rows; 50 paired reference resamples x59 groups.',
  'All models: evaluation N50000, deaths4452, warnings0; reference age<=40 and n500.',
  paste('Total loop seconds:',sum(tim$elapsed_sec)),
  paste('All C mean:',mean(all_c),'median:',median(all_c),'reference SD:',sd(all_c)),
  'PRIVATE evaluation/reference IDs remain server-side, not downloaded.',
  'PDF rendered and visually inspected; original 8x2.5 horizontal boxplot design retained.',
  'Boxplot spread/paired deltas reflect reference-subset variability, not cohort-level confidence intervals.',
  'Current protocol is conditional; historical seed/upstream-loop recovery not claimed.'),file.path(out,'verification_local.txt'))
cat('PASS: full keyed results, source identity, aggregate contracts and paired reference summaries\n')
