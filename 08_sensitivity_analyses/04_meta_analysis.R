#!/usr/bin/env Rscript
# Single source: current Table 1, five independent cohorts, final analysis-sample SD.
suppressPackageStartupMessages({library(digest);library(jsonlite)})
entry <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo <- dirname(dirname(normalizePath(entry,mustWork=TRUE)))
source(file.path(repo,'approved/meta_functions.R'))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L || args[3]!='--allow-provisional')
  stop('Usage: Rscript --vanilla 04_meta_analysis.R TABLE1_CSV NEW_OUTPUT --allow-provisional')
input <- normalizePath(args[1],mustWork=TRUE); out <- args[2]
if(file.exists(out))stop('Output exists; choose a new directory')
before <- digest(file=input,algo='sha256')
contract <- fromJSON(file.path(repo,'approved/table1/data/contract.json'))
if(before!=contract$approved_source_sha256)
  stop('Unapproved Table1/endpoint snapshot: review and explicitly version the source contract before pooling')
d <- table1_meta_input(read.csv(input,stringsAsFactors=FALSE,check.names=FALSE))
pool <- loo <- list()
for (metric in c('DISCO','DM')) {
  x <- d[d$metric==metric,]
  pool[[metric]] <- dl_pool(x$beta,x$SE,metric)
  loo[[metric]] <- do.call(rbind,lapply(seq_len(nrow(x)),function(i)
    dl_pool(x$beta[-i],x$SE[-i],metric,x$cohort[i])))
}
stopifnot(identical(before,digest(file=input,algo='sha256')))
dir.create(out,recursive=TRUE)
write.csv(d,file.path(out,'R04_input_estimates.csv'),row.names=FALSE)
write.csv(do.call(rbind,pool),file.path(out,'R04_pooled_estimates.csv'),row.names=FALSE)
write.csv(do.call(rbind,loo),file.path(out,'R04_leave_one_cohort_out.csv'),row.names=FALSE)
write_json(list(status='NUMERIC_RECOMPUTED_CONDITIONAL_PROVISIONAL_ENDPOINTS',
  input_basename=basename(input),input_sha256=before,source_unchanged=TRUE,rows=nrow(d),
  model='cohort-specific Table1 core adjustment; race additionally NHANES/UKB',
  SD='sample SD (N-1), final common DISCO/DM analysis sample; DM log-transformed',
  endpoint_version=unique(d$endpoint_policy_version),source_contract=contract$version,
  limits=c('Five independent cohorts only; UKB modalities are not independent cohorts',
    'NHANES/UKB include race; other cohorts do not. Adjustment is not identical',
    'No pooling of crude, age+sex-only, per-unit, old SD or logistic OR estimates',
    'Numerical meta-analysis does not certify raw survival definitions',
    'Five-cohort random-effects Wald inference is imprecise; algorithm unchanged')),
  file.path(out,'audit_status.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(do.call(rbind,pool)[c('metric','k','HR','lower95','upper95','I2')],row.names=FALSE)
