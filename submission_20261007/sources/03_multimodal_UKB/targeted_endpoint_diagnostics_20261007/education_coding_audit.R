#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Aggregate comparison only: first recorded qualification versus stored derived value.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
setDTthreads(1)
out<-normalizePath(commandArgs(TRUE)[1]);base<-.disco_path("legacy","")
paths<-c(clinical=file.path(base,'UKB_Data/extended_covariates_clinical.rds'),
         proteomics=file.path(base,'UKB_Data/extended_covariates_proteomics.rds'),
         clinicalprot=file.path(base,'Result2/UKB.ClinicalProt.DISCO.rds'))
h<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
rows<-list();cross<-list();prot<-as.data.table(readRDS(paths['proteomics']))
for(mod in names(paths)){
 d<-if(mod=='clinicalprot')prot[match(readRDS(paths['clinicalprot'])$ID,prot$ID)] else as.data.table(readRDS(paths[mod]))
 q<-as.matrix(d[,paste0('edu_q',0:5),with=FALSE]);stopifnot(all(is.na(q)|q%in%c(-7,-3,1:6)))
 none<-rowSums(q==-7,na.rm=TRUE)>0;positive<-rowSums(q>0,na.rm=TRUE)>0
 current<-ifelse(d$edu_q0==-7,7,ifelse(d$edu_q0<0,NA,d$edu_q0))
 derived<-as.numeric(d$edu_high);derived[none]<-7
 equal<-(is.na(current)&is.na(derived))|(!is.na(current)&!is.na(derived)&current==derived)
 equal[is.na(equal)]<-FALSE
 min_positive<-apply(q,1,function(x){x<-x[!is.na(x)&x>0];if(length(x))min(x) else NA_real_})
 max_positive<-apply(q,1,function(x){x<-x[!is.na(x)&x>0];if(length(x))max(x) else NA_real_})
 observed<-!is.na(d$edu_high)
 rows[[mod]]<-data.table(modality=mod,N=nrow(d),all_qualifications_none=sum(none),none_and_positive=sum(none&positive),
  old_first_missing=sum(is.na(current)),derived_missing=sum(is.na(derived)),
  coding_different=sum(!equal),both_available_different=sum(!is.na(current)&!is.na(derived)&current!=derived),
  first_missing_derived_available=sum(is.na(current)&!is.na(derived)),
  first_available_derived_missing=sum(!is.na(current)&is.na(derived)),
  edu_high_matches_min_positive=sum(d$edu_high[observed]==min_positive[observed],na.rm=TRUE),
  edu_high_matches_max_positive=sum(d$edu_high[observed]==max_positive[observed],na.rm=TRUE),
  edu_high_available=sum(observed))
 cross[[mod]]<-data.table(current,derived)[,.(N=.N),by=.(current,derived)][,modality:=mod]
}
fwrite(rbindlist(rows),file.path(out,'education_coding_summary.csv'))
fwrite(rbindlist(cross),file.path(out,'education_coding_crosstab.csv'))
stopifnot(identical(h,vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))))
write_json(list(status='EXECUTED',source_hashes=as.list(h),individual_rows_returned=FALSE,sources_unchanged=TRUE,
 current='edu_q0, -7 recoded 7, other negative missing',
 comparator='saved edu_high positive categories; -7 in any six qualification arrays recoded 7',
 model_change='None; this is a coding comparison, not a covariate-policy decision'),
 file.path(out,'education_coding_QA.json'),pretty=TRUE,auto_unbox=TRUE)
print(rbindlist(rows))
