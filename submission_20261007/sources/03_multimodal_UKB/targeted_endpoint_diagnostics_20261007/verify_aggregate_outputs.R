#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Independent checks of downloaded, deidentified aggregate statistical outputs.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
setDTthreads(1)
root<-.disco_path("results","");out<-file.path(root,'03_multimodal_UKB/targeted_endpoint_diagnostics_20261007')
d<-fread(file.path(out,'three_modalities_raw_unit_models.csv'))
z<-fread(file.path(out,'three_modalities_perSD_models.csv'))
stopifnot(nrow(d)==216L,nrow(z)==144L,all(is.finite(d$beta)),all(d$N>=d$events),
          all(d$low<=d$HR & d$HR<=d$high))
cc<-fread(file.path(out,'three_modalities_new_fit_coefficients.csv'))
vv<-fread(file.path(out,'three_modalities_new_fit_full_vcov.csv'))
new<-d[grepl('new targeted fit',source)]
beta_error<-se_error<-0
for(i in seq_len(nrow(new))){
 r<-new[i];keys<-as.list(r[,.(modality,disease,score_column,model)])
 b<-cc[modality==keys$modality & disease==keys$disease & metric==keys$score_column & model==keys$model & term==keys$score_column]
 v<-vv[modality==keys$modality & disease==keys$disease & metric==keys$score_column & model==keys$model & term1==keys$score_column & term2==keys$score_column]
 stopifnot(nrow(b)==1,nrow(v)==1)
 beta_error<-max(beta_error,abs(b$beta-r$beta));se_error<-max(se_error,abs(sqrt(v$covariance)-r$se))
}
converted<-d[grepl('exact linear',source)];conv_error<-0
for(mod in c('clinical','proteomics')){
 prior<-fread(file.path(root,'03_multimodal_UKB/month_disease_prs_execution_20261007',paste0(mod,'_8/models.csv')))
 for(i in which(converted$modality==mod)){
  r<-converted[i];score_term<-if(r$metric=='DISCO')'zD' else 'zDM'
  q<-prior[disease==r$disease & model==r$model & term==score_term & policy=='death_first' & version=='raw_dates']
  stopifnot(nrow(q)==1,r$N==q$N,r$events==q$events)
  conv_error<-max(conv_error,abs(r$beta*r$SD-q$beta),abs(r$se*r$SD-q$se),abs(r$C-q$C))
 }
}
cif<-fread(file.path(out,'ClinicalProt_CIF_aggregate.csv'));aj<-fread(file.path(out,'ClinicalProt_CIF_QA.csv'))
stopifnot(nrow(cif)==1360L,nrow(aj)==80L,all(cif$CIF>=-1e-12 & cif$CIF<=1),all(cif$CIF<=cif$KM+1e-12),
          all(cif$CIF_low<=cif$CIF & cif$CIF<=cif$CIF_high))
j<-fread(file.path(out,'ovarian_strict_current_models.csv'))[term=='Joint_11_vs_00']
ov<-fread(file.path(out,'ovarian_strict_current_full_vcov.csv'))
co<-fread(file.path(out,'ovarian_diagnostic_coefficients.csv'))[run=='tight_11']
ts<-c('zD','zPRS','zD:zPRS');bj<-sum(co[term%in%ts]$beta)
vj<-sum(ov[term1%in%ts & term2%in%ts]$covariance)
stopifnot(abs(exp(bj)-j$HR)<1e-12,abs(sqrt(vj)-j$se)<1e-12)
edu<-fread(file.path(out,'education_coding_summary.csv'));stopifnot(all(edu$coding_different==0),all(edu$old_first_missing==edu$derived_missing))
old_root<-file.path(root,'03_multimodal_UKB/month_disease_prs_execution_20261007')
hashes<-fread(file.path(old_root,'artifact_hashes.csv'))
hashes[,current_sha256:=vapply(file.path(old_root,file),function(p)digest(file=p,algo='sha256'),character(1))]
hashes[,status:=ifelse(current_sha256==sha256,'unchanged',ifelse(file=='PROTOCOL.md','expected documented protocol-text repair','unexpected change'))]
stopifnot(!any(hashes$status=='unexpected change'))
fwrite(hashes,file.path(out,'original_execution_artifact_recheck.csv'))
result<-list(status='PASS',scope='Independent aggregate checks; not raw-data re-adjudication or full PH review',
             raw_unit_rows=nrow(d),perSD_rows=nrow(z),new_fit_rows=nrow(new),converted_rows=nrow(converted),
             new_coefficient_max_error=beta_error,new_SE_full_vcov_max_error=se_error,linear_conversion_max_error=conv_error,
             ClinicalProt_CIF_rows=nrow(cif),ClinicalProt_AJ_checks=nrow(aj),ClinicalProt_AJ_max_error=max(aj$max_AJ_error),
             ovarian_Joint_full_covariance_independent_PASS=TRUE,
             education_postrun_finding='No coding or missingness difference across all three modality populations; no covariate choice remains for this issue',
             prior_execution_artifact_files=nrow(hashes),prior_statistical_outputs_unchanged=TRUE,
             expected_text_repair='Only PROTOCOL.md differs from prior execution manifest; original manifest remains untouched')
write_json(result,file.path(out,'independent_aggregate_QA.json'),pretty=TRUE,auto_unbox=TRUE)
print(result)
