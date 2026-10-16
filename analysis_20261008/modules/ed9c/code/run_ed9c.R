#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(pheatmap);library(RColorBrewer);library(jsonlite);library(digest)})
args<-commandArgs(TRUE);module<-if(length(args))normalizePath(args[1])else normalizePath('05_organ_network/ED9c_partial_correlation_20261008')
root<-Sys.getenv('DISCO_INPUT_ROOT');org<-file.path(root,'Result3/Organ')
paths<-c(raw_clinical=file.path(root,'UKB_Data/UKB.clinical.biomarker.rds'),protein_order=file.path(org,'UKB.prot.imputed.rds'),organ_scores=file.path(org,'UKB.tissue.mc2.rds'),metadata=file.path(root,'UKB_Data/UKB.Mortality.csv'),organ_order=file.path(org,'Tissue.order.rds'),field_map=file.path(root,'Results_final/05_organ_network/ED9c_correlation_audit_20261006/field_mapping.csv'))
sha<-function(p)digest(file=p,algo='sha256');hashes<-vapply(paths,sha,'')
raw<-readRDS(paths['raw_clinical']);protein<-readRDS(paths['protein_order']);ids<-as.character(protein[,'eid']);rm(protein);meta<-read.csv(paths['metadata']);score<-as.matrix(readRDS(paths['organ_scores']));ord<-readRDS(paths['organ_order']);map<-read.csv(paths['field_map'],check.names=FALSE)
stopifnot(nrow(score)==length(ids),ncol(score)==18,!anyDuplicated(ids),!anyDuplicated(raw$id),!anyDuplicated(meta$id),setequal(ord,colnames(score)))
ri<-match(ids,as.character(raw$id));mi<-match(ids,as.character(meta$id));stopifnot(!anyNA(mi));age<-meta$Age[mi];sex<-meta$Gender[mi]
selected<-which(map$selected_ED9c);stopifnot(length(selected)==25,all(map$raw_field[selected] %in% colnames(raw)))
x<-as.matrix(raw[ri,map$raw_field[selected],drop=FALSE]);rm(raw,meta);gc()
rows<-list();inverse_errors<-numeric();p_errors<-numeric()
for(i in seq_len(ncol(x))) {
  ok<-complete.cases(x[,i],age,sex,score);N<-sum(ok);design<-model.matrix(~age+factor(sex),data=data.frame(age=age[ok],sex=sex[ok]));q<-qr(design);stopifnot(q$rank==3L)
  res<-qr.resid(q,cbind(x[ok,i],score[ok,,drop=FALSE]));rr<-as.numeric(cor(res[,1],res[,-1,drop=FALSE]));df<-N-4L;logP<-log(2)+pt(-abs(rr)*sqrt(df/(1-rr^2)),df,log.p=TRUE);P<-exp(logP)
  # Independent covariance-inversion partial correlation and F-test P.
  independent_r<-vapply(seq_len(ncol(score)),function(j){z<-cbind(x[ok,i],score[ok,j],age[ok],as.numeric(factor(sex[ok])));inv<-solve(cor(z));-inv[1,2]/sqrt(inv[1,1]*inv[2,2])},numeric(1))
  independent_logP<-pf(rr^2*df/(1-rr^2),1,df,lower.tail=FALSE,log.p=TRUE)
  inverse_errors<-c(inverse_errors,rr-independent_r);p_errors<-c(p_errors,logP-independent_logP)
  marker<-map$saved_label[selected[i]];display<-switch(marker,BUN='Urea',Chol='CHOL',Lp='Lp(a)',marker)
  rows[[i]]<-data.frame(biomarker=marker,display_label=display,raw_field=map$raw_field[selected[i]],ukb_field_id=map$ukb_field_id[selected[i]],organ=colnames(score),N=N,df=df,signed_r=rr,absolute_r=abs(rr),P=P,logP=logP)
}
tab<-do.call(rbind,rows);tab$q_BH_450<-p.adjust(tab$P,'BH');tab$nominal_star<-ifelse(tab$P<.01,'**',ifelse(tab$P<.05,'*',''))
stopifnot(nrow(tab)==450,!anyDuplicated(tab[,c('biomarker','organ')]),max(abs(inverse_errors))<1e-10,max(abs(p_errors))<1e-7)
labels<-vapply(rows,function(d)d$display_label[1],'');r<-p<-matrix(NA_real_,25,18,dimnames=list(labels,colnames(score)))
for(i in seq_along(rows)){r[i,]<-rows[[i]]$absolute_r;p[i,]<-rows[[i]]$P}
colour<-t(r)[ord,,drop=FALSE];numbers<-t(ifelse(p<.01,'**',ifelse(p<.05,'*','')))[ord,,drop=FALSE];stopifnot(identical(dimnames(colour),dimnames(numbers)))
groups<-c(rep('Immune',4),rep('Oxygen',2),rep('Liver',6),rep('Kidney',4),rep('Lipid',7),rep('Glucose',2));ann<-data.frame(Biomarkers=groups,row.names=labels);cols<-brewer.pal(6,'Set2');names(cols)<-unique(groups)
# Preserve the actual original producer's pdf() device. This host's optional
# cairo device is unavailable (X11 runtime libraries are not installed).
pdf(file.path(module,'figures/ED9c_AGE_SEX_PARTIAL_CORRELATION_NEW_CALCULATION.pdf'),height=3.5,width=6)
pheatmap(colour,annotation_col=ann,annotation_colors=list(Biomarkers=cols),cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=numbers,number_color='black',fontsize_number=8,color=colorRampPalette(brewer.pal(7,'Reds'))(50),treeheight_row=15)
dev.off()
stopifnot(file.exists(file.path(module,'figures/ED9c_AGE_SEX_PARTIAL_CORRELATION_NEW_CALCULATION.pdf')))
write.csv(tab,file.path(module,'tables/ED9c_age_sex_450_cells.csv'),row.names=FALSE);write.csv(r,file.path(module,'tables/ED9c_absolute_partial_r.csv'));write.csv(p,file.path(module,'tables/ED9c_nominal_P.csv'))
saveRDS(list(raw=x,organ_scores=score,age=age,sex=sex),file.path(module,'input/private/ED9c_prepared_input.rds'));Sys.chmod(file.path(module,'input/private/ED9c_prepared_input.rds'),'0600')
oldtab<-read.csv(file.path(root,'Results_final/05_organ_network/ED9c_correlation_audit_20261006/pass_1/age_sex_candidate_450_cells.csv'))
cmp<-merge(tab,oldtab,by=c('biomarker','organ'),suffixes=c('_fresh','_accepted'));stopifnot(nrow(cmp)==450)
qa<-list(status='PASS',historical_r_P_source_recovered=FALSE,scope='Fresh original raw clinical input + saved organ-score input to explicitly specified age+sex partial r/P; not historical r/P or complete organ-score training recovery',cells=450,adjustments='Both variables linearly residualized against age + sex using identical cellwise complete samples',imputation='None',winsorization='None',log_transform='None',df='N-4',N_range=range(tab$N),signed_r_range=range(tab$signed_r),negative_correlations=sum(tab$signed_r<0),nominal_P_lt05=sum(tab$P<.05),nominal_P_lt01=sum(tab$P<.01),independent_inverse_correlation_max_error=max(abs(inverse_errors)),independent_F_test_logP_max_error=max(abs(p_errors)),accepted_20261006_new_candidate_r_max_error=max(abs(cmp$signed_r_fresh-cmp$signed_r_accepted)),accepted_20261006_new_candidate_logP_max_error=max(abs(cmp$logP_fresh-cmp$logP_accepted)),label_correction='Raw U, UKB field 30670, is Urea; not BUN',fill='Absolute partial Pearson r; signed r retained in machine-readable table',stars='Nominal P: **<.01; *<.05. BH q in table, not star layer.',input_sha256=as.list(hashes),input_unchanged=identical(hashes,vapply(paths,sha,'')))
stopifnot(qa$input_unchanged);write_json(qa,file.path(module,'results/numerical_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17);writeLines(capture.output(sessionInfo()),file.path(module,'logs/sessionInfo.txt'));print(qa)
