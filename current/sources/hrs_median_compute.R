#!/usr/bin/env Rscript
base<-Sys.getenv('DISCO_RUN_OUTPUT')
parent<-dirname(base)
source(file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/hrs_functions.R'))
dir.create(file.path(base,'numeric'),recursive=TRUE,showWarnings=FALSE)
p<-readRDS(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/04_external_validation/HRS_imputed_reference/20261001/numeric/prepared_INTERNAL.rds'))
# User-requested simple deterministic imputation. NO MICE outputs are read.
median_fill<-function(X,train){
 med<-vapply(X,function(v)median(v[train],na.rm=TRUE),0.)
 stopifnot(all(is.finite(med)))
 a<-as.matrix(X);obs<-is.finite(a)
 for(j in seq_len(ncol(a)))a[!obs[,j],j]<-med[j]
 stopifnot(all(is.finite(a)),identical(unname(a[obs]),unname(as.matrix(X)[obs])))
 list(X=a,medians=med)
}
train<-!p$id%in%p$cid
im<-median_fill(p$X,train)
r<-make_scores(im$X,p,train,mode='full',verify=TRUE);d<-r$data;maps<-metric_map(p)
old<-readRDS(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/04_external_validation/HRS_tracker_reanalysis/20261001/numeric/endpoints_INTERNAL.rds'))
oi<-match(old$data$person_id,d$person_id)
stopifnot(length(oi)==3923,!anyNA(oi),max(abs(d$NHANES18_DISCO[oi]-old$data$DISCO))<1e-9,
 max(abs(d$NHANES18_DM[oi]-old$data$DM))<1e-9)
rows<-list()
for(m in names(maps))for(mod in c('crude','adj','accel','age_sex')){
 if(m=='Age'&&mod!='crude')next
 a<-fit_c(d,maps[[m]],mod)
 rows[[length(rows)+1L]]<-data.frame(metric=m,model=mod,N=nrow(d),deaths=sum(d$status),t(a),check.names=FALSE)
}
tab<-do.call(rbind,rows);tab$SE<-sqrt(tab$Cvar)
write.csv(tab,file.path(base,'numeric/Cindex_all_models.csv'),row.names=FALSE)
write.csv(r$checks,file.path(base,'numeric/score_verification.csv'),row.names=FALSE)
write.csv(data.frame(field=names(im$medians),median_analysis_scale=im$medians),file.path(base,'numeric/imputation_medians.csv'),row.names=FALSE)
saveRDS(list(data=d,table=tab,maps=maps,specs=r$specs,checks=r$checks,medians=im$medians),
 file.path(base,'numeric/models_INTERNAL.rds'))

# Complete-case and missingness diagnostics preserve the same imputed reference
# and scores; they change ONLY the evaluated sample, not feature selection.
miss<-rowMeans(is.na(p$X[p$eval_index,p$features]))
hasA1c<-!is.na(p$X$HbA1c[p$eval_index])
subsets<-list(All=rep(TRUE,nrow(d)),Previous35_complete=miss==0,
 At_most20pct_missing=miss<=.2,At_most50pct_missing=miss<=.5,HbA1c_observed=hasA1c)
sens<-list()
keymetrics<-c('NHANES18_DISCO','NHANES18_DM','HRS18_DISCO','HRS18_DM','HRS35_DISCO','HRS35_DM','HRS36_DISCO','HRS36_DM','GrimAge')
for(s in names(subsets))for(m in keymetrics){v<-d[subsets[[s]],];a<-fit_c(v,maps[[m]],'adj')
 sens[[length(sens)+1L]]<-data.frame(subset=s,metric=m,N=nrow(v),deaths=sum(v$status),C=a['C'])}
write.csv(do.call(rbind,sens),file.path(base,'numeric/missingness_sensitivity.csv'),row.names=FALSE)
flow<-p$flow
flow$imputation<-'Single deterministic per-marker median on analysis scale; no MICE, no age/sex/outcome regression'
flow$median_fit_pool_N<-sum(train);flow$full_fit_imputation_rule<-'Medians from 5855 non-methylation HRS participants; shared imputation for NHANES18 and HRS-reference comparisons'
flow$eval_missing_HbA1c<-sum(!hasA1c);flow$eval_missing_HbA1c_percent<-mean(!hasA1c)*100
flow$eval_more_than50pct_missing_primary<-sum(miss>.5)
flow$reference_persons_1954_or1953<-'1954 reference persons; leave the scored person out when present (1953 for overlapping persons)'
flow$unread_encrypted_supplements_not_imputed<-TRUE
jsonlite::write_json(flow,file.path(base,'numeric/sample_flow.json'),pretty=TRUE,auto_unbox=TRUE)
write.csv(read.csv(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/04_external_validation/HRS_imputed_reference/20261001/numeric/missingness_inventory.csv')),file.path(base,'numeric/missingness_inventory.csv'),row.names=FALSE)
print(tab[tab$model=='adj'&tab$metric%in%keymetrics,c('metric','N','deaths','C','SE')])

# 5-fold household-held-out validation. Fit medians anew in each fold, then
# recompute HRS reference, age weights, residuals and Cox models without test HHs.
cv<-list();cvcheck<-list()
for(k in 1:5){
 training<-is.na(p$bloodfold)|p$bloodfold!=k
 stopifnot(!any(training[p$eval_index[p$fold==k]]),all(training[p$eval_index[p$fold!=k]]))
 ix<-median_fill(p$X,training)
 rr<-make_scores(ix$X,p,training,mode=paste0('cv',k),verify=TRUE)
 v<-rr$data;rr$checks$fold<-k;cvcheck[[k]]<-rr$checks
 for(m in c(keymetrics,'Age'))for(mod in c('crude','adj','accel')){
  if(m=='Age'&&mod!='crude')next
  tr<-v[p$fold!=k,];te<-v[p$fold==k,];tr$s<-tr[[maps[[m]]]];te$s<-te[[maps[[m]]]]
  if(mod=='accel'){lmage<-lm(s~age,data=tr);tr$s<-resid(lmage);te$s<-te$s-predict(lmage,newdata=te)}
  mu<-mean(tr$s);sg<-sd(tr$s);tr$s<-(tr$s-mu)/sg;te$s<-(te$s-mu)/sg
  fm<-if(mod=='adj')Surv(time,status)~s+age else Surv(time,status)~s
  fit<-coxph(fm,data=tr,ties='efron');te$lp<-predict(fit,newdata=te,type='lp')
  cc<-concordance(Surv(time,status)~lp,data=te,reverse=TRUE)
  cv[[length(cv)+1L]]<-data.frame(fold=k,metric=m,model=mod,N=nrow(te),deaths=sum(te$status),C=cc$concordance,
   numerator=unname(cc$count['concordant']+.5*cc$count['tied.x']),denominator=sum(cc$count[c('concordant','discordant','tied.x')]))
 }
 cat('Median-imputation full-pipeline CV fold',k,'done\n')
}
cvtab<-do.call(rbind,cv);agg<-aggregate(cbind(numerator,denominator)~metric+model,data=cvtab,sum);agg$CV_C<-agg$numerator/agg$denominator
write.csv(cvtab,file.path(base,'numeric/CV_each_fold.csv'),row.names=FALSE)
write.csv(agg,file.path(base,'numeric/CV_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,cvcheck),file.path(base,'numeric/CV_score_verification.csv'),row.names=FALSE)
cat('Median-only computation and validation complete\n');print(agg[agg$model=='adj',c('metric','model','CV_C')])
