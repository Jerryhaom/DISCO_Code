#!/usr/bin/env Rscript
# Protein-set sensitivity of the actual ED10d models. Raw inputs are immutable.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE)))
root<-Sys.getenv('DISCO_DATA_ROOT');base<-file.path(root,'Results_final')
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(Rcpp);library(digest)})
setDTthreads(1);options(warn=2);Sys.umask('0077')
for(z in c('PRIVATE','numeric','qa','logs','output/pdf'))dir.create(file.path(out,z),recursive=TRUE,showWarnings=FALSE)
Sys.chmod(file.path(out,'PRIVATE'),'0700')
logmsg<-function(...){cat(format(Sys.time()),...,'\n');flush.console()}
old<-file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20261002_d_original');rebuild<-file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20261002_local_rebuild')
exact<-file.path(base,'05_organ_network/corrected_organ/20260926')
map<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
oh<-fread(file.path(rebuild,'results/Oh_local_coefficients.csv'));sets<-list();coefrows<-list()
for(o in map$clock_organ)for(m in c('Oh','Goeminne')){
 v<-if(m=='Oh')unlist(oh[organ==o,-c('organ','intercept')]) else unlist(fread(file.path(rebuild,'results',paste0(o,'_coefs_GTEx_4x_FC.csv')))[1,-1])
 v<-v[is.finite(v)&v!=0];stopifnot(!anyDuplicated(names(v)))
 sets[[paste(o,m,sep='_')]]<-names(v)
 coefrows[[paste(o,m)]]<-data.frame(organ=o,method=m,protein=names(v),coefficient=as.numeric(v),fold=if(m=='Oh')'local single LASSO' else 'local outer fold 1')
}
fwrite(rbindlist(coefrows),file.path(out,'numeric/selected_proteins.csv'))
des<-readRDS(file.path(exact,'locked_design.rds'));ids<-as.character(des$ids)
meta<-readRDS(file.path(rebuild,'PRIVATE/meta.rds'));meta<-meta[match(ids,meta$eid),]
stopifnot(identical(ids,as.character(meta$eid)))
ref<-which(ids %in% as.character(des$refids));stopifnot(length(ref)==615,identical(ref,which(meta$age<=40)))
orig<-readRDS(file.path(root,'Result3/Organ/UKB.prot.imputed.rds'))
stopifnot(identical(as.numeric(orig[,1]),as.numeric(des$ids)));rownames(orig)<-ids
cla<-readRDS(file.path(root,'Result3/Organ/UKB.Prot.class.rds'))
orgsets<-setNames(lapply(map$disco_organ,function(o)unique(cla$Gene[cla$specific.tissue==o])),map$clock_organ)
needed<-unique(c(unlist(sets),unlist(orgsets)))
logmsg('Read raw Olink, preserving original analysis IDs and reference')
e<-new.env();load(file.path(root,'UKB_Data/20240321_prot_instance0.RData'),envir=e)
p<-e$prot.dat;rm(e);coding<-fread(file.path(root,'Result2/Proteins/coding143.tsv'))
cn<-sub(';.*','',coding$meaning[match(colnames(p)[-(1:2)],as.character(coding$coding))])
stopifnot(!anyNA(cn),!anyDuplicated(cn),all(needed %in% cn),!anyDuplicated(p$eid),all(ids %in% as.character(p$eid)))
colnames(p)[-(1:2)]<-cn
X<-as.matrix(p[match(ids,as.character(p$eid)),needed,drop=FALSE]);rownames(X)<-ids;rm(p);gc(FALSE)
bad<-which(is.finite(X)&abs(X)>1e300,arr.ind=TRUE)
if(nrow(bad)){stopifnot(all(colnames(X)[bad[,2]]=='VWC2L'));X[bad]<-NA_real_}
stopifnot(!any(is.infinite(X)))
missing<-colSums(is.na(X));med<-apply(X,2,median,na.rm=TRUE);stopifnot(all(is.finite(med)))
for(j in seq_len(ncol(X)))X[is.na(X[,j]),j]<-med[j]
common<-intersect(needed,colnames(orig));err<-max(abs(X[,common,drop=FALSE]-orig[,common,drop=FALSE]))
stopifnot(err<1e-12,all(is.finite(X)))
logmsg('Original median-filled values reproduced for',length(common),'proteins; new proteins',length(setdiff(needed,common)))
fwrite(data.frame(protein=needed,imputed_N=unname(missing),median=unname(med),in_original_1553=needed %in% common),file.path(out,'numeric/preprocessing.csv'))
rm(orig);gc(FALSE)
sourceCpp(file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/exact_score.cpp'),cacheDir=file.path(out,'PRIVATE/cpp_cache'))
agecor<-as.vector(cor(X,meta$age));names(agecor)<-colnames(X)
mu<-colMeans(X[ref,,drop=FALSE]);ss<-apply(X[ref,,drop=FALSE],2,sd)
stopifnot(all(ss>0),all(is.finite(agecor)))
probe<-unique(c(1,ref[c(1,2,615)],seq(100,nrow(X),length.out=16)%/%1))
checks<-coverage<-list();scores<-list();original_scores<-readRDS(file.path(exact,'pass_1/organ36_exact_scores.rds'))
score<-function(pp,key){
 xx<-X[,pp,drop=FALSE];rr<-xx[ref,,drop=FALSE];cc<-cor(rr)
 w<-abs(outer(agecor[pp],agecor[pp]));diag(w)<-0;stopifnot(sum(w)>0);w<-w/sum(w)
 y<-organ_score(xx,mu[pp],ss[pp],cc,w,length(ref),FALSE)
 brute<-vapply(probe,function(i)log(length(ref)^2*sum(w*(cor(rbind(rr,xx[i,]))-cc)^2)),0.)
 maxerr<-max(abs(y[probe]-brute));stopifnot(all(is.finite(y)),maxerr<1e-9)
 checks[[key]]<<-data.frame(score=key,p=length(pp),direct_checks=length(probe),max_abs_error=maxerr)
 y
}
for(o in map$clock_organ){
 scores[[paste0(o,'_Original')]]<-score(orgsets[[o]],paste0(o,'_Original'))
 stopifnot(max(abs(scores[[paste0(o,'_Original')]]-original_scores[,map$disco_organ[match(o,map$clock_organ)]]))<1e-10)
 for(m in c('Oh','Goeminne')){
  key<-paste(o,m,sep='_');pp<-sets[[key]]
  coverage[[key]]<-data.frame(organ=o,method=m,p=length(pp),reference_N=length(ref),weight_N=nrow(X),available=length(pp)>=2,
    original_panel_p=length(orgsets[[o]]),overlap_original=length(intersect(pp,orgsets[[o]])),reason=if(length(pp)<2)'One selected protein; no off-diagonal correlation pairs' else '')
  if(length(pp)>=2)scores[[key]]<-score(pp,key)
 }
 logmsg('Exact scores verified',o)
}
fwrite(rbindlist(checks),file.path(out,'qa/direct_score_checks.csv'))
fwrite(rbindlist(coverage),file.path(out,'numeric/protein_set_counts.csv'))
scores<-do.call(cbind,scores);rownames(scores)<-ids
saveRDS(list(scores=scores,sets=sets,agecor=agecor,mu=mu,sd=ss,medians=med,ref_ids=ids[ref]),file.path(out,'PRIVATE/scores.rds'),compress=FALSE)
models<-readRDS(file.path(rebuild,'PRIVATE/common_prediction_models.rds'))
previous<-fread(file.path(old,'numeric/panel_d_Cindex.csv'));rows<-pairs<-association<-list();fits_private<-list();ir<-ip<-ih<-0L
for(o in map$clock_organ){
 item<-models[[o]];te<-item$test;fd<-item$fits$DISCO$model
 stopifnot(identical(rownames(fd),rownames(item$fits$Oh$model)),identical(rownames(fd),rownames(item$fits$Goeminne$model)))
 tr<-data.frame(time=fd[[1]][,1],status=fd[[1]][,2],age=fd$age,sex=fd$sex,D=fd$D,O=item$fits$Oh$model$O,G=item$fits$Goeminne$model$G,row.names=rownames(fd))
 stopifnot(nrow(tr)==18648,nrow(te)==4294,sum(te$status)==460,!length(intersect(rownames(tr),te$eid)))
 for(m in c('Oh','Goeminne'))if(paste(o,m,sep='_') %in% colnames(scores)){
  k<-if(m=='Oh')'DO' else 'DG';tr[[k]]<-scores[rownames(tr),paste(o,m,sep='_')];te[[k]]<-scores[as.character(te$eid),paste(o,m,sep='_')]
 }
 terms<-c(DISCO='D',DISCO_Oh='DO',DISCO_Goeminne='DG',Oh='O',Goeminne='G');terms<-terms[terms %in% names(tr)]
 fs<-list();fit<-function(f)coxph(f,data=tr,ties='efron',x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-10,iter.max=50))
 for(m in names(terms)){
  k<-terms[[m]];fs[[m]]<-fit(reformulate(c(if(grepl('^DISCO',m))'age',k),'Surv(time,status)'))
  cal<-lm(reformulate('age',k),tr);den<-sd(resid(cal));an<-paste0(k,'_accel')
  tr[[an]]<-resid(cal)/den;te[[an]]<-(te[[k]]-predict(cal,newdata=te))/den
  stopifnot(abs(cor(tr[[an]],tr$age))<1e-10)
 }
 for(m in names(terms))fs[[paste0(m,'_accel')]]<-fit(reformulate(paste0(terms[[m]],'_accel'),'Surv(time,status)'))
 cc<-do.call(concordance,c(unname(fs),list(newdata=te,timewt='n',influence=2)))
 cv<-cc$concordance;V<-cc$var;stopifnot(min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
 for(j in seq_along(fs)){
  meth<-sub('_accel$','',names(fs)[j]);grp<-if(grepl('_accel$',names(fs)[j]))'Accel' else 'Original'
  lp<-predict(fs[[j]],newdata=te,type='lp');direct<-concordance(Surv(te$time,te$status)~lp,reverse=TRUE)
  stopifnot(abs(cv[j]-direct$concordance)<1e-12,abs(V[j,j]-direct$var)<1e-12)
  if(meth %in% c('DISCO','Oh','Goeminne')){
   q<-previous[organ==o&group==grp&score==meth];stopifnot(nrow(q)==1,abs(cv[j]-q$C)<1e-12,abs(sqrt(V[j,j])-q$SE)<1e-12)
  }
  se<-sqrt(V[j,j]);ir<-ir+1L;rows[[ir]]<-data.frame(organ=o,display=map$display[match(o,map$clock_organ)],group=grp,score=meth,N=nrow(te),deaths=sum(te$status),C=unname(cv[j]),SE=se,lower95=max(0,cv[j]-qnorm(.975)*se),upper95=min(1,cv[j]+qnorm(.975)*se),formula=paste(deparse(formula(fs[[j]])),collapse=''),available=TRUE)
 }
 if(!'DISCO_Oh' %in% names(terms))for(grp in c('Original','Accel')){
  ir<-ir+1L;rows[[ir]]<-data.frame(organ=o,display=map$display[match(o,map$clock_organ)],group=grp,score='DISCO_Oh',N=nrow(te),deaths=sum(te$status),C=NA_real_,SE=NA_real_,lower95=NA_real_,upper95=NA_real_,formula='Undefined with one selected protein',available=FALSE)
 }
 for(grp in c('Original','Accel'))for(m in c('Oh','Goeminne')){
  suffix<-if(grp=='Accel')'_accel' else '';a<-paste0('DISCO_',m,suffix)
  if(!a %in% names(fs))next
  for(b in c(paste0(m,suffix),paste0('DISCO',suffix))){
   ai<-match(a,names(fs));bi<-match(b,names(fs));delta<-cv[ai]-cv[bi];se<-sqrt(max(0,V[ai,ai]+V[bi,bi]-2*V[ai,bi]));ip<-ip+1L
   pairs[[ip]]<-data.frame(organ=o,group=grp,contrast=paste(sub('_accel$','',a),'minus',sub('_accel$','',b)),delta=unname(delta),SE=se,lower=unname(delta)-qnorm(.975)*se,upper=unname(delta)+qnorm(.975)*se,P=2*pnorm(-abs(delta/se)))
  }
 }
 fits_private[[o]]<-list(fits=fs,test=te,covariance=V)
 logmsg('Cox and direct concordance checks passed',o)
}
result<-rbindlist(rows);paired<-rbindlist(pairs);paired[,q_BH:=p.adjust(P,'BH'),by=.(group,contrast)]
fwrite(result,file.path(out,'numeric/panel_d_matched_Cindex.csv'));fwrite(paired,file.path(out,'numeric/paired_differences.csv'))
saveRDS(fits_private,file.path(out,'PRIVATE/models.rds'))
write_json(list(original_d_all_54_C_and_SE_identical=TRUE,exact_score_kernel_checks=length(checks),one_protein_not_defined='Oh Heart',raw_imputation_max_error=err,original_existing_proteins_verified=length(common),additional_proteins=length(setdiff(needed,common)),invalid_VWC2L_cells_set_missing=nrow(bad),N=length(ids),reference_N=length(ref),reference_test_overlap=sum(as.character(des$refids) %in% as.character(models[[1]]$test$eid)),test_N=4294,test_deaths=460,train_N=18648,train_deaths=1950,primary_models_sex_adjusted=FALSE,age_weight_source='Same 51922-person original cohort; recompute and normalize per protein set',imputation='Per-protein median in the fixed original 51922-person cohort; no other transformations',strict_full_pipeline_independent_test=FALSE,original_source_files_untouched=TRUE),file.path(out,'qa/checks.json'),auto_unbox=TRUE,pretty=TRUE)
sources<-c(file.path(root,'UKB_Data/20240321_prot_instance0.RData'),file.path(root,'Result3/Organ/UKB.prot.imputed.rds'),file.path(rebuild,'results/Oh_local_coefficients.csv'),file.path(rebuild,'results',paste0(map$clock_organ,'_coefs_GTEx_4x_FC.csv')),file.path(rebuild,'PRIVATE/common_prediction_models.rds'),file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/exact_score.cpp'))
fwrite(data.frame(path=sources,sha256=vapply(sources,digest,'',file=TRUE,algo='sha256')),file.path(out,'qa/input_sha256.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'logs/sessionInfo.txt'))
Sys.chmod(list.files(file.path(out,'PRIVATE'),full.names=TRUE),'0600')
Sys.chmod(file.path(out,'PRIVATE/cpp_cache'),'0700')
logmsg('COMPLETE',sum(result$available),'estimated rows,',sum(!result$available),'explicitly undefined rows')
