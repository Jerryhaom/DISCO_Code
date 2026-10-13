#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(survival);library(glmnet);library(Rcpp);library(igraph);library(ggplot2);library(ggrepel);library(ggpubr);library(jsonlite);library(digest)})
args<-commandArgs(TRUE); module<-if(length(args))normalizePath(args[1])else normalizePath('05_organ_network/Fig4c_raw_reconstruction_20261008');root<-Sys.getenv('DISCO_INPUT_ROOT')
paths<-c(raw_proteomics=file.path(root,'UKB_Data/20240321_prot_instance0.RData'),raw_clinical=file.path(root,'UKB_Data/UKB.clinical.biomarker.rds'),metadata=file.path(root,'UKB_Data/UKB.Mortality.csv'),saved_clean_clinical=file.path(root,'Result1/Clean.UKB.Clinical.rds'),saved_selected_protein=file.path(root,'Result2/Selected.UKB.Prot.rds'),saved_protein_fit=file.path(root,'Result2/Lasso.model.prot.rds'),saved_joint_fit=file.path(root,'Result2/Lasso.model.ClinicalProt.rds'),saved_joint_score=file.path(root,'Result2/UKB.ClinicalProt.DISCO.rds'),coding=file.path(root,'Result2/Proteins/coding143.tsv'))
sha<-function(p)digest(file=p,algo='sha256');hashes<-vapply(paths,sha,'');tic<-proc.time()[3];cat('Reading actual raw inputs\n');flush.console()
meta<-read.csv(paths['metadata']);stopifnot(!anyDuplicated(meta$id))
rawc<-readRDS(paths['raw_clinical']);keep<-rowSums(is.na(rawc))<.5*ncol(rawc);cl<-as.data.frame(rawc[keep,,drop=FALSE]);rm(rawc)
median_impute<-function(x){if(!is.numeric(x))stop('Unexpected nonnumeric measurement');x[is.na(x)]<-median(x,na.rm=TRUE);x}
cl[]<-lapply(cl,median_impute);names(cl)<-make.names(names(cl),unique=TRUE)
# Actual value-matched producer is Raw/Result1/Cal_UKB_DISCO.R: no scaling
# or clipping before this saved clean input. The newer Result1 same-name file
# has an additional HD preprocessing branch which does NOT match the object.
names(cl)[1]<-'sampleID';savedc<-readRDS(paths['saved_clean_clinical']);stopifnot(identical(names(cl),names(savedc)),identical(as.character(cl$sampleID),as.character(savedc$sampleID)))
clinical_delta<-as.matrix(cl[,-1])-as.matrix(savedc[,-1]);clinical_diff<-max(abs(clinical_delta));write.csv(data.frame(marker=names(cl)[-1],max_abs_difference=apply(abs(clinical_delta),2,max),mean_abs_difference=colMeans(abs(clinical_delta))),file.path(module,'results/raw_clinical_QC_saved_comparison.csv'),row.names=FALSE);stopifnot(clinical_diff<1e-10);rm(savedc,clinical_delta);gc();cat('Clinical QC/median-imputation matches value-confirmed Raw producer:',clinical_diff,'\n');flush.console()
e<-new.env();load(paths['raw_proteomics'],envir=e);prot<-e$prot.dat;rm(e)
ns<-rowSums(is.na(prot));nc<-colSums(is.na(prot));rowkeep<-ns<.5*ncol(prot);colkeep<-nc<.1*nrow(prot);prot<-as.data.frame(prot[rowkeep,colkeep,drop=FALSE]);prot[]<-lapply(prot,median_impute);names(prot)<-make.names(names(prot),unique=TRUE)
mi<-match(prot$eid,meta$id);stopifnot(!anyNA(mi));info<-meta[mi,,drop=FALSE];pos<-which(info$time>0);info<-info[pos,,drop=FALSE];prot<-prot[match(info$id,prot$eid),,drop=FALSE]
stopifnot(!anyDuplicated(prot$eid),all(is.finite(as.matrix(prot))),identical(as.character(prot$eid),as.character(info$id)));cat('Protein raw QC dimensions:',nrow(prot),'x',ncol(prot)-1L,'\n');flush.console()
oldpfit<-readRDS(paths['saved_protein_fit']);oldjfit<-readRDS(paths['saved_joint_fit']);oldp<-readRDS(paths['saved_selected_protein']);oldselected<-rownames(as.matrix(coef(oldpfit,s='lambda.min')))[as.matrix(coef(oldpfit,s='lambda.min'))[,1]!=0]
savedp_idx<-match(oldp$eid,prot$eid);stopifnot(!anyNA(savedp_idx),all(names(oldp)[-1] %in% names(prot)));protein_diff<-max(abs(as.matrix(oldp[,-1])-as.matrix(prot[savedp_idx,names(oldp)[-1],drop=FALSE])));stopifnot(protein_diff<1e-10)
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(134);foldp<-sample(rep(seq_len(5),length.out=nrow(prot)))
cat('Starting transparent fixed-seed protein LASSO CV (historical protein fold IDs are not stored)\n');flush.console()
fitp<-cv.glmnet(as.matrix(prot[,-1]),Surv(info$time,info$status),type.measure='C',family='cox',nfolds=5,foldid=foldp,cox.ties='breslow')
bp<-rownames(as.matrix(coef(fitp,s='lambda.min')))[as.matrix(coef(fitp,s='lambda.min'))[,1]!=0];cat('New selected proteins:',length(bp),'lambda',fitp$lambda.min,'\n');flush.console()
common<-cl$sampleID[cl$sampleID %in% prot$eid];ci<-match(common,cl$sampleID);pi<-match(common,prot$eid);im<-match(common,meta$id);stopifnot(!anyNA(ci),!anyNA(pi),!anyNA(im));infojoint<-meta[im,,drop=FALSE]
x<-as.matrix(data.frame(cl[ci,-1,drop=FALSE],prot[pi,bp,drop=FALSE]));stopifnot(all(infojoint$time>0),all(is.finite(x)))
set.seed(134);foldj<-sample(rep(seq_len(5),length.out=nrow(x)));cat('Starting joint clinical/protein LASSO CV:',nrow(x),'x',ncol(x),'\n');flush.console()
fitj<-cv.glmnet(x,Surv(infojoint$time,infojoint$status),type.measure='C',family='cox',nfolds=5,foldid=foldj,cox.ties='breslow')
bj<-rownames(as.matrix(coef(fitj,s='lambda.min')))[as.matrix(coef(fitj,s='lambda.min'))[,1]!=0];selected_x<-x[,bj,drop=FALSE];age<-infojoint$Age;sex<-infojoint$Gender;ri<-which(age<=40)
cc<-apply(selected_x,2,function(z)cor(age,z));w<-abs(cc%o%cc);diag(w)<-0;w<-w/sum(w);ref<-selected_x[ri,,drop=FALSE]
# Exact pair update already independently validated for the same R algorithm.
sourceCpp(file.path(Sys.getenv('DISCO_RELEASE_ROOT'),'modules/ed4/code/score_kernel.cpp'),cacheDir=file.path(module,'compile_cache'))
sc<-score_two_weights(selected_x,colMeans(ref),cov(ref),cor(ref),w,w,nrow(ref))[,1]
literal<-vapply(seq_len(5),function(i)log(nrow(ref)^2*sum(w*(cor(rbind(ref,selected_x[i,]))-cor(ref))^2)),numeric(1));score_error<-max(abs(literal-sc[1:5]));stopifnot(score_error<1e-8,all(is.finite(sc)))
design<-model.matrix(~age+sex);res<-qr.resid(qr(design),cbind(sc,selected_x));partial<-as.numeric(cor(res[,1],res[,-1,drop=FALSE]));pval<-vapply(seq_along(bj),function(i)cor.test(res[,1],res[,i+1])$p.value,numeric(1))
pn<-read.delim(paths['coding'],check.names=FALSE);pnl<-unlist(strsplit(pn$meaning,';'));pn$prot<-pnl[seq(1,length(pnl),2)];labels<-bj;isprot<-bj %in% bp;labels[isprot]<-pn$prot[match(as.numeric(sub('^X','',bj[isprot])),pn$coding)];stopifnot(!anyNA(labels),!anyDuplicated(labels))
cm<-abs(cor(selected_x));dimnames(cm)<-list(labels,labels);ed<-which(cm>.2 & cm<1,arr.ind=TRUE);edges<-data.frame(from=rownames(cm)[ed[,1]],to=colnames(cm)[ed[,2]],weight=cm[ed]);g<-graph_from_data_frame(edges,directed=TRUE);de<-degree(g,mode='all');tab<-data.frame(feature=bj,label=labels,partial_r=partial,P_residual_cor_test_Nminus2=pval,degree=unname(de[labels]));ok<-is.finite(tab$degree)
test<-cor.test(abs(tab$partial_r[ok]),tab$degree[ok],method='spearman');top<-tab$label[order(tab$degree,decreasing=TRUE)][1:min(10,sum(ok))];tab$display_label<-ifelse(tab$label %in% top,tab$label,'')
plotdata<-data.frame(dg=tab$degree,cc=tab$partial_r,lb1=tab$display_label)
pdf(file.path(module,'figures/Fig4c_TRANSPARENT_RAW_QC_LASSO_RECONSTRUCTION.pdf'),height=4,width=5)
print(ggplot(plotdata,aes(dg,abs(cc)))+geom_point(color='grey',size=2)+geom_smooth(method='lm',color='black')+ylim(0,.65)+xlim(0,350)+xlab('Degree in Correlation Network')+ylab('Partial correlation with DISCO')+stat_cor(method='spearman',cor.coef.name='rho',label.y.npc=.95)+geom_text_repel(aes(label=lb1),size=3,max.overlaps=15,seed=134)+labs(title='Relationship between biomarker correlations with\n\n       DISCO and biomarker network connectivity')+theme_classic()+theme(axis.title=element_text(size=14),plot.title=element_text(size=14,lineheight=.5)))
dev.off()
saveRDS(list(protein_fit=fitp,joint_fit=fitj,eids=common,selected_features=bj,score=sc,protein_fold_id=foldp,joint_fold_id=foldj,seed=134,RNGkind=RNGkind()),file.path(module,'input/private/rebuilt_models_and_scores.rds'));Sys.chmod(file.path(module,'input/private/rebuilt_models_and_scores.rds'),'0600')
write.csv(tab,file.path(module,'tables/Fig4c_new_statistics.csv'),row.names=FALSE);write.csv(edges,file.path(module,'results/new_graph_edges.csv'),row.names=FALSE)
oldb<-rownames(as.matrix(coef(oldjfit,s='lambda.min')))[as.matrix(coef(oldjfit,s='lambda.min'))[,1]!=0]
qa<-list(status='PASS',historical_full_pipeline_recovered=FALSE,branch='Transparent complete raw QC/imputation -> fixed-seed 5-fold protein LASSO -> fixed-seed joint LASSO -> age-weighted DISCO -> age+sex partial correlations -> newly recomputed graph/scatter',historical_model_objects_not_modified=TRUE,protein_original_fold_assignment='Not stored or recovered; seed134 chosen transparently consistent with joint original code',joint_original_seed=134,clinical_raw_QC_max_error=clinical_diff,protein_raw_QC_saved_selected_max_error=protein_diff,protein_N=nrow(prot),protein_predictors=ncol(prot)-1L,selected_proteins=length(bp),historical_selected_proteins=length(oldselected),protein_lambda_min_new=fitp$lambda.min,protein_lambda_min_historical=oldpfit$lambda.min,joint_N=nrow(x),selected_joint_features=length(bj),historical_joint_features=length(oldb),joint_feature_set_identical=setequal(bj,oldb),reference_n=nrow(ref),literal_score_max_error=score_error,spearman_rho=unname(test$estimate),P=test$p.value,graph_nodes=vcount(g),graph_directed_edges=ecount(g),input_sha256=as.list(hashes),input_unchanged=identical(hashes,vapply(paths,sha,'')),elapsed_sec=unname(proc.time()[3]-tic),R_version=R.version.string,glmnet=as.character(packageVersion('glmnet')))
stopifnot(qa$input_unchanged);write_json(qa,file.path(module,'results/numerical_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17);print(qa);writeLines(capture.output(sessionInfo()),file.path(module,'logs/sessionInfo.txt'))
