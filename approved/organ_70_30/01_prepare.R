#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(entry))
root<-Sys.getenv('DISCO_DATA_ROOT');base<-file.path(root,'Results_final');old<-file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison'),'20261002_local_rebuild')
.libPaths(c(file.path(base,'09_figures_tables/Rlibs'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
setDTthreads(1);Sys.umask('0077');options(warn=1)
for(z in c('PRIVATE','qa','logs','numeric','panels','vendor'))dir.create(file.path(out,z),recursive=TRUE,showWarnings=FALSE)
Sys.chmod(file.path(out,'PRIVATE'),'0700');sink(file.path(out,'logs/01_prepare.log'),split=TRUE)
stopifnot(file.exists(file.path(out,'PROTOCOL.txt')))
orgs<-c('Adipose','Brain','Heart','Intestine','Kidney','Liver','Lung','Muscle','Pancreas')
meta_all<-readRDS(file.path(old,'PRIVATE/meta.rds'));meta_all$eid<-as.character(meta_all$eid);stopifnot(!anyDuplicated(meta_all$eid))
e<-new.env(parent=baseenv());load(file.path(root,'UKB_Data/20240321_prot_instance0.RData'),envir=e);p<-e$prot.dat;rm(e)
rawids<-as.character(p$eid);stopifnot(!anyDuplicated(rawids),all(p$ins_index==0),all(rawids%in%meta_all$eid))
meta_all<-meta_all[match(rawids,meta_all$eid),];stopifnot(identical(meta_all$eid,rawids))
coding<-fread(file.path(root,'Result2/Proteins/coding143.tsv'));cn<-sub(';.*','',coding$meaning[match(names(p)[-(1:2)],as.character(coding$coding))]);stopifnot(!anyNA(cn),!anyDuplicated(cn))
Xall<-as.matrix(p[,-(1:2),drop=FALSE]);dimnames(Xall)<-list(rawids,cn);rm(p);gc(FALSE)
bad<-which(is.finite(Xall)&abs(Xall)>1e300,arr.ind=TRUE);stopifnot(nrow(bad)==1,all(colnames(Xall)[bad[,2]]=='VWC2L'));Xall[bad]<-NA_real_;stopifnot(!any(is.infinite(Xall)))
miss_n<-rowSums(is.na(Xall));oh_pass<-miss_n<=1000;go_pass<-rowMeans(is.na(Xall))<.49948682860075
legal<-is.finite(meta_all$age)&meta_all$age>0&!is.na(meta_all$sex)&meta_all$sex%in%0:1&is.finite(meta_all$time)&meta_all$time>0&!is.na(meta_all$status)&meta_all$status%in%0:1
stage_masks<-list(raw_baseline=rep(TRUE,nrow(Xall)),Oh_participant_QC=oh_pass,both_participant_QC=oh_pass&go_pass,common_with_legal_endpoint=oh_pass&go_pass&legal)
flow<-rbindlist(lapply(seq_along(stage_masks),function(j){k<-stage_masks[[j]];prior<-if(j==1)k else stage_masks[[j-1]];data.frame(stage=names(stage_masks)[j],N=sum(k),deaths=sum(meta_all$status[k]==1,na.rm=TRUE),removed_from_previous=sum(prior&!k),removed_deaths=sum(meta_all$status[prior&!k]==1,na.rm=TRUE))}))
fwrite(flow,file.path(out,'numeric/eligibility_flow.csv'))
reasons<-data.frame(rule=c('Oh_missing_proteins_gt1000','Go_missing_fraction_at_least_cutoff','nonpositive_or_missing_age','invalid_or_missing_sex','nonpositive_or_missing_time','invalid_or_missing_status'),N=c(sum(!oh_pass),sum(!go_pass),sum(!is.finite(meta_all$age)|meta_all$age<=0),sum(is.na(meta_all$sex)|!meta_all$sex%in%0:1),sum(!is.finite(meta_all$time)|meta_all$time<=0),sum(is.na(meta_all$status)|!meta_all$status%in%0:1)))
fwrite(reasons,file.path(out,'numeric/eligibility_rule_counts_nonexclusive.csv'))
elig<-which(stage_masks$common_with_legal_endpoint);elig<-elig[order(as.numeric(rawids[elig]),rawids[elig])]
m<-meta_all[elig,];stopifnot(nrow(m)>1000,!anyNA(m),all(m$status%in%0:1),all(m$sex%in%0:1))
# Pool sparse age bins by the protocol's deterministic group-wide widths.
keys<-character(nrow(m));widths<-numeric(nrow(m));pool_records<-list()
for(s in 0:1)for(v in 0:1){
 ii<-which(m$sex==s&m$status==v);if(!length(ii))next
 width<-10;bins<-floor(m$age[ii]/width)*width
 if(length(ii)<4){width<-Inf;bins<-rep(0,length(ii))}else{
  while(any(table(bins)<4)){width<-width*2;bins<-floor(m$age[ii]/width)*width;stopifnot(width<100000)}
 }
 keys[ii]<-paste0('sex',s,'_status',v,'_width',if(is.finite(width))sprintf('%05d',width)else'all','_age',sprintf('%05d',bins));widths[ii]<-width
 pool_records[[length(pool_records)+1]]<-data.frame(sex=s,status=v,N=length(ii),initial_nonempty_age_bins=length(table(floor(m$age[ii]/10)*10)),effective_width=if(is.finite(width))as.character(width)else'all',effective_nonempty_age_bins=length(table(bins)),sparse_handling_applied=width!=10)
}
stopifnot(all(nzchar(keys)));strata<-sort(unique(keys));tab<-data.table(stratum=strata,N=as.integer(table(factor(keys,levels=strata))))
tab[,minimum_train:=1L];tab[,maximum_train:=ifelse(N==1L,1L,N-1L)];tab[,fractional_remainder:=.70*N-floor(.70*N)];tab[,train_N:=pmax(minimum_train,pmin(maximum_train,as.integer(floor(.70*N))))]
target<-as.integer(floor(.70*nrow(m)));diff<-target-sum(tab$train_N)
if(diff>0){oo<-order(-tab$fractional_remainder,tab$stratum);while(diff>0)for(j in oo)if(diff>0&&tab$train_N[j]<tab$maximum_train[j]){tab$train_N[j]<-tab$train_N[j]+1L;diff<-diff-1L}}
if(diff<0){oo<-order(tab$fractional_remainder,tab$stratum);while(diff<0)for(j in oo)if(diff<0&&tab$train_N[j]>tab$minimum_train[j]){tab$train_N[j]<-tab$train_N[j]-1L;diff<-diff+1L}}
stopifnot(sum(tab$train_N)==target,all(tab$train_N>=tab$minimum_train),all(tab$train_N<=tab$maximum_train))
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20261006)
tridx<-integer();for(j in seq_len(nrow(tab))){ii<-which(keys==tab$stratum[j]);tridx<-c(tridx,ii[sample.int(length(ii),size=tab$train_N[j],replace=FALSE)])}
tridx<-sort(tridx);teidx<-setdiff(seq_len(nrow(m)),tridx);stopifnot(length(tridx)==target,!length(intersect(tridx,teidx)),length(tridx)+length(teidx)==nrow(m))
trids<-m$eid[tridx];teids<-m$eid[teidx];ids<-c(trids,teids);allidx<-elig[c(tridx,teidx)]
meta<-m[c(tridx,teidx),];rownames(meta)<-NULL;meta$split<-rep(c('train','test'),c(length(trids),length(teids)));meta$stratum<-keys[c(tridx,teidx)]
X<-Xall[allidx,,drop=FALSE];stopifnot(identical(rownames(X),ids),identical(meta$eid,ids));rm(Xall);gc(FALSE)
train<-seq_along(trids);test<-length(trids)+seq_along(teids)
fwrite(meta,file.path(out,'PRIVATE/meta.csv'));saveRDS(meta,file.path(out,'PRIVATE/meta.rds'))
saveRDS(list(train=trids,test=teids,seed=20261006,target_train=target,strata=tab,eligibility=flow,protocol_sha256=digest(file.path(out,'PROTOCOL.txt'),algo='sha256',file=TRUE)),file.path(out,'PRIVATE/common_split.rds'))
tab[,test_N:=N-train_N];fwrite(tab,file.path(out,'numeric/split_strata.csv'));fwrite(rbindlist(pool_records),file.path(out,'numeric/sparse_strata_handling.csv'))
balance<-rbindlist(lapply(list(train=train,test=test),function(ii)data.frame(N=length(ii),deaths=sum(meta$status[ii]),male=sum(meta$sex[ii]==1),age_mean=mean(meta$age[ii]),age_sd=sd(meta$age[ii]),age_min=min(meta$age[ii]),age_q25=as.numeric(quantile(meta$age[ii],.25)),age_median=median(meta$age[ii]),age_q75=as.numeric(quantile(meta$age[ii],.75)),age_max=max(meta$age[ii]))),idcol='split');fwrite(balance,file.path(out,'numeric/split_balance.csv'))
oldcoh<-readRDS(file.path(old,'PRIVATE/common_cohorts.rds'))[[1]]
membership_origin<-rbindlist(lapply(c('old_common_train','old_common_test','outside_old_common'),function(role){k<-if(role=='old_common_train')meta$eid%in%as.character(oldcoh$train)else if(role=='old_common_test')meta$eid%in%as.character(oldcoh$test)else !meta$eid%in%as.character(c(oldcoh$train,oldcoh$test));data.frame(old_role=role,N=sum(k),new_train=sum(k&meta$split=='train'),new_test=sum(k&meta$split=='test'))}))
fwrite(membership_origin,file.path(out,'numeric/prior_cohort_membership_counts.csv'))
missing_tr<-colMeans(is.na(X[train,,drop=FALSE]));ohcols<-cn[missing_tr<=.1];gocols<-cn[missing_tr<.1]
fs<-fromJSON(file.path(old,'results/feature_sets.json'));cla<-readRDS(file.path(root,'Result3/Organ/UKB.Prot.class.rds'))
mapping<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE];stopifnot(setequal(mapping$clock_organ,orgs))
atlas<-setNames(lapply(mapping$disco_organ,function(o)unique(cla$Gene[cla$specific.tissue==o])),mapping$clock_organ)
fullfs<-fs;fs$Oh<-lapply(fs$Oh,intersect,ohcols);fs$Goeminne<-lapply(fs$Goeminne,intersect,gocols)
stopifnot(all(unlist(atlas)%in%cn),all(vapply(fs$Oh,length,0L)>0),all(vapply(fs$Goeminne,length,0L)>0))
write_json(list(Oh=fs$Oh,Goeminne=fs$Goeminne,DISCO=atlas),file.path(out,'numeric/feature_sets.json'),auto_unbox=FALSE,pretty=TRUE)
coverage<-rbindlist(lapply(orgs,function(o)data.frame(organ=o,Oh_candidates_original=length(fullfs$Oh[[o]]),Oh_candidates_train_QC=length(fs$Oh[[o]]),Goeminne_candidates_original=length(fullfs$Goeminne[[o]]),Goeminne_candidates_train_QC=length(fs$Goeminne[[o]]),DISCO_fixed_atlas=length(atlas[[o]]))))
fwrite(coverage,file.path(out,'numeric/feature_coverage.csv'));fwrite(data.frame(protein=cn,training_missing_fraction=missing_tr,Oh_QC=cn%in%ohcols,Goeminne_QC=cn%in%gocols),file.path(out,'numeric/training_QC.csv'))
writeBin(as.double(X),file.path(out,'PRIVATE/raw_locked.bin'),size=8,endian='little');write_json(list(shape=dim(X),columns=cn,order='F',dtype='<f8',train_N=length(train),test_N=length(test)),file.path(out,'PRIVATE/raw_locked.json'),auto_unbox=TRUE)
saveRDS(list(X=X,train=train,test=test,Oh_columns=ohcols,Goeminne_columns=gocols,meta=meta,atlas=atlas),file.path(out,'PRIVATE/prepared.rds'),compress=FALSE)
ref<-train[meta$age[train]<=40];stopifnot(length(ref)>=3,!length(intersect(ids[ref],teids)))
summary<-list(raw_N=nrow(meta_all),eligible_N=length(ids),train_N=length(train),train_deaths=sum(meta$status[train]),test_N=length(test),test_deaths=sum(meta$status[test]),target_ratio=.7,actual_train_fraction=length(train)/length(ids),split_seed=20261006,reference_N=length(ref),train_test_overlap=0,reference_test_overlap=0,raw_P=ncol(X),Oh_QC_P=length(ohcols),Goeminne_QC_P=length(gocols),invalid_cells_raw=nrow(bad),all_parameter_fits_train_only=TRUE,original_inputs_modified=FALSE,original_date_source_audit_passed=FALSE,endpoint='frozen_previous_local_allcause_pending_date_chain_audit',test_previously_seen_in_prior_analyses=TRUE,training_Goeminne_knn_official=TRUE,test_Goeminne_knn_modification='Frozen raw training donors, inductive 10-nearest-neighbour sklearn transform')
write_json(summary,file.path(out,'qa/preparation_ready_before_Go_knn.json'),auto_unbox=TRUE,pretty=TRUE)
print(flow);print(balance);print(coverage);cat('Shared 70:30 members frozen; train/test/ref:',length(train),length(test),length(ref),'\n');flush.console()
source(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison'),'20261002_official_reproduction/vendor/Goeminne/scripts/R/Olink_Processing.R'))
M<-X[train,gocols,drop=FALSE];cat('TRAIN ONLY official Goeminne impute.knn',nrow(M),'x',ncol(M),'\n');flush.console()
z<-as.matrix(impute_olink(M));stopifnot(all(is.finite(z)),identical(dimnames(z),dimnames(M)),max(abs(z[!is.na(M)]-M[!is.na(M)]))==0)
selected<-sort(unique(unlist(fs$Goeminne)));writeBin(as.double(z[,selected,drop=FALSE]),file.path(out,'PRIVATE/Goeminne_train_selected.bin'),size=8,endian='little');write_json(list(shape=c(nrow(z),length(selected)),columns=selected,order='F',dtype='<f8'),file.path(out,'PRIVATE/Goeminne_train_selected.json'),auto_unbox=TRUE)
sources<-c(file.path(root,'UKB_Data/20240321_prot_instance0.RData'),file.path(root,'Result2/Proteins/coding143.tsv'),file.path(root,'Result3/Organ/UKB.Prot.class.rds'),file.path(old,'PRIVATE/meta.rds'),file.path(old,'PRIVATE/common_cohorts.rds'),file.path(old,'results/feature_sets.json'),file.path(out,'PROTOCOL.txt'),file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison'),'20261002_official_reproduction/vendor/Goeminne/scripts/R/Olink_Processing.R'))
fwrite(data.frame(path=normalizePath(sources),sha256=vapply(sources,digest,'',file=TRUE,algo='sha256')),file.path(out,'qa/input_sha256.csv'))
write_json(summary,file.path(out,'qa/preparation.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'logs/R_sessionInfo.txt'));pf<-list.files(file.path(out,'PRIVATE'),full.names=TRUE,recursive=FALSE);Sys.chmod(pf[!file.info(pf)$isdir],'0600');Sys.chmod(pf[file.info(pf)$isdir],'0700');cat('Preparation and official TRAIN ONLY Goeminne KNN complete\n');sink()
