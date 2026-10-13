# Minimal corrected propagation of the original Result6 score and plotting chain.
# Original files are never overwritten. Independently executed pass_1/pass_2.
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
library(Hmisc);library(survival);library(survminer);library(ggplot2)
B <- file.path(Sys.getenv('DISCO_DATA_ROOT'),'Result6')
stage <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
BASE <- file.path(stage,'downstream')
pass <- commandArgs(TRUE)[1];stopifnot(pass%in%c('pass_1','pass_2'))
O <- file.path(BASE,pass);dir.create(O,recursive=TRUE,showWarnings=FALSE)
logcon <- file(file.path(O,'run.log'),'wt');sink(logcon,split=TRUE)
inputs<-c(file.path(B,c('Cal_Rugao_omics.R','Plot_surv.R','Plot_subsystems2.R','microbe/Species.txt','Rugao.Met.DISCO.rds','Rugao.Met.DISCO2.rds','Rugao.Meta.DISCO2.rds','Rugao.selected.microbe.rds','Rugao.subsystems.rds','Rugao.subsystems2.rds','Rugao.2019.mortality.rds','Analysis.2019.x.xlsx','Rugao/2019Rugao.data.y1.csv')),file.path(stage,'nmr_corrected/pass_1/Rugao.Met.DISCO2_CORRECTED.rds'))
hashes<-vapply(inputs,function(f)digest::digest(file=f,algo='sha256'),character(1))
write.csv(data.frame(path=inputs,sha256=hashes),file.path(O,'input_sha256.csv'),row.names=FALSE)
cat('Started',format(Sys.time()),'\n')
readsrc <- function(f)readRDS(file.path(B,f))
micro <- read.table(file.path(B,'microbe/Species.txt'),header=TRUE)
colnames(micro) <- as.numeric(gsub('RG','',colnames(micro)))
rownames(micro) <- micro[,1];micro <- micro[,-1]
micro <- micro[rowSums(micro==0)<.95*ncol(micro),]
samples <- as.numeric(colnames(micro));meta <- apply(micro,1,Hmisc::impute)
stopifnot(!anyDuplicated(samples),all(is.finite(meta)))
oldmet0 <- readsrc('Rugao.Met.DISCO2.rds');oldmet <- oldmet0[match(samples,oldmet0$ID),]
newmet0 <- readRDS(file.path(stage,'nmr_corrected/pass_1/Rugao.Met.DISCO2_CORRECTED.rds'));newmet <- newmet0[match(samples,newmet0$ID),]
sel <- function(m,field='mc') {
 cc <- apply(meta,2,function(x)cor(m[[field]],x,method='spearman',use='pairwise.complete.obs'))
 pp <- apply(meta,2,function(x)suppressWarnings(cor.test(m[[field]],x,method='spearman')$p.value))
 bb <- names(cc)[pp<.05];ri <- which(m[[field]]<quantile(m[[field]],na.rm=TRUE)[2]);ww <- abs(outer(cc[bb],cc[bb]));diag(ww) <- 0;ww <- ww/sum(ww)
 list(cc=cc,p=pp,bb=bb,refidx=ri,weight=ww)
}
old <- sel(oldmet);new <- sel(newmet);hddef <- sel(oldmet,'hd')
saved <- readsrc('Rugao.selected.microbe.rds')
stopifnot(identical(old$bb,names(saved)[-1]),identical(samples,saved$ID),max(abs(meta[,old$bb]-as.matrix(saved[,-1])))==0)
# pass_2 independently recomputes all ranks in every augmented sample and then
# computes their Pearson matrix; this is mathematically true Spearman correlation.
scor <- if(pass=='pass_1')function(x)cor(x,method='spearman') else function(x)cor(apply(x,2,rank,ties.method='average'))
disco <- function(d,ref,w,label) {
 stopifnot(ncol(d)==ncol(ref),identical(dim(w),rep(ncol(d),2)),all(is.finite(d)),all(is.finite(ref)))
 r0<-scor(ref);n<-nrow(ref)
 v<-vapply(seq_len(nrow(d)),function(i){
  if(i%%250==0)cat(label,i,'/',nrow(d),format(Sys.time()),'\n')
  r1<-scor(rbind(ref,d[i,]));log(sum((r1-r0)^2*w)*n^2)
 },numeric(1));stopifnot(all(is.finite(v)));v
}
oldscore<-disco(meta[,old$bb],meta[old$refidx,old$bb],old$weight,'old microbe')
newscore<-disco(meta[,new$bb],meta[new$refidx,new$bb],new$weight,'new microbe')
oldpd<-readsrc('Rugao.Meta.DISCO2.rds');oldpd<-oldpd[match(samples,oldpd$ID),]
old_score_error<-max(abs(oldscore-oldpd$mc));cat('Old microbe saved mc match',old_score_error,'\n')
# The HD branch has a different selection/reference and is unaffected by Met mc.
# Verify direct Mahalanobis definitions but preserve saved hd: installed BioAge is
# a local compatibility shim and is NOT claimed to reproduce the historical API.
xhd<-meta[,hddef$bb];refhd<-xhd[hddef$refidx,];md2<-mahalanobis(xhd,colMeans(refhd),cov(refhd))
hdcheck<-data.frame(definition=c('squared Mahalanobis','sqrt Mahalanobis','log squared','log sqrt'),max_abs_error=vapply(list(md2,sqrt(md2),log(md2),log(sqrt(md2))),function(x)max(abs(x-oldpd$hd)),numeric(1)))
write.csv(hdcheck,file.path(O,'unaffected_DM_definition_probe.csv'),row.names=FALSE)
newpd<-oldpd;newpd$mc<-newscore
saveRDS(newpd,file.path(O,'Rugao.Meta.DISCO2_CORRECTED.rds'))
saveRDS(data.frame(ID=samples,meta[,new$bb]),file.path(O,'Rugao.selected.microbe_CORRECTED.rds'))
saveRDS(list(old=old,new=new,HD=hddef,ID=samples,oldscore=oldscore,newscore=newscore,reference_membership_changed=sum(xor(samples%in%samples[old$refidx],samples%in%samples[new$refidx]))),file.path(O,'score_definition_and_verification.rds'))
write.csv(data.frame(ID=samples,old_DISCO_reference=seq_along(samples)%in%old$refidx,corrected_DISCO_reference=seq_along(samples)%in%new$refidx,DM_reference=seq_along(samples)%in%hddef$refidx),file.path(O,'reference_membership_INTERNAL.csv'),row.names=FALSE)
write.csv(data.frame(ID=samples,old_mc=oldpd$mc,replayed_old_mc=oldscore,corrected_mc=newscore),file.path(O,'microbe_corrected_scores_INTERNAL.csv'),row.names=FALSE)
write.csv(data.frame(species=colnames(meta),old_rho=old$cc,old_P=old$p,new_rho=new$cc,new_P=new$p,old_selected=colnames(meta)%in%old$bb,new_selected=colnames(meta)%in%new$bb),file.path(O,'microbe_selection_comparison.csv'),row.names=FALSE)
# S15: report both DISCO and DM reference summaries explicitly. The old source
# reuses the HD reference under one generic reference header (Cal_Rugao:375–426).
ct<-function(x,y){z<-suppressWarnings(cor.test(x,y,method='spearman'));c(rho=unname(z$estimate),P=z$p.value)}
tab<-do.call(rbind,lapply(new$bb,function(b){x<-meta[,b];data.frame(Biomarker=b,N_all=length(x),mean_all=mean(x),sd_all=sd(x),N_reference_DISCO=length(new$refidx),mean_reference_DISCO=mean(x[new$refidx]),sd_reference_DISCO=sd(x[new$refidx]),P_KS_DISCO=suppressWarnings(ks.test(x,x[new$refidx])$p.value),N_reference_DM=length(hddef$refidx),mean_reference_DM=mean(x[hddef$refidx]),sd_reference_DM=sd(x[hddef$refidx]),P_KS_DM=suppressWarnings(ks.test(x,x[hddef$refidx])$p.value),rho_microbe_DISCO=ct(x,newscore)[1],P_microbe_DISCO=ct(x,newscore)[2],rho_microbe_DM=ct(x,oldpd$hd)[1],P_microbe_DM=ct(x,oldpd$hd)[2],rho_met_DISCO=ct(x,newmet$mc)[1],P_met_DISCO=ct(x,newmet$mc)[2],rho_met_DM=ct(x,oldmet$hd)[1],P_met_DM=ct(x,oldmet$hd)[2])}))
write.csv(tab,file.path(O,'S15_corrected_30_species.csv'),row.names=FALSE)
death<-readsrc('Rugao.2019.mortality.rds');stopifnot(!anyDuplicated(death$x2019ID))
survdata<-function(pd){m<-death[match(pd$ID,death$x2019ID),];m$time[is.na(m$status)]<-24;m$status[is.na(m$status)]<-0;pd$time<-m$time;pd$status<-m$status;pd}
newms<-survdata(newmet0);news<-survdata(newpd);olds<-survdata(oldpd)
saveRDS(newms,file.path(O,'Rugao.Met.surv.rds'));saveRDS(news,file.path(O,'Rugao.Meta.surv.rds'))
modelrows<-function(d,label){do.call(rbind,lapply(c('raw','per_SD'),function(scale_type){if(scale_type=='per_SD')d$mc<-as.numeric(scale(d$mc));do.call(rbind,lapply(c('score_only','age_sex'),function(adj){fit<-coxph(as.formula(if(adj=='score_only')'Surv(time,status)~mc' else 'Surv(time,status)~age+sex+mc'),data=d);s<-summary(fit);k<-match('mc',rownames(s$coefficients));data.frame(dataset=label,scale=scale_type,model=adj,N=fit$n,deaths=fit$nevent,HR=s$conf.int[k,'exp(coef)'],lo=s$conf.int[k,'lower .95'],hi=s$conf.int[k,'upper .95'],P=s$coefficients[k,'Pr(>|z|)'],C=unname(s$concordance[1]))}))}))}
clinical<-readsrc('Rugao.subsystems.rds');clinicalhd<-readsrc('Rugao.subsystems2.rds')
rawphe<-suppressWarnings(readxl::read_xlsx(file.path(B,'Analysis.2019.x.xlsx')))
aggregate<-function(met,mi){d<-clinical;d$Met<-met$mc[match(d$ID,met$ID)];d$Meta<-mi$mc[match(d$ID,mi$ID)];stopifnot(identical(names(d)[3:12],c(names(clinical)[3:10],'Met','Meta')));d$mc<-rowMeans(d[,3:12]);d$age<-as.numeric(rawphe$Age[match(d$ID,rawphe$X2019ID)]);d$sex<-rawphe$Sex[match(d$ID,rawphe$X2019ID)];survdata(d)}
oldagg<-aggregate(oldmet0,oldpd);newagg<-aggregate(newmet0,newpd)
saveRDS(newagg,file.path(O,'Fig2f_aggregate_corrected.rds'))
models<-rbind(modelrows(olds,'microbe_old_saved'),modelrows(news,'microbe_corrected'),modelrows(oldagg,'aggregate_old_source_mean'),modelrows(newagg,'aggregate_corrected_source_mean'))
write.csv(models,file.path(O,'mortality_model_comparison.csv'),row.names=FALSE)
counts<-data.frame(metric=c('microbe_people','microbe_species_retained','old_selected','new_selected','old_reference','new_reference','HD_reference','reference_symdiff','old_microbe_score_max_abs_error','microbe_score_Spearman','microbe_score_max_change','aggregate_complete_scores','aggregate_score_Spearman'),value=c(nrow(meta),ncol(meta),length(old$bb),length(new$bb),length(old$refidx),length(new$refidx),length(hddef$refidx),sum(xor(seq_along(samples)%in%old$refidx,seq_along(samples)%in%new$refidx)),old_score_error,cor(oldpd$mc,newscore,method='spearman'),max(abs(newscore-oldpd$mc)),sum(is.finite(newagg$mc)),cor(oldagg$mc,newagg$mc,method='spearman',use='complete.obs')))
write.csv(counts,file.path(O,'key_results.csv'),row.names=FALSE);print(counts);print(models)
# Replay exact original drawing paragraphs with only controlled path/read dispatch
# and explicit printing required in noninteractive execution.
replay<-function(file,lines,env,tag){txt<-readLines(file.path(B,file),warn=FALSE);txt<-txt[intersect(lines,seq_along(txt))];txt<-txt[!grepl('^ggsurvplot\\(fit1,xlim=',txt)];writeLines(txt,file.path(O,paste0('executed_original_',tag,'.R')));env$pdf<-function(file,...)grDevices::pdf(file.path(O,file),...);env$readRDS<-function(file){if(file%in%c('Rugao.Met.surv.rds','Rugao.Meta.surv.rds'))return(base::readRDS(file.path(O,file)));if(file=='Rugao.Meta.DISCO2.rds')return(newpd);readsrc(file)};env$saveRDS<-function(object,file,...)base::saveRDS(object,file.path(O,file),...);env$read.csv<-function(file,...)utils::read.csv(file.path(B,file),...);for(e in parse(text=txt)){for(n in intersect(c('disco.sub1','hd.sub1','pd.vd1','dd4'),ls(env)))assign(n,get(n,env),globalenv());v<-withVisible(eval(e,env));if(v$visible&&(inherits(v$value,'ggsurvplot')||inherits(v$value,'ggplot')))print(v$value)};env}
pdf(file=NULL)
env<-new.env(parent=globalenv());env<-replay('Plot_surv.R',238:352,env,'Plot_surv_238_352')
saveRDS(env$fit1,file.path(O,'Fig2e_age_timescale_survfit.rds'));saveRDS(env$p1,file.path(O,'Fig2e_original_ggplot_corrected.rds'));write.csv(as.data.frame(table(env$pd.vd1$group1)),file.path(O,'Fig2e_extreme_groups.csv'),row.names=FALSE)
env<-new.env(parent=globalenv());env$rugao.phe<-rawphe;env$mortality<-death
env<-replay('Plot_subsystems2.R',181:370,env,'Plot_subsystems2_181_370')
saveRDS(env$fit1,file.path(O,'Fig2f_age_timescale_survfit.rds'));saveRDS(env$p1,file.path(O,'Fig2f_original_ggplot_corrected.rds'))
write.csv(as.data.frame(table(env$pd.vd1$risk1)),file.path(O,'Fig2f_extreme_groups.csv'),row.names=FALSE)
while(dev.cur()>1)dev.off()
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'))
stopifnot(identical(hashes,vapply(inputs,function(f)digest::digest(file=f,algo='sha256'),character(1))))
cat('Complete',format(Sys.time()),'\n');sink();close(logcon)
