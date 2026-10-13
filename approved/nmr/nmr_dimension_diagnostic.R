# Read-only diagnostic of the exact original NMR construction block.
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
base <- file.path(Sys.getenv('DISCO_DATA_ROOT'),'Result6')
out <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
inputs <- file.path(base,c('Cal_Rugao_omics.R','Rugao/NMR.csv','Rugao.Phe.DISCO2.rds','Rugao.PCs.met.rds','Rugao.Met.DISCO.rds','Rugao.Met.DISCO2.rds'))
before <- tools::md5sum(inputs)
met <- read.csv(inputs[2])[,1:317]
pca <- FactoMineR::PCA(met[,5:ncol(met)],graph=FALSE,ncp=50)
met1 <- pca$ind$coord
dis <- readRDS(inputs[3]);dis <- dis[match(met$X2019ID,dis$ID),]
stopifnot(!anyNA(dis$ID))
cc <- apply(met1,2,function(x)cor(dis$mc,x,use='pairwise.complete.obs'))
pp <- apply(met1,2,function(x)cor.test(dis$mc,x,use='pairwise.complete.obs')$p.value)
bb <- names(cc)[pp<.05];cc <- cc[pp<.05]
ref <- met1[which(dis$mc<quantile(dis$mc)[2]),]
weight <- abs(cc %*% t(cc));diag(weight)<-0;weight<-weight/sum(as.numeric(weight))
saved_pcs <- readRDS(inputs[4])
stopifnot(identical(colnames(saved_pcs)[-1],bb))
# Load only the reviewed dmc function definition, never execute the source's writes.
for(expr in parse(inputs[1]))if(is.call(expr) && identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name('dmc')))eval(expr)
warnings <- character()
legacy <- withCallingHandlers(dmc(met1,ref,weight),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
dimension_matched <- dmc(met1[,bb],ref[,bb],weight)
comparison <- data.frame(ID=met$X2019ID,legacy_recycled=legacy,matched_23PC=dimension_matched)
for(f in inputs[5:6]){
 x<-readRDS(f);x<-x[match(comparison$ID,x$ID),]
 comparison[[basename(f)]]<-x$mc
}
write.csv(head(comparison,10),file.path(out,'NMR_first10_score_diagnostic.csv'),row.names=FALSE)
write.csv(comparison,file.path(out,'NMR_all2120_score_diagnostic.csv'),row.names=FALSE)
# Independent covariance rank-one update: does not call the legacy dmc/cor loop.
m<-nrow(ref);r23<-ref[,bb];mu<-colMeans(r23)
centered<-sweep(r23,2,mu);scatter<-crossprod(centered)
rcor<-cov2cor(scatter)
independent<-apply(met1[,bb],1,function(x){
 delta<-x-mu;aug<-scatter+m/(m+1)*tcrossprod(delta)
 ccnew<-cov2cor(aug)
 log(sum((rcor-ccnew)^2*weight)*m*m)
})
stopifnot(max(abs(dimension_matched-independent))<1e-8)
metrics<-data.frame(comparison=c('legacy vs saved1','legacy vs saved2','selected23 vs saved1','selected23 vs independent'),
                   max_abs_difference=c(max(abs(legacy-comparison[[4]])),max(abs(legacy-comparison[[5]])),
                                        max(abs(dimension_matched-comparison[[4]])),max(abs(dimension_matched-independent))),
                   Pearson_r=c(cor(legacy,comparison[[4]]),cor(legacy,comparison[[5]]),cor(dimension_matched,comparison[[4]]),cor(dimension_matched,independent)),
                   Spearman_r=c(cor(legacy,comparison[[4]],method='spearman'),cor(legacy,comparison[[5]],method='spearman'),
                                cor(dimension_matched,comparison[[4]],method='spearman'),cor(dimension_matched,independent,method='spearman')))
write.csv(metrics,file.path(out,'NMR_full_verification.csv'),row.names=FALSE)
saveRDS(list(ID=met$X2019ID,score=dimension_matched,selected_PCs=bb,reference_ID=met$X2019ID[dis$mc<quantile(dis$mc)[2]],
             weight=weight,source='Cal_Rugao_omics.R one-line dimension repair: dmc(met1[,bb],ref[,bb],weight)',
             independent_max_error=max(abs(dimension_matched-independent)),
             status='Dimension-matched candidate, not silently promoted as final manuscript analysis'),
        file.path(out,'NMR_selected23_dimension_matched_candidate.rds'))
table14<-do.call(rbind,lapply(seq_len(50),function(i){
 a<-cor.test(met1[,i],dimension_matched,method='spearman');b<-cor.test(met1[,i],dis$mc,method='spearman')
 data.frame(PC=paste0('PC',i),cumulative_variance=pca$eig[i,3],
            all_mean_sd=sprintf('%.2f (%.2f)',mean(met1[,i]),sd(met1[,i])),
            reference_mean_sd=sprintf('%.2f (%.2f)',mean(ref[,i]),sd(ref[,i])),
            ks_p=ks.test(met1[,i],ref[,i])$p.value,met_disco_rho=unname(a$estimate),met_disco_p=a$p.value,
            clinical_rho=unname(b$estimate),clinical_p=b$p.value)
}))
write.csv(table14,file.path(out,'S14_NMR_statistics_corrected.csv'),row.names=FALSE)
write.csv(data.frame(n_all=2120,n_reference=nrow(ref),n_PCs_table=ncol(met1),n_PCs_selected=length(bb),
                    covariance_vector_length=ncol(met1)^2,weight_vector_length=length(weight),warnings=length(warnings),
                    saved_selected_PCs_identical=TRUE,source_inputs_unchanged=identical(before,tools::md5sum(inputs))),
          file.path(out,'NMR_dimension_diagnostic.csv'),row.names=FALSE)
writeLines(unique(warnings),file.path(out,'NMR_dimension_warnings.txt'))
print(metrics);print(unique(warnings))
