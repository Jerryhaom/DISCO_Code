#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE)))
base<-file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final');.libPaths(c(file.path(base,'09_figures_tables/Rlibs'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival);library(Rcpp);library(digest)})
setDTthreads(1);Sys.umask('0077');options(warn=1);sink(file.path(out,'logs/04_compare.log'),split=TRUE)
input<-Sys.getenv('DISCO_FROZEN_SCORE_SNAPSHOT',out)
ready<-fromJSON(file.path(input,'qa/clocks_complete.json'));stopifnot(isTRUE(ready$all_selected_converged),ready$test_participants_used_in_fitting==0)
ds<-readRDS(file.path(input,'PRIVATE/DISCO_original.rds'));ids<-as.character(ds$ids);trix<-ds$train;teix<-ds$test
pred<-fread(file.path(input,'PRIVATE/clock_predictions.csv'),colClasses=c(eid='character'),data.table=FALSE)
stopifnot(!anyDuplicated(pred[,c('eid','organ','method')]),nrow(pred)==length(ids)*18)
meta<-ds$meta;meta$sex<-factor(meta$sex,levels=0:1,labels=c('Female','Male'));orgs<-colnames(ds$scores)
mapping<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
selected<-fromJSON(file.path(input,'numeric/new_train_selected_proteins.json'));matched<-list();directchecks<-list();coverage<-list()
sourceCpp(file.path(Sys.getenv('DISCO_APPROVED_CODE_ROOT'),'resources/exact_score.cpp'),cacheDir=file.path(out,'PRIVATE/cpp_cache'))
probes<-unique(c(1,ds$reference[c(1,2,length(ds$reference))],teix[c(1,2,length(teix))]))
for(o in orgs)for(m in c('Oh','Goeminne')){
 key<-paste(m,o,sep='_');pp<-selected[[key]];coverage[[key]]<-data.frame(organ=o,method=m,p=length(pp),available=length(pp)>=2,reference_N=length(ds$reference),reference_test_overlap=0,reason=if(length(pp)<2)'Less than two train-selected proteins; no off-diagonal correlation pairs' else '')
 if(length(pp)<2)next
 stopifnot(all(pp%in%colnames(ds$X)));cc<-cor(ds$X[ds$reference,pp,drop=FALSE]);w<-abs(outer(ds$agecor[pp],ds$agecor[pp]));diag(w)<-0;stopifnot(sum(w)>0);w<-w/sum(w)
 s<-organ_score(ds$X[,pp,drop=FALSE],ds$mu[pp],ds$sd[pp],cc,w,length(ds$reference),FALSE)
 direct<-vapply(probes,function(i)log(length(ds$reference)^2*sum(w*(cor(rbind(ds$X[ds$reference,pp,drop=FALSE],ds$X[i,pp,drop=FALSE]))-cc)^2)),0.)
 err<-max(abs(s[probes]-direct));stopifnot(all(is.finite(s)),err<1e-9);directchecks[[key]]<-data.frame(method=m,organ=o,p=length(pp),probes=length(probes),max_abs_error=err);matched[[key]]<-s
}
fwrite(rbindlist(coverage),file.path(out,'numeric/matched_protein_set_counts.csv'));fwrite(rbindlist(directchecks),file.path(out,'qa/matched_DISCO_direct_checks.csv'))
allCI<-allPairs<-allHR<-allPH<-allCor<-allCal<-drows<-mdrows<-allCoh<-calchecks<-list();fitted_private<-list();warnings_out<-list()
fitcox<-function(formula,tr,label){withCallingHandlers(coxph(formula,data=tr,ties='efron',x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-10,iter.max=50)),warning=function(w){warnings_out[[length(warnings_out)+1]]<<-data.frame(model=label,warning=conditionMessage(w));invokeRestart('muffleWarning')})}
hrrow<-function(f,term,o,m,n){b<-unname(coef(f)[term]);se<-sqrt(vcov(f)[term,term]);data.frame(organ=o,method=m,model=n,term=term,N=f$n,deaths=f$nevent,beta=b,SE=se,HR=exp(b),lower95=exp(b-qnorm(.975)*se),upper95=exp(b+qnorm(.975)*se),P=2*pnorm(-abs(b/se)),formula=paste(deparse(formula(f)),collapse=''),fit_sample='common_train_only',scale='per one training linear-age-residual SD; age+sex adjustment',endpoint_status='PROVISIONAL frozen local 2022-12 endpoint pending raw censor/death-date audit')}
for(o in orgs){
 dd<-meta;dd$D<-ds$scores[,o]
 for(m in c('Oh','Goeminne')){pp<-pred[pred$organ==o&pred$method==m,];stopifnot(setequal(pp$eid,ids));k<-if(m=='Oh')'O' else 'G';dd[[k]]<-pp$score[match(ids,pp$eid)]}
 for(m in c('Oh','Goeminne')){key<-paste(m,o,sep='_');if(key %in%names(matched))dd[[if(m=='Oh')'DO' else 'DG']]<-matched[[key]]}
 stopifnot(nrow(dd)==length(ids),all(complete.cases(dd)),all(dd$time>0));tr<-dd[trix,];te<-dd[teix,];stopifnot(nrow(tr)==length(trix),nrow(te)==length(teix),!length(intersect(tr$eid,te$eid)))
 allCoh[[o]]<-data.frame(organ=o,train_N=nrow(tr),train_deaths=sum(tr$status),test_N=nrow(te),test_deaths=sum(te$status),reference_N=length(ds$reference),reference_test_overlap=0,endpoint='PROVISIONAL frozen local all-cause mortality; raw dates/censoring pending')
 terms<-c(DISCO='D',Oh='O',Goeminne='G',DISCO_Oh='DO',DISCO_Goeminne='DG');terms<-terms[terms%in%names(tr)]
 for(k in terms){
   f<-lm(reformulate('age',k),tr);den<-sd(resid(f));stopifnot(den>0);tr[[paste0(k,'a')]]<-resid(f)/den;te[[paste0(k,'a')]]<-(te[[k]]-predict(f,newdata=te))/den
   stopifnot(abs(cor(tr[[paste0(k,'a')]],tr$age))<1e-10,abs(sd(tr[[paste0(k,'a')]])-1)<1e-10)
   allCal[[paste(o,k)]]<-data.frame(organ=o,score=k,calibration='linear_age',intercept=unname(coef(f)[1]),age_slope=unname(coef(f)[2]),training_residual_SD=den,training_N=nrow(tr),test_N_used=0)
   if(k %in%c('O','G')){
     m<-if(k=='O')'Oh' else 'Goeminne';pp<-pred[pred$organ==o&pred$method==m,];err<-max(abs(dd[[k]]-predict(f,newdata=dd)-pp$gap_train_calibrated[match(ids,pp$eid)]));zerr<-max(abs(c(tr[[paste0(k,'a')]],te[[paste0(k,'a')]])-pp$gap_z_trainingSD[match(ids,pp$eid)]));stopifnot(err<1e-8,zerr<1e-8)
     calchecks[[paste(o,k)]]<-data.frame(organ=o,method=m,Python_vs_R_gap_max_error=err,Python_vs_R_training_SD_gap_max_error=zerr)
   }
 }
 # Fixed train-nuisance adjusted test correlation, not a test-fitted partial-r model.
 for(m in c('Oh','Goeminne')){
   k<-if(m=='Oh')'O' else 'G';fd<-lm(D~age+sex,tr);fp<-lm(reformulate(c('age','sex'),k),tr)
   rd<-te$D-predict(fd,newdata=te);rp<-te[[k]]-predict(fp,newdata=te);r<-cor(rd,rp);n<-nrow(te);CI<-tanh(atanh(r)+c(-1,1)*qnorm(.975)/sqrt(n-3))
   allCor[[paste(o,m)]]<-data.frame(organ=o,method=m,N=n,r=r,lower95=CI[1],upper95=CI[2],P=2*pt(-abs(r*sqrt((n-2)/(1-r^2))),n-2),adjustment='apply train-fitted age+sex coefficients to test scores',estimator='Pearson r of train-calibrated test residuals; not test-fitted partial r')
   hforms<-list(D_only=Surv(time,status)~age+sex+Da,P_only=reformulate(c('age','sex',paste0(k,'a')),'Surv(time,status)'),joint=reformulate(c('age','sex','Da',paste0(k,'a')),'Surv(time,status)'))
   for(hn in names(hforms)){
     f<-fitcox(hforms[[hn]],tr,paste('HR',o,m,hn));stopifnot(all(is.finite(coef(f))))
     for(term in intersect(c('Da',paste0(k,'a')),names(coef(f))))allHR[[paste(o,m,hn,term)]]<-hrrow(f,term,o,m,hn)
     if(hn=='joint'){ph<-cox.zph(f,transform='km')$table;allPH[[paste(o,m)]]<-data.frame(organ=o,method=m,term=rownames(ph),chisq=ph[,'chisq'],df=ph[,'df'],P=ph[,'p'],fit_sample='common_train_only')}
   }
 }
 forms<-list(Fair_Age_sex=Surv(time,status)~age+sex,Fair_DISCO=Surv(time,status)~age+sex+D,Fair_Oh=Surv(time,status)~age+sex+O,Fair_Goeminne=Surv(time,status)~age+sex+G,Fair_DISCO_Oh=Surv(time,status)~age+sex+D+O,Fair_DISCO_Goeminne=Surv(time,status)~age+sex+D+G,
   Display_DISCO=Surv(time,status)~age+D,Display_Oh=Surv(time,status)~O,Display_Goeminne=Surv(time,status)~G,Accel_DISCO=Surv(time,status)~Da,Accel_Oh=Surv(time,status)~Oa,Accel_Goeminne=Surv(time,status)~Ga)
 for(m in c('Oh','Goeminne')){
   k<-if(m=='Oh')'DO' else 'DG';if(!k%in%names(tr))next
   forms[[paste0('Matched_Display_DISCO_',m)]]<-reformulate(c('age',k),'Surv(time,status)');forms[[paste0('Matched_Accel_DISCO_',m)]]<-reformulate(paste0(k,'a'),'Surv(time,status)')
   forms[[paste0('Matched_Fair_DISCO_',m)]]<-reformulate(c('age','sex',k),'Surv(time,status)')
 }
 fits<-lapply(names(forms),function(n)fitcox(forms[[n]],tr,paste('C',o,n)));names(fits)<-names(forms);stopifnot(all(vapply(fits,function(f)f$n==length(trix)&&all(is.finite(coef(f))),TRUE)))
 cc<-do.call(concordance,c(unname(fits),list(newdata=te,timewt='n',influence=2)));cv<-unname(cc$concordance);V<-cc$var;stopifnot(length(cv)==length(fits),all(is.finite(V)),min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
 lpall<-sapply(fits,predict,newdata=te,type='lp');stopifnot(all(is.finite(lpall)))
 for(j in seq_along(fits)){
   direct<-concordance(Surv(te$time,te$status)~lpall[,j],reverse=TRUE);stopifnot(abs(cv[j]-direct$concordance)<1e-12,abs(V[j,j]-direct$var)<1e-12)
   model<-names(fits)[j];se<-sqrt(V[j,j]);rr<-data.frame(organ=o,model=model,N=nrow(te),deaths=sum(te$status),Cindex=cv[j],SE=se,lower95=max(0,cv[j]-qnorm(.975)*se),upper95=min(1,cv[j]+qnorm(.975)*se),formula=paste(deparse(formula(fits[[j]])),collapse=''),Cox_fit_sample='common_train_only',evaluation='locked_common_test',preprocessing='strict training-only frozen parameters',endpoint_status='PROVISIONAL frozen local 2022-12 endpoint pending raw censor/death-date audit')
   allCI[[paste(o,model)]]<-rr
   if(grepl('^(Display_|Accel_)',model)){
     score<-sub('^(Display_|Accel_)','',model);group<-if(grepl('^Display_',model))'Original' else 'Accel';drows[[paste(o,model)]]<-data.frame(organ=o,display=mapping$display[match(o,mapping$clock_organ)],group=group,score=score,N=nrow(te),deaths=sum(te$status),C=cv[j],SE=se,lower95=rr$lower95,upper95=rr$upper95,formula=rr$formula)
   }
   if(grepl('^(Matched_Display_|Matched_Accel_)',model)){
     score<-sub('^(Matched_Display_|Matched_Accel_)','',model);group<-if(grepl('^Matched_Display_',model))'Original' else 'Accel';mdrows[[paste(o,model)]]<-data.frame(organ=o,display=mapping$display[match(o,mapping$clock_organ)],group=group,score=score,N=nrow(te),deaths=sum(te$status),C=cv[j],SE=se,lower95=rr$lower95,upper95=rr$upper95,formula=rr$formula,available=TRUE)
   }
 }
 contrasts<-list(c('Fair_DISCO','Fair_Oh'),c('Fair_DISCO','Fair_Goeminne'),c('Fair_DISCO_Oh','Fair_Oh'),c('Fair_DISCO_Goeminne','Fair_Goeminne'),c('Fair_DISCO_Oh','Fair_DISCO'),c('Fair_DISCO_Goeminne','Fair_DISCO'),c('Display_DISCO','Display_Oh'),c('Display_DISCO','Display_Goeminne'),c('Accel_DISCO','Accel_Oh'),c('Accel_DISCO','Accel_Goeminne'))
 for(m in c('Oh','Goeminne'))if(paste0('Matched_Fair_DISCO_',m)%in%names(fits)){
   contrasts<-c(contrasts,list(c(paste0('Matched_Fair_DISCO_',m),paste0('Fair_',m)),c(paste0('Matched_Fair_DISCO_',m),'Fair_DISCO'),c(paste0('Matched_Display_DISCO_',m),paste0('Display_',m)),c(paste0('Matched_Accel_DISCO_',m),paste0('Accel_',m))))
 }
 for(pair in contrasts){a<-match(pair[1],names(fits));b<-match(pair[2],names(fits));delta<-cv[a]-cv[b];se<-sqrt(max(0,V[a,a]+V[b,b]-2*V[a,b]));allPairs[[paste(o,paste(pair,collapse='_'))]]<-data.frame(organ=o,contrast=paste(pair,collapse=' minus '),N=nrow(te),deaths=sum(te$status),delta_Cindex=delta,SE=se,lower95=delta-qnorm(.975)*se,upper95=delta+qnorm(.975)*se,P=if(se>0)2*pnorm(-abs(delta/se))else NA_real_,CI='paired IJ conditional on frozen trained models')}
 for(m in c('Oh','Goeminne'))if(!paste0(if(m=='Oh')'DO'else'DG')%in%names(tr))for(g in c('Original','Accel'))mdrows[[paste(o,m,g)]]<-data.frame(organ=o,display=mapping$display[match(o,mapping$clock_organ)],group=g,score=paste0('DISCO_',m),N=nrow(te),deaths=sum(te$status),C=NA_real_,SE=NA_real_,lower95=NA_real_,upper95=NA_real_,formula='Undefined with fewer than two train-selected proteins',available=FALSE)
 fitted_private[[o]]<-list(fits=fits,test=te,train_ids=tr$eid,risks=lpall,covariance=V,train_calibration_records=allCal[grepl(paste0('^',o,' '),names(allCal))]);cat('Verified all held-out C and paired covariance',o,'models',length(fits),'\n');flush.console()
}
ci<-rbindlist(allCI);pairs<-rbindlist(allPairs);pairs[,q_BH_9:=p.adjust(P,'BH'),by=contrast];hr<-rbindlist(allHR);hr[,q_BH_18:=p.adjust(P,'BH'),by=.(model,term)];ph<-rbindlist(allPH);ph[,q_BH_18:=p.adjust(P,'BH'),by=term]
fwrite(ci,file.path(out,'numeric/all_test_Cindex.csv'));fwrite(ci[grepl('^Fair_',model)],file.path(out,'numeric/fair_test_Cindex.csv'));fwrite(pairs,file.path(out,'numeric/paired_Cindex_differences.csv'));fwrite(hr,file.path(out,'numeric/training_HR.csv'));fwrite(ph,file.path(out,'numeric/training_joint_PH.csv'));fwrite(rbindlist(allCor),file.path(out,'numeric/test_train_adjusted_correlations.csv'));fwrite(rbindlist(allCal),file.path(out,'numeric/training_age_calibration.csv'));fwrite(rbindlist(calchecks),file.path(out,'qa/Python_R_gap_checks.csv'));fwrite(rbindlist(allCoh),file.path(out,'numeric/cohorts.csv'));fwrite(rbindlist(drows),file.path(out,'numeric/panel_d_Cindex.csv'));fwrite(rbindlist(mdrows),file.path(out,'numeric/panel_d_matched_Cindex.csv'))
if(length(warnings_out))fwrite(rbindlist(warnings_out),file.path(out,'qa/Cox_warnings.csv'))
saveRDS(fitted_private,file.path(out,'PRIVATE/comparison_models.rds'));saveRDS(list(matched=matched,selected=selected),file.path(out,'PRIVATE/matched_scores.rds'))
write_json(list(strict_new_train_only_pipeline=TRUE,shared_stratified_split_seed=20261006,target_train_fraction=.7,original_publication_numeric_reproduction=FALSE,original_date_source_audit_passed=FALSE,training_Cox_only=TRUE,heldout_Cindex_models=nrow(ci),fair_models=54,d_display_models=54,matched_available=sum(rbindlist(mdrows)$available),matched_undefined=sum(!rbindlist(mdrows)$available),paired_contrasts=nrow(pairs),training_HR_terms=nrow(hr),train_N=length(trix),train_deaths=sum(meta$status[trix]),test_N=length(teix),test_deaths=sum(meta$status[teix]),reference_N=length(ds$reference),reference_test_overlap=0,all_direct_C_and_SE_checks_pass=TRUE,all_direct_matched_score_checks_pass=TRUE,Python_R_age_gap_checks_pass=TRUE,Cox_warning_count=length(warnings_out),CI_excludes_retraining_uncertainty=TRUE,current_or_original_figures_replaced=FALSE),file.path(out,'qa/comparison_complete.json'),pretty=TRUE,auto_unbox=TRUE)
pf<-list.files(file.path(out,'PRIVATE'),full.names=TRUE,recursive=FALSE);Sys.chmod(pf[!file.info(pf)$isdir],'0600');Sys.chmod(pf[file.info(pf)$isdir],'0700');writeLines(capture.output(sessionInfo()),file.path(out,'logs/comparison_sessionInfo.txt'));sink()
