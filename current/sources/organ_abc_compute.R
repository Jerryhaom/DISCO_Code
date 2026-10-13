#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(entry))
suppressPackageStartupMessages({library(survival);library(data.table);library(jsonlite)})
setDTthreads(1);options(warn=2)
for(x in c('numeric','qa','logs','panels','output/pdf'))dir.create(file.path(out,x),recursive=TRUE,showWarnings=FALSE)
src<-file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20261002_local_rebuild')
stopifnot(fromJSON(file.path(src,'qa/final_audit.json'))$statistics_second_pass_identical)
C<-fread(file.path(src,'results/common_partial_correlations.csv'))
H<-fread(file.path(src,'results/common_HR.csv'))
oldC<-fread(file.path(src,'results/common_test_Cindex.csv'))
models<-readRDS(file.path(src,'PRIVATE/common_prediction_models.rds'))
map<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
organs<-map$clock_organ;stopifnot(setequal(organs,names(models)))
C[,`:=`(partial_r=r,main_figure=TRUE,adjustment='linear_age')]
hh<-H[model%in%c('D_only','P_only','joint')]
hh[,`:=`(term=ifelse(term=='Dz','D','P'),main_figure=TRUE,adjustment='linear_age')]
fwrite(C,file.path(out,'numeric/panel_a_correlations.csv'))
fwrite(hh[term=='D'],file.path(out,'numeric/panel_b_HR.csv'))
fwrite(hh[term=='P'],file.path(out,'numeric/panel_c_HR.csv'))
results<-pairs<-calibration<-independent<-list();ni<-nc<-np<-nr<-0L
for(o in organs){
 item<-models[[o]];te<-item$test
 fd<-item$fits$DISCO$model
 stopifnot(identical(rownames(fd),rownames(item$fits$Oh$model)),identical(rownames(fd),rownames(item$fits$Goeminne$model)))
 tr<-data.frame(time=fd[[1]][,1],status=fd[[1]][,2],age=fd$age,sex=fd$sex,D=fd$D,
  O=item$fits$Oh$model$O,G=item$fits$Goeminne$model$G,row.names=rownames(fd))
 stopifnot(nrow(tr)==18648,nrow(te)==4294,sum(te$status)==460,!length(intersect(rownames(tr),te$eid)))
 fits<-item$fits[c('DISCO','Oh','Goeminne')]
 for(j in seq_along(c('D','O','G'))){
  key<-c('D','O','G')[j];method<-c('DISCO','Oh','Goeminne')[j];an<-paste0(key,'_accel')
  cal<-lm(reformulate('age',key),tr)
  res<-resid(cal);den<-sd(res)
  tr[[an]]<-res/den;te[[an]]<-(te[[key]]-predict(cal,newdata=te))/den
  stopifnot(abs(cor(tr[[an]],tr$age))<1e-10)
  form<-reformulate(an,'Surv(time,status)')
  fits[[paste0(method,'_accel')]]<-coxph(form,data=tr,ties='efron',x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-10,iter.max=50))
  nc<-nc+1L;calibration[[nc]]<-data.frame(organ=o,score=method,train_N=nrow(tr),test_N=nrow(te),
   intercept=unname(coef(cal)[1]),age_slope=unname(coef(cal)[2]),training_residual_SD=den,
   training_residual_age_r=cor(tr[[an]],tr$age),test_gap_age_r=cor(te[[an]],te$age),
   definition='OLS score~age fitted only in common training; fixed coefficients applied to test; gap-only Cox')
 }
 cc<-do.call(concordance,c(unname(fits),list(newdata=te,timewt='n',influence=2)))
 cv<-cc$concordance;V<-cc$var
 stopifnot(all(eigen(V,symmetric=TRUE,only.values=TRUE)$values > -1e-10))
 for(j in seq_along(fits)){
  method<-sub('_accel$','',names(fits)[j]);group<-if(j<=3)'Age + sex' else 'Accel'
  direct<-concordance(Surv(te$time,te$status)~predict(fits[[j]],newdata=te,type='lp'),reverse=TRUE)
  stopifnot(abs(cv[j]-direct$concordance)<1e-12)
  if(j<=3){q<-oldC[organ==o&model==method];stopifnot(nrow(q)==1,abs(q$Cindex-cv[j])<1e-12,abs(q$SE-sqrt(V[j,j]))<1e-12)}
  nr<-nr+1L;se<-sqrt(V[j,j]);results[[nr]]<-data.frame(organ=o,display=map$display[match(o,organs)],group=group,score=method,
    N=nrow(te),deaths=sum(te$status),C=unname(cv[j]),SE=se,lower95=max(0,cv[j]-qnorm(.975)*se),upper95=min(1,cv[j]+qnorm(.975)*se),
    Cox_fit='common_train',evaluation='common_test',formula=if(j<=3)'age + sex + score' else 'training-calibrated age gap only')
 }
 for(g in 0:1)for(j in 2:3){a<-1+3*g;b<-j+3*g;se<-sqrt(V[a,a]+V[b,b]-2*V[a,b]);delta<-cv[a]-cv[b]
  np<-np+1L;pairs[[np]]<-data.frame(organ=o,group=c('Age + sex','Accel')[g+1],contrast=paste('DISCO minus',c('DISCO','Oh','Goeminne')[j]),
    delta=delta,SE=se,lower=delta-qnorm(.975)*se,upper=delta+qnorm(.975)*se,P=2*pnorm(-abs(delta/se)))
 }
 # Recompute all displayed correlations and HRs from the actual common test sample.
 for(method in c('Oh','Goeminne')){
  meth_loop<-method
  key<-if(method=='Oh')'O' else 'G'
  rr<-cor(resid(lm(D~age+sex,te)),resid(lm(reformulate(c('age','sex'),key),te)))
  stopifnot(abs(rr-C[organ==o&method==meth_loop,partial_r])<1e-12)
  h<-te;h$Dz<-as.numeric(scale(resid(lm(D~age,h))));h$Pz<-as.numeric(scale(resid(lm(reformulate('age',key),h))))
  for(model in c('D_only','P_only','joint')){
   mod_loop<-model
   term<-switch(model,D_only='Dz',P_only='Pz',joint=c('Dz','Pz'))
   fit<-coxph(reformulate(c('age','sex',term),'Surv(time,status)'),h,control=coxph.control(eps=1e-10,iter.max=50))
   for(t in term){q<-H[organ==o&method==meth_loop&model==mod_loop&term==t];stopifnot(nrow(q)==1)
    err<-abs(coef(fit)[t]-q$beta);stopifnot(err<1e-10)
    ni<-ni+1L;independent[[ni]]<-data.frame(organ=o,method=method,model=model,term=t,beta_error=unname(err))
   }
  }
 }
}
fwrite(rbindlist(results),file.path(out,'numeric/panel_d_Cindex.csv'))
pr<-rbindlist(pairs);pr[,q_BH_9:=p.adjust(P,'BH'),by=.(group,contrast)]
fwrite(pr,file.path(out,'numeric/panel_d_paired_differences.csv'))
fwrite(rbindlist(calibration),file.path(out,'qa/age_gap_calibration.csv'))
fwrite(rbindlist(independent),file.path(out,'qa/independent_HR_checks.csv'))
write_json(list(heatmap_values=18,HR_points=72,Cindex_points=54,common_test_N=4294,common_test_deaths=460,
  all_age_sex_Cindex_equal_verified_rebuild=TRUE,all_Cindex_directions_checked=TRUE,
  age_gap_calibration_training_only=TRUE,independent_HR_terms_checked=ni,old_outputs_overwritten=FALSE),
  file.path(out,'qa/numeric_checks.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: real data independently checked; panel data ready.\n')
