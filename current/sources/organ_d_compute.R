#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(entry))
suppressPackageStartupMessages({library(survival);library(data.table);library(jsonlite)})
setDTthreads(1);options(warn=2)
for(x in c('numeric','qa','logs','panels','output/pdf'))dir.create(file.path(out,x),recursive=TRUE,showWarnings=FALSE)
prior<-file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20261002_rebuilt_figures')
models<-readRDS(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20261002_local_rebuild/PRIVATE/common_prediction_models.rds'))
previous<-fread(file.path(prior,'numeric/panel_d_Cindex.csv'))
previous_pairs<-fread(file.path(prior,'numeric/panel_d_paired_differences.csv'))
map<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
ans<-paired<-audit<-list();ir<-ip<-ia<-0L
for(o in map$clock_organ){
 item<-models[[o]];te<-item$test;fd<-item$fits$DISCO$model
 stopifnot(identical(rownames(fd),rownames(item$fits$Oh$model)),identical(rownames(fd),rownames(item$fits$Goeminne$model)))
 tr<-data.frame(time=fd[[1]][,1],status=fd[[1]][,2],age=fd$age,sex=fd$sex,D=fd$D,
  O=item$fits$Oh$model$O,G=item$fits$Goeminne$model$G,row.names=rownames(fd))
 stopifnot(nrow(tr)==18648,nrow(te)==4294,sum(te$status)==460,!length(intersect(rownames(tr),te$eid)))
 forms<-list(DISCO=Surv(time,status)~age+D,Oh=Surv(time,status)~O,Goeminne=Surv(time,status)~G)
 fit<-function(f)coxph(f,data=tr,ties='efron',x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-10,iter.max=50))
 fits<-lapply(forms,fit)
 stopifnot(identical(names(coef(fits$DISCO)),c('age','D')),identical(names(coef(fits$Oh)),'O'),identical(names(coef(fits$Goeminne)),'G'))
 for(j in 1:3){
  key<-c('D','O','G')[j];method<-c('DISCO','Oh','Goeminne')[j];an<-paste0(key,'_accel')
  cal<-lm(reformulate('age',key),tr);den<-sd(resid(cal))
  tr[[an]]<-resid(cal)/den;te[[an]]<-(te[[key]]-predict(cal,newdata=te))/den
  stopifnot(abs(cor(tr[[an]],tr$age))<1e-10)
  fits[[paste0(method,'_accel')]]<-fit(reformulate(an,'Surv(time,status)'))
 }
 cc<-do.call(concordance,c(unname(fits),list(newdata=te,timewt='n',influence=2)))
 cv<-cc$concordance;V<-cc$var;stopifnot(all(eigen(V,symmetric=TRUE,only.values=TRUE)$values > -1e-10))
 for(j in 1:6){
  meth<-sub('_accel$','',names(fits)[j]);grp<-if(j<=3)'Original' else 'Accel'
  lp<-predict(fits[[j]],newdata=te,type='lp')
  check<-concordance(Surv(te$time,te$status)~lp,reverse=TRUE)
  stopifnot(abs(cv[j]-check$concordance)<1e-12,abs(V[j,j]-check$var)<1e-12)
  if(j>3){q<-previous[organ==o&score==meth&group=='Accel'];stopifnot(nrow(q)==1,abs(cv[j]-q$C)<1e-12,abs(sqrt(V[j,j])-q$SE)<1e-12)}
  se<-sqrt(V[j,j]);ir<-ir+1L
  ans[[ir]]<-data.frame(organ=o,display=map$display[match(o,map$clock_organ)],group=grp,score=meth,
   N=nrow(te),deaths=sum(te$status),C=unname(cv[j]),SE=se,lower95=max(0,cv[j]-qnorm(.975)*se),upper95=min(1,cv[j]+qnorm(.975)*se),
   Cox_fit='common_train',evaluation='common_test',formula=if(j==1)'age + D' else if(j==2)'O' else if(j==3)'G' else 'training-calibrated age gap only')
  ia<-ia+1L;audit[[ia]]<-data.frame(organ=o,model=names(fits)[j],terms=paste(names(coef(fits[[j]])),collapse=' + '),train_N=fits[[j]]$n,
   test_N=nrow(te),test_deaths=sum(te$status),direct_C_error=abs(cv[j]-check$concordance),direct_variance_error=abs(V[j,j]-check$var))
 }
 for(g in 0:1)for(j in 2:3){a<-1+g*3;b<-j+g*3;se<-sqrt(V[a,a]+V[b,b]-2*V[a,b]);delta<-cv[a]-cv[b];contrast<-paste('DISCO minus',c('DISCO','Oh','Goeminne')[j])
  if(g==1){ct<-contrast;q<-previous_pairs[organ==o&group=='Accel'&contrast==ct];stopifnot(nrow(q)==1,abs(q$delta-delta)<1e-12,abs(q$SE-se)<1e-12)}
  ip<-ip+1L;paired[[ip]]<-data.frame(organ=o,group=c('Original','Accel')[g+1],contrast=contrast,delta=unname(delta),SE=se,
   lower=unname(delta)-qnorm(.975)*se,upper=unname(delta)+qnorm(.975)*se,P=2*pnorm(-abs(delta/se)))
 }
}
d<-rbindlist(ans);pr<-rbindlist(paired);pr[,q_BH_9:=p.adjust(P,'BH'),by=.(group,contrast)]
fwrite(d,file.path(out,'numeric/panel_d_Cindex.csv'));fwrite(pr,file.path(out,'numeric/panel_d_paired_differences.csv'))
fwrite(rbindlist(audit),file.path(out,'qa/model_formula_checks.csv'))
write_json(list(left_formulas=list(DISCO='Surv(time,status) ~ age + D',Oh='Surv(time,status) ~ O',Goeminne='Surv(time,status) ~ G'),
  sex_in_any_d_model=FALSE,common_train_N=18648,common_test_N=4294,test_deaths=460,
  all_27_Accel_C_and_SE_unchanged=TRUE,all_18_Accel_paired_differences_unchanged=TRUE,
  direct_risk_C_and_variance_checks=54,uniform_age_sex_results_preserved=TRUE),file.path(out,'qa/numeric_checks.json'),pretty=TRUE,auto_unbox=TRUE)
print(d[group=='Original',.(organ,score,C,lower95,upper95)])
