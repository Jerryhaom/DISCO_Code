options(stringsAsFactors=FALSE, warn=1)
if (nzchar(Sys.getenv('DISCO_R_LIBS'))) .libPaths(c(Sys.getenv('DISCO_R_LIBS'), .libPaths()))
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1)
root <- normalizePath(Sys.getenv('DISCO_RESULTS_ROOT'), mustWork=TRUE)
mod <- normalizePath(Sys.getenv('DISCO_RUN_OUTPUT'), mustWork=TRUE)
src <- c(scores=file.path(root,'03_multimodal_UKB/inputs_private_20261007/af73f791c853e41e01225259450bb82a8c72de5e069c65b50c378cd795d95642_scores_7030.csv.gz'),
         matched=file.path(root,'03_multimodal_UKB/inputs_private_20261007/private_materialized/data/hm/UKB/organ_endpoint_repairs_20261007_run01/matched_scores.csv'),
         endpoint=file.path(root,'03_multimodal_UKB/inputs_private_20261007/018ac88b4321c6cef8a8259367ab1fc0717122fa31dbe7a6978cf6f5d949208e_private_organ_endpoint.csv.gz'),
         covariates=file.path(root,'03_multimodal_UKB/inputs_private_20261007/6a881e52a35fb9d132657f93d787bed108a8e5733af1413ed50576e10480c139_input_organ.csv.gz'),
         calibration=file.path(root,'05_organ_network/endpoint_repairs_20261007/output/clock/training_calibration.csv'),
         prior_HR=file.path(root,'05_organ_network/ED10bc_mortality_fully_adjusted_20261010/results/training_HR_fully_adjusted_72.csv'))
hashes <- vapply(src,function(p)digest(file=p,algo='sha256'),'')
s <- fread(src['scores'], colClasses=c(eid='character'))
mt <- fread(src['matched'], colClasses=c(eid='character'))
stopifnot(nrow(s)==44831,nrow(mt)==44831,!anyDuplicated(s$eid),!anyDuplicated(mt$eid),sum(s$is_train)==31381,sum(s$is_test)==13450)
mi <- match(s$eid,mt$eid);stopifnot(!anyNA(mi))
cov <- fread(src['covariates'],select=c('eid','age','sex','bmi_raw','race','edu_q0','smoke_status','alcohol_freq'),colClasses=c(eid='character'))
ci <- match(s$eid,cov$eid)
stopifnot(!anyNA(ci),max(abs(s$age-cov$age[ci]))<1e-9,all(as.character(s$sex)==ifelse(cov$sex[ci]=='Female','0','1')))
s[,`:=`(sex=factor(cov$sex[ci]),bmi_raw=cov$bmi_raw[ci],race_raw=cov$race[ci],edu_q0=cov$edu_q0[ci],smoke_status=cov$smoke_status[ci],alcohol_freq=cov$alcohol_freq[ci])]
s[,race:=factor(ifelse(is.na(race_raw)|race_raw=='','White',race_raw),levels=c('White','Non-White'))]
s[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
s[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
s[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
s[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
covars <- c('age','sex','bmi','race','edu','smoke','drink')
e <- fread(src['endpoint'],select=c('eid','time','event','reason'),colClasses=c(eid='character'))
ei <- match(s$eid,e$eid);stopifnot(!anyNA(ei),!anyDuplicated(e$eid))
s[,`:=`(time=e$time[ei],event=as.integer(e$event[ei]))]
stopifnot(all(s$time>0),all(s$event%in%0:1),sum(s$event[s$is_train==1])==3318,sum(s$event[s$is_test==1])==1423)
cal <- fread(src['calibration']);prior <- fread(src['prior_HR'])
organs <- c('Adipose','Brain','Heart','Intestine','Kidney','Liver','Lung','Muscle','Pancreas')
cal_checks <- list()
for(o in organs)for(m in c('Oh','Goeminne')){
  pc <- if(m=='Oh')'O'else'G'
  for(code in c(paste0('D',pc),pc)){
    rawcol <- if(code==pc)paste(m,o,sep='_')else paste('M',m,o,sep='_')
    values <- if(code==pc)s[[rawcol]]else mt[[rawcol]][mi]
    cc <- cal[organ==o & score==code];stopifnot(nrow(cc)==1,cc$train_N==31381,cc$test_N_used==0,all(is.finite(values)))
    f <- lm(values[s$is_train==1]~s$age[s$is_train==1])
    errors <- c(intercept=abs(coef(f)[1]-cc$age_intercept),slope=abs(coef(f)[2]-cc$age_slope),SD=abs(sd(resid(f))-cc$age_residual_SD))
    stopifnot(max(errors)<1e-9)
    s[[rawcol]] <- values
    s[[paste0('z_',rawcol)]] <- (values-cc$age_intercept-cc$age_slope*s$age)/cc$age_residual_SD
    cal_checks[[length(cal_checks)+1]] <- data.table(organ=o,method=m,score=code,max_error=max(errors))
  }
}
needs <- c('time','event',covars,grep('^z_',names(s),value=TRUE))
tr <- s[is_train==1 & complete.cases(s[,..needs])]
stopifnot(nrow(tr)==30590,sum(tr$event)==3180)
HR <- PH <- coefrows <- scalechecks <- clockchecks <- warnings <- list()
fit <- function(form,dat,key){withCallingHandlers(coxph(form,data=dat,ties='efron',x=TRUE,y=TRUE,model=TRUE,control=coxph.control(eps=1e-9,iter.max=40)),warning=function(w){warnings[[length(warnings)+1]]<<-data.table(model=key,warning=conditionMessage(w));invokeRestart('muffleWarning')})}
for(o in organs){
 for(m in c('Oh','Goeminne')){
  pc <- if(m=='Oh')'O'else'G';dc <- paste0('D',pc)
  dat <- copy(tr[,c('time','event',covars),with=FALSE])
  dat$D <- tr[[paste('M',m,o,sep='_')]];dat$P <- tr[[paste(m,o,sep='_')]]
  dat$Da <- tr[[paste0('z_M_',m,'_',o)]];dat$Pa <- tr[[paste0('z_',m,'_',o)]]
  term_list <- list(D_only='Da',P_only='Pa',joint=c('Da','Pa'))
  for(model in names(term_list)){
   terms <- term_list[[model]];key <- paste(o,m,model,sep='|')
   form <- reformulate(c(covars,terms),'Surv(time,event)')
   f <- fit(form,dat,key);stopifnot(all(is.finite(coef(f))))
   rf <- fit(reformulate(c(covars,sub('a$','',terms)),'Surv(time,event)'),dat,paste(key,'raw_scale',sep='|'))
   for(term in terms){
    b <- unname(coef(f)[term]);se <- sqrt(vcov(f)[term,term]);raw <- sub('a$','',term)
    cc <- cal[organ==o & score==if(raw=='D')dc else pc]
    err <- c(beta=abs(b-unname(coef(rf)[raw])*cc$age_residual_SD),SE=abs(se-sqrt(vcov(rf)[raw,raw])*cc$age_residual_SD),loglik=abs(tail(f$loglik,1)-tail(rf$loglik,1)))
    stopifnot(max(err[c('beta','SE')])<1e-8,err['loglik']<1e-6)
    scalechecks[[length(scalechecks)+1]] <- data.table(organ=o,method=m,model=model,term=term,beta_error=err['beta'],SE_error=err['SE'],loglik_error=err['loglik'])
    HR[[length(HR)+1]] <- data.table(organ=o,method=m,model=model,term=if(term=='Da')paste0(dc,'a')else paste0(pc,'a'),score=if(term=='Da')'Matched DISCO'else m,N=f$n,deaths=f$nevent,beta=b,SE=se,HR=exp(b),lower95=exp(b-qnorm(.975)*se),upper95=exp(b+qnorm(.975)*se),P=2*pnorm(-abs(b/se)),formula=paste(deparse(form),collapse=''),scale='training_age_residual_SD',covariates='age + sex + BMI + ethnicity + education + smoking + drinking')
    if(model=='P_only'){
     model_key <- model;old <- prior[organ==o & method==m & model==model_key]
     stopifnot(nrow(old)==1)
     error <- max(abs(c(b-old$beta,se-old$SE)));stopifnot(error<1e-8)
     clockchecks[[length(clockchecks)+1]] <- data.table(organ=o,method=m,max_beta_SE_error=error)
    }
   }
   cs <- summary(f)$coefficients
   coefrows[[length(coefrows)+1]] <- data.table(organ=o,method=m,model=model,term=rownames(cs),beta=cs[,'coef'],SE=cs[,'se(coef)'],P=cs[,'Pr(>|z|)'])
   ph <- cox.zph(f,transform='km')$table
   PH[[length(PH)+1]] <- data.table(organ=o,method=m,model=model,term=rownames(ph),chisq=ph[,'chisq'],df=ph[,'df'],P=ph[,'p'],N=f$n,deaths=f$nevent)
  }
 }
 cat('Completed',o,'matched-protein mortality models\n');flush.console()
}
h <- rbindlist(HR);ph <- rbindlist(PH)
stopifnot(nrow(h)==72,!anyDuplicated(h[,.(organ,method,model,term)]))
h[,q_BH_72:=p.adjust(P,'BH')]
fwrite(h,file.path(mod,'results/matched_mortality_HR_72.csv'))
fwrite(ph,file.path(mod,'results/matched_proportional_hazards.csv'))
fwrite(rbindlist(coefrows),file.path(mod,'results/matched_full_coefficients.csv'))
fwrite(rbindlist(scalechecks),file.path(mod,'qa/raw_scale_equivalence.csv'))
fwrite(rbindlist(clockchecks),file.path(mod,'qa/clock_single_model_identity.csv'))
fwrite(rbindlist(cal_checks),file.path(mod,'qa/calibration_identity.csv'))
if(length(warnings))fwrite(rbindlist(warnings),file.path(mod,'results/model_warnings.csv'))
stopifnot(identical(hashes,vapply(src,function(p)digest(file=p,algo='sha256'),'')))
j <- h[model=='joint' & score=='Matched DISCO']
write_json(list(status='PASS',N=30590,deaths=3180,training_calibration_N=31381,test_N_used_for_calibration=0,reference_N=344,models=54,HR_terms=72,covariates=covars,score_basis='per training age-residual SD',max_clock_replay_error=max(rbindlist(clockchecks)$max_beta_SE_error),max_calibration_error=max(rbindlist(cal_checks)$max_error),warning_count=length(warnings),joint_matched_DISCO_HR_range=range(j$HR),joint_matched_DISCO_CI_above_one=sum(j$lower95>1),joint_matched_DISCO_FDR_significant=sum(j$q_BH_72<.05),source_hashes=as.list(hashes)),file.path(mod,'results/matched_model_verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(mod,'results/R_sessionInfo.txt'))
print(j[,.(organ,method,HR,lower95,upper95,P,q_BH_72)])
