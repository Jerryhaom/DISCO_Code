#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Targeted diagnostic only. Preserve original score scaling, covariates, endpoint and variables.
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1)
args <- commandArgs(TRUE); out <- args[1]
base <- .disco_path("private","month_disease_prs_execution_20261007")
paths <- file.path(base,c('input_proteomics.csv','private_endpoint_28.csv','PRS_input.csv','fit_models.R'))
before <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
m <- fread(paths[1]);stopifnot(!anyDuplicated(m$eid))
m[,zD:=as.numeric(scale(mc))];m[,sex:=factor(sex)]
m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
d <- copy(m[sex=='Female'])
prs <- fread(paths[3],select=c('eid','26232-0.0'));stopifnot(!anyDuplicated(prs$eid))
p <- prs[['26232-0.0']][match(d$eid,prs$eid)];d[,zPRS:=as.numeric(scale(p))]
# Read through Python once to retain only ovarian endpoint rows in RAM.
cmd <- sprintf("python3 -c %s", shQuote(paste0(
 "import csv,sys; f=open(", dQuote(paths[2], q=FALSE),
 "); r=csv.DictReader(f); w=csv.DictWriter(sys.stdout,fieldnames=r.fieldnames); w.writeheader(); ",
 "[w.writerow(x) for x in r if x['disease']=='ovarian_ca']")))
ep <- fread(cmd=cmd);stopifnot(!anyDuplicated(ep$eid))
ix <- match(d$eid,ep$eid);stopifnot(!anyNA(ix))
for(v in c('time','event','tie','reason'))set(d,j=v,value=ep[[v]][ix])
d <- d[is.na(reason)|reason==''];d[,status:=as.integer(event==1)]
rhs <- 'zD*zPRS+age+smoke+drink+edu+race'
fml <- as.formula(paste('Surv(time,status)~',rhs))
need <- all.vars(fml);d <- d[complete.cases(d[,..need]) & is.finite(time) & time>=0]
stopifnot(nrow(d)==25139,sum(d$status)==133)
old <- fread(file.path(base,'proteomics_28/coefficients.csv'))[disease=='ovarian_ca' & policy=='death_first' & model=='f1_interaction']
old_joint <- fread(file.path(base,'proteomics_28/models.csv'))[disease=='ovarian_ca' & policy=='death_first' & model=='f1_interaction' & term=='Joint_11_vs_00']
configs <- list(default=list(),tight_11=list(eps=1e-11,iter.max=60),tight_12=list(eps=1e-12,iter.max=80,toler.chol=1e-14),
                old_initial=list(eps=1e-11,iter.max=60),days=list(eps=1e-11,iter.max=60))
rows <- list();coefrows <- list();warnrows <- list();gradientrows <- list();fits <- list()
for(label in names(configs)) {
  data <- copy(d);if(label=='days')data[,time:=time*(365.25/12)]
  warnings <- character(); ctrl <- do.call(coxph.control,configs[[label]])
  call_args <- list(formula=fml,data=data,ties='efron',singular.ok=FALSE,x=TRUE,y=TRUE,control=ctrl)
  if(label=='old_initial')call_args$init <- old$beta
  f <- withCallingHandlers(do.call(coxph,call_args),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  b <- coef(f);V <- vcov(f);det <- coxph.detail(f);score <- colSums(det$score)
  delta <- as.numeric(V%*%score);names(delta)<-names(b)
  eigenvalues <- eigen(V,symmetric=TRUE,only.values=TRUE)$values
  L <- setNames(rep(0,length(b)),names(b));L[c('zD','zPRS','zD:zPRS')]<-1
  joint <- sum(L*b);jointse <- sqrt(as.numeric(t(L)%*%V%*%L))
  rows[[label]] <- data.table(run=label,N=f$n,events=f$nevent,iterations=f$iter,warning_count=length(warnings),
    warnings=paste(warnings,collapse=' | '),loglik=tail(f$loglik,1),max_abs_score=max(abs(score)),
    max_abs_Newton_step=max(abs(delta)),max_Newton_step_in_SE=max(abs(delta/sqrt(diag(V)))),
    covariance_min_eigenvalue=min(eigenvalues),covariance_condition=max(eigenvalues)/min(eigenvalues),
    max_beta_diff_from_saved=max(abs(b-old$beta)),interaction_beta=unname(b['zD:zPRS']),
    interaction_SE=unname(sqrt(diag(V))['zD:zPRS']),Joint_HR=exp(joint),Joint_low=exp(joint-qnorm(.975)*jointse),
    Joint_high=exp(joint+qnorm(.975)*jointse),Joint_P=2*pnorm(-abs(joint/jointse)),
    Joint_HR_diff_from_saved=exp(joint)-old_joint$HR)
  coefrows[[label]] <- data.table(run=label,term=names(b),beta=unname(b),SE=sqrt(diag(V)))
  gradientrows[[label]] <- data.table(run=label,term=names(score),score=unname(score),Newton_step=delta)
  warnrows[[label]] <- if(length(warnings))data.table(run=label,message=warnings) else NULL
  fits[[label]] <- f
}
tab <- rbindlist(rows);fwrite(tab,file.path(out,'ovarian_diagnostic_summary.csv'))
fwrite(rbindlist(coefrows),file.path(out,'ovarian_diagnostic_coefficients.csv'))
fwrite(rbindlist(gradientrows),file.path(out,'ovarian_diagnostic_gradients.csv'))
if(any(lengths(warnrows)))fwrite(rbindlist(warnrows),file.path(out,'ovarian_diagnostic_warnings.csv'))
ph <- tryCatch(as.data.table(cox.zph(fits[['tight_11']],transform='km')$table,keep.rownames='term'),error=function(e)data.table(error=conditionMessage(e)))
fwrite(ph,file.path(out,'ovarian_targeted_PH.csv'))
current <- fits[['tight_11']];b<-coef(current);V<-vcov(current);s<-summary(current)
stopifnot(tab[run=='tight_11',warning_count]==0L,tab[run=='tight_11',max_beta_diff_from_saved]<1e-6,
          tab[run=='tight_11',covariance_min_eigenvalue]>0)
current_rows <- lapply(c('zD','zPRS','zD:zPRS','Joint_11_vs_00'),function(term){
 L<-setNames(rep(0,length(b)),names(b));if(term=='Joint_11_vs_00')L[c('zD','zPRS','zD:zPRS')]<-1 else L[term]<-1
 estimate<-sum(L*b);se<-sqrt(as.numeric(t(L)%*%V%*%L))
 data.table(modality='proteomics',task='28',disease='ovarian_ca',policy='death_first',version='raw_dates',
            model='f1_interaction',term=term,N=current$n,events=current$nevent,beta=estimate,se=se,
            HR=exp(estimate),low=exp(estimate-qnorm(.975)*se),high=exp(estimate+qnorm(.975)*se),P=2*pnorm(-abs(estimate/se)),
            C=unname(s$concordance[1]),C_se=unname(s$concordance[2]),covariates=rhs,
            numerical_control='eps=1e-11; iter.max=60; otherwise unchanged coxph defaults',
            diagnostic_flag='default premature-convergence warning reproduced; finite stable solution confirmed',
            historical_results_replaced=FALSE)
})
fwrite(rbindlist(current_rows),file.path(out,'ovarian_strict_current_models.csv'))
vv<-as.data.table(as.table(V));setnames(vv,c('term1','term2','covariance'))
vv[,`:=`(modality='proteomics',task='28',disease='ovarian_ca',policy='death_first',version='raw_dates',model='f1_interaction')]
fwrite(vv,file.path(out,'ovarian_strict_current_full_vcov.csv'))
after <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1));stopifnot(identical(before,after))
write_json(list(scope='targeted ovarian numerical diagnostic; not replacement of main results',
    original_inputs_unchanged=TRUE,source_hashes=as.list(before),endpoint='saved raw_dates death_first',
    variables_unchanged=TRUE,score_scaling_unchanged=TRUE,model_formula=paste(deparse(fml),collapse=''),
    cox_code_sha256=digest(file=normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),algo='sha256')),
    file.path(out,'ovarian_diagnostic_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
print(tab)
