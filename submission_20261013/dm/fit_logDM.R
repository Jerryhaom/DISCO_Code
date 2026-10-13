#!/usr/bin/env Rscript
# Harmonize DM HRs to natural-log DM per SD. Retain all existing risk sets,
# covariates and frozen modality scaling populations used in Fig3c/Table16.
options(warn=1)
if (nzchar(Sys.getenv('DISCO_R_LIBS'))) .libPaths(c(Sys.getenv('DISCO_R_LIBS'), .libPaths()))
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1)
root <- normalizePath(Sys.getenv('DISCO_RESULTS_ROOT'), mustWork=TRUE)
module <- normalizePath(Sys.getenv('DISCO_RUN_OUTPUT'), mustWork=TRUE)
for(d in c('input','results','tables','logs','qa','output/pdf','staging'))
  dir.create(file.path(module,d),recursive=TRUE,showWarnings=FALSE)
private <- file.path(root,'03_multimodal_UKB/inputs_private_20261007/private_materialized/data/hm/UKB')
paths <- c(clinical=file.path(private,'month_disease_prs_execution_20261007/input_clinical.csv'),
  proteomics=file.path(private,'month_disease_prs_execution_20261007/input_proteomics.csv'),
  clinicalprot=file.path(private,'targeted_endpoint_diagnostics_20261007.JaNZnJ/input_clinicalprot.csv'),
  endpoint=file.path(private,'month_disease_prs_execution_20261007/private_endpoint_8.csv'))
expected <- read_json(file.path(root,'03_multimodal_UKB/targeted_endpoint_diagnostics_20261007/three_modalities_execution_QA.json'))$source_hashes
hashes <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
stopifnot(all(hashes==unlist(expected[names(paths)])))
oldpath <- file.path(root,'03_multimodal_UKB/Fig3c_PA_adjusted_20261010/tables/Fig3c_8cov_perSD_HR_48.csv')
old <- fread(oldpath); stopifnot(nrow(old)==48L)
flow <- fread(file.path(root,'03_multimodal_UKB/Fig3c_PA_adjusted_20261010/tables/Fig3c_8cov_sample_flow_24.csv'))
ep <- fread(paths['endpoint'],select=c('eid','disease','time','days','event','tie','reason'))
diseases <- c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')
labels <- c(clinical='Clinical',proteomics='Prot',clinicalprot='ClinicalProt')
res <- scales <- comparisons <- coefficients <- ph_results <- warnings <- list()
fig <- copy(old)
for(mod in names(labels)) {
  m <- fread(paths[mod]); stopifnot(!anyDuplicated(m$eid),all(is.finite(m$hd)),all(m$hd>0))
  ldm <- log(m$hd); mu <- mean(ldm); sdev <- sd(ldm)
  m[,z_logDM := (ldm-mu)/sdev]
  # Exact algebraic equivalence of logged distance and logged squared distance.
  equiv <- max(abs(as.numeric(scale(ldm))-as.numeric(scale(2*ldm))))
  stopifnot(equiv<1e-12)
  scales[[mod]] <- data.table(modality=mod,metric='DM',transform='natural log',N=nrow(m),mean=mu,SD=sdev,
    SD_population='full saved modality scoring population before endpoint/covariate exclusions',
    squared_distance_equivalence_error=equiv)
  m[,sex:=factor(sex)]
  m[,smoke:=factor(ifelse(smoke_status<0,NA,smoke_status),levels=0:2)]
  m[,drink:=factor(ifelse(alcohol_freq<0,NA,alcohol_freq),levels=1:6,ordered=TRUE)]
  m[,edu:=factor(ifelse(edu_q0==-7,7,ifelse(edu_q0<0,NA,edu_q0)),levels=1:7)]
  m[,race:=factor(ifelse(is.na(race)|race=='','White',race),levels=c('White','Non-White'))]
  m[,bmi:=ifelse(bmi_raw<0,NA,bmi_raw)]
  stopifnot(all(na.omit(m$ipaq_group)%in%0:2))
  m[,activity:=factor(ipaq_group,levels=0:2,labels=c('Low','Moderate','High'))]
  for(dis in diseases) {
    e <- ep[disease==dis]; ix <- match(m$eid,e$eid);stopifnot(!anyNA(ix),!anyDuplicated(e$eid))
    dd <- copy(m)
    for(v in c('time','days','event','tie','reason'))set(dd,j=v,value=e[[v]][ix])
    dd <- dd[is.na(reason)|reason==''];dd[,status:=as.integer(event==1)]
    covs <- c('age','sex','bmi','smoke','drink','edu','race','activity')
    adj <- dd[complete.cases(dd[,c('time','status','mc','hd',covs),with=FALSE]) & is.finite(time) & time>0]
    prior <- old[modality==mod & disease==dis & metric=='DM']
    fl <- flow[modality==mod & disease==dis]
    stopifnot(nrow(prior)==1L,nrow(fl)==1L,
      nrow(adj)==prior$N,sum(adj$status)==prior$events,
      nrow(dd)==prior$N_crude,sum(dd$status)==prior$events_crude,
      digest(adj$eid,algo='sha256')==fl$current_risk_set_id_sha256)
    for(model in c('Crude','Adjusted')) {
      d <- if(model=='Crude')dd else adj
      rhs <- c('z_logDM',if(model=='Adjusted')covs)
      fml <- as.formula(paste('Surv(time,status)~',paste(rhs,collapse='+')))
      fit <- withCallingHandlers(coxph(fml,data=d,ties='efron',singular.ok=FALSE,x=TRUE,y=TRUE,
        control=coxph.control(eps=1e-11,iter.max=60)),warning=function(w) {
          warnings[[length(warnings)+1L]] <<- data.table(modality=mod,disease=dis,model=model,message=conditionMessage(w))
          invokeRestart('muffleWarning')
        })
      sm <- summary(fit);b <- unname(coef(fit)['z_logDM']);se <- unname(sqrt(vcov(fit)['z_logDM','z_logDM']))
      cc <- unname(sm$concordance);p <- 2*pnorm(-abs(b/se))
      stopifnot(fit$n==nrow(d),fit$nevent==sum(d$status),is.finite(b),is.finite(se),se>0)
      row <- data.table(modality=mod,display_modality=unname(labels[mod]),disease=dis,metric='DM',model=model,
        N=fit$n,events=fit$nevent,beta=b,se=se,HR=exp(b),low=exp(b-qnorm(.975)*se),high=exp(b+qnorm(.975)*se),P=p,
        log10_P=(log(2)+pnorm(-abs(b/se),log.p=TRUE))/log(10),
        C=cc[1],C_se=cc[2],C_low=max(0,cc[1]-qnorm(.975)*cc[2]),C_high=min(1,cc[1]+qnorm(.975)*cc[2]),
        SD=sdev,score_mean=mu,SD_input_N=nrow(m),score_column='z_logDM',unit='per 1 SD of natural-log DM',
        SD_population='full saved modality scoring population before endpoint/covariate exclusions',
        formula=paste(deparse(fml),collapse=''),risk_set_sha256=digest(d$eid,algo='sha256'))
      res[[length(res)+1L]] <- row
      coefficients[[length(coefficients)+1L]] <- data.table(modality=mod,disease=dis,model=model,
        term=names(coef(fit)),beta=unname(coef(fit)),se=sqrt(diag(vcov(fit))))
      if(model=='Crude') {
        # A strictly monotone score transform leaves the score-only ranking unchanged.
        stopifnot(abs(row$C-prior$C_crude)<1e-10)
      } else {
        idx <- which(fig$modality==mod & fig$disease==dis & fig$metric=='DM')
        for(v in c('beta','se','HR','low','high','P','SD','score_column'))set(fig,i=idx,j=v,value=row[[v]])
        set(fig,i=idx,j='C_adjusted',value=row$C);set(fig,i=idx,j='C_adjusted_se',value=row$C_se)
        set(fig,i=idx,j='covariates',value=paste(rhs,collapse='+'))
        ph <- cox.zph(fit,transform='km',terms=TRUE,global=TRUE)$table
        ph_results[[length(ph_results)+1L]] <- data.table(modality=mod,disease=dis,metric='DM',term=rownames(ph),
          chisq=ph[,'chisq'],df=ph[,'df'],P=ph[,'p'])
        comparisons[[length(comparisons)+1L]] <- data.table(modality=mod,disease=dis,N=fit$n,events=fit$nevent,
          old_HR=prior$HR,new_HR=row$HR,old_low=prior$low,new_low=row$low,old_high=prior$high,new_high=row$high,
          old_P=prior$P,new_P=row$P,old_C_adjusted=prior$C_adjusted,new_C_adjusted=row$C,
          crude_C_retained=prior$C_crude)
      }
      rm(fit);invisible(gc(FALSE))
    }
    cat('Completed',mod,dis,'N',nrow(adj),'events',sum(adj$status),'\n');flush.console()
    fwrite(rbindlist(res),file.path(module,'tables/DM_logHR_models_48.csv'))
  }
  rm(m);invisible(gc(FALSE))
}
final <- rbindlist(res); stopifnot(nrow(final)==48L)
fig[,score_transform:=ifelse(metric=='DM','natural log of saved positive DM','unchanged DISCO')]
stopifnot(identical(fig[metric=='DISCO',names(old),with=FALSE],old[metric=='DISCO']))
stopifnot(identical(fig$C_crude,old$C_crude),identical(fig$N,old$N),identical(fig$events,old$events))
fwrite(fig,file.path(module,'tables/Fig3c_8cov_perSD_HR_48.csv'))
fwrite(rbindlist(scales),file.path(module,'tables/logDM_scaling_3.csv'))
fwrite(rbindlist(comparisons),file.path(module,'tables/DM_HR_before_after_24.csv'))
fwrite(rbindlist(coefficients),file.path(module,'tables/DM_logHR_coefficients.csv'))
ph <- rbindlist(ph_results);ph[,q_BH:=p.adjust(P,'BH'),by=term];fwrite(ph,file.path(module,'tables/DM_logHR_PH_tests.csv'))
if(length(warnings))fwrite(rbindlist(warnings),file.path(module,'logs/model_warnings.csv'))
stopifnot(identical(hashes,vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))))
write_json(list(status='DM_LOG_TRANSFORM_HR_UNIFIED',models=48L,adjusted_models=24L,crude_models=24L,
  covariates=c('age','sex','BMI','smoking','drinking','education','ethnicity','physical activity'),
  endpoint_and_risk_sets_unchanged=TRUE,DISCO_rows_unchanged=TRUE,crude_C_index_unchanged=TRUE,
  SD_population='previously specified full modality scoring population retained',
  DM_transform='natural logarithm before SD standardization',HRS_log_squared_equivalent_after_standardization=TRUE,
  input_hashes=as.list(hashes),inputs_unchanged=TRUE,warning_count=length(warnings)),
  file.path(module,'results/fit_verification.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(module,'logs/sessionInfo.txt'))
cat('COMPLETE: 48 log-DM models; all original risk sets and crude C-indices verified.\n')
