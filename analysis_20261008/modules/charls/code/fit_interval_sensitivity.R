#!/usr/bin/env Rscript
# True interval-censored PH sensitivity. No assumed exact death dates.
Sys.umask('0077')
args <- commandArgs(TRUE)
run <- if(length(args)) args[1] else 'primary'
B <- if(length(args)>1) as.integer(args[2]) else 200L
scenario <- if(length(args)>2) args[3] else 'observed_wave'
endpoint_description <- if(scenario=='retrospective_strict')
 'Observed alive/dead wave intervals plus 10 unambiguous solar-calendar deaths before2018 recovered from2020 exit; 1 conflicting report excluded' else
 '(last confirmed alive wave, first confirmed dead wave]; last-alive right censoring'
root <- Sys.getenv('DISCO_RESULTS_INPUT')
out <- Sys.getenv('DISCO_RUN_MODULE')
.libPaths(c(file.path(out,'rlib'),.libPaths()))
suppressPackageStartupMessages({library(icenReg);library(survival);library(jsonlite);library(digest);library(doParallel);library(doRNG)})
script_path <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
script_sha <- digest(file=script_path,algo='sha256')
for(p in c('results','logs','private/models')) dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
Sys.chmod(file.path(out,'private/models'),'0700')
source_path <- file.path(out,'private/input',if(scenario=='retrospective_strict')'retrospective_strict_interval_frame.rds' else 'wave_interval_endpoint_frame.rds')
source_sha <- digest(file=source_path,algo='sha256')
full <- readRDS(source_path)
fixed_sd <- c(DISCO=sd(full$mc),DM=sd(full$log_hd))
fixed_mean <- c(DISCO=mean(full$mc),DM=mean(full$log_hd))
excluded <- if('analysis_excluded'%in%names(full)) full$analysis_excluded else full$no_postbaseline_followup
d <- full[!excluded,]
d$age <- d$age - mean(full$age) # numerical centering only; same linear-age covariate.
stopifnot(nrow(d)>11000L,sum(d$status)>=1096L,!any(d$endpoint_conflict))
warning_rows <- list(); results <- list(); timing <- list()
capture <- function(expr,metric,method) withCallingHandlers(expr,warning=function(w){
 warning_rows[[length(warning_rows)+1L]] <<- data.frame(metric,method,warning=conditionMessage(w));invokeRestart('muffleWarning')
})
if(B>0){
 cluster <- parallel::makeCluster(4L)
 parallel::clusterCall(cluster,function(lib){.libPaths(c(lib,.libPaths()));suppressPackageStartupMessages(library(icenReg));TRUE},file.path(out,'rlib'))
 registerDoParallel(cluster)
}
for(metric in c('DISCO','DM')){
 d$score_z <- (if(metric=='DISCO')d$mc else d$log_hd)-fixed_mean[metric]
 d$score_z <- d$score_z/fixed_sd[metric]
 model_seed <- 20261008L + match(metric,c('DISCO','DM'))
 set.seed(model_seed)
 if(B>0){registerDoParallel(cluster);registerDoRNG(model_seed)}
 t0 <- proc.time()[3]
 fit <- capture(ic_sp(Surv(lower_nominal,upper_nominal,type='interval2') ~ score_z+age+sex+smoke+drink+education,
                      data=d,model='ph',bs_samples=B,useMCores=B>0,B=c(0,1)),metric,'semiparametric_interval_PH')
 timing[[metric]] <- unname(proc.time()[3]-t0)
 beta <- coef(fit)['score_z']
 se <- if(B>0)sqrt(diag(vcov(fit)))['score_z'] else NA_real_
 results[[length(results)+1L]] <- data.frame(run,metric,method='semiparametric_interval_PH',N=nrow(d),deaths=sum(d$status),
    primary_SD=fixed_sd[metric],beta=beta,SE=se,HR=exp(beta),lower95=exp(beta-1.95996398454005*se),upper95=exp(beta+1.95996398454005*se),
    P=if(B>0)2*pnorm(-abs(beta/se)) else NA_real_,CI_method=if(B>0)paste0('Wald normal; covariance from ',B,' individual bootstrap resamples') else 'Point fit only',
    bootstrap_samples=B,ordinary_Harrell_C=NA_real_,covariates='age + sex + smoke + drink + education',
    endpoint=endpoint_description,stringsAsFactors=FALSE)
 saveRDS(fit,file.path(out,'private/models',paste0(run,'_',metric,'_interval_PH.rds')))
 # Parametric checks use the nominal wave spacing, not individual elapsed dates.
 for(dist in c('weibull','exponential')){
   af <- capture(survreg(Surv(lower_nominal,upper_nominal,type='interval2')~score_z+age+sex+smoke+drink+education,data=d,dist=dist),metric,paste0(dist,'_interval_PH'))
   cb <- coef(af)['score_z'];sc <- af$scale;v <- vcov(af)
   gradient <- rep(0,nrow(v));names(gradient)<-rownames(v);gradient['score_z'] <- -1/sc
   if('Log(scale)'%in%names(gradient)) gradient['Log(scale)'] <- cb/sc
   ph_beta <- -cb/sc;ph_se <- sqrt(drop(t(gradient)%*%v%*%gradient))
   results[[length(results)+1L]] <- data.frame(run,metric,method=paste0(dist,'_interval_PH_nominal_spacing'),N=nrow(d),deaths=sum(d$status),
      primary_SD=fixed_sd[metric],beta=ph_beta,SE=ph_se,HR=exp(ph_beta),lower95=exp(ph_beta-1.95996398454005*ph_se),upper95=exp(ph_beta+1.95996398454005*ph_se),
      P=2*pnorm(-abs(ph_beta/ph_se)),CI_method='Delta-method Wald CI including estimated Weibull log-scale when applicable',
      bootstrap_samples=0,ordinary_Harrell_C=NA_real_,covariates='age + sex + smoke + drink + education',
      endpoint=paste(endpoint_description,'; parametric sensitivity uses nominal 0/2/4/7 spacing, not precise calendar time'),stringsAsFactors=FALSE)
   saveRDS(af,file.path(out,'private/models',paste0(run,'_',metric,'_',dist,'_interval_PH.rds')))
 }
 # Descriptive concordance only among provably ordered pairs; no imputed death time.
 score <- if(metric=='DISCO')d$mc else d$log_hd
 numerator <- denominator <- 0
 for(i in which(!is.na(d$upper_nominal))){
   j <- which(d$lower_nominal >= d$upper_nominal[i])
   denominator <- denominator+length(j)
   numerator <- numerator+sum(score[i]>score[j])+0.5*sum(score[i]==score[j])
 }
 write.csv(data.frame(metric,N=nrow(d),deaths=sum(d$status),ordered_pair_N=denominator,
   concordant_plus_half_tied=numerator,identifiable_pair_concordance=numerator/denominator,
   definition='Event i is definitely earlier only when upper_i <= lower_j; score ties receive 0.5. Not ordinary Harrell C or an interval-concordance CI.'),
   file.path(out,'results',paste0(run,'_',metric,'_identifiable_pair_concordance.csv')),row.names=FALSE)
}
if(B>0)parallel::stopCluster(cluster)
result <- do.call(rbind,results);rownames(result)<-NULL
write.csv(result,file.path(out,'results',paste0(run,'_interval_PH_results.csv')),row.names=FALSE)
write.csv(if(length(warning_rows))do.call(rbind,warning_rows) else data.frame(metric=character(),method=character(),warning=character()),
          file.path(out,'results',paste0(run,'_model_warnings.csv')),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'logs',paste0(run,'_sessionInfo.txt')))
write_json(list(run=run,source_sha256=source_sha,original_source_unchanged=identical(source_sha,digest(file=source_path,algo='sha256')),
  executed_code_sha256=script_sha,executed_code_unchanged_during_run=identical(script_sha,digest(file=script_path,algo='sha256')),
  full_primary_N=nrow(full),evaluation_N=nrow(d),deaths=sum(d$status),scenario=scenario,
  excluded_baseline_only_N=sum(full$no_postbaseline_followup),excluded_total_N=sum(excluded),
  fixed_primary_SD=as.list(fixed_sd),fixed_primary_mean=as.list(fixed_mean),B=B,seed_DISCO=20261009L,seed_DM=20261010L,
  package_versions=lapply(c('icenReg','MLEcens','survival','Rcpp','RcppEigen','doParallel','doRNG'),function(p)list(package=p,version=as.character(packageVersion(p)))),
  reproducible_parallel_bootstrap='doRNG per-iteration streams, metric-specific fixed seed; not worker scheduling RNG',
  elapsed_seconds=timing,warning_count=length(warning_rows),source_contains_precise_individual_baseline_dates=FALSE,
  primary_Table1_Meta_ED3_ED8_not_modified=TRUE),file.path(out,'results',paste0(run,'execution.json')),auto_unbox=TRUE,pretty=TRUE,digits=16)
print(result[,c('metric','method','N','deaths','HR','lower95','upper95','P')],row.names=FALSE)
