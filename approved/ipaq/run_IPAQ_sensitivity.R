#!/usr/bin/env Rscript
# UKB clinical mortality: current Table 1 core adjustment with/without IPAQ.
suppressPackageStartupMessages({library(survival); library(jsonlite); library(digest)})
options(stringsAsFactors=FALSE, digits=17)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: Rscript --vanilla run_IPAQ_sensitivity.R DATA_ROOT NEW_OUTPUT_DIRECTORY")
data_root <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
if (file.exists(out) || dir.exists(out)) stop("Output already exists; use a new directory")
entry <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)), mustWork=TRUE)
here <- dirname(entry)
reference_path <- file.path(Sys.getenv("DISCO_APPROVED_CODE_ROOT"),"table1/data/Table1_analysis_sample_SD.csv")
inputs <- c(file.path(data_root,"UKB_Data",c("extended_covariates_clinical.rds","covariates_raw_full.rds","ukb_race_21000.csv")),
            reference_path,entry,file.path(here,"protocol_before_run.txt"))
stopifnot(all(file.exists(inputs)))
hash_file <- function(p) digest(file=p,algo="sha256")
hash_before <- vapply(inputs,hash_file,character(1))
u <- readRDS(inputs[1]); raw <- readRDS(inputs[2]); race <- read.csv(inputs[3],colClasses="character")
stopifnot(!anyNA(u$ID),!anyDuplicated(u$ID),!anyNA(raw$eid),!anyDuplicated(raw$eid),
          !anyNA(race$eid),!anyDuplicated(race$eid))
ri <- match(u$ID,raw$eid); ei <- match(u$ID,race$eid)
stopifnot(!anyNA(ri),!anyNA(ei),all(u$sex %in% c("Male","Female")),
          isTRUE(all.equal(u$ipaq_group,raw$ipaq_group[ri],tolerance=0,check.attributes=FALSE)))
mode_impute <- function(x) {
  counts <- table(x)
  if (!length(counts)) stop("No observed values for mode imputation")
  modes <- names(counts)[counts==max(counts)]
  if (length(modes)!=1L) stop("Tied modes require explicit review")
  x[is.na(x)] <- modes[1]
  x
}
race_raw <- race$race[ei]
stopifnot(all(is.na(race_raw) | race_raw %in% c("","White","Non-White")))
race_missing <- is.na(race_raw) | race_raw==""
race_filled <- race_raw; race_filled[race_missing] <- "White"
edu <- as.numeric(raw$edu_q0[ri])
d <- data.frame(mc=as.numeric(u$mc),log_hd=log(as.numeric(u$hd)),time=as.numeric(u$time),
                status=as.numeric(u$status),age=as.numeric(u$age),sex=ifelse(u$sex=="Male",1,0))
d$smoke <- mode_impute(ifelse(u$smoke_status %in% c(1,2),1,ifelse(u$smoke_status==0,0,NA)))
d$drink <- mode_impute(ifelse(u$alcohol_freq==6,0,ifelse(u$alcohol_freq %in% 1:5,1,NA)))
d$education <- mode_impute(ifelse(edu %in% c(1,2),1,ifelse(edu %in% c(-7,3,4,5,6),0,NA)))
d$race <- ifelse(race_filled=="White","White","Non-White")
for (v in c("sex","smoke","drink","education","race")) d[[v]] <- factor(d[[v]])
ipa_raw <- as.numeric(u$ipaq_group)
stopifnot(all(is.na(ipa_raw) | ipa_raw %in% 0:2),all(d$status %in% c(0,1)),
          all(is.finite(d$time)),all(d$time>=0),all(is.finite(d$mc)),all(is.finite(d$log_hd)))
d$ipaq <- factor(ipa_raw,levels=0:2,labels=c("Low","Moderate","High"))
core <- c("age","sex","smoke","drink","education","race")
base_keep <- complete.cases(d[c("time","status","mc","log_hd",core)])
same_keep <- base_keep & !is.na(d$ipaq)
stopifnot(all(base_keep),sum(same_keep)>0)
cohort_hash <- function(keep) digest(paste(sort(as.character(u$ID[keep]),method="radix"),collapse="\n"),algo="sha256",serialize=FALSE)
counts <- data.frame(population=c("Table1_full","IPAQ_observed_common","IPAQ_missing_excluded"),
                     n=c(sum(base_keep),sum(same_keep),sum(base_keep & !same_keep)),
                     deaths=c(sum(d$status[base_keep]),sum(d$status[same_keep]),sum(d$status[base_keep & !same_keep])))
spec <- list(Table1_base_full_sample=list(keep=base_keep,covariates=core),
             Core_same_IPAQ_sample=list(keep=same_keep,covariates=core),
             Core_plus_IPAQ=list(keep=same_keep,covariates=c(core,"ipaq")))
reference <- read.csv(reference_path,stringsAsFactors=FALSE)
rows <- coefficient_rows <- list()
for (metric in c("DISCO","DM")) {
  variable <- if(metric=="DISCO") "mc" else "log_hd"
  for (model in names(spec)) {
    keep <- spec[[model]]$keep
    x <- d[keep,,drop=FALSE]
    score_sd <- sd(x[[variable]])
    score_mean <- mean(x[[variable]])
    x$score_z <- (x[[variable]]-score_mean)/score_sd
    form <- reformulate(c("score_z",spec[[model]]$covariates),response="Surv(time,status)")
    warnings <- character()
    fit <- withCallingHandlers(coxph(form,data=x,ties="efron",na.action=na.fail,
                                    x=FALSE,y=FALSE,model=FALSE),
      warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
    sm <- summary(fit)
    b <- unname(coef(fit)["score_z"]); se <- unname(sqrt(diag(vcov(fit)))["score_z"])
    ci <- unname(confint(fit)["score_z",])
    p <- unname(sm$coefficients["score_z","Pr(>|z|)"])
    log_p <- log(2)+pnorm(-abs(b/se),log.p=TRUE)
    stopifnot(fit$n==sum(keep),fit$nevent==sum(x$status),all(is.finite(coef(fit))),
              is.finite(score_sd),score_sd>0,abs(sd(x$score_z)-1)<1e-12)
    key <- paste(metric,model,sep="__")
    rows[[key]] <- data.frame(cohort="UKB clinical",metric=metric,model=model,n=fit$n,deaths=fit$nevent,
      HR_per_SD=exp(b),CI95_lower=exp(ci[1]),CI95_upper=exp(ci[2]),beta_per_SD=b,SE_per_SD=se,
      P=p,log_P=log_p,neg_log10_P=-log_p/log(10),P_report=if(log_p<log(1e-300)) "<1e-300" else sprintf("%.3e",p),
      HR_report=sprintf("%.3f (%.3f–%.3f)",exp(b),exp(ci[1]),exp(ci[2])),
      score_SD=score_sd,score_mean=score_mean,score_variable=variable,
      formula=paste(deparse(form),collapse=""),covariates=paste(spec[[model]]$covariates,collapse="+"),
      adjusted_model_C_index=unname(sm$concordance[1]),adjusted_model_C_index_SE=unname(sm$concordance[2]),
      cohort_sha256=cohort_hash(keep),warning_count=length(warnings),warnings=paste(unique(warnings),collapse=" | "),
      stringsAsFactors=FALSE)
    cc <- sm$coefficients
    coefficient_rows[[key]] <- data.frame(metric=metric,model=model,term=rownames(cc),
       beta=cc[,"coef"],se=cc[,"se(coef)"],P=cc[,"Pr(>|z|)"],row.names=NULL)
    if(model=="Table1_base_full_sample") {
      old <- reference[reference$cohort=="UKB clinical" & reference$metric==metric,,drop=FALSE]
      stopifnot(nrow(old)==1L,fit$n==old$n,fit$nevent==old$events,
        abs(exp(b)-old$adjusted_HR_per_analysis_SD)<1e-10,
        abs(exp(ci[1])-old$adjusted_HR_CI95_lower)<1e-10,
        abs(exp(ci[2])-old$adjusted_HR_CI95_upper)<1e-10,
        abs(score_sd-old$score_SD_analysis_sample)<1e-12)
    }
    cat(sprintf("%s %s N=%d deaths=%d HR=%.9f [%.9f, %.9f] P=%s\n",metric,model,fit$n,fit$nevent,
                exp(b),exp(ci[1]),exp(ci[2]),rows[[key]]$P_report))
  }
}
models <- do.call(rbind,rows);rownames(models)<-NULL
coefficients <- do.call(rbind,coefficient_rows);rownames(coefficients)<-NULL
stopifnot(nrow(models)==6L,!anyDuplicated(paste(models$metric,models$model)))
comparison <- do.call(rbind,lapply(c("DISCO","DM"),function(metric) {
  a<-models[models$metric==metric & models$model=="Core_same_IPAQ_sample",]
  b<-models[models$metric==metric & models$model=="Core_plus_IPAQ",]
  stopifnot(a$n==b$n,a$deaths==b$deaths,identical(a$cohort_sha256,b$cohort_sha256),a$score_SD==b$score_SD)
  data.frame(metric=metric,n=b$n,deaths=b$deaths,HR_without_IPAQ=a$HR_per_SD,HR_with_IPAQ=b$HR_per_SD,
             percent_HR_change=100*(b$HR_per_SD/a$HR_per_SD-1),same_participants=TRUE,same_score_SD=TRUE)
}))
after <- vapply(inputs,hash_file,character(1));stopifnot(identical(hash_before,after))
dir.create(out,recursive=TRUE)
write.csv(models,file.path(out,"UKB_IPAQ_models.csv"),row.names=FALSE,na="NA")
write.csv(models[c("metric","model","n","deaths","HR_report","P_report")],
          file.path(out,"UKB_IPAQ_display.csv"),row.names=FALSE)
write.csv(coefficients,file.path(out,"coefficients.csv"),row.names=FALSE)
write.csv(comparison,file.path(out,"same_sample_comparison.csv"),row.names=FALSE)
write.csv(counts,file.path(out,"sample_flow.csv"),row.names=FALSE)
write.csv(data.frame(category=c("Low","Moderate","High"),code=0:2,
           n=as.integer(table(d$ipaq[same_keep])),deaths=as.integer(tapply(d$status[same_keep],d$ipaq[same_keep],sum))),
          file.path(out,"IPAQ_categories.csv"),row.names=FALSE)
write.csv(data.frame(path=inputs,sha256=hash_before),file.path(out,"source_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
write_json(list(status="COMPUTED_PENDING_INDEPENDENT_VERIFICATION",input_hashes_unchanged=TRUE,
  full_population_Table1_replication_passed=TRUE,IPAQ_matches_raw=TRUE,
  common_n=sum(same_keep),common_deaths=sum(d$status[same_keep]),IPAQ_missing_excluded=sum(base_keep & !same_keep),
  race_missing_imputed_White=sum(race_missing),same_sample_comparisons=TRUE,warning_count=sum(models$warning_count),
  source_field="UKB IPAQ activity group; 22032; coding100700",IPAQ_coding="0 Low / 1 Moderate / 2 High",
  main_Table1_unchanged=TRUE,PH_test_performed=FALSE,
  limitations=c("Complete-case analysis for IPAQ may be selective","Uses existing scores and mortality endpoints; does not rebuild their full raw-data chain",
                "Simple core-covariate imputation does not model imputation uncertainty","Statistical significance is not evidence of absence of residual confounding")),
  file.path(out,"audit.json"),pretty=TRUE,auto_unbox=TRUE)
