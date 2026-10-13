#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(survival);library(jsonlite);library(digest)})
args <- commandArgs(trailingOnly=TRUE);out <- args[1]
paths <- file.path(out,c('analysis_input.csv','private_endpoint.csv'))
hash <- function() vapply(paths,function(p) digest(file=p,algo='sha256'),character(1))
before <- hash()
d <- read.csv(paths[1],colClasses=c(eid='character'));e <- read.csv(paths[2],colClasses=c(eid='character'))
stopifnot(!anyDuplicated(d$eid),!anyDuplicated(e$eid),identical(d$eid,e$eid))
score_guard <- digest(d[c('eid','mc','hd')],algo='sha256')
# Independent endpoint calculation in R using day-resolution dates.
for(v in c('visit','blood','loss','death0','death1','death','end')) e[[v]] <- as.Date(e[[v]])
cutoff <- as.Date('2022-12-31')
death <- pmin(e$death0,e$death1,na.rm=TRUE);death[is.infinite(as.numeric(death))] <- NA
stopifnot(identical(as.numeric(death),as.numeric(e$death)))
eligible <- e$reason==''
expected_status <- as.integer(!is.na(death)&death<=cutoff)
expected_end <- pmin(e$loss,cutoff,na.rm=TRUE)
expected_end[expected_status==1] <- death[expected_status==1]
expected_days <- as.numeric(expected_end-e$blood)
stopifnot(all(e$status[eligible]==expected_status[eligible]),
 all(e$end[eligible]==expected_end[eligible]),all(e$days[eligible]==expected_days[eligible]),
 all(abs(e$time[eligible]-expected_days[eligible]/(365.25/12))<1e-10),
 all(e$days[eligible]>0),all(is.na(death[eligible])|death[eligible]>e$blood[eligible]),
 all(is.na(e$blood[e$reason=='no_valid_blood_date'])),
 all(death[e$reason=='baseline_or_prebaseline_death']<=e$blood[e$reason=='baseline_or_prebaseline_death']),
 all(expected_days[e$reason=='nonpositive_followup']<=0))
for(v in c('mc','hd','time','status','age','sex')) d[[v]] <- as.numeric(unclass(d[[v]]))
d$status[!d$status %in% c(0,1)] <- NA
d$log_hd <- log(d$hd)
covs <- c('age','sex','smoke','drink','education','race')
for(v in covs[-1]) d[[v]] <- factor(d[[v]])
old <- d
keep_old <- complete.cases(old[c('mc','time','status',covs)])
stopifnot(identical(keep_old,complete.cases(old[c('log_hd','time','status',covs)])))
d$time <- e$time;d$status <- e$status
complete <- complete.cases(d[c('mc','time','status',covs)])
stopifnot(identical(complete,complete.cases(d[c('log_hd','time','status',covs)])))
keep <- eligible&complete
stopifnot(all(is.finite(d$mc[keep])),all(is.finite(d$log_hd[keep])),all(is.finite(d$time[keep])))
new <- d[keep,];oldfit <- old[keep_old,]
fit <- function(data,label,version) {
 metric <- if(label=='DISCO') 'mc' else 'log_hd'
 s <- sd(data[[metric]]);data$score_z <- as.numeric(scale(data[[metric]]))
 warnings <- character()
 run <- function(f) withCallingHandlers(coxph(f,data=data,ties='efron'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 fml <- as.formula(paste('Surv(time,status)~',paste(c('score_z',covs),collapse='+')))
 f <- run(fml);crude <- run(Surv(time,status)~score_z)
 b <- unname(coef(f)['score_z']);se <- unname(sqrt(diag(vcov(f)))['score_z']);ci <- exp(b+c(-1,1)*qnorm(.975)*se)
 cc <- unname(summary(crude)$concordance)
 direct <- concordance(Surv(time,status)~score_z,data=data,reverse=coef(crude)[1]>0)
 stopifnot(abs(cc[1]-direct$concordance)<1e-12,abs(cc[2]-sqrt(direct$var))<1e-12)
 # A second fit in raw days verifies month rescaling preserves all event ordering/ties.
 if(version=='new') {
   data$time <- data$time*(365.25/12);dayfit <- run(fml)
   stopifnot(max(abs(coef(dayfit)-coef(f)))<1e-10,max(abs(vcov(dayfit)-vcov(f)))<1e-10)
 }
 data.frame(cohort='UKB clinical',metric=label,version=version,n=f$n,events=f$nevent,model='A',
 adjusted_HR_per_analysis_SD=exp(b),adjusted_HR_CI95_lower=ci[1],adjusted_HR_CI95_upper=ci[2],
 beta=b,SE=se,adjusted_P=2*pnorm(-abs(b/se)),crude_C_index=cc[1],crude_C_index_SE=cc[2],
 crude_C_index_CI95_lower=max(0,cc[1]-qnorm(.975)*cc[2]),crude_C_index_CI95_upper=min(1,cc[1]+qnorm(.975)*cc[2]),
 score_SD_analysis_sample=s,score_mean_analysis_sample=mean(data[[metric]]),
 standardized_formula=paste(deparse(fml),collapse=''),
 scaling='One within-cohort sample SD in the effective Table 1 Cox analysis population; DM log-transformed before scaling',
 warning_count=length(warnings),warnings=paste(unique(warnings),collapse=' | '))
}
tab <- do.call(rbind,lapply(c('DISCO','DM'),function(m) rbind(fit(oldfit,m,'old'),fit(new,m,'new'))))
write.csv(tab,file.path(out,'UKB_models_old_new.csv'),row.names=FALSE)
reason <- e$reason;reason[reason==''] <- 'endpoint_eligible'
flow <- as.data.frame(table(endpoint_reason=reason,old_model_included=keep_old,new_model_included=keep))
flow <- flow[flow$Freq>0,]
write.csv(flow,file.path(out,'N_change_decomposition.csv'),row.names=FALSE)
transition <- as.data.frame(table(old_included=keep_old,new_included=keep,old_event=old$status,new_event=e$status,useNA='ifany'))
write.csv(transition[transition$Freq>0,],file.path(out,'event_change_decomposition.csv'),row.names=FALSE)
stopifnot(identical(before,hash()),identical(score_guard,digest(d[c('eid','mc','hd')],algo='sha256')))
json <- list(independent_R_endpoint_QA_passed=TRUE,day_vs_month_cox_equal=TRUE,crude_C_direct_check_passed=TRUE,
 scores_unchanged=TRUE,remote_input_hashes=as.list(before),score_ID_object_sha256=score_guard,
 input_n=nrow(d),old_n=sum(keep_old),old_events=sum(old$status[keep_old]),new_n=sum(keep),new_events=sum(d$status[keep]),
 retained_n=sum(keep&keep_old),added_n=sum(keep&!keep_old),removed_n=sum(!keep&keep_old),
 endpoint_eligible_n=sum(eligible),endpoint_eligible_covariate_or_score_missing_n=sum(eligible&!complete),
 common_retained_event_0_to_1=sum(keep&keep_old&old$status==0&d$status==1,na.rm=TRUE),
 common_retained_event_1_to_0=sum(keep&keep_old&old$status==1&d$status==0,na.rm=TRUE),
 followup_month_quantiles=as.list(quantile(new$time,c(0,.25,.5,.75,1))),model_warning_count=sum(tab$warning_count))
write_json(json,file.path(out,'model_QA.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'remote_R_sessionInfo.txt'))
print(tab[c('metric','version','n','events','adjusted_HR_per_analysis_SD','adjusted_HR_CI95_lower','adjusted_HR_CI95_upper','crude_C_index')],row.names=FALSE)
