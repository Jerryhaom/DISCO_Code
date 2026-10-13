#!/usr/bin/env Rscript
# Reconcile actual participant-level evidence. Export summaries only.
Sys.umask('0077')
suppressPackageStartupMessages({library(readxl);library(survival);library(jsonlite);library(digest)})
root <- Sys.getenv('DISCO_RESULTS_INPUT')
out <- Sys.getenv('DISCO_RUN_MODULE')
for(p in c('input','results','tables','logs','private/input'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
Sys.chmod(file.path(out,'private'),'0700');Sys.chmod(file.path(out,'private/input'),'0700')
prepared <- file.path(root,'04_external_validation/CHARLS_followup_sensitivity_20261007/private/input/primary_prepared_CHARLS_frames.rds')
sources <- c(prepared,
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/CHARLS.2011.disco.rds'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/2015/Health_Status_and_Functioning_2013.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/2015/Health_Status_and_Functioning_2015.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/Function/Health_Status_and_Functioning_2018.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/Mortality/Death_2013.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/Mortality/Death_2015.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/Mortality/Death_2018.xlsx'),
 file.path(Sys.getenv('DISCO_INPUT_ROOT'), 'Result5/Frailty/Mortality/Death_2020.xlsx'))
before <- vapply(sources,function(p)digest(file=p,algo='sha256'),character(1))
read <- function(p) suppressMessages(as.data.frame(read_xlsx(p,col_types='text',.name_repair='unique')))
key <- function(x) trimws(as.character(x))
d <- readRDS(prepared)$prepared_table1_common_frame
e <- readRDS(sources[2]);m <- match(d$id,key(e$id));stopifnot(!anyNA(m),identical(d$id1,key(e$id1[m])))
h13 <- read(sources[3]);h15 <- read(sources[4]);h18 <- read(sources[5])
x13 <- read(sources[6]);s15 <- read(sources[7]);s18 <- read(sources[8]);x20 <- read(sources[9])
dups <- do.call(rbind,lapply(seq_along(list(h13,h15,h18,x13,s15,s18,x20)),function(i){
 t<-list(h13,h15,h18,x13,s15,s18,x20)[[i]];k<-key(t$ID)
 data.frame(source=basename(sources[i+2]),N=nrow(t),missing_ID=sum(is.na(k)|k==''),duplicate_ID_rows=sum(duplicated(k)),duplicate_nonmissing_ID_rows=sum(duplicated(k)&!is.na(k)&k!=''))
}))
write.csv(dups,file.path(out,'results/source_ID_duplicate_audit.csv'),row.names=FALSE)
for(t in list(x13,s15,s18,x20))stopifnot(!anyDuplicated(key(t$ID)))
status15 <- s15$died[match(d$id1,key(s15$ID))];status18 <- s18$died[match(d$id1,key(s18$ID))]
dead13 <- d$id1%in%key(x13$ID);dead15 <- !is.na(status15)&status15=='1 Died';dead18 <- !is.na(status18)&status18=='1 Died'
alive13 <- d$id1%in%key(h13$ID);alive15 <- !is.na(status15)&status15=='0 Alive';alive18 <- !is.na(status18)&status18=='0 Alive'
health15 <- d$id1%in%key(h15$ID);health18 <- d$id1%in%key(h18$ID)
conflicts <- data.frame(rule=c('dead2013_and_health2013','dead2013_and_alive2015','dead2013_and_alive2018','dead2015_and_alive2018','dead2015_and_health2015','dead2018_and_health2018'),
 N=c(sum(dead13&alive13),sum(dead13&alive15),sum(dead13&alive18),sum(dead15&alive18),sum(dead15&health15),sum(dead18&health18)))
write.csv(conflicts,file.path(out,'results/status_evidence_conflicts.csv'),row.names=FALSE)
baseline_alive <- rep(TRUE,nrow(d))
alive_matrix <- cbind(baseline2011=baseline_alive,observed2013=alive13,observed2015=alive15,observed2018=alive18)
dead_matrix <- cbind(dead2013=dead13,dead2015=dead15,dead2018=dead18)
wave_grid <- c(0,2,4,7) # explicitly nominal observation-wave elapsed scale, not fabricated individual dates.
last_alive <- apply(alive_matrix,1,function(z)max(wave_grid[z]))
first_dead <- apply(dead_matrix,1,function(z)if(any(z))min(c(2,4,7)[z]) else Inf)
endpoint_conflict <- is.finite(first_dead)&last_alive>=first_dead
lower <- last_alive;upper <- first_dead
no_followup <- !is.finite(upper)&lower==0
# Keep original historical event definitions for reconciliation. No change made to old objects.
stopifnot(identical(as.numeric(is.finite(first_dead)),as.numeric(d$status)),sum(endpoint_conflict)==0L)
wave <- data.frame(wave=c('baseline2011','2013','2015','2018'),alive_evidence_N=colSums(alive_matrix),
 died_evidence_N=c(0,colSums(dead_matrix)),nominal_time_years=wave_grid)
write.csv(wave,file.path(out,'results/observed_wave_status_counts.csv'),row.names=FALSE)
by_interval <- as.data.frame(table(lower_nominal=lower,upper_nominal=ifelse(is.infinite(upper),'right_censored',as.character(upper))),stringsAsFactors=FALSE)
by_interval <- by_interval[by_interval$Freq>0,]
write.csv(by_interval,file.path(out,'results/interval_endpoint_counts.csv'),row.names=FALSE)
excluded1088 <- d$status==0 & !alive18
rescue <- as.data.frame(table(old_2018_unverified=excluded1088,last_observed_alive_nominal=last_alive),stringsAsFactors=FALSE)
rescue <- rescue[rescue$Freq>0,]
write.csv(rescue,file.path(out,'results/unverified_2018_last_alive_recovery.csv'),row.names=FALSE)
validnum <- function(x,lo,hi){y<-suppressWarnings(as.numeric(x));y[is.na(y)|y<lo|y>hi]<-NA_real_;y}
calendar13 <- x13$exb002[match(d$id1,key(x13$ID))]
year13 <- validnum(x13$exb001_1[match(d$id1,key(x13$ID))],2011,2014)
month13 <- validnum(x13$exb001_2[match(d$id1,key(x13$ID))],1,12)
calendar20 <- x20$exb002[match(d$id1,key(x20$ID))]
year20 <- validnum(x20$exb001_1[match(d$id1,key(x20$ID))],2011,2020)
month20 <- validnum(x20$exb001_2[match(d$id1,key(x20$ID))],1,12)
day20 <- validnum(x20$exb001_3[match(d$id1,key(x20$ID))],1,31)
dt_counts <- data.frame(source=c('2013_exit','2020_exit'),cohort_matched_N=c(sum(dead13),sum(d$id1%in%key(x20$ID))),
 valid_death_year_N=c(sum(!is.na(year13)),sum(!is.na(year20))),valid_death_year_month_N=c(sum(!is.na(year13)&!is.na(month13)),sum(!is.na(year20)&!is.na(month20))),
 valid_death_year_month_day_N=c(0,sum(!is.na(year20)&!is.na(month20)&!is.na(day20))))
write.csv(dt_counts,file.path(out,'results/raw_death_date_availability.csv'),row.names=FALSE)
calcounts <- function(x,label){t<-as.data.frame(table(calendar=x,useNA='ifany'));t$source<-label;t[,c('source','calendar','Freq')]}
cal_counts <- rbind(calcounts(calendar13,'2013_exit'),calcounts(calendar20,'2020_exit'))
cal_counts <- cal_counts[cal_counts$Freq>0,]
write.csv(cal_counts,file.path(out,'results/death_calendar_counts.csv'),row.names=FALSE)
date15valid <- !is.na(validnum(s15$iyear[match(d$id1,key(s15$ID))],2015,2016))&!is.na(validnum(s15$imonth[match(d$id1,key(s15$ID))],1,12))
date18valid <- !is.na(validnum(s18$iyear[match(d$id1,key(s18$ID))],2018,2018))&!is.na(validnum(s18$imonth[match(d$id1,key(s18$ID))],1,12))
dd <- d;dd$lower_nominal <- lower;dd$upper_nominal <- ifelse(is.infinite(upper),NA_real_,upper);dd$endpoint_conflict <- endpoint_conflict;dd$no_postbaseline_followup <- no_followup
dd$alive13 <- alive13;dd$alive15 <- alive15;dd$alive18 <- alive18;dd$calendar13 <- calendar13;dd$year13 <- year13;dd$month13 <- month13
dd$calendar20 <- calendar20;dd$year20 <- year20;dd$month20 <- month20;dd$day20 <- day20
saveRDS(dd,file.path(out,'private/input/wave_interval_endpoint_frame.rds'));Sys.chmod(file.path(out,'private/input/wave_interval_endpoint_frame.rds'),'0600')
saveRDS(list(d15=s15,d18=s18,x13=x13,x20=x20),file.path(out,'private/input/actual_date_sources.rds'));Sys.chmod(file.path(out,'private/input/actual_date_sources.rds'),'0600')
write.csv(data.frame(path=sources,sha256=before),file.path(out,'input/endpoint_source_manifest.csv'),row.names=FALSE)
stopifnot(identical(before,vapply(sources,function(p)digest(file=p,algo='sha256'),character(1))))
write_json(list(status='OBSERVED_WAVE_INTERVAL_ENDPOINTS_BUILT_NOT_EXACT_PERSONAL_DATES',N=nrow(d),events=sum(is.finite(first_dead)),
 interval_endpoint_conflict_N=sum(endpoint_conflict),no_postbaseline_alive_or_death_evidence_N=sum(no_followup),
 interval_eligible_N=sum(!no_followup),censored_N=sum(!is.finite(first_dead)&!no_followup),
 old_2018_unverified_N=sum(excluded1088),recovered2015_alive_N=sum(excluded1088&last_alive==4),recovered2013_alive_N=sum(excluded1088&last_alive==2),
 evidence2013='Actual Health_Status_and_Functioning_2013 individual questionnaire IDs, not death-only unmatched presumed alive',
 evidence2015='Explicit died=0 in full2015 mortality table; full health questionnaire checked for conflicts',
 evidence2018='Explicit died=0 in full2018 mortality table; full health questionnaire checked for conflicts',
 nominal_scale=c('baseline2011=0','2013=2','2015=4','2018=7'),
 missing_2011_individual_blood_or_visit_dates=TRUE,
 valid2015_personal_interview_year_month_N=sum(date15valid),valid2018_personal_interview_year_month_N=sum(date18valid),
 death_dates_not_assumed_equal_interview_dates=TRUE,raw_death_calendar_must_be_verified_before_Gregorian_conversion=TRUE,
 old_source_statuses_match_new_first_death_evidence=TRUE,original_files_unchanged=TRUE,participant_values_exported_publicly=FALSE),
 file.path(out,'results/endpoint_construction.json'),auto_unbox=TRUE,pretty=TRUE)
cat('N',nrow(d),'events',sum(is.finite(first_dead)),'eligible',sum(!no_followup),'no_followup',sum(no_followup),'conflicts',sum(endpoint_conflict),'\n')
cat('Recovered among1088: 2015alive',sum(excluded1088&last_alive==4),'2013alive',sum(excluded1088&last_alive==2),'baselineonly',sum(excluded1088&last_alive==0),'\n')
