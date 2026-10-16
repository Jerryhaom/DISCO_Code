#!/usr/bin/env Rscript
# Retrospective deaths before all 2018 visits: solar-calendar-only, no alive evidence contradiction.
Sys.umask('0077')
suppressPackageStartupMessages({library(jsonlite);library(digest)})
root <- Sys.getenv('DISCO_RUN_MODULE')
d <- readRDS(file.path(root,'private/input/wave_interval_endpoint_frame.rds'))
s <- readRDS(file.path(root,'private/input/actual_date_sources.rds'))
map <- function(table,col) table[[col]][match(d$id1,as.character(table$ID))]
ym <- function(y,m){v<-suppressWarnings(as.numeric(y))*12+suppressWarnings(as.numeric(m));v[is.na(y)|is.na(m)]<-NA_real_;v}
ym15 <- ym(map(s$d15,'iyear'),map(s$d15,'imonth'))
ym18 <- ym(map(s$d18,'iyear'),map(s$d18,'imonth'))
solar <- !is.na(d$calendar20)&d$calendar20=='1 Solar'&!is.na(d$year20)&!is.na(d$month20)&!is.na(d$day20)
death_date <- as.Date(rep(NA_character_,nrow(d)))
death_date[solar] <- as.Date(sprintf('%04d-%02d-%02d',d$year20[solar],d$month20[solar],d$day20[solar]))
stopifnot(all(!is.na(death_date[solar])))
reported_ym <- ym(d$year20,d$month20)
before13 <- solar & d$alive13 & d$year20<2013
before15 <- solar & d$alive15 & !is.na(ym15) & reported_ym<ym15
before18 <- solar & d$alive18 & !is.na(ym18) & reported_ym<ym18
before13[is.na(before13)]<-FALSE;before15[is.na(before15)]<-FALSE;before18[is.na(before18)]<-FALSE
conflict <- before13|before15|before18
same15 <- solar&d$alive15&!is.na(ym15)&reported_ym==ym15
same18 <- solar&d$alive18&!is.na(ym18)&reported_ym==ym18
same15[is.na(same15)]<-FALSE;same18[is.na(same18)]<-FALSE
same_month_ambiguous <- same15|same18
# Use only dates after the latest possible baseline campaign year and before any 2018 visit.
# The model uses a conservative observation interval ending at the 2018 wave, not an imputed death day.
eligible_death <- solar & d$status==0 & d$year20>=2013 & d$year20<=2017 & !conflict & !same_month_ambiguous
eligible_death[is.na(eligible_death)]<-FALSE
stopifnot(all(d$lower_nominal[eligible_death]<7))
new <- d;new$status_legacy <- d$status
new$status[eligible_death] <- 1
new$upper_nominal[eligible_death] <- 7
new$retrospective_verified_solar_death <- eligible_death
new$retrospective_date_alive_conflict <- conflict
new$same_month_order_unresolved <- same_month_ambiguous
new$no_postbaseline_followup <- is.na(new$upper_nominal)&new$lower_nominal==0
new$analysis_excluded <- new$no_postbaseline_followup|conflict
stopifnot(identical(d$mc,new$mc),identical(d$log_hd,new$log_hd),all(new$lower_nominal[!new$analysis_excluded]<ifelse(is.na(new$upper_nominal[!new$analysis_excluded]),Inf,new$upper_nominal[!new$analysis_excluded])))
saveRDS(new,file.path(root,'private/input/retrospective_strict_interval_frame.rds'))
Sys.chmod(file.path(root,'private/input/retrospective_strict_interval_frame.rds'),'0600')
counts <- data.frame(rule=c('full2020_exit_records','full2020_nonmissing_year_month_day','full2020_solar_calendar','full2020_lunar_calendar',
 'cohort_matched2020_exit','cohort_valid_solar_exact_report','cohort_valid_lunar_exact_report',
 'solar_report_year_before2018','death_before_later_alive_evidence','same_interview_month_order_unresolved',
 'retrospective_pre2018_nonconflicting_deaths_used','2011_report_only_baseline_not_used','baseline_only_rescued_by_retrospective_death',
 'new_interval_analysis_N','new_interval_analysis_deaths','remaining_baseline_only_excluded','conflict_excluded'),
 N=c(nrow(s$x20),sum(!is.na(s$x20$exb001_1)&!is.na(s$x20$exb001_2)&!is.na(s$x20$exb001_3)),
 sum(s$x20$exb002=='1 Solar',na.rm=TRUE),sum(s$x20$exb002=='2 Lunar',na.rm=TRUE),
 sum(d$id1%in%as.character(s$x20$ID)),sum(solar),sum(d$calendar20=='2 Lunar'&!is.na(d$year20)&!is.na(d$month20)&!is.na(d$day20),na.rm=TRUE),
 sum(solar&d$year20<2018,na.rm=TRUE),sum(conflict),sum(same_month_ambiguous),sum(eligible_death),
 sum(solar&d$year20==2011&d$no_postbaseline_followup,na.rm=TRUE),sum(eligible_death&d$no_postbaseline_followup),
 sum(!new$analysis_excluded),sum(new$status[!new$analysis_excluded]),sum(new$no_postbaseline_followup),sum(conflict)))
write.csv(counts,file.path(root,'results/retrospective_death_date_audit.csv'),row.names=FALSE)
changes<-as.data.frame(table(old_last_alive=d$lower_nominal[eligible_death],old_event=d$status[eligible_death],new_upper=new$upper_nominal[eligible_death]));changes<-changes[changes$Freq>0,]
write.csv(changes,file.path(root,'results/retrospective_endpoint_changes_aggregate.csv'),row.names=FALSE)
write_json(list(status='STRICT_SOLAR_RETROSPECTIVE_INTERVAL_BRANCH_BUILT',
  death_count=sum(eligible_death),analysis_N=sum(!new$analysis_excluded),analysis_deaths=sum(new$status[!new$analysis_excluded]),
  conflicting_death_date_alive_evidence_N=sum(conflict),same_month_order_uncertain_N=sum(same_month_ambiguous),
  public_exact_death_dates_are_reported_not_administratively_certified=TRUE,
  no_lunar_conversion_or_midmonth_imputation_performed=TRUE,
  death_year2018_not_added_because_original_individual2018_visit_window_uncertain=TRUE,
  death_year2011_baseline_only_not_added_because_before_after_blood_not_provable=TRUE,
  final_interval='Known last alive wave to before2018 wave; exact report not converted to pretend exact elapsed time',
  score_primary_SD_unchanged=TRUE,primary_results_not_modified=TRUE),file.path(root,'results/retrospective_branch_construction.json'),auto_unbox=TRUE,pretty=TRUE)
print(counts,row.names=FALSE)
