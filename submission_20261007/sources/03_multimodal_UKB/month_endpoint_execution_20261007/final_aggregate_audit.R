suppressPackageStartupMessages(library(jsonlite))
out <- commandArgs(trailingOnly=TRUE)[1]
d <- read.csv(file.path(out,'analysis_input.csv'));e <- read.csv(file.path(out,'private_endpoint.csv'))
stopifnot(identical(d$eid,e$eid));keep <- e$reason==''
reason_names <- c('missing_raw_ID','missing_blood_ID','inconsistent_tube_dates','no_valid_blood_date','baseline_or_prebaseline_death','nonpositive_followup')
counts <- setNames(lapply(reason_names,function(r)sum(e$reason==r)),reason_names)
write_json(list(exclusion_counts=counts,
 included_actual_blood_differs_from_visit=sum(keep&e$blood_visit_different==1),
 included_old_nonpositive_followup=sum(keep&d$time<=0),
 included_old_zero_followup_deaths=sum(keep&d$time==0&d$status==1),
 included_followup_changed_vs_old_years_times_12=sum(keep&abs(e$time-d$time*12)>1e-9,na.rm=TRUE),
 excluded_old_events=sum(d$status[!keep]==1),
 new_months_median=median(e$time[keep]),
 policy='1900-01-01 tube values ignored; use sole distinct valid tube date; isolate multiple valid dates; no assessment-date imputation'),
 file.path(out,'additional_aggregate_QA.json'),pretty=TRUE,auto_unbox=TRUE)
