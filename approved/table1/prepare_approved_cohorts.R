# Table 1 missing-data rules approved by the author on 2026-10-06.
# Score construction, observed covariate coding and model covariates are unchanged.
prepare_approved_cohorts <- function(results_root) {
  legacy_path <- file.path(Sys.getenv("DISCO_APPROVED_CODE_ROOT"),"table1/source_prepare.R")
  legacy <- new.env(parent=globalenv())
  preparation_warnings <- character()
  stopped <- FALSE
  for (expr in parse(legacy_path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("for"))) {
      stopped <- TRUE
      break
    }
    withCallingHandlers(eval(expr, legacy), warning=function(w) {
      preparation_warnings <<- c(preparation_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }
  stopifnot(stopped, length(legacy$cohorts)==5L)
  cohorts <- legacy$cohorts

  # Compute the mode in observed, binary-coded education among the scored
  # CHARLS cohort (not the reference group). Preserve missingness until filling.
  education_raw <- suppressWarnings(as.numeric(legacy$ch$edu))
  education_observed <- ifelse(is.na(education_raw), NA_real_,
    ifelse(education_raw %in% c(1,2,3), 0, 1))
  education_counts <- table(education_observed)
  stopifnot(length(education_counts)>0L)
  education_modes <- names(education_counts)[education_counts==max(education_counts)]
  # Do not silently choose between tied modes.
  stopifnot(length(education_modes)==1L)
  cohorts$CHARLS$d$education <- legacy$impute_mode(education_observed)

  # Keep an explicit missing flag: imputed White is not observed White.
  race_raw <- legacy$race$race[match(legacy$ukb$ID, legacy$race$eid)]
  race_missing <- is.na(race_raw) | race_raw==""
  race_filled <- race_raw
  race_filled[race_missing] <- "White"
  cohorts[["UKB clinical"]]$d$race <- ifelse(race_filled=="White", "White", "Non-White")

  # For the locked inputs the explicit approved rules are numerically identical
  # to the previous implementation. A changed input requires renewed review.
  unchanged <- vapply(names(cohorts), function(nm)
    identical(cohorts[[nm]]$d, legacy$cohorts[[nm]]$d), logical(1))
  stopifnot(all(unchanged))
  audit <- list(
    author_approved_date="2026-10-06",
    CHARLS=list(input_n=nrow(legacy$ch), missing_education_n=sum(is.na(education_observed)),
      mode_population="Observed binary-coded education in the full scored CHARLS cohort, before imputation",
      observed_low_n=unname(education_counts["0"]), observed_high_n=unname(education_counts["1"]),
      mode=unname(education_modes), rule="Cohort-mode single imputation"),
    UKB=list(input_n=nrow(legacy$ukb), missing_race_n=sum(race_missing),
      rule="Missing or empty race is singly imputed as White; observed categories are unchanged"),
    all_prepared_cohort_frames_identical_to_previous=all(unchanged),
    cohort_checks=as.list(unchanged),
    covariate_selection_approved=FALSE,
    limitation="Author-approved simple imputation; no claim that imputed values are observed or that imputation uncertainty is modeled")
  list(cohorts=cohorts, preparation_warnings=unique(preparation_warnings), audit=audit)
}
