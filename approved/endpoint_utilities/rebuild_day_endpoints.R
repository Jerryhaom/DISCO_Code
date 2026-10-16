# Day-level UKB endpoint reconstruction. Explicit policy choices are mandatory.
# Generic vectors only: this file contains no participant data or identifiers.

endpoint_date <- function(x, label) {
  if (inherits(x, "Date")) return(x)
  if (!is.character(x)) stop(label, " must be Date or ISO YYYY-MM-DD text")
  blank <- is.na(x) | x == ""
  if (any(!blank & !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)))
    stop(label, " contains non-ISO or special-coded dates; resolve them first")
  y <- as.Date(x)
  if (any(!blank & (is.na(y) | format(y, "%Y-%m-%d") != x)))
    stop(label, " contains an invalid date")
  y
}

earliest_recorded_date <- function(dates) {
  if (!is.data.frame(dates) && !is.matrix(dates)) stop("dates must be a date-column table")
  if (!ncol(dates)) stop("At least one diagnosis-date instance is required")
  z <- vapply(seq_len(ncol(dates)), function(j)
    as.numeric(endpoint_date(dates[, j], paste0("diagnosis instance ", j))),
    numeric(nrow(dates)))
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  y <- apply(z, 1, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE))
  as.Date(y, origin = "1970-01-01")
}

build_day_endpoints <- function(id, index_date, death_date, loss_date,
                                diagnosis_date, admin_date, loss_rule,
                                baseline_same_day, disease_death_tie) {
  # No administrative cutoff, loss policy, baseline or tie rule is defaulted.
  if (anyDuplicated(id) || anyNA(id)) stop("Unique, nonmissing participant IDs required")
  n <- length(id)
  if (!n || any(vapply(list(index_date, death_date, loss_date, diagnosis_date), length, 0L) != n))
    stop("Date vectors must have the same nonzero length as id")
  if (!loss_rule %in% c("censor", "ignore")) stop("Explicit loss_rule: censor or ignore")
  if (!baseline_same_day %in% c("prevalent", "incident"))
    stop("Explicit baseline_same_day: prevalent or incident")
  if (!disease_death_tie %in% c("exclude", "disease_first", "death_first"))
    stop("Explicit disease_death_tie: exclude, disease_first or death_first")
  baseline <- endpoint_date(index_date, "index_date")
  death <- endpoint_date(death_date, "death_date")
  loss <- endpoint_date(loss_date, "loss_date")
  disease <- endpoint_date(diagnosis_date, "diagnosis_date")
  admin <- endpoint_date(admin_date, "admin_date")
  if (length(admin) == 1L) admin <- rep(admin, n)
  if (length(admin) != n || anyNA(admin)) stop("Explicit, nonmissing administrative cutoff required")
  if (anyNA(baseline)) stop("Unverified/missing baseline dates cannot be imputed here")
  if (any(admin < baseline)) stop("Administrative cutoff precedes baseline")
  end <- pmin(admin, death, if (loss_rule == "censor") loss else rep(as.Date(NA), n), na.rm = TRUE)
  eligible_mortality <- end > baseline
  mortality_event <- !is.na(death) & death == end & eligible_mortality
  before <- !is.na(disease) & (disease < baseline |
    (baseline_same_day == "prevalent" & disease == baseline))
  tie <- !is.na(disease) & !is.na(death) & disease == death & disease <= end
  eligible_incidence <- eligible_mortality & !before & !(disease_death_tie == "exclude" & tie)
  after_baseline <- !is.na(disease) & (disease > baseline |
    (baseline_same_day == "incident" & disease == baseline))
  event <- eligible_incidence & after_baseline & disease <= end &
    !(disease_death_tie == "death_first" & tie)
  # For cause-specific Cox, competing death is censored at its actual date.
  # For CIF, the same record is explicitly coded 2 rather than a disease event.
  incidence_end <- end
  incidence_end[event] <- disease[event]
  cif_event <- ifelse(event, 1L, ifelse(eligible_incidence & mortality_event, 2L, 0L))
  days_mortality <- as.numeric(end - baseline)
  days_incidence <- as.numeric(incidence_end - baseline)
  data.frame(id = id, index_date = baseline, administrative_cutoff = admin,
    observation_end = end, mortality_days = days_mortality,
    mortality_years = days_mortality / 365.25, mortality_event = as.integer(mortality_event),
    eligible_mortality = eligible_mortality, prevalent_disease = before,
    eligible_incidence = eligible_incidence, incidence_days = days_incidence,
    incidence_years = days_incidence / 365.25, incidence_event = as.integer(event),
    cif_event = cif_event, diagnosis_after_end = !is.na(disease) & disease > end,
    same_day_disease_death = tie, death_before_or_on_baseline = !is.na(death) & death <= baseline,
    loss_before_or_on_baseline = !is.na(loss) & loss <= baseline,
    stringsAsFactors = FALSE)
}
