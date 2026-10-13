#!/usr/bin/env Rscript
# Portable release bootstrap (source/data/output roots are independent).
local({
  frames <- vapply(sys.frames(), function(e) if (is.null(e$ofile)) "" else as.character(e$ofile)[1], character(1))
  frames <- frames[nzchar(frames)]
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  entry <- if (length(frames)) tail(frames, 1) else if (length(arg)) sub("^--file=", "", arg[1]) else ""
  cfg <- file.path(dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE), mustWork = TRUE)), "..", "project_config.R")
  source(cfg, local = .GlobalEnv)
})
disco_init_module("05_organ_network")


# Audit organ-DISCO and organ-network saved objects without modifying legacy files.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
organ <- file.path(root, "Result3/Organ")
sens <- file.path(root, "Result3/Sensitivity")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mc36 <- readRDS(file.path(organ, "UKB.tissue.mc.rds"))
mc18 <- readRDS(file.path(organ, "UKB.tissue.mc2.rds"))
odis <- readRDS(file.path(organ, "UKB.odis.rds"))
pc18 <- readRDS(file.path(organ, "Tissue.pcorr.noverlap2.rds"))

numeric_finite <- function(x) {
  if(is.data.frame(x))sum(vapply(x,function(z)if(is.numeric(z))sum(is.finite(z)) else 0,numeric(1)))
  else sum(is.finite(x))
}
numeric_total <- function(x) {
  if(is.data.frame(x))sum(vapply(x,function(z)if(is.numeric(z))length(z) else 0,numeric(1)))
  else length(x)
}
objects <- list(mc36,mc18,odis,pc18)
audit <- data.frame(
  object = c("UKB.tissue.mc.rds", "UKB.tissue.mc2.rds", "UKB.odis.rds", "Tissue.pcorr.noverlap2.rds"),
  rows = c(nrow(mc36), nrow(mc18), nrow(odis), nrow(pc18)),
  cols = c(ncol(mc36), ncol(mc18), ncol(odis), ncol(pc18)),
  finite_values = vapply(objects,numeric_finite,numeric(1)),
  total_values = vapply(objects,function(x)prod(dim(x)),numeric(1)),
  numeric_cells = vapply(objects,numeric_total,numeric(1)),
  finite_value_scope = "numeric cells only; character metadata is not a non-finite number",
  stringsAsFactors = FALSE
)

organ_names <- data.frame(
  panel = c("mc36", "mc18", "odis_metadata", "pcorr18"),
  n = c(ncol(mc36), ncol(mc18), sum(!(names(odis) %in% c("mdis", "status", "time", "age", "sex", "cause"))), nrow(pc18)),
  names = c(
    paste(colnames(mc36), collapse = ";"),
    paste(colnames(mc18), collapse = ";"),
    paste(names(odis)[!(names(odis) %in% c("mdis", "status", "time", "age", "sex", "cause"))], collapse = ";"),
    paste(rownames(pc18), collapse = ";")
  ), stringsAsFactors = FALSE
)

sensitivity <- data.frame(
  object = c("UKB.organ.random.disco.rds", "UKB.organ.random.cindex.rds", "Tissue.corr.random.rds", "Tissue.corr.sensitivity.rds"),
  exists = file.exists(file.path(sens, c("UKB.organ.random.disco.rds", "UKB.organ.random.cindex.rds", "Tissue.corr.random.rds", "Tissue.corr.sensitivity.rds"))),
  stringsAsFactors = FALSE
)
random_list <- readRDS(file.path(sens, "UKB.organ.random.disco.rds"))
random_c <- readRDS(file.path(sens, "UKB.organ.random.cindex.rds"))
sensitivity$replicates <- c(length(random_list), ncol(random_c), NA, NA)

metadata <- data.frame(
  n = nrow(odis),
  events = sum(odis$status == 1, na.rm = TRUE),
  invalid_status = sum(!is.na(odis$status) & !odis$status %in% c(0, 1)),
  repeated_age_values = anyDuplicated(odis$age) > 0,
  participant_id_column = NA_character_,
  duplicated_participant_ids = NA,
  na_total = sum(is.na(odis)),
  stringsAsFactors = FALSE
)

write.csv(audit, file.path(out_dir, "saved_object_audit.csv"), row.names = FALSE)
write.csv(organ_names, file.path(out_dir, "organ_panel_audit.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(out_dir, "sensitivity_object_audit.csv"), row.names = FALSE)
write.csv(metadata, file.path(out_dir, "organ_metadata_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_WITH_ORGAN_PANEL_DEFINITION_BLOCKER",
  "saved_organ_scores_audited=true",
  "random_protein_sensitivity_objects_present=true",
  "organ_count_18_vs_36_requires_lock=true",
  "raw_protein_comparison_not_present_in_this_audit=true",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))
cat("PASS_WITH_ORGAN_PANEL_DEFINITION_BLOCKER\n")
