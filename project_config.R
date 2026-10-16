# One configuration for source code, restricted inputs, and generated results.
.disco_config_file <- local({
  files <- vapply(sys.frames(), function(e) if (is.null(e$ofile)) "" else as.character(e$ofile)[1], character(1))
  files <- files[nzchar(files)]
  if (!length(files)) stop("Source project_config.R, or use Rscript run.R <module/script.R>")
  normalizePath(tail(files, 1), mustWork = TRUE)
})
.disco_code_root <- dirname(.disco_config_file)
.disco_initial_wd <- getwd()
disco_code_path <- function(...) file.path(.disco_code_root, ...)
.disco_env_path <- function(name, fallback) {
  p <- path.expand(Sys.getenv(name, unset = fallback))
  if (!grepl("^(/|[A-Za-z]:|\\\\)", p)) p <- file.path(.disco_initial_wd, p)
  normalizePath(p, mustWork = FALSE)
}
.disco_data_root <- .disco_env_path("DISCO_DATA_ROOT", disco_code_path("data"))
disco_project_root <- function() .disco_data_root
.disco_output_root <- .disco_env_path("DISCO_OUTPUT_ROOT", disco_code_path("results"))
disco_results_root <- function() .disco_output_root
# Child Rscript processes inherit the resolved roots, even after setwd().
Sys.setenv(DISCO_DATA_ROOT=.disco_data_root, DISCO_OUTPUT_ROOT=.disco_output_root)
for (.key in c("DISCO_PRS_FILE", "DISCO_RULAS_WAVE2_FILE", "DISCO_NHANES4_FILE")) {
  if (nzchar(Sys.getenv(.key))) do.call(Sys.setenv, setNames(list(.disco_env_path(.key, "")), .key))
}
disco_path <- function(...) file.path(disco_project_root(), ...)
disco_results_path <- function(...) file.path(disco_results_root(), ...)
disco_external_path <- function(key, fallback) .disco_env_path(key, disco_path(fallback))
disco_init_module <- function(module) {
  path <- disco_results_path(module)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  setwd(path)
  invisible(path)
}

# Cox nobs() counts events, not analyzed participants. Model tables use fit$n.
# Empty or non-finite fitted results must not be reported as a successful audit.
disco_validate_cox_rows <- function(rows) {
  required <- c("n", "events", "c_index", "se")
  if (!is.data.frame(rows) || !nrow(rows) || !all(required %in% names(rows)))
    stop("No complete Cox result rows; audit cannot pass", call. = FALSE)
  vals <- as.matrix(rows[, required, drop = FALSE])
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(rows$n <= 0) ||
      any(rows$events <= 0) || any(rows$events > rows$n) ||
      any(rows$c_index < 0 | rows$c_index > 1) || any(rows$se < 0))
    stop("Invalid Cox result metadata or concordance; audit cannot pass", call. = FALSE)
  invisible(TRUE)
}
