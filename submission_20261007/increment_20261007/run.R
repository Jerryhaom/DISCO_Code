#!/usr/bin/env Rscript
# New standalone entrypoint. Historical repository outputs are not required.
args <- commandArgs(TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(dirname(file_arg))
source(file.path(root, "bootstrap.R")); initialize_increment(root)
mode <- if (length(args)) args[1] else "smoke"
out <- new_output(if (length(args) >= 2L) args[2] else file.path(getwd(),
  sprintf("DISCO_%s_%s", mode, format(Sys.time(), "%Y%m%dT%H%M%S"))))
stopifnot(mode %in% c("ed7a", "smoke", "full", "generate_full", "reference_figures"))
if (mode == "ed7a") {
  source(file.path(root, "plot_ED7a.R"), local = FALSE)
  plot_ED7a(root, out)
} else {
  source(file.path(root, "generate_simulation.R"))
  source(file.path(root, "fit_simulation.R"))
  source(file.path(root, "plot_simulation.R"))
  if (mode == "reference_figures") {
    reference <- fread(file.path(root, "input/reference_models_wide_96.csv"))
    plot_simulation(root, out, reference, reference_only = TRUE)
  } else {
    manifest <- generate_simulation(root, out, smoke = mode == "smoke")
    if (mode != "generate_full") {
      fitted <- fit_simulation(root, out, manifest)
      plot_simulation(root, out, fitted, smoke = mode == "smoke")
    }
  }
}
sources <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
input_files <- switch(mode,
  ed7a = file.path(root, c("input/ED7a_event_time_curves.csv", "input/ED7a_C_annotations.csv")),
  reference_figures = file.path(root, "input/reference_models_wide_96.csv"),
  generate_full = character(),
  c(list.files(file.path(out, "generated_data"), "_squared_data.csv$", full.names = TRUE),
    file.path(root, "input/reference_models_384.csv")))
stopifnot(all(startsWith(normalizePath(input_files), paste0(root, "/")) |
  startsWith(normalizePath(input_files), paste0(out, "/"))))
write_json(list(mode = mode, output = out, R = R.version.string,
  source_root = root, working_directory = getwd(), package_library_paths = .libPaths(),
  scientific_input_sha256 = as.list(setNames(vapply(input_files,
    function(p) digest(file = p, algo = "sha256"), character(1)), input_files)),
  scientific_inputs_only_inside_extract_or_new_output = TRUE,
  source_sha256 = as.list(setNames(vapply(sources, function(p) digest(file = p, algo = "sha256"), character(1)),
    substring(sources, nchar(root) + 2L))),
  participants_bundled = FALSE, restricted_root_environment_required = FALSE,
  scope = if (mode == "smoke") "Four new synthetic datasets and16 new models; not all96/full S6-S7" else mode),
  file.path(out, "runtime_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("COMPLETE", mode, "OUTPUT", out, "\n")
