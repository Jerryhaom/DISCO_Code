# Local package library only; no participant-data environment is needed.
initialize_increment <- function(root) {
  lib <- Sys.getenv("DISCO_INCREMENT_RLIB", unset = "")
  if (nzchar(lib)) {
    if (!dir.exists(lib)) stop("DISCO_INCREMENT_RLIB does not exist")
    .libPaths(c(normalizePath(lib), .libPaths()))
  }
  packages <- read.csv(file.path(root, "package_versions.csv"), stringsAsFactors = FALSE)
  actual <- vapply(packages$package, function(p) {
    if (!requireNamespace(p, quietly = TRUE)) return("MISSING")
    as.character(packageVersion(p))
  }, character(1))
  mismatch <- actual != packages$version
  if (any(mismatch)) stop(paste("Dependency version mismatch:",
    paste(packages$package[mismatch], actual[mismatch], "expected", packages$version[mismatch], collapse = "; ")))
  if (getRversion() != "4.5.3") warning("Tested R version is4.5.3; numerical replay must be checked on this R version")
  options(stringsAsFactors = FALSE, warn = 1, device = function(...) grDevices::pdf(file = NULL, ...))
  suppressPackageStartupMessages({
    library(data.table); library(jsonlite); library(digest); library(survival)
    library(coxed); library(simsurv); library(mvtnorm)
    library(ggplot2); library(ggpubr); library(ggthemes); library(ggsci); library(reshape2)
  })
  setDTthreads(1)
  invisible(packages)
}

new_output <- function(out) {
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
    stop("Choose an empty output directory; existing files are not overwritten")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out)
}
