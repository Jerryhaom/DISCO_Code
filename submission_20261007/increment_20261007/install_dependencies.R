#!/usr/bin/env Rscript
# Official CRAN source installation in a user-chosen LOCAL library only.
args <- commandArgs(TRUE); stopifnot(length(args) >= 1L)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(dirname(file_arg))
lib <- args[1]
system_libs <- normalizePath(unique(c(.Library, .Library.site, .libPaths())), mustWork = FALSE)
candidate <- normalizePath(lib, mustWork = FALSE)
if (candidate %in% system_libs) stop("Choose a task-local library, not an existing global library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE); lib <- normalizePath(lib)
.libPaths(c(lib, .libPaths()))
lock <- read.csv(file.path(root, "package_versions.csv"), stringsAsFactors = FALSE)
cache <- if (length(args) >= 2L) args[2] else tempfile("cran_sources_")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
pins <- list(
  PermAlgo = list(url = "https://cran.r-project.org/src/contrib/PermAlgo_1.2.tar.gz",
    sha = "aa2c774d6c6dcfeec882c1936e8723ef49bd36030fb10c17ca60bb9d4a519443"),
  simsurv = list(url = "https://cran.r-project.org/src/contrib/Archive/simsurv/simsurv_1.0.0.tar.gz",
    sha = "799a37538080391dde194a2d73645933c1e278b2566dc09cf0f08090d7506299"),
  coxed = list(url = "https://cran.r-project.org/src/contrib/Archive/coxed/coxed_0.3.3.tar.gz",
    sha = "d0d6cb8fea9516b3c63b34d0d81f3804c18a07f97a83e51555575c8ed4c75626")
)
order <- c(setdiff(lock$package, c("PermAlgo", "simsurv", "coxed")), "PermAlgo", "simsurv", "coxed")
for (p in order) {
  v <- lock$version[match(p, lock$package)]
  installed <- if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING"
  if (installed == v) {cat("EXACT_VERSION_AVAILABLE", p, v, "\n"); next}
  tar <- file.path(cache, paste0(p, "_", v, ".tar.gz"))
  if (!file.exists(tar)) {
    urls <- if (!is.null(pins[[p]])) pins[[p]]$url else c(
      paste0("https://cran.r-project.org/src/contrib/", p, "_", v, ".tar.gz"),
      paste0("https://cran.r-project.org/src/contrib/Archive/", p, "/", p, "_", v, ".tar.gz"))
    ok <- FALSE
    for (u in urls) {
      ok <- tryCatch({download.file(u, tar, mode = "wb", method = "libcurl"); TRUE}, error = function(e) FALSE)
      if (ok) break
    }
    if (!ok) stop(paste("Official pinned source unavailable", p, v))
  }
  if (!is.null(pins[[p]])) {
    if (!requireNamespace("digest", quietly = TRUE)) stop("Install locked digest before pinned generator dependencies")
    if (digest::digest(file = tar, algo = "sha256") != pins[[p]]$sha) stop(paste("Source SHA256 mismatch", p))
  }
  install.packages(tar, repos = NULL, type = "source", lib = lib)
  actual <- read.dcf(file.path(lib, p, "DESCRIPTION"))[1, "Version"]
  if (actual != v) stop(paste("Installed version does not match lock", p))
}
cat("Local dependency installation complete:", lib, "\n")
cat("Set DISCO_INCREMENT_RLIB to that directory before run.R\n")
