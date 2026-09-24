#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE), mustWork = TRUE))
if (!length(args)) stop("Usage: Rscript run.R <module/script.R> [script arguments]")
script <- normalizePath(file.path(repo, args[1]), mustWork = TRUE)
if (!startsWith(script, paste0(repo, "/")) || !grepl("[.]R$", script)) stop("Choose an R script inside this repository")
# Invoke a child so commandArgs(trailingOnly=TRUE) remains meaningful to replay scripts.
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(args[-1])))
quit(status = status)
