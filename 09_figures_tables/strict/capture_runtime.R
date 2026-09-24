# Loaded only by an explicitly isolated strict replay invocation.
# Scientific calculations and source plotting expressions are not modified.
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
.strict_inputs <- list()
.strict_active <- FALSE
strict_record_input <- function(path,kind) {
  if (!.strict_active) return(invisible(NULL))
  if (!is.character(path) || length(path)!=1 || is.na(path) || !file.exists(path)) return(invisible(NULL))
  resolved <- normalizePath(path,mustWork=TRUE)
  roots <- c(disco_project_root(), disco_results_root(), .disco_code_root)
  if (!any(vapply(roots, function(p) startsWith(resolved, paste0(p, "/")), logical(1))) || dir.exists(resolved)) return(invisible(NULL))
  key <- paste(kind,resolved,sep="|")
  if (!is.null(.strict_inputs[[key]])) return(invisible(NULL))
  fi <- file.info(resolved)
  .strict_inputs[[key]] <<- list(path=resolved,kind=kind,bytes=unname(fi$size),
    mtime=as.character(fi$mtime),md5_at_read=unname(tools::md5sum(resolved)),cwd=getwd())
  invisible(NULL)
}
strict_trace_reader <- function(name,namespace,arg) {
  ns <- asNamespace(namespace)
  if (!exists(name,envir=ns,inherits=FALSE)) return(invisible(NULL))
  tracer <- substitute(strict_record_input(ARG,KIND),list(ARG=as.name(arg),KIND=paste0(namespace,"::",name)))
  trace(name,tracer=tracer,where=ns,print=FALSE)
}
strict_trace_reader("readRDS","base","file")
strict_trace_reader("load","base","file")
strict_trace_reader("read.table","utils","file")
for (fun in c("read_csv","read_tsv","read_delim")) {
  if (requireNamespace("readr",quietly=TRUE)) strict_trace_reader(fun,"readr","file")
}
for (fun in c("read_xlsx","read_excel","read_xls")) {
  if (requireNamespace("readxl",quietly=TRUE)) strict_trace_reader(fun,"readxl","path")
}
if (requireNamespace("BioAge",quietly=TRUE)) {
  description <- packageDescription("BioAge")
  if (grepl("stub|replacement",paste(description$Title,description$Description),ignore.case=TRUE)) {
    trace("hd_calc",where=asNamespace("BioAge"),print=FALSE,
          tracer=quote(stop("STRICT GUARD: local BioAge stub cannot stand in for historical hd_calc")))
  }
}
strict_finish <- function(tag,status) {
  dir.create(file.path(out,"runtime"),showWarnings=FALSE)
  blocks <- Filter(function(x) identical(x$tag,tag),provenance)
  # Current installed versions are an observed environment, NOT the historical one.
  pkgs <- lapply(sort(loadedNamespaces()),function(n) {
    d <- packageDescription(n)
    list(package=n,version=as.character(packageVersion(n)),path=find.package(n),
         description_md5=unname(tools::md5sum(file.path(find.package(n),"DESCRIPTION"))),
         built=d$Built %||% NA_character_,remote_sha=d$RemoteSha %||% NA_character_)
  })
  record <- list(tag=tag,status=status,R=R.version.string,platform=R.version$platform,
    command=commandArgs(),cwd=getwd(),libPaths=.libPaths(),locale=Sys.getlocale(),
    timezone=Sys.timezone(),RNGkind=RNGkind(),rng_present=exists(".Random.seed",envir=.GlobalEnv),
    packages=pkgs,inputs=unname(.strict_inputs),source_blocks=blocks,
    warnings="Observed current runtime only. No historical-environment equivalence claimed.")
  jsonlite::write_json(record,file.path(out,"runtime",paste0(tag,".json")),auto_unbox=TRUE,pretty=TRUE,na="null")
  writeLines(capture.output(sessionInfo()),file.path(out,"runtime",paste0(tag,"_sessionInfo.txt")))
}
`%||%` <- function(x,y) if(is.null(x)) y else x
