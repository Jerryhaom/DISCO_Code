#!/usr/bin/env Rscript
# Inventory source schemas and candidate date fields; no participant values exported.
Sys.umask('0077')
suppressPackageStartupMessages({library(readxl);library(jsonlite);library(digest)})
out <- Sys.getenv('DISCO_RUN_MODULE')
for(p in c('input','results','figures','tables','logs','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
Sys.chmod(file.path(out,'private'),'0700')
all <- system2(Sys.which('rg'),c('--files',Sys.getenv('DISCO_INPUT_ROOT')),stdout=TRUE)
paths <- all[grepl('/CHARLS/|/Frailty/',all,ignore.case=TRUE)&grepl('\\.(xlsx|csv|txt|rds)$',all,ignore.case=TRUE)&!grepl('/Results_final/',all,fixed=TRUE)]
hashes <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
write.csv(data.frame(path=paths,sha256=hashes),file.path(out,'input/all_CHARLS_local_source_inventory.csv'),row.names=FALSE)
unique_paths <- paths[!duplicated(hashes)]
rows <- list();schema <- list()
for(p in unique_paths){
 extension <- tolower(tools::file_ext(p))
 sheets <- if(extension=='xlsx') excel_sheets(p) else 'single'
 for(s in sheets){
  d <- tryCatch(if(extension=='xlsx')as.data.frame(read_xlsx(p,sheet=s)) else if(extension=='rds')readRDS(p) else if(extension=='csv')read.csv(p,check.names=FALSE,stringsAsFactors=FALSE) else read.delim(p,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
  if(!is.data.frame(d))next
  names <- names(d)
  candidate <- names[grepl('date|year|month|day|iyear|imonth|died|exit|exb001|exb002|blood|time|interview|alive|death',names,ignore.case=TRUE)]
  schema[[paste(p,s,sep='::')]]<-list(N=nrow(d),P=ncol(d),columns=names,candidate_date_status_fields=candidate)
  for(v in candidate){
   x<-d[[v]]; vals <- unique(as.character(x[!is.na(x)]))
   # Only low-cardinality status and date components, not identifiers or health values.
   public_date_component <- grepl('^(iyear|imonth|iday|exb001_|died|interviewyear|interviewmonth)',v,ignore.case=TRUE)
   rows[[length(rows)+1L]]<-data.frame(path=p,sheet=s,N=nrow(d),variable=v,nonmissing=sum(!is.na(x)),unique_n=length(vals),
     distinct_values=if(public_date_component&&length(vals)<=70)paste(sort(vals),collapse=' | ') else 'not exported',stringsAsFactors=FALSE)
  }
 }
}
write_json(schema,file.path(out,'input/local_unique_schema_inventory.json'),auto_unbox=TRUE,pretty=TRUE)
write.csv(if(length(rows))do.call(rbind,rows) else data.frame(),file.path(out,'results/candidate_date_status_fields.csv'),row.names=FALSE)
cat('Paths',length(paths),'unique hashes',length(unique_paths),'schema tables',length(schema),'candidate fields',length(rows),'\n')
