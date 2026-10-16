#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
options(stringsAsFactors=FALSE)
.libPaths(c(.disco_path("results","09_figures_tables/Rlibs"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite);library(digest)})
if(requireNamespace('ggthemes',quietly=TRUE))theme_few<-ggthemes::theme_few else {
  # Exact local ggthemes::theme_few definition; server lacks that package.
  theme_few<-function(base_size=12,base_family=''){
    gray<-'#4D4D4D';black<-'#000000'
    theme_bw(base_size=base_size,base_family=base_family)+
      theme(line=element_line(colour=gray),rect=element_rect(fill='white',colour=NA),text=element_text(colour=black),
        axis.ticks=element_line(colour=gray),legend.key=element_rect(colour=NA),panel.border=element_rect(colour=gray),
        panel.grid=element_blank(),strip.background=element_rect(fill='white',colour=NA))
  }
}
entry<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]));out<-dirname(entry)
data_root<-if(dir.exists(file.path(out,'output/organ')))file.path(out,'output/organ')else file.path(out,'organ')
dir.create(file.path(out,'panels'),showWarnings=FALSE);dir.create(file.path(out,'qa'),showWarnings=FALSE)
r<-fread(file.path(data_root,'organ_mortality_models.csv'))[branch=='age_sex'&riskset!='legacy_delete']
di<-fread(file.path(data_root,'organ_disease_models.csv'))[branch=='age_sex'&policy=='death_first']
stopifnot(nrow(r)==126,nrow(di)==144,!anyDuplicated(r[,.(organ,outcome)]),!anyDuplicated(di[,.(organ,outcome)]),
  !any(r$needs_model_review),!any(di$needs_model_review))
orig<-Sys.getenv('DISCO_ORGAN_ORIGINAL_PLOT',.disco_path("legacy","Result3/Organ/Plot_Fig4E.R"));src<-readLines(orig,warn=FALSE)
stopifnot(grepl('ggplot(rr,aes(Var1,Var2))',src[343],fixed=TRUE),grepl('facet_wrap(tt~.',src[349],fixed=TRUE))
writeLines(src[341:361],file.path(out,'qa/Fig4a_actual_original_plot_block.txt'))
order<-readRDS(Sys.getenv('DISCO_ORGAN_ORDER',.disco_path("legacy","Result3/Organ/Tissue.order.rds")))
ends<-rev(unique(c('All-cause','Heart disease','Respiratory disease','Cancer','Cerebrovascular disease',"Alzheimer's disease",'Others',
  'Obesity','Hypertension','Diabetes','Cancer','MI','Stroke','COPD','Dementia')))
for(scale in c('unit','SD')){
  assemble<-function(x,tt)data.table(Var1=gsub('[.]',' ',x$organ),Var2=x$outcome,value=x[[paste0('HR_',scale)]],
    p=ifelse(x$P<.01,'**',ifelse(x$P<.05,' *','')),tt=tt)
  rr<-rbind(assemble(r,'Mortality'),assemble(di,'Chronic disease incidence'))
  rr$Var1<-factor(rr$Var1,levels=order);rr$Var2<-factor(rr$Var2,levels=ends)
  rr$tt<-factor(rr$tt,levels=c('Mortality','Chronic disease incidence'));stopifnot(!anyNA(rr$Var1),!anyNA(rr$Var2))
  p<-ggplot(rr,aes(Var1,Var2))+geom_tile(aes(fill=value),color='white')+
    geom_text(aes(label=sprintf('%0.2f',value)),size=3.5,hjust=.7)+geom_text(aes(label=p),hjust=-.8,vjust=.7)+
    facet_wrap(tt~.,scales='free_y',nrow=2)+scale_fill_distiller(palette='Reds',direction=1)+theme_few()+
    labs(fill=if(scale=='SD')'HR per SD'else'HR per unit')+xlab('Organ and tissues')+ylab('')+
    theme(panel.spacing.x=unit(.5,'lines'),axis.text.x=element_text(size=14,color='black',angle=30,hjust=1),
      axis.text.y=element_text(size=12,color='black'),strip.text=element_text(size=16,color='black'),axis.title=element_text(size=14))
  path<-file.path(out,'panels',paste0('Fig4a_age_sex_NEW_ENDPOINT_per_',scale,'.pdf'));stopifnot(!file.exists(path))
  # Match the original cairo_pdf device as well as the 12x6-inch geometry.
  withCallingHandlers(cairo_pdf(path,width=12,height=6,family='Helvetica'),
    warning=function(w)stop('Cairo device unavailable: ',conditionMessage(w)))
  stopifnot(file.exists(path));print(p);invisible(dev.off());stopifnot(file.info(path)$size>0)
  bb<-ggplot_build(p);stopifnot(nrow(bb$data[[1]])==270,identical(bb$data[[2]]$label,sprintf('%0.2f',rr$value)))
  fwrite(rr,file.path(out,'qa',paste0('Fig4a_per_',scale,'_plot_data.csv')))
}
write_json(list(status='VALUE_CHECKS_PASSED',cells=270,mortality_cells=126,disease_cells=144,original_source=orig,
  original_source_sha256=digest(file=orig,algo='sha256'),original_block='341-361',original_dimensions=c(width=12,height=6),
  PDF_device='cairo_pdf, as in confirmed original plotting block',
  primary_branch='age+sex',disease_primary='death-first; baseline prevalent excluded; competing death censored',
  SD='branch common complete population before disease-specific prevalence exclusions; fixed across outcomes',
  choices='per-unit preserves original scale; per-SD is separately named candidate, not silently substituted',
  layout='same aes(Var1,Var2), two facet_wrap rows and original font sizes'),file.path(out,'qa/Fig4a_plot_QA.json'),pretty=TRUE,auto_unbox=TRUE)
cat('DONE Fig4a both rows updated; 270 exact cells per scale\n')
