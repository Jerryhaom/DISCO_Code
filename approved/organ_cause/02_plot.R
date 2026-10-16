#!/usr/bin/env Rscript
# Fig4a: original Plot_Fig4E.R aes/facet/theme, corrected mortality only.
options(stringsAsFactors=FALSE)
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggthemes);library(reshape2);library(jsonlite);library(digest)})
entry<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)));base<-dirname(entry)
root<-Sys.getenv("DISCO_DATA_ROOT")
args<-commandArgs(TRUE);pass<-if(length(args))args[1] else "pass_1";stopifnot(pass%in%c("pass_1","pass_2"))
out<-file.path(base,"figures",pass);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
r<-fread(file.path(base,pass,"corrected_model_results.csv"))
order<-readRDS(file.path(root,"Result3/Organ/Tissue.order.rds"))
original<-file.path(root,"Result3/Organ/Plot_Fig4E.R")
source_lines<-readLines(original,warn=FALSE)
stopifnot(grepl("ggplot(rr,aes(Var1,Var2))",source_lines[343],fixed=TRUE),
  grepl("facet_wrap(tt~., scales = \"free_y\",nrow=2)",source_lines[349],fixed=TRUE))
writeLines(source_lines[341:361],file.path(out,"original_Plot_Fig4E_lines341_361.txt"))
stars<-function(p)ifelse(is.na(p),"",ifelse(p<.01,"**",ifelse(p<.05,"*","")))
disease_hr<-readRDS(file.path(root,"Result3/Organ/mortality/res.rds"))
disease_stars<-readRDS(file.path(root,"Result3/Organ/mortality/pl.rds"))
di<-melt(disease_hr);names(di)<-c("Var1","Var2","value");di$p<-as.character(melt(disease_stars)$value)
di$tt<-"Chronic disease incidence"
expfile<-file.path(root,"Results_final/05_organ_network/outputs/organ_expanded_covariates/organ_chronic_disease_expanded.csv")
ex<-fread(expfile);dex<-data.frame(Var1=ex$organ,Var2=ex$disease,value=ex$Model2_HR,p=stars(ex$Model2_P),tt="Chronic disease incidence")
plot_input_paths<-c(original,file.path(base,pass,"corrected_model_results.csv"),
  file.path(root,"Result3/Organ/Tissue.order.rds"),file.path(root,"Result3/Organ/mortality/res.rds"),
  file.path(root,"Result3/Organ/mortality/pl.rds"),expfile,entry)
plot_input_hashes<-vapply(plot_input_paths,function(x)digest(file=x,algo="sha256"),character(1))
mort_order<-c("All-cause","Heart disease","Respiratory disease","Cancer","Cerebrovascular disease","Alzheimer's disease","Others")
end_order<-rev(unique(c(mort_order,colnames(disease_hr))))
device_log<-list();outputs<-list()
draw<-function(d,path,unit,combined){
  d$Var1<-factor(d$Var1,levels=order);d$Var2<-factor(d$Var2,levels=end_order)
  d$tt<-factor(d$tt,levels=c("Mortality","Chronic disease incidence"));d$p[d$p=="*"]<-" *"
  stopifnot(!anyNA(d$Var1),!anyNA(d$Var2))
  # Preserve original panel coordinates, number positioning, two-row wrapping,
  # organ order, colors, font sizes, and 12x6-inch combined canvas.
  p<-ggplot(d,aes(Var1,Var2))+
    geom_tile(aes(fill=value),color="white")+
    geom_text(aes(label=ifelse(is.na(value),"NA",sprintf("%0.2f",value))),size=3.5,hjust=.7)+
    geom_text(aes(label=p),hjust=-.8,vjust=.7)+
    facet_wrap(tt~.,scales="free_y",nrow=2)+
    scale_fill_distiller(palette="Reds",direction=1,na.value="grey90")+
    theme_few()+labs(fill=if(unit=="per_SD")"HR per SD" else "HR per unit")+
    xlab("Organ and tissues")+ylab("")+
    theme(panel.spacing.x=grid::unit(.5,"lines"),
      axis.text.x=element_text(size=14,color="black",angle=30,hjust=1),
      axis.text.y=element_text(size=12,color="black"),strip.text=element_text(size=16,color="black"),
      axis.title=element_text(size=14))
  h<-if(combined)6 else 3.65
  device<-"cairo_pdf";device_warnings<-character()
  before_device<-dev.cur()
  ok<-tryCatch(withCallingHandlers({grDevices::cairo_pdf(path,height=h,width=12);dev.cur()!=before_device},
    warning=function(w){device_warnings<<-c(device_warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)FALSE)
  if(!ok){device<-"pdf_Helvetica_fallback";grDevices::pdf(path,height=h,width=12,family="Helvetica",useDingbats=FALSE)}
  print(p);dev.off()
  stopifnot(file.exists(path),file.info(path)$size>0)
  device_log[[length(device_log)+1]]<<-list(file=basename(path),device=device,warnings=device_warnings,width=12,height=h)
  fwrite(d,file.path(out,paste0(sub(".pdf$","",basename(path)),"_plot_data.csv")),na="NA")
  outputs[[length(outputs)+1]]<<-data.frame(file=path,sha256=digest(file=path,algo="sha256"),rows=nrow(d),NA_cells=sum(is.na(d$value)),unit=unit,
    mortality_recalculated=TRUE,disease_recalculated=FALSE)
}
for(model_branch in c("age_sex","expanded_legacy")){
  d<-r[branch==model_branch]
  stopifnot(nrow(d)==126)
  for(unit in c("per_unit","per_SD")){
    hh<-paste0("HR_",unit)
    md<-data.frame(Var1=gsub("[.]"," ",d$organ),Var2=d$outcome,value=d[[hh]],p=stars(d$P),tt="Mortality")
    draw(md,file.path(out,paste0("Fig4a_mortality_",model_branch,"_",unit,".pdf")),unit,FALSE)
    if(unit=="per_unit"){
      dd<-if(model_branch=="age_sex")di else dex
      draw(rbind(md,dd),file.path(out,paste0("Fig4a_",model_branch,"_per_unit_corrected_mortality_PRESERVED_disease.pdf")),unit,TRUE)
    }
  }
}
fwrite(rbindlist(outputs),file.path(out,"figure_manifest.csv"))
stopifnot(identical(plot_input_hashes,vapply(plot_input_paths,function(x)digest(file=x,algo="sha256"),character(1))))
fwrite(data.table(file=plot_input_paths,sha256=plot_input_hashes),file.path(out,"plot_input_manifest.csv"))
write_json(list(original_plot=original,original_plot_sha256=digest(file=original,algo="sha256"),
  drawing_lines="341-361",numeric_only_mortality_change=TRUE,
  adaptations=c("corrected mortality HR/P","explicit per-unit/per-SD legend","NA warning cells","mortality-only 3.65-inch output"),
  original_aes="Var1,Var2",original_facets="facet_wrap(tt~.,scales=free_y,nrow=2)",
  original_combined_dimensions=c(width=12,height=6),font_sizes=list(x_text=14,y_text=12,strip=16,axis_title=14,value_mm=3.5),
  disease_input=list(main=file.path(root,"Result3/Organ/mortality/res.rds"),expanded=expfile,raw_diagnostic_dates_verified=FALSE),
  combined_units="per-unit for both rows; not yet a whole-figure per-SD release",devices=device_log),
  file.path(out,"plot_provenance.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Completed six PDFs",out,"\n")
