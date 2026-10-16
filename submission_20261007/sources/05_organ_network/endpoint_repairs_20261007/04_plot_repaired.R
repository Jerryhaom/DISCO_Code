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
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table);library(jsonlite);library(grid);library(digest)})
setDTthreads(1)
entry<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]));out<-dirname(entry)
orig<-file.path(out,'../organ_clock_comparison');fig<-file.path(out,'panels');dir.create(fig,showWarnings=FALSE)
map<-fread(file.path(orig,'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
organs<-map$clock_organ;labs<-setNames(map$display,map$clock_organ)
C<-fread(file.path(out,'output/clock/test_train_adjusted_correlations.csv'))
C[,`:=`(main_figure=TRUE,adjustment='linear_age',partial_r=r)]
H<-fread(file.path(out,'output/clock/training_HR.csv'))
stopifnot(nrow(H)==72,!anyDuplicated(H[,.(organ,method,model,term)]),all(H$scale=='age_residual_SD'))
H[,term:=ifelse(term=='Da','D','P')];H[,`:=`(main_figure=TRUE,adjustment='linear_age',lower=lower95,upper=upper95)]
# Reuse only the actual original plotting definitions, not script filename guesses.
src<-parse(file.path(orig,'20260928/03_plot_panels.R'))
for(e in src)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('base','make_panels'))eval(e)
base<-base+theme(plot.title=element_text(size=9,face='plain'))
p<-make_panels('linear_age')
p$a$layers[[1]]<-geom_tile(width=1,height=1,color='white',linewidth=.22)
p$a<-p$a+labs(tag=NULL)+theme(axis.text.x=element_text(size=9.5,angle=30,hjust=1,color='black'),axis.text.y=element_text(size=9.5,color='black'))
for(n in c('b','c'))p[[n]]<-p[[n]]+labs(tag=NULL,x='HR per training age-residual SD (95% CI)')
savepdf<-function(g,name,w,h){
  path<-file.path(fig,name);stopifnot(!file.exists(path));pdf(path,width=w,height=h,family='Helvetica',useDingbats=FALSE)
  print(g);invisible(dev.off());stopifnot(file.info(path)$size>0)
}
savepdf(p$a,'ED10a_correlations_new_endpoint.pdf',9.6,2.2)
savepdf(p$b,'ED10b_DISCO_before_after_new_endpoint.pdf',4.8,4.1)
savepdf(p$c,'ED10c_ProtAge_before_after_new_endpoint.pdf',4.8,4.1)
for(n in c('b','c')){
  bb<-ggplot_build(p[[n]]);ss<-p[[n]]$data;stopifnot(nrow(ss)==36,
    max(abs(bb$data[[2]]$x-log10(ss$lower)))<1e-12,max(abs(bb$data[[2]]$xend-log10(ss$upper)))<1e-12,
    max(abs(bb$data[[3]]$x-log10(ss$HR)))<1e-12)
}
aa<-ggplot_build(p$a);stopifnot(nrow(p$a$data)==18,identical(aa$data[[2]]$label,sprintf('%.2f',p$a$data$partial_r)))
# Original paired Original / Accel nine-organ 3x3 design.
font<-theme_classic(base_family='Helvetica',base_size=7)+theme(axis.title=element_text(size=8),axis.text=element_text(size=7,color='black'),
  plot.title=element_text(size=8,hjust=.5,face='plain'),plot.margin=margin(2,8,2,8),axis.line=element_line(linewidth=.35),
  axis.ticks=element_line(linewidth=.3),axis.ticks.length=unit(2,'pt'))
ranges<-builtdata<-list()
d<-fread(file.path(out,'output/clock/panel_d_Cindex.csv'),data.table=FALSE);stopifnot(nrow(d)==54)
d$score<-factor(d$score,levels=c('Goeminne','Oh','DISCO'));d$display<-labs[d$organ]
left<-d[d$group=='Original',];original_limits<-c(min(.58,floor((min(left$lower95)-.008)*100)/100),max(.82,ceiling((max(left$upper95)+.008)*100)/100))
original_breaks<-c(.60,.80);original_breaks<-original_breaks[original_breaks>=original_limits[1]&original_breaks<=original_limits[2]]
source(file.path(orig,'20261006_strict_70_30/logs/strict_adapted_make_one.R'))
original_make_one<-make_one
make_one<-function(organ,grp){
  q<-original_make_one(organ,grp)
  if(grp=='Accel')q<-q+scale_y_discrete(labels=c(DISCO='DISCO Accel',Oh='Oh gap',Goeminne='Goeminne gap'),
    expand=expansion(add=.48))+theme(axis.text.y=element_text(size=6.5,color='black'))
  if(grp=='Accel'){
    tick<-ranges[[paste(organ,grp)]]$limits
    q<-q+scale_x_continuous(breaks=tick,labels=function(z)sprintf('%.2f',z),expand=expansion(mult=0))
    ranges[[paste(organ,grp)]]$breaks<<-tick
  }
  q
}
cells<-lapply(organs,function(o){q<-make_one(o,'Original')+make_one(o,'Accel')+plot_layout(widths=c(1,1))+
  plot_annotation(title=labs[o],theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=.5),plot.margin=margin(6,7,3,7)));wrap_elements(full=q)})
p$d<-wrap_plots(cells,ncol=3,nrow=3)+plot_annotation(title='Mortality discrimination',caption='C-index (95% CI)',
  theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=0,margin=margin(l=54)),plot.caption=element_text(family='Helvetica',size=9,hjust=.5),plot.margin=margin(6,7,6,7)))
savepdf(p$d,'ED10d_Original_Accel_3x3_new_endpoint.pdf',9.6,4.6)
full<-((p$a+labs(tag='a'))/((p$b+labs(tag='b'))|(p$c+labs(tag='c')))/(wrap_elements(full=p$d)+labs(tag='d')))+
  plot_layout(heights=c(1,2,2.25))
full<-full&theme(plot.tag=element_text(size=14,face='bold'),plot.tag.position=c(0,.99))
savepdf(full,'Organ_clock_comparison_NEW_ENDPOINT_REPAIRED.pdf',9.6,11.1)
original_ranges<-ranges
font<-font+theme(axis.text=element_text(size=8,color='black'),plot.title=element_text(size=9))
d<-fread(file.path(out,'output/clock/panel_d_matched_combined_Cindex.csv'),data.table=FALSE)
stopifnot(nrow(d)==90,all(is.finite(d$C)),all(d$available),!anyDuplicated(d[,c('organ','group','score')]))
d$score<-factor(d$score,levels=c('Goeminne','Oh','DISCO_Goeminne','DISCO_Oh','DISCO'));d$display<-labs[d$organ]
left<-d[d$group=='Original',];original_limits<-c(min(.54,floor((min(left$lower95)-.008)*100)/100),max(.82,ceiling((max(left$upper95)+.008)*100)/100))
original_breaks<-c(.55,.80);original_breaks<-original_breaks[original_breaks>=original_limits[1]&original_breaks<=original_limits[2]]
ranges<-builtdata<-list();source(file.path(orig,'20261006_strict_70_30/logs/strict_matched_adapted_make_one.R'))
original_matched_make_one<-make_one
make_one<-function(organ,grp){
  q<-original_matched_make_one(organ,grp)
  if(grp=='Accel')q<-q+scale_y_discrete(drop=FALSE,labels=c(DISCO='DISCO Accel',DISCO_Oh='DISCO-Oh Accel',
    DISCO_Goeminne='DISCO-Go Accel',Oh='Oh gap',Goeminne='Goeminne gap'),expand=expansion(add=.48))+
      theme(axis.text.y=element_text(size=7,color='black'))
  if(grp=='Accel'){
    tick<-ranges[[paste(organ,grp)]]$limits
    q<-q+scale_x_continuous(breaks=tick,labels=function(z)sprintf('%.2f',z),expand=expansion(mult=0))
    ranges[[paste(organ,grp)]]$breaks<<-tick
  }
  q
}
cells<-lapply(organs,function(o){q<-make_one(o,'Original')+make_one(o,'Accel')+plot_layout(widths=c(1,1))+
  plot_annotation(title=labs[o],theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=.5),plot.margin=margin(6,7,3,7)));wrap_elements(full=q)})
matched_plot<-wrap_plots(cells,ncol=3,nrow=3)+plot_annotation(caption='C-index (95% CI)',theme=theme(plot.caption=element_text(family='Helvetica',size=10,hjust=.5),plot.margin=margin(6,6,6,6)))
savepdf(matched_plot,'S9_matched_proteins_90cells_new_endpoint.pdf',11.8,6.8)
write_json(list(status='VALUE_CHECKS_PASSED',HR_b_points=36,HR_c_points=36,heatmap_cells=18,d_cells=54,S9_cells=90,
  actual_abc_source='20260928/03_plot_panels.R make_panels',actual_d_source='20261006_strict_70_30/logs/strict_adapted_make_one.R',
  actual_S9_source='20261006_strict_70_30/logs/strict_matched_adapted_make_one.R',HR_scale='training age-residual SD; previous primary specification retained',
  before_after_models='age+sex; then mutually adjusted for same-organ alternative score',d_ranges=original_ranges,S9_ranges=ranges,
  caveats=c('HR is training association summary; PH failures remain','d Original has unequal age information; infer comparative age-independent prediction from Fair models',
    'CI excludes repeated split, retraining and imputation uncertainty','new endpoint full scores evaluated, not merely assumed unchanged')),
  file.path(out,'qa/plot_value_checks.json'),pretty=TRUE,auto_unbox=TRUE)
fwrite(H,file.path(out,'qa/HR_plot_point_table.csv'))
cat('DONE repaired clock figure panels: a18 b36 c36 d54 S9 90\n')
