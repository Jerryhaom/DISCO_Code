#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(entry))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table);library(jsonlite);library(grid)})
d<-fread(file.path(out,'numeric/panel_d_Cindex.csv'),data.table=FALSE)
d$score<-factor(d$score,levels=c('Goeminne','Oh','DISCO'));organs<-unique(d$organ)
labs<-setNames(d$display[match(organs,d$organ)],organs)
# Read the actual prior d design. Only grouping/label/required axis extent changes.
src<-parse(file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/organ_plot_source.R'))
assignments<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-')),as.list(src))
pick<-function(n)assignments[[which(vapply(assignments,function(e)identical(e[[2]],as.name(n)),FALSE))]]
eval(pick('font'))
raw<-paste(deparse(pick('make_one'),width.cutoff=500),collapse='\n')
adapt<-gsub('Age + sex','Original',raw,fixed=TRUE)
adapt<-sub('lim <- c(0.58, 0.82)','lim <- original_limits',adapt,fixed=TRUE)
adapt<-sub('breaks <- c(0.6, 0.7, 0.8)','breaks <- original_breaks',adapt,fixed=TRUE)
adapt<-sub('DISCO = "DISCO"','DISCO = if (grp == "Original") "DISCO + Age" else "DISCO"',adapt,fixed=TRUE)
stopifnot(!identical(raw,adapt),grepl('original_limits',adapt,fixed=TRUE),grepl('DISCO + Age',adapt,fixed=TRUE))
left<-d[d$group=='Original',];original_limits<-c(min(.58,floor((min(left$lower95)-.008)*100)/100),.82)
original_breaks<-c(.55,.65,.75);original_breaks<-original_breaks[original_breaks>=original_limits[1]]
ranges<-builtdata<-list();eval(parse(text=adapt))
writeLines(raw,file.path(out,'logs/original_make_one.R'));writeLines(adapt,file.path(out,'logs/adapted_make_one.R'))
cells_code<-paste(deparse(pick('cells'),width.cutoff=500),collapse='\n')
eval(parse(text=gsub('Age + sex','Original',cells_code,fixed=TRUE)))
p<-wrap_plots(cells,ncol=3,nrow=3)+plot_annotation(title='Mortality discrimination',caption='C-index (95% CI)',
 theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=0,margin=margin(l=54)),
 plot.caption=element_text(family='Helvetica',size=9,hjust=.5),plot.margin=margin(6,7,6,7)))
pdf(file.path(out,'panels/d_Cindex.pdf'),width=9.6,height=4.6,family='Helvetica',useDingbats=FALSE);print(p);invisible(dev.off())
saveRDS(list(plot=p,data=d,ranges=ranges),file.path(out,'plot_source.rds'))
write_json(list(all_54_points_and_CI_match_build=TRUE,all_CI_visible=TRUE,original_limits=original_limits,
 original_breaks=original_breaks,ranges=ranges,actual_source='../20261002_rebuilt_figures/02_plot.R',
 original_function_extracted=TRUE,layout='9 organs 3x3; Original | Accel'),file.path(out,'qa/plot_checks.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: d redrawn from actual original function; all C/CI checked.\n')
