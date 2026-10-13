#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);out<-dirname(normalizePath(entry))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table);library(jsonlite);library(grid)})
setDTthreads(1)
map<-fread(file.path(file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison'),'20260928/results/organ_mapping.csv'))[main_figure==TRUE]
organs<-map$clock_organ;labs<-setNames(map$display,map$clock_organ)
C<-fread(file.path(out,'numeric/panel_a_correlations.csv'))
H<-rbind(fread(file.path(out,'numeric/panel_b_HR.csv')),fread(file.path(out,'numeric/panel_c_HR.csv')))
# Reuse the actual original plotting function and theme, not a filename guess.
src<-parse(file.path(Sys.getenv('DISCO_CODE_ROOT'),'sources/organ_abc_plot_source.R'))
for(e in src)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('base','make_panels'))eval(e)
base<-base+theme(plot.title=element_text(size=9,face='plain'))
p<-make_panels('linear_age')
p$a$layers[[1]]<-geom_tile(width=1,height=1,color='white',linewidth=.22)
p$a<-p$a+ggplot2::labs(tag=NULL)+theme(axis.text.x=element_text(size=9.5,angle=30,hjust=1,color='black'),axis.text.y=element_text(size=9.5,color='black'))
for(n in c('b','c'))p[[n]]<-p[[n]]+ggplot2::labs(tag=NULL,x='HR per age-residual SD (95% CI)')
savepdf<-function(g,name,w,h){pdf(file.path(out,'panels',name),width=w,height=h,family='Helvetica',useDingbats=FALSE);print(g);invisible(dev.off())}
savepdf(p$a,'a_correlations.pdf',9.6,2.2)
savepdf(p$b,'b_DISCO_HR.pdf',4.8,4.1)
savepdf(p$c,'c_ProtAge_HR.pdf',4.8,4.1)

# Retain the actual v2 d design: 3x3 organs, two columns inside each organ.
d<-fread(file.path(out,'numeric/panel_d_Cindex.csv'),data.table=FALSE)
d$score<-factor(d$score,levels=c('Goeminne','Oh','DISCO'))
font<-theme_classic(base_family='Helvetica',base_size=7)+theme(
 axis.title=element_text(size=8),axis.text=element_text(size=7,color='black'),
 plot.title=element_text(size=8,hjust=.5,face='plain'),
 plot.margin=margin(2,8,2,8),axis.line=element_line(linewidth=.35),
 axis.ticks=element_line(linewidth=.3),axis.ticks.length=grid::unit(2,'pt'))
ranges<-builtdata<-list()
make_one<-function(organ,grp){
 sub<-d[d$organ==organ&d$group==grp,]
 if(grp=='Age + sex'){lim<-c(.58,.82);breaks<-c(.60,.70,.80)} else {
  lim<-c(floor((min(sub$lower95)-.008)*100)/100,ceiling((max(sub$upper95)+.008)*100)/100)
  breaks<-c(lim[1],round(mean(lim),2),lim[2])
 }
 stopifnot(all(sub$lower95>lim[1]),all(sub$upper95<lim[2]))
 ranges[[paste(organ,grp)]]<<-list(organ=organ,group=grp,limits=lim,breaks=breaks)
 yl<-c(DISCO='DISCO',Oh='Oh ProtAge',Goeminne='Goeminne ProtAge')
 q<-ggplot(sub,aes(C,score,fill=score))+geom_col(width=.65)+
  geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.24,linewidth=.32)+
  scale_fill_manual(values=c(DISCO='#E64B35',Oh='grey70',Goeminne='grey70'),guide='none')+
  scale_y_discrete(labels=yl,expand=expansion(add=.48))+
  scale_x_continuous(breaks=breaks,labels=function(z)sprintf('%.2f',z),expand=expansion(mult=0))+
  coord_cartesian(xlim=lim)+ggplot2::labs(x=NULL,y=NULL,title=grp)+font
 if(grp=='Accel')q<-q+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
 bb<-ggplot_build(q);stopifnot(max(abs(bb$data[[1]]$x-sub$C))<1e-12,
  max(abs(bb$data[[2]]$xmin-sub$lower95))<1e-12,max(abs(bb$data[[2]]$xmax-sub$upper95))<1e-12)
 builtdata[[paste(organ,grp)]]<<-bb$data
 q
}
cells<-lapply(organs,function(organ){
 q<-make_one(organ,'Age + sex')+make_one(organ,'Accel')+plot_layout(widths=c(1,1))+
  plot_annotation(title=labs[organ],theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=.5),plot.margin=margin(6,7,3,7)))
 wrap_elements(full=q)
})
p$d<-wrap_plots(cells,ncol=3,nrow=3)+plot_annotation(title='Mortality discrimination',caption='C-index (95% CI)',
 theme=theme(plot.title=element_text(family='Helvetica',size=9,face='plain',hjust=0,margin=margin(l=54)),
 plot.caption=element_text(family='Helvetica',size=9,hjust=.5),plot.margin=margin(6,7,6,7)))
savepdf(p$d,'d_Cindex.pdf',9.6,4.6)
# Test every plotted value, label and confidence interval against its input.
aa<-ggplot_build(p$a);stopifnot(identical(aa$data[[2]]$label,sprintf('%.2f',p$a$data$partial_r)),
 all(abs(aa$data[[1]]$xmax-aa$data[[1]]$xmin-1)<1e-10),all(abs(aa$data[[1]]$ymax-aa$data[[1]]$ymin-1)<1e-10))
for(n in c('b','c')){
 bb<-ggplot_build(p[[n]]);ss<-p[[n]]$data
 stopifnot(nrow(ss)==36,max(abs(bb$data[[2]]$x-log10(ss$lower)))<1e-12,
 max(abs(bb$data[[2]]$xend-log10(ss$upper)))<1e-12,max(abs(bb$data[[3]]$x-log10(ss$HR)))<1e-12)
}
saveRDS(list(plots=p,C=C,H=H,d=d,ranges=ranges),file.path(out,'plot_source.rds'))
write_json(list(all_18_heatmap_labels_checked=TRUE,all_72_HR_points_and_CI_checked=TRUE,
 all_54_Cindex_points_and_CI_checked=TRUE,all_d_intervals_visible=TRUE,heatmap_cells_touch=TRUE,
 font='Helvetica',body_pt=7,title_pt=9,heatmap_axis_pt=9.5,layout='a; b|c; d 3x3',ranges=ranges,
 actual_plot_source_abc='../20260928/03_plot_panels.R',actual_plot_source_d='../20261002_layout_v2/01_plot.R'),
 file.path(out,'qa/plot_value_checks.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: 4 vector panels; all point estimates, labels and CI endpoints checked.\n')
