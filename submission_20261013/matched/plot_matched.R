if (nzchar(Sys.getenv('DISCO_R_LIBS'))) .libPaths(c(Sys.getenv('DISCO_R_LIBS'), .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(jsonlite)})
mod <- normalizePath(Sys.getenv('DISCO_RUN_OUTPUT'), mustWork=TRUE)
options(device=function(...)grDevices::cairo_pdf(file=file.path(mod,'qa/measurement.pdf'),family='Arial',...))
stopifnot(read_json(file.path(mod,'results/matched_model_verification.json'))$status=='PASS')
h <- fread(file.path(mod,'results/matched_mortality_HR_72.csv'))
organs <- c('Adipose','Brain','Heart','Intestine','Kidney','Liver','Lung','Muscle','Pancreas')
display <- c('Adipose tissue','Brain','Heart muscle','Intestine','Kidney','Liver','Lung','Skeletal muscle','Pancreas')
limits <- c(min(.93,min(h$lower95)/1.025),max(1.8,max(h$upper95)*1.025))
breaks <- c(.75,1,1.25,1.5,1.75,2)[c(.75,1,1.25,1.5,1.75,2)>=limits[1] & c(.75,1,1.25,1.5,1.75,2)<=limits[2]]
plots <- list()
for(m in c('Oh','Goeminne')){
 d <- copy(h[method==m])
 d[,score_label:=factor(ifelse(score=='Matched DISCO','Matched DISCO',paste(m,'organ clock')),levels=c('Matched DISCO',paste(m,'organ clock')))]
 d[,stage:=factor(ifelse(model=='joint','Mutually adjusted','Fully adjusted'),levels=c('Fully adjusted','Mutually adjusted'))]
 d[,y:=10-match(organ,organs)+ifelse(model=='joint',-.105,.105)]
 p <- ggplot(d,aes(HR,y,color=stage,shape=stage))+
  geom_vline(xintercept=1,linetype='dashed',color='#888888',linewidth=.35)+
  geom_segment(aes(x=lower95,xend=upper95,yend=y),linewidth=.45)+geom_point(size=1.8)+
  facet_grid(.~score_label)+scale_x_log10(limits=limits,breaks=breaks,labels=function(x)format(x,nsmall=2,trim=TRUE))+
  scale_y_continuous(breaks=9:1,labels=display,limits=c(.5,9.5),expand=c(0,0))+
  scale_color_manual(values=c('Fully adjusted'='#9A9A9A','Mutually adjusted'='#23799C'),name=NULL)+
  scale_shape_manual(values=c('Fully adjusted'=16,'Mutually adjusted'=17),name=NULL)+
  labs(title=paste0('Matched DISCO vs ',m,': mortality'),x='HR per 1 training age-residual SD (95% CI)',y=NULL)+
  theme_classic(base_family='Arial',base_size=8)+theme(strip.background=element_blank(),strip.text=element_text(size=8.5),plot.title=element_text(size=9,hjust=.5),axis.text=element_text(color='black',size=8),axis.text.y=element_text(size=8),axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.title.x=element_text(size=8),legend.position='bottom',legend.text=element_text(size=8),panel.spacing=grid::unit(2.5,'mm'),plot.margin=margin(3,3,3,3))
 b <- ggplot_build(p)
 stopifnot(nrow(b$data[[2]])==36,nrow(b$data[[3]])==36,all(is.finite(b$data[[3]]$x)),max(abs(b$data[[2]]$x-log10(d$lower95)))<1e-12,max(abs(b$data[[2]]$xend-log10(d$upper95)))<1e-12,max(abs(b$data[[3]]$x-log10(d$HR)))<1e-12)
 plots[[m]] <- p
}
p <- wrap_plots(plots,ncol=2,guides='collect')+plot_annotation(tag_levels='a') & theme(legend.position='bottom',plot.tag=element_text(size=10,face='bold'))
grDevices::cairo_pdf(file=file.path(mod,'output/Supplementary_Fig_S8_matched_DISCO_clock_HR.pdf'),width=180/25.4,height=90/25.4,family='Arial')
print(p)
invisible(grDevices::dev.off())
ggsave(file.path(mod,'output/Supplementary_Fig_S8_matched_DISCO_clock_HR.png'),p,width=180,height=90,units='mm',dpi=400,bg='white')
write_json(list(status='PASS',points=72,width_mm=180,height_mm=90,all_HR_and_CI_coordinates_verified=TRUE,x_limits=limits,legend=c('Fully adjusted','Mutually adjusted')),file.path(mod,'qa/figure_coordinate_verification.json'),pretty=TRUE,auto_unbox=TRUE)
