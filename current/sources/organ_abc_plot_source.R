#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table);library(grid)})
setDTthreads(1)
out <- file.path(Sys.getenv('DISCO_DATA_ROOT'), 'Results_final/05_organ_network/organ_clock_comparison/20260928')
res <- file.path(out,'results'); figs <- file.path(out,'output/pdf'); qa <- file.path(out,'qa')
dir.create(figs,recursive=TRUE,showWarnings=FALSE)
C <- fread(file.path(res,'partial_correlations.csv'))
H <- fread(file.path(res,'Cox_HR.csv'))
map <- fread(file.path(res,'organ_mapping.csv'))[main_figure==TRUE]
co <- fread(file.path(res,'analysis_cohorts.csv'))[main_figure==TRUE]
organs <- map$clock_organ
labs <- setNames(map$display,map$clock_organ)
base <- theme_classic(base_size=8,base_family='Helvetica') + theme(
 plot.title=element_text(size=9,face='bold',hjust=0),
 axis.title=element_text(size=8),axis.text=element_text(size=7,color='black'),
 legend.title=element_text(size=7),legend.text=element_text(size=7),
 strip.background=element_blank(),strip.text=element_text(size=8,face='bold'),
 plot.tag=element_text(size=12,face='bold'),plot.tag.position=c(0,1),
 plot.margin=margin(7,7,5,7),panel.spacing=unit(8,'pt'))

make_panels <- function(adj) {
 cc <- copy(C[main_figure==TRUE & adjustment==adj])
 hh <- copy(H[main_figure==TRUE & adjustment==adj])
 cc[,organ:=factor(organ,levels=organs)]
 cc[,method:=factor(method,levels=c('Goeminne','Oh'))]
 cc[,label:=sprintf('%.2f',partial_r)]
 cc[,ink:=ifelse(abs(partial_r)>.62,'white','black')]
 pa <- ggplot(cc,aes(organ,method,fill=partial_r)) +
   geom_tile(color='white',linewidth=.6) + geom_text(aes(label=label,color=ink),size=2.65) +
   scale_color_identity() + scale_fill_gradient2(low='#3265A8',mid='white',high='#C74742',midpoint=0,limits=c(-1,1),
     breaks=c(-1,0,1),name='Partial r') +
   guides(fill=guide_colorbar(display='rectangles')) +
   scale_x_discrete(labels=labs) + scale_y_discrete(labels=c(Oh='Oh',Goeminne='Goeminne')) +
   labs(x=NULL,y=NULL,title='Age- and sex-adjusted correlation',tag='a') + base +
   theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=30,hjust=1),
     legend.position='right',legend.key.height=unit(12,'pt'),legend.key.width=unit(6,'pt'))
 forest <- function(which_term,title,tag) {
   ss <- copy(hh[term==which_term])
   ss[,stage:=factor(ifelse(model=='joint','Mutually adjusted','Age + sex'),levels=c('Age + sex','Mutually adjusted'))]
   ss[,method:=factor(method,levels=c('Oh','Goeminne'))]
   ss[,yy:=length(organs)+1-match(organ,organs)+ifelse(stage=='Age + sex',.13,-.13)]
   xr <- range(c(ss$lower,ss$upper,1)); lim <- exp(log(xr)+c(-.055,.065))
   breaks <- c(.5,.75,1,1.25,1.5,2,3,4,5)
   breaks <- breaks[breaks>=lim[1] & breaks<=lim[2]]
   if(length(breaks)>5) breaks <- breaks[unique(round(seq(1,length(breaks),length.out=5)))]
   ggplot(ss,aes(x=HR,y=yy,color=stage,shape=stage)) +
     geom_vline(xintercept=1,linetype='dashed',color='#888888',linewidth=.35) +
     geom_segment(aes(x=lower,xend=upper,yend=yy),linewidth=.45) +
     geom_point(size=1.6,stroke=.35) +
     facet_wrap(~method,nrow=1) +
     scale_color_manual(values=c('Age + sex'='#9A9A9A','Mutually adjusted'='#23799C'),name=NULL) +
     scale_shape_manual(values=c('Age + sex'=16,'Mutually adjusted'=17),name=NULL) +
     scale_x_log10(limits=lim,breaks=breaks,labels=function(z)format(z,trim=TRUE,nsmall=2)) +
     scale_y_continuous(breaks=seq_along(organs),labels=rev(unname(labs)),limits=c(.5,length(organs)+.5),expand=c(0,0)) +
     labs(x='HR per 1 SD (95% CI)',y=NULL,title=title,tag=tag) + base +
     theme(legend.position='bottom',legend.key.width=unit(11,'pt'),legend.spacing.x=unit(1,'pt'),
       legend.margin=margin(0,0,0,0),axis.line.y=element_blank(),axis.ticks.y=element_blank()) +
     guides(color=guide_legend(nrow=1),shape=guide_legend(nrow=1))
 }
 pb <- forest('D','Organ DISCO and mortality','b')
 pc <- forest('P','Organ protein age and mortality','c')
 list(a=pa,b=pb,c=pc)
}
save_pdf <- function(p,name,w,h) {
 pdf(file.path(figs,name),width=w,height=h,family='Helvetica',useDingbats=FALSE)
 print(p); dev.off()
}
p <- make_panels('linear_age')
save_pdf(p$a,'a_adjusted_correlations.pdf',7.4,2.2)
save_pdf(p$b,'b_DISCO_HR.pdf',4.8,4.1)
save_pdf(p$c,'c_protein_age_HR.pdf',4.8,4.1)
caption <- paste0('Exploratory UKB associations; published age models; common complete cases within each organ. Training overlap unresolved.\n',
  'All models adjust for age and sex; mutually adjusted models include both scores. HRs summarize follow-up; proportional hazards checks failed.')
combined <- (p$a / (p$b | p$c)) + plot_layout(heights=c(1,2.1)) +
 plot_annotation(caption=caption,theme=theme(plot.caption=element_text(size=7,hjust=0,family='Helvetica'),plot.margin=margin(4,5,4,5)))
save_pdf(combined,'Organ_age_DISCO_abc_EXPLORATORY.pdf',9.6,6.8)
ps <- make_panels('spline_age_df3')
sens <- (ps$a / (ps$b | ps$c)) + plot_layout(heights=c(1,2.1)) +
 plot_annotation(caption='Sensitivity: natural-spline age (df = 3) and sex adjustment. Exploratory UKB associations; training overlap unresolved.',
  theme=theme(plot.caption=element_text(size=7,hjust=0,family='Helvetica'),plot.margin=margin(4,5,4,5)))
save_pdf(sens,'Organ_age_DISCO_abc_SPLINE_SENSITIVITY.pdf',9.6,6.8)
saveRDS(p,file.path(res,'panel_plot_objects.rds'))
cat('Created 5 vector PDFs\n')
