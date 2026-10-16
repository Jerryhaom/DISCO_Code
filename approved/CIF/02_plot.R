#!/usr/bin/env Rscript
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);here<-dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE)))
out<-file.path(here,'output/pdf');dir.create(out,recursive=TRUE,showWarnings=FALSE)
a<-fread(file.path(here,'pass_1/curves_aggregate.csv'))
labels<-c('Obesity','Hypertension','Diabetes','Cancer','Myocardial infarction','Stroke','COPD','Dementia')
groups<-c('DISCO Q1','DISCO Q4','DM Q1','DM Q4')
cols<-c('DISCO Q1'='#E64B35FF','DISCO Q4'='#E64B35FF','DM Q1'='#4DBBD5FF','DM Q4'='#4DBBD5FF')
ltys<-c('DISCO Q1'='solid','DISCO Q4'='dashed','DM Q1'='solid','DM Q4'='dashed')
# Same original ED7 encodings and 2x4 organization; no fitted C-index annotation.
make<-function(scale,method){
 q<-a[policy=='exclude_ambiguous' & time_scale==scale & n_risk>=20]
 q[,label:=factor(label,levels=labels)];q[,group:=factor(group,levels=groups)]
 q[,`:=`(estimate=get(method),lo=get(paste0(method,'_low')),hi=get(paste0(method,'_high')))]
 g<-ggplot(q,aes(time,estimate,colour=group,fill=group,linetype=group))+
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.10,colour=NA,show.legend=FALSE)+
  geom_step(linewidth=.65)+facet_wrap(~label,ncol=4,scales='free_y')+
  scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+scale_linetype_manual(values=ltys)+
  scale_y_continuous(labels=function(x)sprintf('%.0f%%',100*x),expand=expansion(mult=c(.01,.05)))+
  labs(x=if(scale=='followup')'Follow-up (years)' else 'Age (years)',
   y=if(method=='CIF')'Disease cumulative incidence' else '1 - Kaplan-Meier',title='UKB')+
  theme_classic(base_size=13,base_family='Helvetica')+
  theme(plot.title=element_text(size=17,face='bold',hjust=.5),strip.text=element_text(size=12,face='bold'),
   strip.background=element_blank(),panel.spacing=unit(1.3,'lines'),
   axis.title=element_text(size=15),axis.text=element_text(size=11,colour='black'),
   legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=13),
   legend.key.width=unit(1.6,'cm'),legend.margin=margin(t=5,b=0),plot.margin=margin(12,14,10,10))
 if(scale=='followup')g<-g+scale_x_continuous(breaks=c(0,5,10,15),limits=c(0,15.7))
 g
}
for(scale in c('followup','attained_age'))for(method in c('CIF','KM')){
 name<-paste0('ED7_',method,'_sensitivity_',scale,'.pdf')
 ggsave(file.path(out,name),make(scale,method),width=11,height=6,
  device=function(filename,...)grDevices::pdf(file=filename,...,useDingbats=FALSE))
}
write_json(list(primary='ED7_CIF_sensitivity_followup.pdf',companion='ED7_CIF_sensitivity_attained_age.pdf',
 comparators='KM PDFs use identical first-event risk sets; not replacements for the existing ED7 artwork',
 curve_tail_rule='Plot while at least 20 participants remain at risk in the disease-free state',
 confidence='95% confidence intervals from survival::survfit',policy='exclude_ambiguous',
 original_dates_verified=FALSE,adjustment='Unadjusted quartile curves',
 device='R vector pdf, Helvetica, useDingbats=FALSE; Cairo unavailable without X11'),
 file.path(here,'plot_definitions.json'),pretty=TRUE,auto_unbox=TRUE)
