#!/usr/bin/env Rscript
# New sensitivity displays. Not historical KM/attained-age figure replacements.
Sys.umask('0077')
root <- Sys.getenv('DISCO_RESULTS_INPUT')
out <- Sys.getenv('DISCO_RUN_MODULE')
suppressPackageStartupMessages({library(ggplot2);library(survival);library(ggsci);library(patchwork);library(jsonlite)})
comparison <- read.csv(file.path(out,'tables/CHARLS_endpoint_sensitivity_comparison.csv'))
keep <- c('primary_original','semiparametric_interval_PH','retrospective_semiparametric_interval_PH','weibull_interval_PH_nominal_spacing','exponential_interval_PH_nominal_spacing')
labels <- c('Original Cox','Interval PH','Retrospective interval PH','Weibull PH','Exponential PH')
d <- comparison[comparison$analysis%in%keep,]
d$method <- factor(labels[match(d$analysis,keep)],levels=rev(labels))
d$metric <- factor(d$metric,levels=c('DISCO','DM'))
p <- ggplot(d,aes(HR,method,color=metric))+
 geom_vline(xintercept=1,linewidth=.35,linetype=2,color='grey60')+
 geom_errorbar(aes(xmin=lower95,xmax=upper95),width=.13,orientation='y',linewidth=.8)+geom_point(size=3.3)+
 facet_wrap(~metric,nrow=1)+scale_color_manual(values=c(DISCO='#E64B35',DM='#4DBBD5'))+
 scale_x_continuous(limits=c(1,1.65),breaks=c(1,1.2,1.4,1.6))+
 labs(title='CHARLS',x='Mortality HR per original 1 SD',y=NULL)+
 theme_classic(base_size=14)+theme(legend.position='none',plot.title=element_text(size=18,face='bold'),
  axis.title.x=element_text(size=16),axis.text=element_text(size=13),strip.text=element_text(size=16),
  strip.background=element_blank(),panel.spacing=grid::unit(1.4,'lines'))
pdf(file.path(out,'figures/CHARLS_interval_endpoint_HR_sensitivity.pdf'),width=9.5,height=4.4,family='Helvetica',useDingbats=FALSE);print(p);dev.off()
frames <- readRDS(file.path(out,'private/input/panel_interval_group_frames.rds'))
curve_rows <- list();model_summaries <- list()
grid <- c(0,2,4,7)
for(panel in names(frames)){
 d <- frames[[panel]];d <- d[!d$no_postbaseline_followup,]
 groups <- if(panel=='ED3b')c('DISCO Low Risk','DISCO High Risk','DM Low Risk','DM High Risk') else levels(d$group)
 for(g in groups){
  z <- d[as.character(d$group)==g,]
  fit <- survfit(Surv(lower_nominal,upper_nominal,type='interval2')~1,data=z)
  s <- summary(fit,times=grid,extend=TRUE)
  stopifnot(length(s$surv)==4,all(diff(s$surv)<=1e-9),all(s$surv>=0&s$surv<=1))
  lab <- if(panel=='ED3b')sub('High Risk','Q4',sub('Low Risk','Q1',g)) else
    paste(sub(' [14]$','',g),paste0('DISCO Q',sub('^.* ','',g)),sep=', ')
  curve_rows[[length(curve_rows)+1L]] <- data.frame(panel,group=g,label=lab,wave_index=seq_along(grid),
   nominal_coordinate=grid,display_wave=c('2011','2013','2015','2018'),survival=s$surv,
   definition='Turnbull NPMLE evaluated at observation-wave boundaries; connecting lines are guides, not exact event-time KM')
  model_summaries[[length(model_summaries)+1L]]<-data.frame(panel,group=g,N=nrow(z),deaths=sum(z$status),
   followup='Interval censored (last confirmed alive, first dead); right censored last alive',estimator='survival::survfit interval2 Turnbull NPMLE')
 }
}
curves <- do.call(rbind,curve_rows)
write.csv(curves,file.path(out,'results/interval_survival_at_observation_waves.csv'),row.names=FALSE)
write.csv(do.call(rbind,model_summaries),file.path(out,'results/interval_survival_curve_group_QA.csv'),row.names=FALSE)
make_panel <- function(panel,title){
 q <- curves[curves$panel==panel,]
 q$label <- factor(q$label,levels=unique(q$label))
 ggplot(q,aes(wave_index,survival,color=label,group=label))+
  geom_line(linewidth=.8)+geom_point(size=2)+scale_color_npg()+
  scale_x_continuous(breaks=1:4,labels=c('2011','2013','2015','2018'),expand=expansion(mult=c(.03,.03)))+
  scale_y_continuous(limits=c(.5,1),breaks=seq(.5,1,.1))+
  labs(title=title,x='Observed survey wave',y='Interval-estimated survival',color=NULL)+
  theme_classic(base_size=14)+theme(plot.title=element_text(size=17,face='bold'),axis.title=element_text(size=16),
   axis.text=element_text(size=13),legend.position='bottom',legend.text=element_text(size=12),legend.key.width=grid::unit(1,'lines'))+
  guides(color=guide_legend(ncol=2,byrow=TRUE))
}
curves_plot <- make_panel('ED3b','CHARLS: DISCO / DM') + make_panel('ED8b','CHARLS: DISCO + frailty')+
 plot_layout(widths=c(1,1))
pdf(file.path(out,'figures/CHARLS_observed_wave_interval_survival_sensitivity.pdf'),width=11,height=5,family='Helvetica',useDingbats=FALSE);print(curves_plot);dev.off()
write_json(list(figures=list('CHARLS_interval_endpoint_HR_sensitivity.pdf','CHARLS_observed_wave_interval_survival_sensitivity.pdf'),
  forest_uses_approved_primary_SD=TRUE,interval_PH_CI='200 reproducible bootstrap resamples; not interval-width-as-error-bar',
  survival_curve_estimator='Turnbull NPMLE evaluated at observation-wave boundaries',
  survival_CI_shown=FALSE,survival_lines_are_guides_between_observed_waves=TRUE,
  forest_retrospective_model='10 unambiguous solar dates before2018 recovered from2020 exit, 1 conflicting date excluded',
  survival_curve_branch='Observed-wave1096-death branch; retrospective1106-death branch only in forest/table, not these curves',
  source_table_sha256=digest::digest(file=file.path(out,'tables/CHARLS_endpoint_sensitivity_comparison.csv'),algo='sha256'),
  source_group_frame_sha256=digest::digest(file=file.path(out,'private/input/panel_interval_group_frames.rds'),algo='sha256'),
  plot_code_sha256=digest::digest(file=file.path(out,'code/plot_interval_sensitivity.R'),algo='sha256'),
  original_ED3b_mapping='Result4/Validate_Charls.R:205-220 score quartiles; 257-299 attained-age legacy KM producer',
  original_ED8b_mapping='Result5/Asso_DISCO_FI.R:569-581 frailty/Q1-Q4 selection; 609-640 attained-age legacy KM producer',
  minimal_style_reuse='Original panel grouping and NPG color order, large scientific-figure labels; new interval estimand requires wave axis, not attained-age axis',
  PDF_device='grDevices::pdf Helvetica; cairo shared library unavailable because local X11 libraries are absent; statistical inputs unchanged',
  original_figures_and_CURRENT_not_overwritten=TRUE),file.path(out,'results/figure_scope.json'),auto_unbox=TRUE,pretty=TRUE)
