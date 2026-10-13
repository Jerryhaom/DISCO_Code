#!/usr/bin/env Rscript
# Reference-layout forest: eight-covariate HRs, explicitly labeled crude C.
options(device=function(...)grDevices::pdf(NULL,...),warn=1)
root <- normalizePath(Sys.getenv('DISCO_RESULTS_ROOT'), mustWork=TRUE)
.libPaths(c(file.path(root,'09_figures_tables/Rlibs'),.libPaths()))
suppressPackageStartupMessages({library(forestploter);library(grid);library(gtable);library(jsonlite)})
module <- normalizePath(Sys.getenv('DISCO_RUN_OUTPUT'), mustWork=TRUE)
main <- read.csv(file.path(module,'tables/Fig3c_8cov_perSD_HR_48.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(main)==48L,all(main$model=='extended_BMI_IPAQ'),all(main$unit=='per 1 SD'))
mods <- c('ClinicalProt','Prot','Clinical')
colors <- c('#00A087FF','#3C5488FF')
all_ci <- range(c(main$low,main$high))
axis_limits <- c(.9,2.3)
if(all_ci[1]<axis_limits[1]) axis_limits[1] <- floor(all_ci[1]*10)/10
if(all_ci[2]>axis_limits[2]) axis_limits[2] <- ceiling(all_ci[2]*10)/10
axis_ticks <- c(.9,1.3,1.7,2.3)
if(axis_limits[1]<.9) axis_ticks <- c(axis_limits[1],axis_ticks)
if(axis_limits[2]>2.3) axis_ticks <- c(axis_ticks,axis_limits[2])
half <- function(diseases,titles) {
  ds <- list();ests <- lows <- ups <- list(DISCO=numeric(),DM=numeric())
  for(i in seq_along(diseases)) {
    ds[[length(ds)+1L]] <- c(titles[i],rep('',6))
    for(m in c('DISCO','DM')) {
      ests[[m]] <- c(ests[[m]],NA);lows[[m]] <- c(lows[[m]],NA);ups[[m]] <- c(ups[[m]],NA)
    }
    for(mod in mods) {
      row <- main[main$disease==diseases[i]&main$display_modality==mod,]
      a <- row[row$metric=='DISCO',];b <- row[row$metric=='DM',]
      stopifnot(nrow(a)==1L,nrow(b)==1L,a$N==b$N,a$events==b$events)
      ds[[length(ds)+1L]] <- c(paste0('    ',mod),sprintf('%d (%.2f%%)',a$events,a$events_pct),
        sprintf('%.2f',a$HR),sprintf('%.3f',a$C_crude),sprintf('%.2f',b$HR),sprintf('%.3f',b$C_crude),strrep(' ',27))
      for(m in c('DISCO','DM')) {
        t <- row[row$metric==m,]
        ests[[m]] <- c(ests[[m]],t$HR);lows[[m]] <- c(lows[[m]],t$low);ups[[m]] <- c(ups[[m]],t$high)
      }
    }
  }
  tab <- as.data.frame(do.call(rbind,ds),stringsAsFactors=FALSE)
  names(tab) <- c('Diseases/Metrics','Events (%)','DISCO\nHR','DISCO\nC-index','DM\nHR','DM\nC-index','')
  tm <- forest_theme(base_size=15,base_family='Helvetica',ci_pch=15,ci_fill=colors,ci_col=colors,
    ci_lwd=2.5,ci_lty=1,legend_position='none',legend_value=c('DISCO','DM'),
    refline_gp=gpar(lwd=1,lty='dashed',col='grey20'),
    core=list(bg_params=list(fill='white'),padding=unit(c(2,3),'mm')),
    colhead=list(fg_params=list(fontface='bold'),padding=unit(c(2,3),'mm')))
  p <- forest(tab,est=ests,lower=lows,upper=ups,sizes=.7,ci_column=7,ref_line=1,nudge_y=.4,
    x_trans='log',xlim=axis_limits,ticks_at=axis_ticks,ticks_digits=1,theme=tm)
  p <- add_border(p,row=1:16,col=1:6,where='top',gp=gpar(lwd=.5))
  p <- add_border(p,part='header',col=1:6,where='bottom',gp=gpar(lwd=1.5))
  p <- edit_plot(p,row=c(1,5,9,13),col=1,which='text',gp=gpar(fontface='bold'))
  stopifnot(length(p$widths)==9L)
  p$widths[1] <- unit(3,'mm');p$widths[9] <- unit(3,'mm')
  p$widths[2:8] <- unit(c(1.90,1.65,.70,.85,.60,.85,2.20),'in')
  p$heights[2] <- unit(.60,'in');p$heights[3:18] <- unit(.36,'in')
  legend <- grobTree(
    segmentsGrob(x0=.16,x1=.24,y0=.75,y1=.75,gp=gpar(col=colors[1],lwd=2.5)),
    pointsGrob(x=.20,y=.75,pch=15,size=unit(9,'pt'),gp=gpar(col=colors[1])),
    textGrob('DISCO',x=.26,y=.75,just='left',gp=gpar(fontfamily='Helvetica',fontsize=15)),
    segmentsGrob(x0=.16,x1=.24,y0=.30,y1=.30,gp=gpar(col=colors[2],lwd=2.5)),
    pointsGrob(x=.20,y=.30,pch=15,size=unit(9,'pt'),gp=gpar(col=colors[2])),
    textGrob('DM',x=.26,y=.30,just='left',gp=gpar(fontfamily='Helvetica',fontsize=15)))
  p <- gtable_add_grob(p,legend,t=2,l=8,clip='off',name='two_metric_legend')
  for(k in which(p$layout$name=='colhead-fg' & p$layout$l %in% 4:7)) { p$grobs[[k]]$x <- unit(.5,'npc'); p$grobs[[k]]$hjust <- .5 }
  list(plot=p,table=tab,est=ests,lower=lows,upper=ups)
}
left <- half(c('obesity','hypertension','diabetes','cancer'),c('Obesity','Hypertension','Type 2 diabetes','Cancer'))
right <- half(c('mi','stroke','copd','dementia'),c('MI','Stroke','COPD','Dementia'))
# Preserve every text style from the latest font-plus-one reference.
reference <- readRDS(file.path(root,'03_multimodal_UKB/Fig3c_font_plus1_20261010/results/plot_objects.rds'))
copy_text_style <- function(a,b) {
  if(inherits(a,'text')) { stopifnot(inherits(b,'text')); a$gp <- b$gp }
  if(!is.null(a$children)) {stopifnot(length(a$children)==length(b$children));for(i in seq_along(a$children)) a$children[[i]] <- copy_text_style(a$children[[i]],b$children[[i]])}
  if(!is.null(a$grobs)) {stopifnot(length(a$grobs)==length(b$grobs));for(i in seq_along(a$grobs)) a$grobs[[i]] <- copy_text_style(a$grobs[[i]],b$grobs[[i]])}
  a
}
left$plot <- copy_text_style(left$plot,reference$left$plot)
right$plot <- copy_text_style(right$plot,reference$right$plot)
pdf(file.path(module,'output/pdf/Fig3c_REESTIMATED.pdf'),width=20,height=7,family='Helvetica',
  title='Fig3c | age, sex, BMI, smoking, drinking, education, ethnicity and IPAQ activity | per 1 SD')
grid.newpage()
pushViewport(viewport(x=.267,y=.51,width=.449,height=.94));grid.draw(left$plot);popViewport()
pushViewport(viewport(x=.733,y=.51,width=.449,height=.94));grid.draw(right$plot);popViewport()
dev.off()
saveRDS(list(left=left,right=right,plotted=main),file.path(module,'results/plot_objects.rds'))
write.csv(main,file.path(module,'qa/plotted_48.csv'),row.names=FALSE)
write_json(list(status='DRAWN_PENDING_VISUAL_QA',rows=48,width_pt=1440,height_pt=504,
  reference_body_font_pt=15,axis_limits=axis_limits,axis_ticks=axis_ticks,
  all_48_CI_within_axes=TRUE,unadjusted_C_retained=TRUE,events_from_adjusted_cohort=TRUE),
  file.path(module,'results/plot_geometry.json'),auto_unbox=TRUE,pretty=TRUE)
