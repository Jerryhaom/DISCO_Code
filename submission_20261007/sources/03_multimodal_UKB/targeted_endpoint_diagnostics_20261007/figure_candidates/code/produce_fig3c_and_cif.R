#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Reuse the actual near-source forest drawing block; no recovered forest11 claim.
options(device=function(...) grDevices::pdf(file=NULL,...), warn=1)
.libPaths(c(.disco_path("results","09_figures_tables/Rlibs"), .libPaths()))
suppressPackageStartupMessages({library(data.table);library(forestploter);library(grid)
  library(gtable);library(jsonlite);library(ggplot2);library(digest)})
base <- .disco_path("results","03_multimodal_UKB/targeted_endpoint_diagnostics_20261007")
out <- file.path(base,'figure_candidates')
for (d in c('inputs','output/pdf','output/png','qa')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
source_code <- .disco_path("results","03_multimodal_UKB/reproduction/code/plot_forest.R")
src <- c(file.path(base,'Table19_Fig3c_raw_unit_candidate.csv'),
         file.path(base,'three_modalities_raw_unit_models.csv'),
         file.path(base,'ClinicalProt_CIF_aggregate.csv'),file.path(base,'ClinicalProt_CIF_QA.csv'),source_code)
for (p in src) stopifnot(file.copy(p,file.path(out,'inputs',basename(p)),overwrite=TRUE))
sha <- function(p) digest(file=p,algo='sha256')
source_hashes <- vapply(src,sha,'')
raw <- fread(file.path(out,'inputs/three_modalities_raw_unit_models.csv'))
adjusted <- fread(file.path(out,'inputs/Table19_Fig3c_raw_unit_candidate.csv'))
stopifnot(nrow(adjusted)==96L,all(adjusted$unit=='one saved score unit'),
          all(adjusted$policy=='death_first'),all(adjusted$version=='raw_dates'))
mod_map <- c(clinicalprot='ClinicalProt',proteomics='Prot',clinical='Clinical')
mods <- c('ClinicalProt','Prot','Clinical')
colors <- c('#00A087FF','#3C5488FF')
key <- function(d) paste(d$disease,d$modality,d$metric,sep='|')
crude <- raw[model=='crude' & metric %in% c('DISCO','DM')]
stopifnot(nrow(crude)==48L)
all_plotted <- list();plot_check <- list()
draw_forest <- function(selected) {
  main <- if (selected=='extended_BMI') copy(adjusted[model==selected]) else copy(crude)
  stopifnot(nrow(main)==48L,!anyDuplicated(key(main)))
  j <- match(key(main),key(crude));stopifnot(!anyNA(j))
  main[,`:=`(C_crude=crude$C[j],N_crude=crude$N[j],events_crude=crude$events[j],
              percent_crude=100*crude$events[j]/crude$N[j],
              lower_ci=low,upper_ci=high,original_modality=modality)]
  main[,modality:=unname(mod_map[modality])]
  stopifnot(!anyNA(main$modality),all(main$C_crude>0 & main$C_crude<1))
  axis_lim <- if (selected=='extended_BMI') c(.9,2.3) else c(.9,2.5)
  axis_ticks <- if (selected=='extended_BMI') c(.9,1.3,1.7,2.3) else c(.9,1.3,1.7,2.1,2.5)
  # The half() drawing segment below follows actual plot_forest.R lines29-88.
  # Only mapped input names and auxiliary axis limits differ.
  half <- function(diseases,titles) {
    ds<-list();ests<-lows<-ups<-list(DISCO=numeric(),DM=numeric())
    for(i in seq_along(diseases)) {
      ds[[length(ds)+1L]]<-c(titles[i],rep('',6))
      for(m in c('DISCO','DM')) {
        ests[[m]]<-c(ests[[m]],NA);lows[[m]]<-c(lows[[m]],NA);ups[[m]]<-c(ups[[m]],NA)
      }
      for(mod in mods) {
        row<-main[main$disease==diseases[i]&main$modality==mod,]
        a<-row[row$metric=='DISCO',];b<-row[row$metric=='DM',]
        stopifnot(nrow(a)==1L,nrow(b)==1L,a$N_crude==b$N_crude,a$events_crude==b$events_crude)
        ds[[length(ds)+1L]]<-c(paste0('    ',mod),sprintf('%d (%.2f%%)',a$events_crude,a$percent_crude),
          sprintf('%.2f',a$HR),sprintf('%.3f',a$C_crude),sprintf('%.2f',b$HR),sprintf('%.3f',b$C_crude),strrep(' ',27))
        for(m in c('DISCO','DM')) {
          t<-row[row$metric==m,]
          ests[[m]]<-c(ests[[m]],t$HR);lows[[m]]<-c(lows[[m]],t$lower_ci);ups[[m]]<-c(ups[[m]],t$upper_ci)
        }
      }
    }
    tab<-as.data.frame(do.call(rbind,ds),stringsAsFactors=FALSE)
    names(tab)<-c('Diseases/Metrics','Events (%)\nCrude cohort','DISCO\nHR','DISCO\nC-index','DM\nHR','DM\nC-index','')
    stopifnot(all(unlist(lows)>0,na.rm=TRUE))
    range_ci<-range(c(unlist(lows),unlist(ups)),na.rm=TRUE)
    stopifnot(range_ci[1]>=axis_lim[1],range_ci[2]<=axis_lim[2])
    tm<-forest_theme(base_size=12,base_family='Helvetica',ci_pch=15,ci_fill=colors,ci_col=colors,
      ci_lwd=1.5,ci_lty=1,legend_position='none',legend_value=c('DISCO','DM'),
      refline_gp=gpar(lwd=1,lty='dashed',col='grey20'),
      core=list(bg_params=list(fill='white'),padding=unit(c(2,3),'mm')),
      colhead=list(fg_params=list(fontface='plain'),padding=unit(c(2,3),'mm')))
    p<-forest(tab,est=ests,lower=lows,upper=ups,sizes=.4,ci_column=7,ref_line=1,
      x_trans='log',xlim=axis_lim,ticks_at=axis_ticks,ticks_digits=1,theme=tm)
    p<-add_border(p,row=1:16,col=1:6,where='top',gp=gpar(lwd=.5))
    p<-add_border(p,part='header',col=1:6,where='bottom',gp=gpar(lwd=1.5))
    stopifnot(length(p$widths)==9L)
    p$widths[2:8]<-unit(c(1.60,1.60,.65,.79,.57,.79,2.12),'in')
    p$heights[2]<-unit(.43,'in');p$heights[3:18]<-unit(.218,'in')
    legend<-grobTree(
      segmentsGrob(x0=.10,x1=.18,y0=.75,y1=.75,gp=gpar(col=colors[1],lwd=1.5)),
      pointsGrob(x=.14,y=.75,pch=15,size=unit(7,'pt'),gp=gpar(col=colors[1])),
      textGrob('DISCO',x=.20,y=.75,just='left',gp=gpar(fontfamily='Helvetica',fontsize=12)),
      segmentsGrob(x0=.10,x1=.18,y0=.30,y1=.30,gp=gpar(col=colors[2],lwd=1.5)),
      pointsGrob(x=.14,y=.30,pch=15,size=unit(7,'pt'),gp=gpar(col=colors[2])),
      textGrob('DM',x=.20,y=.30,just='left',gp=gpar(fontfamily='Helvetica',fontsize=12)))
    p<-gtable_add_grob(p,legend,t=2,l=8,clip='off',name='two_metric_legend')
    list(plot=p,table=tab,est=ests,lower=lows,upper=ups)
  }
  left<-half(c('obesity','hypertension','diabetes','cancer'),c('Obesity','Hypertension','Type 2 diabetes','Cancer'))
  right<-half(c('mi','stroke','copd','dementia'),c('MI','Stroke','COPD','Dementia'))
  name<-if(selected=='extended_BMI') 'Fig3c_raw_unit_expanded_BMI_NEAR_SOURCE' else 'Fig3c_raw_unit_crude_NEAR_SOURCE'
  pdf(file.path(out,'output/pdf',paste0(name,'.pdf')),width=20,height=5,family='Helvetica',
      title=paste('Near-source new endpoint Fig3c',selected,'not recovered historical forest11'))
  grid.newpage()
  pushViewport(viewport(x=.255,y=.51,width=.47,height=.92));grid.draw(left$plot);popViewport()
  pushViewport(viewport(x=.745,y=.51,width=.47,height=.92));grid.draw(right$plot);popViewport()
  dev.off()
  fwrite(main,file.path(out,'qa',paste0(name,'_plotted.csv')))
  saveRDS(list(left=left[c('table','est','lower','upper')],right=right[c('table','est','lower','upper')],
               plotted=main,model=selected),file.path(out,'qa',paste0(name,'_objects.rds')))
  all_plotted[[selected]] <<- main
  plot_check[[selected]] <<- list(rows=nrow(main),xlim=axis_lim,
      max_ci=max(main$high),min_ci=min(main$low),clipped_confidence_intervals=0L)
}
draw_forest('extended_BMI');draw_forest('crude')

cif<-fread(file.path(out,'inputs/ClinicalProt_CIF_aggregate.csv'))
cif_qa<-fread(file.path(out,'inputs/ClinicalProt_CIF_QA.csv'))
stopifnot(nrow(cif)==1360L,nrow(cif_qa)==80L,max(cif_qa$max_AJ_error)<1e-12)
disease_order<-c('obesity','hypertension','diabetes','cancer','mi','stroke','copd','dementia')
disease_title<-c('Obesity','Hypertension','Type 2 diabetes','Cancer','MI','Stroke','COPD','Dementia')
z<-cif[policy=='death_first' & group!='All' & time_months<=168]
stopifnot(nrow(z)==480L,all(z$CIF<=z$KM+1e-12),all(z$CIF_low<=z$CIF),all(z$CIF_high>=z$CIF),
          all(z$CIF_low>=0),all(z$CIF_high<=1),all(z$at_risk>=.1*z$N))
z[,`:=`(years=time_months/12,disease=factor(disease,levels=disease_order,labels=disease_title),
         Score=ifelse(grepl('qDM',group),'DM','DISCO'),Quartile=ifelse(grepl('_1$',group),'Q1','Q4'))]
fwrite(z,file.path(out,'qa/ClinicalProt_CIF_plotted.csv'))
counts<-copy(cif_qa[policy=='death_first'])
counts[,disease_title:=disease_title[match(disease,disease_order)]]
fwrite(counts,file.path(out,'qa/ClinicalProt_CIF_counts.csv'))
risk_counts<-z[time_months %in% c(0,60,120,168),.(Score,disease,Quartile,years,N,at_risk)]
fwrite(risk_counts,file.path(out,'qa/ClinicalProt_CIF_risk_table.csv'))
plot_cif<-function(score) {
  d<-copy(z[Score==score]);d[,Quartile:=factor(paste(score,Quartile),levels=paste(score,c('Q1','Q4')))]
  d[,CIF_pct:=100*CIF]
  lines<-rbind(d[,.(disease,years,Quartile,risk=100*CIF,Estimator='CIF')],
               d[,.(disease,years,Quartile,risk=100*KM,Estimator='1 - KM')])
  lines[,Estimator:=factor(Estimator,levels=c('CIF','1 - KM'))]
  ggplot()+
    geom_ribbon(data=d,aes(years,ymin=100*CIF_low,ymax=100*CIF_high,fill=Quartile),alpha=.15,show.legend=FALSE)+
    geom_line(data=lines,aes(years,risk,color=Quartile,linetype=Estimator,group=interaction(Quartile,Estimator)),linewidth=.8)+
    facet_wrap(~disease,ncol=4,scales='free_y')+
    scale_color_manual(values=setNames(c('#2878B5','#D65F28'),paste(score,c('Q1','Q4'))))+
    scale_fill_manual(values=setNames(c('#2878B5','#D65F28'),paste(score,c('Q1','Q4'))))+
    scale_linetype_manual(values=c('CIF'='solid','1 - KM'='22'))+
    guides(color=guide_legend(order=1),linetype=guide_legend(order=2))+
    scale_x_continuous(breaks=c(0,5,10,14),limits=c(0,14),expand=expansion(mult=c(0,.01)))+
    scale_y_continuous(expand=expansion(mult=c(0,.09)))+
    labs(title=paste('Clinical + proteomics:',score),x='Follow-up (years)',
         y='Cumulative disease probability (%)',color=NULL,linetype=NULL)+
    theme_classic(base_size=16,base_family='Helvetica')+
    theme(plot.title=element_text(size=20,face='plain',hjust=0,margin=margin(b=14)),
          axis.title=element_text(size=18),axis.title.x=element_text(margin=margin(t=10)),
          axis.title.y=element_text(margin=margin(r=10)),axis.text=element_text(size=13,color='black'),
          strip.background=element_blank(),strip.text=element_text(size=16,face='plain',margin=margin(b=7)),
          legend.position='top',legend.justification='left',legend.text=element_text(size=15),
          legend.key.width=unit(1.8,'lines'),legend.box='horizontal',panel.spacing=unit(1.15,'lines'),
          axis.line=element_line(linewidth=.5),axis.ticks=element_line(linewidth=.5),
          plot.margin=margin(15,18,13,15))
}
plots<-list(DISCO=plot_cif('DISCO'),DM=plot_cif('DM'))
pdf(file.path(out,'output/pdf/ClinicalProt_CIF_sensitivity.pdf'),width=12.8,height=7.4,family='Helvetica',
    title='Clinical plus proteomics, death competing risk sensitivity, annual grid')
for (p in plots) print(p)
dev.off()
saveRDS(list(plots=plots,data=z,qa=cif_qa),file.path(out,'qa/ClinicalProt_CIF_plot_objects.rds'))
stopifnot(identical(source_hashes,vapply(src,sha,'')))
write_json(list(status='NEW_ENDPOINT_NEAR_SOURCE_CANDIDATE',
  forest_source=source_code,forest_source_sha256=sha(source_code),
  original_forest11_recovered=FALSE,exact_historical_pixels_claimed=FALSE,
  sources=data.frame(path=src,sha256=source_hashes),
  forest=plot_check,forest_inches=c(20,5),forest_family='Helvetica',forest_base_pt=12,
  forest_differences=c('Minimal CSV/RDS field and modality adapter','Actual raw-unit new-endpoint estimates',
    'Type 2 diabetes label reflects E11 endpoint','Events/C-index explicitly crude cohort',
    'Crude axis upper limit2.5 avoids clipping','Removed small provenance footnote; external caption provides details'),
  cif=list(pages=2,facets_per_page=8,policy='death_first',grid='Annual estimates, 0 to14 years; lines interpolate annual values',
    competing_event='Death before disease',pointwise_CI='Aalen-Johansen variance normal approximation, clipped0 to1',
    model='Unadjusted stratified by saved-score quartiles, no new model fits',tie_sensitivity_preserved_not_plotted=TRUE,
    minimum_displayed_risk_fraction=min(z$at_risk/z$N),max_precomputed_AJ_error=max(cif_qa$max_AJ_error))),
  file.path(out,'producer_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(c('Figures are transparent new-endpoint candidates, not restored historical forest11.',
  'Fig3c primary: expanded_BMI Cox HR per one saved DISCO or raw DM unit. DM is not log-transformed in this figure.',
  'Covariates: age, sex, BMI, smoking, drinking, education, race; missing race assigned White by author policy.',
  'Forest table Events and C-index are the crude score-only sample, which can exceed the adjusted Cox sample.',
  'Primary and crude plotted CSVs contain separate adjusted N/events and crude N/events.',
  'ClinicalProt, Prot and Clinical are respectively396 selected combined features,367 proteins and58 clinical markers.',
  'CIF pages1/2 show DISCO/DM Q1 andQ4; death before disease is competing event, solidCIF with95%CI, dashed1-KM.',
  'No age adjustment in the quartile curves. These plots do not prove age-independent or incremental prediction.',
  'Annual summary values are connected by straight lines; these are not full individual event-time step curves.',
  'Display ends14years with at least10% of each stratum still at risk. This is display support, not a changed administrative cutoff.',
  'Full disease-specific cohort, event and competing-death counts are inqa/ClinicalProt_CIF_counts.csv;',
  'deaths counts mean deaths observed before that disease, not all deaths of that cohort.',
  'Death-first/disease-first tie sensitivity values are preserved in the frozen80-rowQA but onlydeath-first is plotted.',
  'Historical provider/region coverage, loss-of-follow-up reasons and229 baseline-date anomalies remain unresolved endpoint evidence.',
  'No individual data downloaded, no scoring or cohort-model full rerun performed in this figure task.'),
  file.path(out,'FIGURE_CAPTIONS_AND_LIMITS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'qa/R_sessionInfo.txt'))
cat('Produced2near-source forest PDFs and1two-pageCIF PDF; aggregate inputs only.\n')
