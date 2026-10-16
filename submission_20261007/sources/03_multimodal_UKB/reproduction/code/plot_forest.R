#!/usr/bin/env Rscript
# File-system roots are explicit; no participant data are bundled.
.disco_path <- function(kind, suffix="") {
  key <- switch(kind, private="DISCO_PRIVATE_ROOT", results="DISCO_RESULTS_ROOT", legacy="DISCO_LEGACY_ROOT")
  root <- Sys.getenv(key, unset="")
  if (!nzchar(root)) stop(paste("Required environment variable:", key))
  if (nzchar(suffix)) file.path(root, suffix) else root
}
# Newly authored Fig3c forest design from near-source forestploter code and
# inspected final layout. Re-estimated coefficients only; never original PDF art.
options(device=function(...)grDevices::pdf(file=NULL,...),warn=1)
suppressPackageStartupMessages({library(forestploter);library(grid);library(gtable);library(jsonlite)})
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
input<-normalizePath(args[1],mustWork=TRUE);out<-args[2];selected<-args[3]
allowed<-.disco_path("results","09_figures_tables/outputs/fig3c_reestimated/")
stopifnot(startsWith(out,allowed),!grepl("(^|/)\\.\\.(/|$)",out),!file.exists(out),
          normalizePath(dirname(out),mustWork=TRUE)==dirname(out),
          selected%in%c("M2_corrected","M2_legacy05","M2_legacyR01"))
dir.create(out);stopifnot(normalizePath(out)==out)
x<-readRDS(input);r<-x$results
r<-r[r$time_rule=="nonnegative",]
main<-r[r$model==selected,];crude<-r[r$model=="M1_crude",]
stopifnot(nrow(main)==48L,nrow(crude)==48L)
key<-function(d)paste(d$disease,d$modality,d$metric,sep="|")
j<-match(key(main),key(crude));stopifnot(!anyNA(j))
main$C_crude<-crude$C_index[j];main$N_crude<-crude$N[j];main$events_crude<-crude$events[j]
main$percent_crude<-crude$event_percent[j]
write.csv(main,file.path(out,"plotted_results.csv"),row.names=FALSE,na="")
mods<-c("ClinicalProt","Prot","Clinical")
colors<-c("#00A087FF","#3C5488FF")
half<-function(diseases,titles){
  ds<-list();ests<-lows<-ups<-list(DISCO=numeric(),DM=numeric())
  for(i in seq_along(diseases)){
    ds[[length(ds)+1L]]<-c(titles[i],rep("",6))
    for(m in c("DISCO","DM")){
      ests[[m]]<-c(ests[[m]],NA);lows[[m]]<-c(lows[[m]],NA);ups[[m]]<-c(ups[[m]],NA)
    }
    for(mod in mods){
      row<-main[main$disease==diseases[i]&main$modality==mod,]
      a<-row[row$metric=="DISCO",];b<-row[row$metric=="DM",]
      stopifnot(nrow(a)==1L,nrow(b)==1L,a$N_crude==b$N_crude,a$events_crude==b$events_crude)
      ds[[length(ds)+1L]]<-c(paste0("    ",mod),sprintf("%d (%.2f%%)",a$events_crude,a$percent_crude),
        sprintf("%.2f",a$HR),sprintf("%.3f",a$C_crude),sprintf("%.2f",b$HR),sprintf("%.3f",b$C_crude),strrep(" ",27))
      for(m in c("DISCO","DM")){
        t<-row[row$metric==m,]
        ests[[m]]<-c(ests[[m]],t$HR);lows[[m]]<-c(lows[[m]],t$lower_ci);ups[[m]]<-c(ups[[m]],t$upper_ci)
      }
    }
  }
  tab<-as.data.frame(do.call(rbind,ds),stringsAsFactors=FALSE)
  names(tab)<-c("Diseases/Metrics","Events (%)\nM1 cohort","DISCO\nHR","DISCO\nC-index","DM\nHR","DM\nC-index","")
  stopifnot(all(unlist(lows)>0,na.rm=TRUE))
  # Preserve final axis if every newly estimated confidence interval fits.
  range_ci<-range(c(unlist(lows),unlist(ups)),na.rm=TRUE)
  stopifnot(range_ci[1]>=.9,range_ci[2]<=2.3)
  tm<-forest_theme(base_size=12,base_family="Helvetica",ci_pch=15,ci_fill=colors,ci_col=colors,
    ci_lwd=1.5,ci_lty=1,legend_position="none",legend_value=c("DISCO","DM"),
    refline_gp=gpar(lwd=1,lty="dashed",col="grey20"),
    core=list(bg_params=list(fill="white"),padding=unit(c(2,3),"mm")),
    colhead=list(fg_params=list(fontface="plain"),padding=unit(c(2,3),"mm")))
  p<-forest(tab,est=ests,lower=lows,upper=ups,sizes=.4,ci_column=7,ref_line=1,
    x_trans="log",xlim=c(.9,2.3),ticks_at=c(.9,1.3,1.7,2.3),ticks_digits=1,theme=tm)
  p<-add_border(p,row=1:16,col=1:6,where="top",gp=gpar(lwd=.5))
  p<-add_border(p,part="header",col=1:6,where="bottom",gp=gpar(lwd=1.5))
  # Explicit page geometry; dual-series forest is the seventh column.
  # forestploter adds one padding row/column on each side of the 7-column table.
  stopifnot(length(p$widths)==9L)
  p$widths[2:8]<-unit(c(1.60,1.60,.65,.79,.57,.79,2.12),"in")
  p$heights[2]<-unit(.43,"in")
  p$heights[3:18]<-unit(.218,"in")
  legend<-grobTree(
    segmentsGrob(x0=.10,x1=.18,y0=.75,y1=.75,gp=gpar(col=colors[1],lwd=1.5)),
    pointsGrob(x=.14,y=.75,pch=15,size=unit(7,"pt"),gp=gpar(col=colors[1])),
    textGrob("DISCO",x=.20,y=.75,just="left",gp=gpar(fontfamily="Helvetica",fontsize=12)),
    segmentsGrob(x0=.10,x1=.18,y0=.30,y1=.30,gp=gpar(col=colors[2],lwd=1.5)),
    pointsGrob(x=.14,y=.30,pch=15,size=unit(7,"pt"),gp=gpar(col=colors[2])),
    textGrob("DM",x=.20,y=.30,just="left",gp=gpar(fontfamily="Helvetica",fontsize=12)))
  p<-gtable_add_grob(p,legend,t=2,l=8,clip="off",name="two_metric_legend")
  list(plot=p,table=tab,est=ests,lower=lows,upper=ups)
}
left<-half(c("obesity","hypertension","diabetes","cancer"),c("Obesity","Hypertension","Diabetes","Cancer"))
right<-half(c("mi","stroke","copd","dementia"),c("MI","Stroke","COPD","Dementia"))
pdf(file.path(out,"Fig3c_REESTIMATED.pdf"),width=20,height=5,family="Helvetica",
    title=paste("REESTIMATED_FROM_UKB | Fig3c |",selected,"| final historical source not recovered"))
grid.newpage()
pushViewport(viewport(x=.255,y=.51,width=.47,height=.92));grid.draw(left$plot);popViewport()
pushViewport(viewport(x=.745,y=.51,width=.47,height=.92));grid.draw(right$plot);popViewport()
coding<-switch(selected,M2_corrected="valid no-qualification category retained; unknown ethnicity missing",
  M2_legacy05="legacy education missing; unknown ethnicity separate category",
  M2_legacyR01="legacy education missing; unknown ethnicity assigned White")
grid.text(paste0("RE-ESTIMATED | HR: expanded adjusted M2, per saved score unit; 95% CI. C-index and Events: unadjusted M1 cohort. ",coding,"."),
  x=.025,y=.014,just=c("left","bottom"),gp=gpar(fontsize=8,fontfamily="Helvetica",col="grey25"))
dev.off()
saveRDS(list(left=left[c("table","est","lower","upper")],right=right[c("table","est","lower","upper")],plotted=main),
  file.path(out,"plot_objects.rds"))
write_json(list(status="REESTIMATED_FROM_UKB_NOT_RECOVERED_SOURCE",selected_model=selected,
  input=input,source_md5=unname(tools::md5sum(input)),figure_inches=c(20,5),font="Helvetica",
  graphics_source="New forestploter composition; near-source forest.R table/border/log-axis design; final inspected two-series layout",
  differences=c("New estimates and complete-case M2 cohorts","Explicit Events M1 cohort header", "New provenance/coding footnote", "Explicit current page geometry"),
  original_source_recovered=FALSE,upstream_endpoint_construction_verified=FALSE,exact_original_pixels_claimed=FALSE),
  file.path(out,"figure_evidence.json"),auto_unbox=TRUE,pretty=TRUE)
cat("COMPLETE new Fig3c forest",selected,"\n")
