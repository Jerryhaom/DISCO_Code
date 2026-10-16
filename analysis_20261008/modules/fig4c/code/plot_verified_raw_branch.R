#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(ggpubr)})
args<-commandArgs(TRUE);m<-if(length(args))normalizePath(args[1])else normalizePath('05_organ_network/Fig4c_raw_reconstruction_20261008');t<-read.csv(file.path(m,'tables/Fig4c_new_statistics.csv'));p<-data.frame(dg=t$degree,cc=t$partial_r,lb1=t$display_label)
pdf(file.path(m,'figures/Fig4c_TRANSPARENT_RAW_QC_LASSO_RECONSTRUCTION.pdf'),height=4,width=5)
print(ggplot(p,aes(dg,abs(cc)))+geom_point(color='grey',size=2)+geom_smooth(method='lm',color='black')+ylim(0,.65)+xlim(0,350)+xlab('Degree in Correlation Network')+ylab('Partial correlation with DISCO')+stat_cor(method='spearman',cor.coef.name='rho',label.y.npc=.95)+geom_text_repel(aes(label=lb1),size=3,max.overlaps=15,seed=134)+labs(title='Relationship between biomarker correlations with\n\n       DISCO and biomarker network connectivity')+theme_classic()+theme(axis.title=element_text(size=14),plot.title=element_text(size=14,lineheight=.5)))
dev.off()
