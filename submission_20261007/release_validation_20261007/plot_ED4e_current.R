#!/usr/bin/env Rscript
# Same identified Test_Distance.R ggplot block, from accepted aggregate cells.
entry <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
here <- dirname(normalizePath(entry,mustWork=TRUE))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE)
input <- if(length(args)>=1L)normalizePath(args[1],mustWork=TRUE) else file.path(here,'data','ED4e_matrix_distances_OFFICIAL_SS.csv')
out <- if(length(args)>=2L)args[2] else file.path(getwd(),'ED4e_current_output')
if(file.exists(out))stop('Choose a new output directory')
res <- read.csv(input,stringsAsFactors=FALSE)
stopifnot(nrow(res)==11L,all(c('ci','se','d','dt')%in%names(res)),sum(res$d=='Procrustes.SS')==1L,all(is.finite(res$ci)),all(is.finite(res$se)),all(res$se>0))
res$dt=factor(res$dt,levels = res$dt[c(1:3,5:10,4,11)])
p1=ggplot(res,aes(ci,dt)) +
  geom_bar(stat = 'identity',width =0.8,alpha=0.9,fill='grey') +
  geom_errorbar(aes(xmin=ci-se,xmax=ci+se),width=0.5) +
  coord_cartesian(xlim=c(0.6,0.78)) +
  theme_classic() + xlab('C-index') + ylab('Distance') +
  labs(title = 'NHANES') +
  theme(axis.text.y = element_text(size=10,color='black'),axis.title = element_text(size=12))
dir.create(out,recursive=TRUE)
pdf(file.path(out,'ED4e_OFFICIAL_SS_current.pdf'),height=3,width=5)
print(p1);dev.off()
write.csv(res,file.path(out,'ED4e_drawing_cells.csv'),row.names=FALSE)
writeLines(c('Actual original drawing: Result1/Robust/Test_Distance.R:146-173; device 175-177.',
  'SS reflects official CovTools 0.5.7 corrected comparator; old 0.5.4 SS is not recommended.',
  'Other nine comparator cells retained; Default follows locked table precision, not a refit.',
  'Bars are C-index +/- one model SE, not 95% CI. This public entry does not refit private NHANES.',
  'Source: https://cran.r-project.org/web/packages/CovTools/news/news.html'),file.path(out,'scope.txt'))
