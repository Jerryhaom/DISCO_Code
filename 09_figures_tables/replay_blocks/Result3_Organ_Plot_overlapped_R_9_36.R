odis=readRDS("UKB.odis.rds")
cl2=which(is.na(odis$cause))
length(cl2)

cc=cor(odis[,1:36])

library(reshape2)
library(ggplot2)

cl=readRDS("Selected.tissue.rds")
names(sel[cl])=gsub("[.]"," ",colnames(cc)[cl])

overlapped=matrix(0,ncol = length(cl),nrow = length(cl))
#saveRDS(cl,"Selected.tissue.rds")

for(i in 1:length(cl)){
  for(j in 1:length(cl)){
    i1=cl[i]
    j1=cl[j]
    ol=sel[[i1]][sel[[i1]] %in% sel[[j1]]]
    overlapped[i,j]=length(ol)
    #if(i==j){
    #  overlapped[i,j]=0  }
  }
}

colnames(overlapped)=gsub("[.]"," ",colnames(cc)[cl])
rownames(overlapped)=colnames(overlapped)
