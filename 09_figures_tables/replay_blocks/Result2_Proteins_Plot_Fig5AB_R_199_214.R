cc=cor(train.x[,bb],method = "pearson")
colnames(cc)=bb1
rownames(cc)=colnames(cc)

cc=abs(cc)
cc=reshape2::melt(cc)
cc=na.omit(cc)
cc=cc[cc$value>0.2,]
cc=cc[cc$value<1,]
dim(cc)

#################################
#ig$layout=layout.fruchterman.reingold(ig)
#plot(ig)
ig=graph_from_data_frame(cc)
V(ig)$Size=degree(ig)
