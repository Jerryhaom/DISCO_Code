prot.cla1=readRDS("UKB.Prot.class.rds")
#write.csv(prot.cla1,"UKB.Prot.class.csv",row.names = F)
tissues=gsub("[.]"," ",colnames(overlapped))
length(tissues)
prot.cla2=prot.cla1[prot.cla1$specific.tissue %in% tissues,]
unique(prot.cla2$specific.tissue)

overlapped2=matrix(0,ncol = length(cl),nrow = length(cl))

table(prot.cla2$specific.tissue)
clt=unique(prot.cla2$specific.tissue)
for(i in 1:length(clt)){
  for(j in 1:length(clt)){
    i1=prot.cla2$ENTREZID[prot.cla2$specific.tissue==clt[i]]
    j1=prot.cla2$ENTREZID[prot.cla2$specific.tissue==clt[j]]
    ol=i1[i1 %in% j1]
    overlapped2[i,j]=length(ol)
    #if(i==j){
    #  overlapped2[i,j]=0
    #}
  }
}

colnames(overlapped2)=clt
rownames(overlapped2)=colnames(overlapped2)
