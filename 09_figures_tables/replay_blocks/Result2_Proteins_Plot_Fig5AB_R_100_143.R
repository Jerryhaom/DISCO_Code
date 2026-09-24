pn=read_tsv("Proteins/coding143.tsv")
pnl=unlist(strsplit(pn$meaning,";"))
pn$prot=pnl[seq(1,length(pnl),2)]

bb1=bb
bb1=gsub("^X","",bb)
bb1[40:length(bb1)]=pn$prot[match(bb1[40:length(bb1)],pn$coding)]

names(cc)=bb1
pt=cc[abs(cc)>0.2]

cc=c()
pc=c()
pl1=data.frame(y2=pd$mc, x=pd$age,s=pd$sex)
l2=lm(y2~.,data=pl1)

for(i in 1:length(bb)){
  pl=data.frame(y=train.x[,bb[i]],x=pd$age,s=pd$sex)
  l1=lm(y~.,data=pl)
  #c=cor.test(pl$y,l2$residuals)
  c=cor.test(l1$residuals,l2$residuals)
  cc=c(cc,c$estimate)
  pc=c(pc,c$p.value)

}
names(cc)=bb1
cc[pc<0.05/length(bb)]
length(cc[abs(cc)>0.2])

pt=cc #[abs(cc)>0.2]
pt1=pt[names(pt) %in% pn$prot]

#write.table(pt1,"Proteins/Prot.DISCO.cor.txt",quote = F,sep="\t")

pd=data.frame(cc,cc1,pc)
pd1=pd #[rownames(pd) %in% pn$prot,]
pd1$lb=rownames(pd1)
dim(pd1)

tx=rownames(head(pd1[order(pd1$cc,decreasing = T),],10))
pd1$lb[!pd1$lb %in% tx]=""
library(ggplot2)
library(ggrepel)
pd1$p=p.adjust(pd1$pc,"fdr")
