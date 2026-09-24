cl=which(colnames(cc) %in% sel.t1)
length(cl)

prot.cla$sig=as.factor(prot.cla$sig)
count=matrix(0,ncol = 2,nrow = length(cl))
for(i in 1:length(cl)){
  i1=cl[i]
  l1=prot.cla$sig[prot.cla$Gene %in% sel[[i1]]]
  
  count[i,]=c(length(l1[l1==levels(l1)[1]]),
              length(l1[l1==levels(l1)[2]]))
}
count=data.frame(count)
colnames(count)=c("Enhanced","Enriched")

count$tissue=colnames(overlapped)
count$tissue=factor(count$tissue,
                    levels = count$tissue[order(count$Enhanced+count$Enriched)])

count=melt(count,id="tissue")

pdf("Prot.tissue.sel.pdf",height = 5,width = 4)
ggplot(count,aes(value,tissue)) +
  geom_bar(aes(fill=variable),stat="identity")  +
  labs(fill="") +
  scale_x_continuous(position = "top",expand = c(0.01,0.01)) +
  theme_classic() + theme(legend.position =c(0.7,0.1)) +
  ylab("organ and tissue") + xlab("Count") +
  theme(axis.text.y = element_text(size=14,color = "black"),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14))

dev.off()
