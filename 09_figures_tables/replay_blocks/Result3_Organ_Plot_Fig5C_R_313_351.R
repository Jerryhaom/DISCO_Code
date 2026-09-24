pcorr2=pcorr1
or21=rev(or2)
pcorr2=pcorr2[match(or21,colnames(pcorr2)),
              match(or21,colnames(pcorr2))]
pcorr2[upper.tri(pcorr2,diag = T)]=NA
pcorr=melt(pcorr2)
pcorr=na.omit(pcorr)
#pcorr=pcorr[pcorr$value<1,]
pcorr=data.frame(pcorr)
head(pcorr)
colnames(pcorr)=c("X1","X2","X3")
pcorr$X3=as.numeric(pcorr$X3)
pcorr$X1=gsub("[.]"," ",pcorr$X1)
pcorr$X2=gsub("[.]"," ",pcorr$X2)
pcorr$X1=factor(pcorr$X1,levels = rev(or21))
pcorr$X2=factor(pcorr$X2,levels = or21)
pcorr3=pcorr
pcorr
summary(as.numeric(pcorr$X3))

pdf("Organ.pcorr.pdf",height = 4.5,width = 5.5)
ggplot(pcorr,aes(X1,X2)) +
  geom_tile(aes(fill=X3),color="black",alpha=0.9) +
  geom_text(aes(label=sprintf("%0.2f",X3)),size=2) +
  scale_fill_distiller(palette = "Reds",direction = 1) +
  theme_classic() + labs(fill="Partial\nCorr") +
  ylab("Organ and tissue DISCO") + xlab("") +
  labs(title = "Age- and sex- adjusted partial correlations\n
       between organ/tissue DISCO scores") +
  theme(legend.position = c(0.8,0.7),
        legend.key.size = unit(1,"lines"),
        legend.key.height = unit(0.8,"lines"),
        plot.title = element_text(size=16,lineheight = 0.5,hjust = 1),
        plot.margin = unit(c(0.1,0.3,0.1,0.1)*2,"lines"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.title = element_text(size=14,color = "black")) +
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1,color = "black"))

dev.off()
