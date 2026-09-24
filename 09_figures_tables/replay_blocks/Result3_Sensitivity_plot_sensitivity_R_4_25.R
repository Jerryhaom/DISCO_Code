cindex=readRDS("UKB.organ.random.cindex.rds")
cindex1=cindex
head(cindex1)
library(reshape2)
library(ggplot2)

cindex1=melt(cindex1)
or2=readRDS("../Organ/Tissue.order.rds")
cindex1$Var1=factor(cindex1$Var1,levels = rev(or2))

pdf("Organ.cindex.sensitivity.pdf",height = 4.5,width = 4)
ggplot(cindex1,aes(value,Var1)) +
  geom_boxplot(aes(fill=Var1)) +
  scale_fill_manual(values = col) +
  theme_bw() +
  theme(legend.position = "None",
        axis.text.y = element_text(size=10),color="black") +
  xlab("C-index") + ylab("Organ and tissue DISCO") +
  labs(title = "Sensitivity analysis of mortality prediction by\n
       proteins included in organ/tissue DISCO scores") +
  theme(plot.title = element_text(lineheight = 0.5,hjust= 1,size=12))
dev.off()
