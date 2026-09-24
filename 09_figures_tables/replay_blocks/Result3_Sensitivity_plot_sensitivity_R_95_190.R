pcor1=readRDS("Tissue.corr.sensitivity.rds")
pcor=pcor1
library(reshape2)
pcor$Var1=factor(pcor$Var1,levels = rev(or2))
pcor$Var2=factor(pcor$Var2,levels = rev(or2))
pcor$value[as.numeric(pcor$Var1)>=as.numeric(pcor$Var2)]=NA
pcor$Var2=factor(pcor$Var2,levels = or2)


pdf("Organ.pcorr.sensitivity.pdf",height = 4,width = 5)

ggplot(pcor[!is.na(pcor$value),],aes(Var2,Var1)) +
  geom_tile(aes(fill=value),color="black",alpha=0.9) +
  geom_text(aes(label=sprintf("%0.2f",value)),size=2) +
  scale_fill_distiller(palette = "Reds",direction = 1) +
  theme_classic() + labs(fill="Partial\nCorr") +
  ylab("Organ and tissue DISCO") + xlab("") +
  labs(title = "Sensitivity analysis of organ/tissue DISCO correlations") +
  theme(legend.position = c(0.8,0.7),
        legend.key.size = unit(1,"lines"),
        legend.key.height = unit(0.8,"lines"),
        plot.title = element_text(size=12,hjust = 1),
        axis.text.y = element_text(size=10,color = "black"),
        axis.title = element_text(size=12,color = "black")) +
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1,color = "black"))

dev.off()

##################################
all.pcorr=readRDS("../Organ/Tissue.pcorr.noverlap.rds")
all.pcorr=data.frame(all.pcorr)
all.pcorr$X1=gsub("[.]"," ",all.pcorr$X1)
all.pcorr$X2=gsub("[.]"," ",all.pcorr$X2)
pcor1$raw=all.pcorr$X3[match(paste(pcor1$Var1,pcor1$Var2),
                             paste(all.pcorr$X1,all.pcorr$X2))]
#pcor1$raw[is.na(pcor1$raw)]=all.pcorr$X3[match(paste(pcor1$Var1,pcor1$Var2),
#                             paste(all.pcorr$X2,all.pcorr$X1))][is.na(pcor1$raw)]

pcor1$raw=as.numeric(pcor1$raw)
summary(pcor1$raw)
summary(pcor1$value)

pdf("Organ.pcorr.cmp.pdf",height = 3,width = 3.5)

ggplot(pcor1[!is.na(pcor1$raw),],aes(value,raw)) +
  geom_point(size=2,color="grey") +
  geom_smooth(method = "lm",color="black") +
  theme_classic() + ylab("Parital correlation based on\nspecific protein panel") +
  xlab("Parital correlation based on\nrandomly selected protein panel") +
  ylim(0.25,0.7) +
  labs(title = "Sensitivity of organ/tissue DISCO\n
                correlation to protein composition") +
  stat_cor(#method = "spearman",
           cor.coef.name = "r",
           label.y = 0.7*0.95) +
  theme(plot.title = element_text(size=12,hjust = 1,lineheight = 0.5))

dev.off()

#############################################
pcor=readRDS("Tissue.corr.sensitivity.rds")
pcorr=dcast(pcor,Var1~Var2,value.var = "value")
rownames(pcorr)=pcorr[,1]
pcorr=pcorr[,2:ncol(pcorr)]
deg=apply(pcorr, 1, sum)-1
sort(deg,decreasing = T)

cindex2=cindex1 %>% dplyr::group_by(Var1) %>%
  dplyr::summarise(m=mean(value))
cindex2

cindex2$de=deg[match(cindex2$Var1,names(deg))]
print(cor.test(cindex2$de,cindex2$m,method = "spearman"))
cindex2$Var1=factor(cindex2$Var1,levels = names(col))

pdf("Organ.cindex.degree.pdf",height = 3.5,width = 4)

ggplot(cindex2,aes(m,de)) +
  geom_point(aes(fill=Var1),size=3,shape=21,
             color="black") +
  geom_smooth(method = "lm",color="navy",se=F) +
  stat_cor(cor.coef.name = "r",#method = "spearman",
           label.y.npc = 0.95) +
  scale_fill_manual(values = col) +
  theme_classic() +
  ylim(6,11) +
  #labs() +
  theme(legend.position = 'None') +
  xlab("Mean C-index for mortality \n(200 replicates)") + 
  ylab("Degree centrailty")  +
  labs(title = "Sensitivity analysis of degree centrality-mortality\n
       associations of organ/tissue DISCOs") +
  theme(plot.title = element_text(size=12,hjust =0.5,lineheight = 0.5),
        axis.title = element_text(size=12))

dev.off()
