library(survminer)
pd.vd1=DNAm1
pd.vd1$group1=cut(pd.vd1$GrimAgeMort,
                  breaks=quantile(pd.vd1$GrimAgeMort,probs =seq(0,1,0.25)),
                  include.lowest = T,right = T)

pd.vd1$group1=as.numeric(pd.vd1$group1)
pd.vd1=pd.vd1[pd.vd1$group1 %in% c(1,4),]
pd.vd1$group1[pd.vd1$group1>1]=2
pd.vd1$group1=pd.vd1$group1+2

pd.vd2=DNAm1
pd.vd2$group1=cut(pd.vd2$HorvathAge,
                  breaks=quantile(pd.vd2$HorvathAge,probs =seq(0,1,0.25)),
                  include.lowest = T,right = T)
table(pd.vd2$group1)

pd.vd2$group1=as.numeric(pd.vd2$group1)
pd.vd2=pd.vd2[pd.vd2$group1 %in% c(1,4),]
pd.vd2$group1[pd.vd2$group1>1]=2
pd.vd2$group1=pd.vd2$group1+4

pd.vd3=DNAm1
pd.vd3$group1=cut(pd.vd3$mc,
                  breaks=quantile(pd.vd3$mc,probs =seq(0,1,0.25)),
                  include.lowest = T,right = T)
table(pd.vd3$group1)

pd.vd3$group1=as.numeric(pd.vd3$group1)
pd.vd3=pd.vd3[pd.vd3$group1 %in% c(1,4),]
pd.vd3$group1[pd.vd3$group1>1]=2

pd.vd1=rbind(pd.vd1,pd.vd2,pd.vd3)

pd.vd1$group1=as.factor(pd.vd1$group1)
levels(pd.vd1$group1)=c(paste("DISCO",c("Low","High"),"Risk"),
                        paste("GrimAge",c("Low","High"),"Risk"),
                        paste("HorvathAge",c("Low","High"),"Risk"))

pd.vd1$group2=gsub(" Low Risk","",pd.vd1$group1)
pd.vd1$group2=gsub(" High Risk","",pd.vd1$group2)
table(pd.vd1$group2)
stat1=pd.vd1 %>% dplyr::group_by(group2,group1) %>%
  dplyr::summarise(m=sprintf("%0.1f",mean(age)),
            s=sprintf("%0.1f",sd(age)),
            m1=min(age),m2=max(age))


library(ggpp)
stat1$yy=rep(c(0.95,0.8),3)
pdf("Age.distribution.pdf",height = 5,width = 3)
ggplot() +
  geom_density(data=pd.vd1,aes(age,color=group1),size=1.2) +
  geom_text_npc(data=stat1,aes(npcx=0.65,npcy=yy,
                               color=group1,
                               label=paste(m,s,sep="±")),
                hjust=0) +
  facet_wrap(group2~.,ncol = 1) +
  scale_fill_npg() +
  scale_color_npg() + xlab("Age (Years)") +
  theme_few() +
  theme(legend.position = "None",
        strip.text = element_text(size=16),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=14),
        strip.background = element_blank())
dev.off()
