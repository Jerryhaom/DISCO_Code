library(ggpointdensity)
pd1=pd
pdf("UKB.DISCO.Age.pdf",height = 3,width = 3)
ggplot(pd1,aes(age,mc)) +
  geom_pointdensity(size=0.5) +
  scale_color_distiller(palette = "Spectral") +
  #geom_point(color="grey",size=0.5) +
  geom_smooth(method = "lm", color="black") +
  stat_cor(method = "spearman",
           cor.coef.name = c("rho"),
           label.y.npc = 0.95) +
  scale_y_continuous(limits = c(-6,13),breaks = c(-5,0,5,10)) +
  labs(title = "UKB") +
  theme_few() + xlab("Age (Years)") + ylab("DISCO") +
  theme(legend.position = "None")

dev.off()

pd2=pd1[which(pd1$mc<13 & pd1$mc > -6),]
pdf("UKB.DISCO.HD.pdf",height = 3,width = 3)
ggplot(pd2,aes(log(hd),mc)) +
  geom_pointdensity(size=0.5) +
  scale_color_distiller(palette = "Spectral") +
  #geom_point(color="grey",size=0.5) +
  geom_smooth(method = "lm", color="black") +
  stat_cor(method = "spearman",
           cor.coef.name = c("rho"),
           label.y.npc = 0.95) +
  scale_y_continuous(limits = c(-6,13),breaks = c(-5,0,5,10)) +
  labs(title = "UKB") +
  theme_few() + xlab("DM (Log-transformed)") + ylab("DISCO") +
  theme(legend.position = "None")
dev.off()

##########################################
setwd(disco_path("Result1/"))
pd=readRDS("UKB.Clinical.DISCO.rds")
dd3=pd
#############Survival plots

dd4=dd3
dd4$Group <- cut(dd4$mc,
                 breaks = quantile(dd4$mc,probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE,
                 labels = FALSE,
                 right = TRUE)
dd4=dd4[as.numeric(dd4$Group) %in% c(1,4),]
dd4$Group=as.numeric(dd4$Group)
dd4$Group[dd4$Group>1]=2

dd5=dd3
dd5$Group <- cut(dd5$hd,
                 breaks = quantile(dd5$hd,probs = seq(0, 1, 0.25), na.rm = TRUE),
                 include.lowest = TRUE,
                 labels = FALSE,
                 right = TRUE)
dd5=dd5[as.numeric(dd5$Group) %in% c(1,4),]
dd5$Group=as.numeric(dd5$Group)
dd5$Group[dd5$Group>1]=2
dd5$Group=dd5$Group+2
table(dd5$Group)
dd4=rbind(dd4,dd5)

fit1 <- survfit(Surv(age,age+time,status) ~ Group,  data = dd4)
fit1

pdf("UKB.DISCO2.all.pdf",height = 3,width = 4,onefile = FALSE)
ggsurvplot(fit1, # 创建的拟合对象
           data = dd4,
           conf.int = T, # 显示置信区间
           censor = F,
           pval = F, # 添加P值
           xlab = "Age (Years)",
           ylab = "Survival probability",
           legend.title = "",
           #legend.labs=paste("Q",seq(1,4),sep = ""),
           legend.labs=c(paste("DISCO",c("Low","High"),"Risk"),
                         paste("DM/Entropy",c("Low","High"),"Risk")),
           title="UKB",
           ylim=c(0.48,1),
           xlim=c(38,max(fit1$time)),
           pval.coord=c(0.1,0.85),
           legend=c(0.3,0.35),
           palette = "npg")
dev.off()
