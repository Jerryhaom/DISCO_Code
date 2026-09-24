setwd(disco_path("Result2/"))
pd=readRDS("UKB.Prot.DISCO.rds")
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

fit1 <- survfit(Surv(age,time+age,status) ~ Group,  data = dd4)
fit1

pdf("UKB.Prot.DISCO.all.pdf",height = 3,width = 4,onefile = FALSE)
p=ggsurvplot(fit1, # 创建的拟合对象
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
           title="UKB, Proteomics",
           ylim=c(0.3,1),
           xlim=c(40,max(fit1$time)),
           pval.coord=c(0.1,0.8),
           legend=c(0.32,0.35),
           palette = "npg")
p1=p$plot+theme(legend.background = element_blank(),
                legend.text = element_text(size=12),
                axis.text = element_text(size=10),
                axis.title.y  = element_text(size=14),
                axis.title.x  = element_text(size=14),
                plot.title = element_text(size=14)) 
p1
dev.off()
