res=readRDS("../UKB.clinical.sim2.rds")
library(reshape2)
library(ggthemes)
res1=melt(res,id="N")

pdf("./Sim.clinical.pdf",height = 3,width = 4)
ggplot(res1,aes(N,value)) +
  geom_point(aes(color=variable),size=1,
             position = position_dodge(2)) +
  geom_smooth(aes(group=variable),color="black") +
  ylab("C-index for mortality") +
  xlab("Number of clinical biomarkers") + 
  labs(color="",title = "UKB") + theme_bw() +
  scale_x_continuous(labels = seq(5,55,10),
                     breaks = seq(5,55,10)) +
  theme(axis.title = element_text(size=12))
dev.off()  

#res=readRDS("Num/UKB.clinicalprot.sim.rds")
#colnames(res)=c("DISCO","MD")
#res=data.frame(res)
#head(res)
#res$N=rep(seq(20,420,20),each=50)
#saveRDS(res,"UKB.clinicalprot.sim2.rds")

res=readRDS("UKB.clinicalprot.sim2.rds")
library(ggthemes)
res1=melt(res,id="N")

pdf("./Sim.clinicalprot.pdf",height = 3,width = 5)
ggplot(res1,aes(N,value)) +
  geom_point(aes(color=variable),size=1,
             position = position_dodge(2)) +
  geom_smooth(aes(group=variable),color="black") +
  ylab("C-index for mortality") +
  xlab("Number of clinical and proteomics biomarkers") + 
  labs(color="",title = "UKB") + theme_bw() +
  scale_x_continuous(labels = seq(20,420,80),
                     breaks = seq(20,420,80)) +
  theme(axis.title = element_text(size=12))
dev.off()  

