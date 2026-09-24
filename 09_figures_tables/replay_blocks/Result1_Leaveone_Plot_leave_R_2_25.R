setwd(disco_path("Result1/Leaveone/"))
library(tidyverse)
library(ggthemes)

res=readRDS("UKB.clinical.leaveoneout.rds")
res=data.frame(res)
res=res[res$phe!="id",]

res1=res %>% group_by(phe) %>%
  summarise(m=median(DISCO))

n=res1$phe[order(res1$m)]
n=unique(c("All",n))
res$phe=factor(res$phe,levels=n)
library(tidyverse)

pdf("Leave_one_out2.pdf",height = 2.5,width = 8)
ggplot(res,aes(phe,DISCO)) +
  geom_boxplot() +
  geom_hline(yintercept = mean(res$DISCO[res$phe!="All"]),
             linetype="dashed",color="red",size=1) +
  ylab("DISCO C-index") + xlab("Biomarkers") +
  theme_classic() + labs(title = "Leave-one-out analysis of DISCO") +
  theme(axis.text.x = element_text(size=8,angle = 60,hjust = 1,color="black"))
