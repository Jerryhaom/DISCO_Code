setwd(disco_path("Result1/Robust/"))
res=readRDS("../MC.phe2.rds")
library(reshape2)
res=melt(res,id="phe")
res$phe=as.numeric(res$phe)
res1=res %>% group_by(phe,variable) %>%
  summarise(m=median(value),se=sd(value)/(sqrt(length(value))))

p2=ggplot(res1,aes(phe,m,color=variable)) +
  geom_point() +
  geom_line(aes(group=variable)) +
  geom_errorbar(aes(ymin=m-se,ymax=m+se),width=0.3,color="black") +
  xlab("Biomarker numbers") + ylab("DISCO distribution\nLog-transformed") + 
  labs(color="Statistics",title = "NHANES") +
  scale_color_npg(labels=c("Min","Q1","Median","Mean","Q3","Max")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))

res=readRDS("../MC.phe3.rds")
library(reshape2)
res=melt(res,id="phe")
res$value=as.numeric(res$value)
#res$value=res$value*res$value
res$phe=as.numeric(res$phe)
res1=res %>% group_by(phe,variable) %>%
  summarise(m=median(value),se=sd(value)/(sqrt(length(value))))

p21=ggplot(res1,aes(phe,m)) +
  geom_point() +
  geom_line(aes(group=variable),color="navy") +
  geom_errorbar(aes(ymin=m-se,ymax=m+se),width=0.3,color="black") +
  xlab("Biomarker numbers") +  ylab("Correlation") + 
  labs(color="",title = "NHANES") +
  #scale_x_continuous(breaks = seq(10,90,20)) +
  #scale_color_npg(labels=c("Min","Q1","Median","Mean","Q3","Max")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))

#pdf("./Robust/DISCO.phe2.pdf",height = 3,width = 4)
#dev.off()
###########

res=readRDS("../MC.refsize3.rds")
library(reshape2)
res=melt(res,id="frq")
#res$frq=as.factor(res$frq*100)
res$frq=as.factor(res$frq) #*100
res1=res %>% group_by(frq,variable) %>%
  summarise(m=mean(value),se=sd(value)/(sqrt(length(value))))

p3=ggplot(res1,aes(frq,m,color=variable)) +
  geom_point() +
  geom_line(aes(group=variable)) +
  geom_errorbar(aes(ymin=m-se,ymax=m+se),width=0.3,color="black") +
  xlab("Reference sample size") + ylab("DISCO distribution\nLog-transformed") + 
  labs(color="Statistics",title = "NHANES") +
  #scale_x_continuous(breaks = seq(10,90,20)) +
  scale_color_npg(labels=c("Min","Q1","Median","Mean","Q3","Max")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))

res=readRDS("../MC.refsize4.rds")
library(reshape2)
res=melt(res,id="frq")
res$value=as.numeric(res$value)
#res$value=res$value*res$value
res$frq=as.factor(res$frq) #*100
res1=res %>% group_by(frq,variable) %>%
  summarise(m=mean(value),se=sd(value)/(sqrt(length(value))))

p4=ggplot(res1,aes(frq,m)) +
  geom_point() +
  geom_line(aes(group=variable),color="navy") +
  geom_errorbar(aes(ymin=m-se,ymax=m+se),width=0.3,color="black") +
  xlab("Reference sample size") + ylab("Correlation") + 
  labs(color="",title = "NHANES") +
  #scale_x_continuous(breaks = seq(10,90,20)) +
  #scale_color_npg(labels=c("Min","Q1","Median","Mean","Q3","Max")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))

#pdf("./Robust/DISCO.ref.pdf",height = 3,width = 5.5)
#dev.off()

library(ggpubr)
pdf("DISCO.robust.pdf",height = 5,width = 9)
p11=ggarrange(NULL,p4,NULL,p21,NULL,widths = c(0.05,1,0.12,1,0.25),
              nrow = 1,ncol = 5)
p22=ggarrange(p3,NULL,p2,widths = c(1,0.05,1),
              legend = "right",
              common.legend = T,nrow = 1,ncol = 3)
ggarrange(p22,p11,nrow = 2,ncol = 1,heights = c(0.7,0.4))

dev.off()
