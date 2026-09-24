setwd(disco_path("Result1/Robust/"))
res=readRDS("../MC.refsel.rds")
library(reshape2)
res$ty=factor(res$ty,levels = c("DISCO","Accel"))

p1=ggplot(res[res$ty=="DISCO",],aes(C,th)) +
  geom_bar(stat = 'identity',width =0.8,alpha=0.9,
           fill="grey") +
  coord_cartesian(xlim=c(0.65,0.78)) +
  geom_errorbar(aes(xmin=C-se,xmax=C+se),
                width=0.3,color="black") +
  xlab("C-index") + 
  ylab("Reference age (years)") + 
  labs(title = "DISCO") +
  theme_classic() +
  theme(axis.text.y = element_text(size=10,color="black"),
        axis.title = element_text(size=12))

p2=ggplot(res[res$ty=="Accel",],aes(C,th)) +
  geom_bar(stat = 'identity',width =0.8,alpha=0.9,
           fill="grey") +
  coord_cartesian(xlim=c(0.5,0.57)) +
  geom_errorbar(aes(xmin=C-se,xmax=C+se),
                width=0.3,color="black") +
  xlab("C-index") + 
  ylab("Reference age (years)") + 
  labs(title = "Accel") +
  theme_classic() +
  theme(axis.text.y = element_text(size=10,color="black"),
        axis.title = element_text(size=12))

pdf("./Robust/DISCO.refsel.pdf",height = 3,width = 5.5)
ggarrange(p1,NULL,p2,ncol = 3,
          widths = c(1,0.05,1))

dev.off()
