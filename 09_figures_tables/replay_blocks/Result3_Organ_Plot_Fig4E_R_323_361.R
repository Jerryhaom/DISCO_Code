r1=t(res2)
p1=t(pl1)

r=readRDS("mortality/res.rds")
p=readRDS("mortality/pl.rds")

rr=rbind(melt(r1),melt(r))
rr$p=rbind(melt(p1),melt(p))[,3]
rr$tt=c(rep("Mortality",nrow(melt(r1))),
        rep("Chronic disease incidence",nrow(melt(r))))
rr$Var2=gsub("All cause","All-cause",rr$Var2)
rr$tt=factor(rr$tt,levels = c("Mortality","Chronic disease incidence"))

rr$Var1=factor(rr$Var1,levels = or2)
rr$Var2=factor(rr$Var2,levels = rev(unique(rr$Var2)))
rr$p[rr$p=="*"]=" *"
################

cairo_pdf("Organ.outcomes2.pdf",height = 6,width = 12)

ggplot(rr,aes(Var1,Var2)) +
  geom_tile(aes(fill=value),color="white") +
  geom_text(aes(label=sprintf("%0.2f",value)),
            size=3.5,hjust=0.7) +
  geom_text(aes(label=p),hjust= -0.8,vjust= 0.7) +
  #facet_wrap(tt~.,scales = "free_x") +
  facet_wrap(tt~., scales = "free_y",nrow=2) +
  scale_fill_distiller(palette = "Reds",direction = 1) +
  theme_few() + labs(fill="HR") + 
  xlab("Organ and tissues") + ylab("") +
  #xlab("Outcomes") +
  theme(panel.spacing.x = unit(0.5,"lines"),
        axis.text.x  = element_text(size=14,color="black",
                                    angle = 30,hjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        strip.text = element_text(size=16,color="black"),
        axis.title = element_text(size=14))

dev.off()
