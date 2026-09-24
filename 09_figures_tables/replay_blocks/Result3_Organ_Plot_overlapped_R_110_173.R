overlapped31=overlapped3
overlapped31[upper.tri(overlapped3,diag = F)]=NA
pd=melt(overlapped31)
pd=na.omit(pd)
pd$Var1=factor(pd$Var1,levels = rev(or2))
pd$Var2=factor(pd$Var2,levels = or2)

p1=ggplot(pd,aes(Var2,Var1)) +
  geom_tile(aes(fill=value),color="black") +
  geom_text(aes(label=value)) +
  theme_classic() +
  scale_fill_distiller(palette = "Reds",
                       limits=c(0,256),
                       breaks=c(0,100,200),
                       direction = 1) +
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10)) +
  xlab("") + #ylab("Organ and tissues") +
  ylab("Before feature selection") +
  labs(fill="Counts") +
  theme(legend.direction = "horizontal",
        legend.key.height = unit(0.8,"lines"),
        legend.key.width  = unit(0.8,"lines"),
        legend.position = c(-0.08,1.1),
        legend.title.position = "top",
        plot.margin = unit(c(7,7,0.5,0.5),"lines"))


overlapped11=overlapped[match(or2,rownames(overlapped)),
                        match(or2,rownames(overlapped))]
overlapped41=overlapped11
overlapped41[lower.tri(overlapped41,diag = F)]=NA
pd=melt(overlapped41)
pd=na.omit(pd)
pd$Var1=factor(pd$Var1,levels = rev(or2))
pd$Var2=factor(pd$Var2,levels = or2)

p2=ggplot(pd,aes(Var2,Var1)) +
  geom_tile(aes(fill=value),color="black") +
  geom_text(aes(label=value)) +
  theme_classic() +
  scale_fill_distiller(palette = "Reds",
                       limits=c(0,256),
                       direction = 1) +
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 0),
        axis.text.y = element_text(size = 10),
        legend.position = "None") +
  xlab("After feature selection") + ylab("") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right")  +
  theme(panel.background = element_rect(fill = "transparent"),  # 面板透明‌:ml-citation{ref="3,5" data="citationList"}
    plot.background = element_rect(fill = "transparent",color="transparent"))

p2

pdf("Tissue.Overlapped.prot.pdf",height = 6,width = 7,
    onefile = F)

text="Shared proteins used to calculate organ/tissue DISCO scores"
g<-ggplotGrob(p2)
p1+annotation_custom(g,xmin=1.5,xmax=25,ymin=2,ymax=26)  +
  labs(title = text) +
  theme(plot.title = element_text(hjust = 0,vjust = 32))
dev.off()
