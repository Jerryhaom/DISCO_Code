cor.test(abs(pd1$cc),pd1$dg,method = "spearman")

tx=rownames(head(pd1[order(pd1$dg,decreasing = T),],10))
tx
rownames(head(pd1[order(pd1$dg,decreasing = T),],15))

pd1$lb1=rownames(pd1)
pd1$lb1[!pd1$lb1 %in% tx]=""

pdf("./Prot.de.pdf",height = 4,width = 5)
ggplot(pd1,aes(dg,abs(cc))) +
  geom_point(color="grey",size=2) +
  geom_smooth(method = "lm",color="black") +
  ylim(0,0.65) + xlim(0,350) +
  xlab("Degree in Correlation Network") +
  ylab("Partial correlation with DISCO") +
  ylab("Partial correlation with DISCO") +
  stat_cor(method = "spearman",
           cor.coef.name = c("rho"),
           label.y.npc = 0.95) + 
  geom_text_repel(aes(label=lb1),size=3,max.overlaps = 15) +
  labs(title = "Relationship between protein correlations with\n
       DISCO and protein network connectivity") +
  #geom_segment(aes(x = -0.2, y = -0.2, xend = 0.6, yend = 0.6),
  #             linetype="dashed")  +
  theme_classic() +
  theme(axis.title = element_text(size=14),
        plot.title = element_text(size=14,lineheight = 0.5))

dev.off()
