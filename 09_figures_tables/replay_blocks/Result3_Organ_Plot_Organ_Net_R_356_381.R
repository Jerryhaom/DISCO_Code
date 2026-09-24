levels(dci1$di)=c("Obesity","Hypertension","Diabetes","Cancer",
                  "MI","Stroke","COPD","Dementia")

pdf("Degree_cindex.disease.pdf",height = 5,width = 8)
ggplot(dci1,aes(ci,de)) +
  geom_point(aes(fill=tissue),size=3,shape=21,
             color="black") +
  geom_smooth(method = "lm",color="navy") +
  stat_cor(cor.coef.name = "r",#method = "spearman",
           label.y.npc = 0.95) +
  facet_wrap(di~.,scales = "free_x",nrow = 2) +
  scale_fill_manual(values = col) +
  theme_few() + labs(fill="Organ and tissues") +
  guides(fill=guide_legend(nrow = 3)) +
  theme(#legend.position = "none",#
        legend.position = 'bottom',
        legend.title.position = "top",
        legend.key.height = unit(0.5,"lines"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"lines"),
        legend.box.margin = unit(c(-1,0.1,0.1,0.1),"lines"),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.spacing.x = unit(1,"lines")) +
  ylim(5,11) +
  xlab("C-index for diseases") + ylab("Degree centrailty") 
dev.off()
