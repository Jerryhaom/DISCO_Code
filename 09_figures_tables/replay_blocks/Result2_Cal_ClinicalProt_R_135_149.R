pdf("UKB.Entropy.HD.pdf",height = 3,width =3)
ggplot(pd,aes(log(hd),en)) +
  geom_pointdensity(size=0.5) +
  scale_color_distiller(palette = "Spectral") +
  #geom_smooth(method = "lm",color="black") +
  #geom_point(color="grey",size=0.5) +
  #scale_y_continuous(limits = c(-0.01,0.1)) +
  labs(title = "UKB, Proteomics+Clinical") +
  stat_cor(method = "spearman",
           cor.coef.name = c("rho"),
           label.y.npc=0.95) +
  theme_few() + xlab("DM (Log-transformed)") +
  ylab("Entropy distance") +
  theme(legend.position = "None")
dev.off()
