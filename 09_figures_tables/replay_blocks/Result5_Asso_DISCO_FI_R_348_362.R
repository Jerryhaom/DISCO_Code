pdf("CHARLS.DISCO.Frailty2.pdf",height = 3,width = 3)
ggplot(pd.vd[!is.na(pd.vd$fl),],aes(fl,mc)) +
  geom_violin(aes(fill=fl)) +
  geom_boxplot(width=0.2,fill="white") +
  labs(title = "CHARLS, 2015 (Wave 3)") +
  stat_compare_means(comparisons = list(c("Robust","Prefrail"),
                                        c("Prefrail","Frail"),
                                        c("Robust","Frail"))) +
  theme_few() + xlab("") + ylab("DISCO") +
  scale_fill_manual(values = rev(pal_npg("nrc")(3))) +
  theme(axis.text.x = element_text(size=12,color="black"),
        legend.position = "None") +
  scale_y_continuous(expand = c(0.1,0.1))

dev.off()
