p3=ggplot(MC,aes(mc)) +
  geom_histogram(aes(y=..density..),color="grey",
                 bins = 50,fill="grey",alpha=0.8) +
  geom_density() + #xlim(0,30)
  theme_classic() + xlab("DISCO (Log-transformed)") 

p2=ggplot(MC[MC$raw<20,],aes(raw)) +
  geom_histogram(aes(y=..density..),color="grey",
                 bins = 50,fill="grey",alpha=0.8) +
  geom_density() + #xlim(0,30)
  theme_classic() + xlab("DISCO (Raw)") + labs(title = "NHANES") +
  annotation_custom(grob=ggplotGrob(p3),
                    ymin = 0.1, ymax = 0.6, xmin=6, xmax=20)


pdf("NHANES.DISCO.Distribution.pdf",height = 3,width = 3.5)
p2
dev.off()
