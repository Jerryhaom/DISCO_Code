g2=fortify(ig)


cairo_pdf("Partial.net.pdf",height = 4,width = 6)
ggplot(data = g2, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(colour = "grey80") +
  geom_point(aes(colour = name,size=size)) +
  geom_nodetext(aes(label = name),size=5)+
  theme_blank() +
  scale_color_manual(values = col) +
  scale_size_continuous(range = c(3,10)) +
  guides(color=F) + xlim(-0.05,1.05) + labs(size="Degree")
dev.off()
