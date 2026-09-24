prot.cla=readxl::read_xlsx("ukb_simple_specific_1553.xlsx")
table(prot.cla$sig)
table(prot.cla$specific.tissue)

prot.cla$specific.tissue=gsub(" 1","",prot.cla$specific.tissue)
prot.cla$specific.tissue=factor(prot.cla$specific.tissue,
                                levels = names(sort(table(prot.cla$specific.tissue))))

pdf("Prot.tissue.pdf",height = 5,width = 4)
ggplot(prot.cla,aes(specific.tissue)) +
  geom_bar(stat="count",aes(fill=sig)) +
  coord_flip() +
  scale_fill_discrete(labels=c("Enhanced","Enriched")) +
  labs(fill="") +
  scale_y_continuous(position = "right",expand = c(0.01,0.01)) +
  theme_classic() + theme(legend.position =c(0.7,0.1)) +
  xlab("Oragn and tissues") + ylab("Counts")
