prot.cla=readxl::read_xlsx("ukb_simple_specific_1553.xlsx")
table(prot.cla$sig)
table(prot.cla$specific.tissue)

prot.cla$specific.tissue=gsub(" 1","",prot.cla$specific.tissue)
prot.cla$specific.tissue=factor(prot.cla$specific.tissue,
                                levels = names(sort(table(prot.cla$specific.tissue))))
