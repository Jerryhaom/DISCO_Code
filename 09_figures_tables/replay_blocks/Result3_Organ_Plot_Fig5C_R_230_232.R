sel.t=table(prot.cla$specific.tissue)[table(prot.cla$specific.tissue)>40]
sel.t=names(sel.t)
sel.t=sel.t[!sel.t %in% c("placenta","ovary","testis")]
