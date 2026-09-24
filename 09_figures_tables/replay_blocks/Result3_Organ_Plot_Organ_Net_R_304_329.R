cc=pcorr1
cs1=readRDS("Tissue.cs1.rds")
deg=apply(cc,1,sum)-1
cs1$de=deg[match(cs1$tissue,names(deg))]
table(cs1$variable)
cs2=cs1[cs1$variable=="C-index",]
or2=cs2$tissue[order(cs2$de,cs2$value,decreasing = T)]
col=readRDS("Tissue.col.rds")
col2=col[match(or2,names(col))]


#cs1$tissue=factor(cs1$tissue,levels = or2)
dci=read_xlsx("Disease_cindex.xlsx")
dcit=na.omit(unique(dci[[1]]))
dcic=na.omit(dci[[8]])
dcic=as.numeric(dcic[2:length(dcic)])
dcit=dcit[2:length(dcit)]
length(dcit)
length(dcic)
dci1=data.frame(tissue=rep(dcit,each=8),
                ci=dcic,di=rep(dci[[2]][3:10],18))
dci1$di=factor(dci1$di,levels = dci[[2]][3:10])
dci1$de=deg[match(dci1$tissue,names(deg))]
dci1$tissue=factor(dci1$tissue,levels = or2)


