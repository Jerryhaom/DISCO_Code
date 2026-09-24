de=centralization.degree(ig)$res
names(de)=names(V(ig))
pd1$dg=de[match(rownames(pd1),names(de))]
