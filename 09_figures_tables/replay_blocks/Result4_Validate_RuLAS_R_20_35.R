info=read.csv("Rugao/Rugao.2016.Mortality.csv")
phe=read.csv("Rugao/2016Rugao.data.x.csv")
func=read.csv("Rugao//2016Rugao.data.y1.csv")
phe=phe[match(info$`X2016编号`,phe$X2016ID),]
func=func[match(info$`X2016编号`,func$X2016ID),]
dim(info)
dim(phe)

biomarker=phe
biomarker=apply(biomarker,2,Hmisc::impute)
biomarker=data.frame(biomarker)


#bb=colnames(biomarker)[5:ncol(biomarker)]
bb=colnames(biomarker)[c(5:13,20:67)]
length(bb)
