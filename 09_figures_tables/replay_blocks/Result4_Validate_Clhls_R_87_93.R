biomarker=read.csv("CLHLS//New/biomarker_dataset_CLHLS_2014-1.tab.csv")
biomarker=apply(biomarker,2,Hmisc::impute)
biomarker=data.frame(biomarker)

summary(biomarker$trueage)
ref=biomarker[which(biomarker$trueage<=70),4:ncol(biomarker)]
dim(ref)
