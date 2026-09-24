biomarker=read.table("./CHARLS/Blood_Biomarkers.txt",sep="\t",header = T)
biomarker=apply(biomarker,2,Hmisc::impute)
biomarker=data.frame(biomarker)

info=readRDS("./CHARLS/Charls.sample.info.rds")
info$age=2011-as.numeric(info$Birth)
info=info[match(biomarker$ID,info$ID),]
head(info)

biomarker=biomarker[which(info$age>=45),]
info=info[match(biomarker$ID,info$ID),]
dim(info)

ref=biomarker[which(info$age<=50),2:ncol(biomarker)]
dim(ref)
