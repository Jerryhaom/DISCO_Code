library(readr)
info=read_csv("../UKB_Data/UKB.Mortality.csv")

clinical=readRDS("../Result1/UKB.Clinical.DISCO.rds")
prot=readRDS("UKB.Prot.DISCO.rds")
met=readRDS("UKB.Met.DISCO.rds")

common=clinical$ID[clinical$ID %in% prot$ID]
length(common)

common1=common[common %in% met$ID]
length(common1)

clinical.data=readRDS("../Result1/Clean.UKB.Clinical.rds")
clinical.data=clinical.data[match(common,clinical$ID),]
prot.data=readRDS("Selected.UKB.Prot.rds")
prot.data=prot.data[match(common,prot.data$eid),]
clinical.prot=data.frame(clinical.data,
                         prot.data[,2:ncol(prot.data)])

info1=info[match(common,info$id),]
############ feature selection
train.x <- as.matrix(clinical.prot[,2:ncol(clinical.prot)])
train.y <- Surv(info1$time,info1$status)
dim(train.x)

pd=readRDS("UKB.ClinicalProt.DISCO.rds")

fit_cv=readRDS("Lasso.model.ClinicalProt.rds")
coef=coef(fit_cv, s = "lambda.min")
#coef=coef(fit_cv, s = "lambda.1se")
coef
bb=rownames(coef)[coef[,1]!=0]
length(bb)
bb

cc=apply(train.x[,bb],2,function(x)
  cor(x,pd$mc,method = "spearman"))
pl=apply(train.x[,bb],2,function(x)
  cor.test(x,pd$mc,method = "spearman")$p.value)

cc1=apply(train.x[,bb],2,function(x)
  cor(x,pd$age))

