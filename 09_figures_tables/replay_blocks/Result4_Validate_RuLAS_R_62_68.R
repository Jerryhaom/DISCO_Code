library("FactoMineR")
library("factoextra")
res.pca <- PCA(biomarker[,bb], graph = FALSE,ncp = length(bb))
ind <- get_pca_ind(res.pca)
ind

biomarker1=data.frame(sampleID=biomarker$X2016ID,ind$coord[,1:10])
