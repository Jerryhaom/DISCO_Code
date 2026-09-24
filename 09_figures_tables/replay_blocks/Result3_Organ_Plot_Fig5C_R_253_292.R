prot1=readRDS("UKB.prot.ogran.rds")
res=readRDS("UKB.tissue.mc.rds")
odis=data.frame(res)
cc=cor(odis)

pcorr=matrix(0,ncol = 4,nrow = 400)
k=1
sel.t1=gsub(" ",".",sel.t)
cl=which(colnames(cc) %in% sel.t1)
pcorr1=matrix(0,ncol = length(cl),nrow = length(cl))

sel=readRDS("Lasso.select.prot.rds")
selp=sel

for(i1 in 1:length(cl)){
  print(paste("###########:",i1))
  for(j1 in 1:length(cl)){
    i=cl[i1]
    j=cl[j1]
    if(j>=i){  pcorr1[i1,j1]=1;    next }
    ol=sel[[i]][sel[[i]] %in% sel[[j]]]
    pd=data.frame(y=odis[,i],
                  x=odis[,j],
                  age=info1$Age,
                  sex=info1$Gender)
    if(length(ol)>0){
      pd=cbind(pd,prot1[,ol])
    }
    l1=lm(y~.,data=pd[,-2])
    l2=lm(x~.,data=pd[,-1])
    c=cor(l1$residuals,l2$residuals)
    p=cor.test(l1$residuals,l2$residuals)$p.value
    pcorr[k,]=c(colnames(odis)[i],
                colnames(odis)[j],c,p)
    pcorr1[i1,j1]=c
    pcorr1[j1,i1]=c
    k=k+1
  }
}
pcorr=pcorr[pcorr[,1]!=0,]
