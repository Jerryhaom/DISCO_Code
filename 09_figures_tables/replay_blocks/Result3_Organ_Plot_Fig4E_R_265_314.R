odis=readRDS("UKB.odis.rds")
cl2=which(is.na(odis$cause))
length(cl2)

cc=cor(odis[,1:36])
sel.t1=gsub(" ",".",sel.t)
cl=which(colnames(cc) %in% sel.t1)
res=matrix(0,ncol=length(cl),nrow = 7)
pl=res
k=1
for(i in cl){
  print(paste("########:",i))
  fl=paste("Surv(time,status)~",colnames(odis)[i],"+age+sex")
  c=coxph(as.formula(fl),data=odis)
  s=summary(c)
  cs=c(s$coefficients[1,2])
  p=c(s$coefficients[1,5])
  for(j in cause){
    cl1=which(odis$cause==j)
    d1=odis[c(cl1,cl2),]
    c=coxph(as.formula(fl),data=d1)
    s=summary(c)
    #print(s)
    #cs=c(cs,s$concordance[1])
    cs=c(cs,s$coefficients[1,2])
    p=c(p,s$coefficients[1,5])
  }
  res[,k]=cs
  pl[,k]=p
  k=k+1
}

cause
rownames(res)=c("All cause","Cancer","Heart disease","Cerebrovascular disease",
                "Respiratory disease","Alzheimer's disease","Others")
colnames(res)=colnames(odis)[cl]
pl1=pl
ft <- function(x){
  x1=cut(as.numeric(x),breaks=c(0,0.05,0.01,1),include.lowest = T)
  levels(x1)=c("**","*","")
  x1
}
pl=apply(pl,2, ft)
res1=res[c(1,3,5,2,4,6,7),]
pl1=pl[c(1,3,5,2,4,6,7),]
colnames(res1)=gsub("[.]"," ",colnames(res1))
r=readRDS("mortality/res.rds")
res2=res1[,match(rownames(r),colnames(res1))]
pl1=pl1[,match(rownames(r),colnames(res1))]

