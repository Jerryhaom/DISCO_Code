suppressPackageStartupMessages(library(survival))
make_spec<-function(ref,weight_data,weight_age){
 ca<-apply(weight_data,2,function(v)cor(v,weight_age));w<-abs(outer(ca,ca));diag(w)<-0;w<-w/sum(w)
 cr<-cor(ref);ev<-eigen(cr,symmetric=TRUE,only.values=TRUE)$values
 stopifnot(all(is.finite(w)),min(ev)>1e-10)
 list(N=nrow(ref),p=ncol(ref),mu=colMeans(ref),scatter=crossprod(sweep(ref,2,colMeans(ref))),
 covariance=cov(ref),correlation=cr,weight=w,reference=ref,min_eigen=min(ev),condition=max(ev)/min(ev))
}
score<-function(x,sp,loo=rep(FALSE,nrow(x))){
 out<-cbind(DISCO=rep(NA_real_,nrow(x)),DM=log(mahalanobis(x,sp$mu,sp$covariance)))
 for(i in seq_len(nrow(x))){
  v<-x[i,];n<-sp$N;mu<-sp$mu;S<-sp$scatter;rr<-sp$correlation
  if(loo[i]){n<-n-1;mu<-(sp$N*sp$mu-v)/n
   S<-sp$scatter-(sp$N/n)*tcrossprod(v-sp$mu);rr<-cov2cor(S)
   out[i,'DM']<-log(mahalanobis(matrix(v,nrow=1),mu,S/(n-1))) }
  delta<-v-mu;grown<-S+(n/(n+1))*tcrossprod(delta)
  out[i,'DISCO']<-log(sum((rr-cov2cor(grown))^2*sp$weight)*n^2)
 }
 stopifnot(all(is.finite(out)));out
}
make_scores<-function(X,p,train_blood,mode='full',verify=FALSE){
 X18<-X[,names(p$oldmap)];colnames(X18)<-unname(p$oldmap)
 X35<-X[,p$features];X36<-X[,c(p$features,'HbA1c')]
 ref<-p$drawage<60
 if(mode!='full')ref<-ref&train_blood
 specs<-list(NHANES18=make_spec(p$nhref,as.matrix(p$nhanes[,-1]),p$nhanes$age),
  HRS18=make_spec(X18[ref,],X18[train_blood,],X[train_blood,'AGE']),
  HRS35=make_spec(X35[ref,],X35[train_blood,],X[train_blood,'AGE']),
  HRS36=make_spec(X36[ref,],X36[train_blood,],X[train_blood,'AGE']))
 mats<-list(NHANES18=X18,HRS18=X18,HRS35=X35,HRS36=X36)
 d<-p$evaluation;checks<-list()
 for(n in names(specs)){
  sp<-specs[[n]];xm<-mats[[n]];loo<-n!='NHANES18'&ref[p$eval_index]
  sc<-score(xm[p$eval_index,],sp,loo)
  for(m in colnames(sc))d[[paste0(n,'_',m)]]<-sc[,m]
  chk<-data.frame(scheme=n,ref_N=sp$N,features=sp$p,LOO_N=sum(loo),condition=sp$condition,min_eigen=sp$min_eigen)
  if(verify){
   probe<-unique(c(head(which(loo),3),head(which(!loo),3)))
   errd<-errm<-numeric()
   for(i in probe){rr<-sp$reference
    if(loo[i]){j<-match(p$eval_index[i],which(ref));rr<-rr[-j,,drop=FALSE]}
    v<-xm[p$eval_index[i],];nref<-nrow(rr)
    direct<-log(sum((cor(rr)-cor(rbind(rr,v)))^2*sp$weight)*nref^2)
    standardized<-drop((v-colMeans(rr))/apply(rr,2,sd))
    dm<-log(drop(t(standardized)%*%solve(cor(rr),standardized)))
    errd<-c(errd,abs(direct-sc[i,'DISCO']));errm<-c(errm,abs(dm-sc[i,'DM']))
   }
   chk$literal_DISCO_error<-max(errd);chk$independent_DM_error<-max(errm)
   stopifnot(max(errd)<1e-7,max(errm)<1e-6)
  }
  checks[[n]]<-chk
 }
 list(data=d,checks=do.call(rbind,checks),specs=specs)
}
metric_map<-function(p)c(p$clock_map,setNames(paste0(rep(c('NHANES18','HRS18','HRS35','HRS36'),each=2),c('_DISCO','_DM')),
 paste0(rep(c('NHANES18','HRS18','HRS35','HRS36'),each=2),c('_DISCO','_DM'))),Age='age')
fit_c<-function(d,field,model){
 v<-d;v$s<-v[[field]]
 if(model=='accel')v$s<-resid(lm(s~age,data=v))
 v$s<-as.numeric(scale(v$s))
 f<-switch(model,crude=Surv(time,status)~s,adj=Surv(time,status)~s+age,
  age_sex=Surv(time,status)~s+age+factor(SEX),accel=Surv(time,status)~s)
 if(field=='age')f<-Surv(time,status)~s
 fit<-coxph(f,data=v,ties='efron');v$lp<-predict(fit,type='lp')
 cc<-concordance(Surv(time,status)~lp,data=v,reverse=TRUE,cluster=v$household)
 sm<-summary(fit)
 c(C=unname(cc$concordance),Cvar=unname(cc$var),beta=coef(fit)[1],beta_var=vcov(fit)[1,1])
}
pool_C<-function(v){
 q<-qlogis(v$C);u<-v$Cvar/(v$C*(1-v$C))^2
 M<-length(q);B<-if(M>1)var(q)else 0;U<-mean(u);T<-U+(1+1/M)*B
 df<-if(B>0)(M-1)*(1+U/((1+1/M)*B))^2 else Inf
 z<-qt(.975,df);a<-mean(q)
 data.frame(C=plogis(a),lower95=plogis(a-z*sqrt(T)),upper95=plogis(a+z*sqrt(T)),
 MI_mean_C=mean(v$C),MI_min_C=min(v$C),MI_max_C=max(v$C),between_var_logit=B,
 within_var_logit=U,M=M,df=df)
}
