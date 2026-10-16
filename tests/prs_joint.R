source('approved/prs_joint.R')
library(survival)
set.seed(441)
d<-data.frame(time=rexp(300),status=rbinom(300,1,.7),mc=rnorm(300),prs=rnorm(300))
f0<-coxph(Surv(time,status)~mc+prs,d)
f1<-coxph(Surv(time,status)~mc*prs,d)
h<-disco_joint_from_f1(f1);b<-coef(f1);v<-vcov(f1)
stopifnot(is.na(disco_joint_from_f1()$HR),abs(h$beta-sum(b))<1e-12,
  abs(h$SE-sqrt(sum(v)))<1e-12,abs(h$HR-exp(sum(b)))<1e-12,
  abs(h$lower95-exp(sum(b)-qnorm(.975)*sqrt(sum(v))))<1e-12,
  abs(h$upper95-exp(sum(b)+qnorm(.975)*sqrt(sum(v))))<1e-12,
  inherits(try(disco_joint_from_f1(f0),silent=TRUE),'try-error'),
  inherits(try(disco_joint_from_f1(list(beta=coef(f0),interaction=coef(f1)[3])),silent=TRUE),'try-error'))
cat('8 same-model joint-contrast/covariance and mixed-source rejection checks passed; synthetic only\n')
