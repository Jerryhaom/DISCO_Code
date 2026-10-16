# Regression tests for model metadata and fail-closed result validation.
source('project_config.R')
suppressPackageStartupMessages(library(survival))
x <- data.frame(time=c(1,2,3,4,5,6,7,8),status=c(1,0,1,0,1,0,1,0),score=c(.4,.1,.8,.2,.3,.9,.5,.7))
fit <- coxph(Surv(time,status)~score,data=x)
ci <- summary(fit)$concordance
r <- data.frame(n=fit$n,events=fit$nevent,c_index=unname(ci[1]),se=unname(ci[2]))
stopifnot(r$n==8,r$events==4,nobs(fit)==4)
disco_validate_cox_rows(r)
reject <- function(v) inherits(try(disco_validate_cox_rows(v),silent=TRUE),'try-error')
stopifnot(reject(NULL),reject(data.frame()),reject(r[FALSE,]),
          reject(transform(r,n=NA_real_)),reject(transform(r,n=3)),
          reject(transform(r,events=0)),reject(transform(r,c_index=1.1)),
          reject(transform(r,se=-.1)))
cat('9 model metadata / rejection checks passed\n')
