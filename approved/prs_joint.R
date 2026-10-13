# A joint contrast must use all coefficients and covariance from one interaction fit.
disco_joint_from_f1 <- function(fit=NULL,terms=c('mc','prs','mc:prs')) {
  empty <- data.frame(HR=NA_real_,lower95=NA_real_,upper95=NA_real_,beta=NA_real_,SE=NA_real_,
    status='MISSING_SAME_F1_BETA_AND_VCOV',
    contrast='DISCO=1, PRS=1 versus DISCO=0, PRS=0 in the fitted score scales',
    note='Do not multiply f0 main-effect HRs by an f1 interaction HR',stringsAsFactors=FALSE)
  if(is.null(fit))return(empty)
  if(!inherits(fit,'coxph') || length(terms)!=3L || anyDuplicated(terms) || !grepl(':',terms[3],fixed=TRUE))
    stop('One fitted Cox interaction model and its two main terms plus interaction required')
  b <- coef(fit); V <- vcov(fit)
  if(!all(terms%in%names(b)) || !all(terms%in%rownames(V)) || !identical(rownames(V),colnames(V)))
    stop('All three beta terms and full same-model covariance required')
  v <- V[terms,terms,drop=FALSE]; a <- b[terms]
  if(any(!is.finite(c(a,v))))stop('Non-finite interaction fit; joint contrast cannot be reported')
  variance <- as.numeric(sum(v)); if(variance<=0)stop('Nonpositive joint-contrast variance')
  beta <- sum(a); se <- sqrt(variance)
  data.frame(HR=exp(beta),lower95=exp(beta-qnorm(.975)*se),upper95=exp(beta+qnorm(.975)*se),
    beta=beta,SE=se,status='SAME_F1_COEFFICIENTS_AND_FULL_COVARIANCE',contrast=empty$contrast,
    note='Joint contrast from one interaction model; scale/model/endpoint approval remains separate',stringsAsFactors=FALSE)
}
