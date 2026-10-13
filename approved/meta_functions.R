# Current Table 1 associations, not historical per-unit estimates.
table1_meta_input <- function(d) {
  cohorts <- c('NHANES','UKB clinical','CHARLS','CLHLS','RuLAS wave 2')
  fields <- c('cohort','metric','n','events','model','adjusted_HR_per_analysis_SD',
    'adjusted_HR_CI95_lower','adjusted_HR_CI95_upper','score_SD_analysis_sample',
    'standardized_formula','scaling','warning_count')
  if (!is.data.frame(d) || nrow(d)!=10L || !all(fields %in% names(d)))
    stop('Require ten current Table 1 rows with model/SD metadata')
  if (anyDuplicated(paste(d$cohort,d$metric)) ||
      !setequal(paste(d$cohort,d$metric),as.vector(outer(cohorts,c('DISCO','DM'),paste))))
    stop('Only five independent Table 1 cohorts and DISCO/log(DM) may be pooled')
  scaling <- 'One within-cohort sample SD in the effective Table 1 Cox analysis population; DM log-transformed before scaling'
  if (anyNA(d$scaling) || any(d$scaling != scaling)) stop('Mixed or undocumented SD/scoring estimands')
  race <- d$cohort %in% c('NHANES','UKB clinical')
  expect <- paste0('Surv(time,status)~score_z+age+sex+smoke+drink+education',ifelse(race,'+race',''))
  observed <- gsub('[[:space:]]','',d$standardized_formula)
  if (anyNA(observed) || any(observed!=expect) || any(d$model!=ifelse(race,'A','B')))
    stop('Mixed adjustment models; current Table 1 core adjustment required')
  numeric_fields <- c('n','events','adjusted_HR_per_analysis_SD','adjusted_HR_CI95_lower',
    'adjusted_HR_CI95_upper','score_SD_analysis_sample','warning_count')
  if (!all(vapply(d[numeric_fields],is.numeric,FALSE)) ||
      any(!is.finite(as.matrix(d[numeric_fields]))) || any(d$n<1 | d$events<1 | d$events>d$n) ||
      any(d$score_SD_analysis_sample<=0) || any(d$warning_count!=0)) stop('Invalid or warned Table 1 estimates')
  h <- d$adjusted_HR_per_analysis_SD; lo <- d$adjusted_HR_CI95_lower; hi <- d$adjusted_HR_CI95_upper
  if (any(lo<=0 | lo>=h | hi<=h)) stop('Invalid HR/CI values')
  for (cohort in cohorts) {
    r <- d[d$cohort==cohort,]
    if (length(unique(r$n))!=1 || length(unique(r$events))!=1) stop('DISCO/DM sample mismatch')
  }
  data.frame(cohort=d$cohort,metric=d$metric,n=d$n,events=d$events,
    HR=h,lower95=lo,upper95=hi,beta=log(h),SE=(log(hi)-log(lo))/(2*qnorm(.975)),
    score_SD=d$score_SD_analysis_sample,model_formula=d$standardized_formula,
    adjustment_family='age+sex+smoking+alcohol+education; race additionally NHANES and UKB',
    SD_population='final common within-cohort Cox analysis sample; sample SD (N-1)',
    score_transform=ifelse(d$metric=='DM','natural log DM','DISCO'),
    endpoint_policy_version='Table1_saved_cohort_endpoints_20261006_provisional',
    endpoint_status='PROVISIONAL: raw mortality/censoring definitions not all finally approved',
    UKB_race_missing='White single imputation approved 2026-10-06',
    CHARLS_education_missing='observed binary education scored-cohort mode approved 2026-10-06',
    stringsAsFactors=FALSE)
}
dl_pool <- function(beta, se, metric, left_out=NA_character_) {
  if (length(beta)<2 || length(beta)!=length(se) || any(!is.finite(c(beta,se))) || any(se<=0))
    stop('At least two finite independent estimates with positive SE required')
  k <- length(beta); w <- 1/se^2; fe <- sum(w*beta)/sum(w)
  Q <- sum(w*(beta-fe)^2); df <- k-1
  tau2 <- max(0,(Q-df)/(sum(w)-sum(w^2)/sum(w)))
  wr <- 1/(se^2+tau2); b <- sum(wr*beta)/sum(wr); s <- sqrt(1/sum(wr))
  data.frame(metric=metric,left_out=left_out,k=k,beta=b,SE=s,
    HR=exp(b),lower95=exp(b-qnorm(.975)*s),upper95=exp(b+qnorm(.975)*s),
    P=2*pnorm(-abs(b/s)),tau2=tau2,Q=Q,df=df,I2=if(Q==0)0 else max(0,(Q-df)/Q)*100,
    heterogeneity_P=pchisq(Q,df,lower.tail=FALSE),method='DerSimonian-Laird random effects; Wald normal 95% CI',
    estimand='adjusted association per one final within-cohort analysis-sample SD',
    adjustment='core age+sex+smoking+alcohol+education; race additionally NHANES and UKB',
    endpoint_policy_version='Table1_saved_cohort_endpoints_20261006_provisional',stringsAsFactors=FALSE)
}
