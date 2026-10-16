  cox_models <- list(
    coxph(Surv(time, status) ~ Mahalanobis, data = distance_data),
    coxph(Surv(time, status) ~ DISCO, data = distance_data),
    coxph(Surv(time, status) ~ Entropy, data = distance_data),
    coxph(Surv(time, status) ~ KL, data = distance_data)  # 加入KL距离模型
  )
  
  # 返回四种模型的C-index和p值
  c(
    "Mahalanobis_C" = summary(cox_models[[1]])$concordance[1],
    "Mahalanobis_se" = summary(cox_models[[1]])$concordance[2],
    "DISCO_C" = summary(cox_models[[2]])$concordance[1],
    "DISCO_se" = summary(cox_models[[2]])$concordance[2],
    "Entropy_C" = summary(cox_models[[3]])$concordance[1],
    "Entropy_se" = summary(cox_models[[3]])$concordance[2],
    "KL_C" = summary(cox_models[[4]])$concordance[1],      # 加入KL距离C-index
    "KL_se" = summary(cox_models[[4]])$concordance[2]       # 加入KL距离标准误
  )
