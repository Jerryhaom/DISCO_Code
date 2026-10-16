config <- list(
  n_standard = 1000,               # 标准样本量
  samples = c(10, 50, 100),       # 特征数量（多个值，重命名为samples避免歧义）
  rhos = c(0.1, 0.3, 0.5, 0.7),    # 相关系数
  noise_sd = 0.2,                  # 标准噪声
  high_noise_sd = 1,               # 高噪声
  seed = 123,                      # 主种子
  distributions = c("gompertz", "weibull", "exponential"),  # 生存分布
  censoring_rate = 0.3,            # 删失率
  maxt = 10,                       # 最大时间
  data_dir = "experiments/Data"    # 数据保存目录
)

generate_cov_matrix <- function(p, rho) {
  if (p < 1) stop("特征数量p必须为正整数")
  sigma <- toeplitz(rho^(0:(p-1)))
  diag(sigma) <- 1  # 确保方差为1
  sigma
}

# 生成基础数据（支持单一p值，动态种子确保可重复）
generate_base_data <- function(n, p, rho, seed) {
  set.seed(seed)  # 基于主种子+参数的独立种子
  mu <- rep(0, p)
  sigma <- generate_cov_matrix(p, rho)
  data <- as.data.table(rmvnorm(n = n, mean = mu, sigma = sigma))
  setnames(data, paste0("X", 1:p))  # 命名特征（X1, X2, ..., Xp）
  data
}

# 生成线性预测器（统一逻辑，减少冗余）
generate_lin_pred <- function(base_data, transform_fn, noise_sd, seed) {
  set.seed(seed)  # 噪声种子独立
  transformed_data <- base_data[, lapply(.SD, transform_fn)]  # 向量化变换
  true_betas <- rep(1, ncol(transformed_data))  # 真实系数
  lin_pred <- as.matrix(transformed_data) %*% true_betas  # 线性组合
  scaled_lin_pred <- scale(lin_pred) + rnorm(nrow(base_data), 0, noise_sd)  # 加噪声
  as.vector(scaled_lin_pred)
}

# --------------------------
# 3. 生存数据生成函数优化
# --------------------------
# 生成参数模型生存数据（结构化参数传递）
generate_surv_data <- function(base_data, transform_fn, dist_type, lambda, gamma, 
                               noise_sd, seed) {
  tryCatch({
    # 生成线性预测器（独立种子）
    lin_pred_seed <- seed + 100  # 与基础数据种子区分
    scaled_lin_pred <- generate_lin_pred(base_data, transform_fn, noise_sd, lin_pred_seed)
    
    # 生存数据参数
    surv_args <- list(
      dist = dist_type,
      lambdas = lambda,
      x = data.frame(y = scaled_lin_pred),
      betas = c(y = 1),
      maxt = config$maxt
    )
    if (dist_type != "exponential") {
      surv_args$gammas <- gamma  # 指数分布无gamma参数
    }
    
    # 生成生存数据
    set.seed(seed + 200)  # 生存数据种子
    surv_data <- do.call(simsurv, surv_args)
    
    # 整合数据
    data.table(
      id = surv_data$id,
      base_data,
      time = surv_data$eventtime,
      status = surv_data$status,
      y = scaled_lin_pred
    )
  }, error = function(e) {
    log_error("生成 {dist_type} 数据失败: {e$message}")
    NULL
  })
}

# 生成Cox模型数据（统一格式）
generate_cox_data <- function(base_data, transform_fn, noise_sd, seed) {
  tryCatch({
    # 生成线性预测器
    lin_pred_seed <- seed + 300
    scaled_lin_pred <- generate_lin_pred(base_data, transform_fn, noise_sd, lin_pred_seed)
    
    # 生成Cox数据
    set.seed(seed + 400)
    simdata <- sim.survdata(
      T = 1000,
      num.data.frames = 1,
      X = data.frame(y = scaled_lin_pred),
      beta = 1,
      censor = config$censoring_rate,
      spline = FALSE
    )
    
    # 整合数据
    data.table(
      id = 1:nrow(simdata$data),
      base_data,
      time = simdata$data$y / 100,  # 时间缩放
      status = as.integer(simdata$data$failed),
      y = scaled_lin_pred
    )
  }, error = function(e) {
    log_error("生成Cox数据失败: {e$message}")
    NULL
  })
}
