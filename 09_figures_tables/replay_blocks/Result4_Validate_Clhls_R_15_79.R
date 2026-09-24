entropy_distance_parallel_fast <- function(data, reference_data, n_cores = NULL) {
  # 数据预处理 - 移除常数列
  ref_data <- as.matrix(reference_data)
  test_data <- as.matrix(data)
  
  # 找出并移除常数列
  col_sd_ref <- apply(ref_data, 2, sd)
  col_sd_test <- apply(test_data, 2, sd)
  keep_cols <- (col_sd_ref > 0) & (col_sd_test > 0)
  
  if (sum(!keep_cols) > 0) {
    cat("Removed", sum(!keep_cols), "constant columns\n")
  }
  
  ref_data <- ref_data[, keep_cols, drop = FALSE]
  test_data <- test_data[, keep_cols, drop = FALSE]
  
  # 检查数据维度
  n_ref <- nrow(ref_data)
  p <- ncol(ref_data)
  
  # 预计算参考统计量
  ref_mean <- colMeans(ref_data)
  cov_ref <- cov(ref_data)
  ref_log_det <- determinant(cov_ref, logarithm = TRUE)$modulus
  
  # 设置并行核心数
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  # 注册并行后端
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 导出必要数据到工作进程
  clusterExport(cl, c("ref_data", "ref_mean", "n_ref", "p", "ref_log_det"), 
                envir = environment())
  
  # 并行计算每个测试样本
  result <- foreach(i = 1:nrow(test_data), .combine = c, .packages = "matrixcalc") %dopar% {
    new_point <- as.numeric(test_data[i, ])
    
    # 使用增量公式更新协方差矩阵（更高效）
    # 增量更新均值
    mu_new <- (n_ref * ref_mean + new_point) / (n_ref + 1)
    
    # 增量更新协方差（Sherman-Morrison公式变体）
    delta <- new_point - ref_mean
    S_ref <- (n_ref - 1) * cov_ref
    S_new <- S_ref + (n_ref / (n_ref + 1)) * outer(delta, delta)
    cov_new <- S_new / n_ref
    
    # 计算对数行列式
    new_log_det <- determinant(cov_new, logarithm = TRUE)$modulus
    
    # 返回熵距离
    new_log_det - ref_log_det
  }
  
  # 关闭并行后端
  stopCluster(cl)
  
  return(result)
}
