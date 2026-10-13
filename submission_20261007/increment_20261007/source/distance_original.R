mahalanobis_distance <- function(new_point, reference_data) {
  center <- colMeans(reference_data)
  cov_matrix <- cov(reference_data)
  
  # 处理协方差矩阵奇异性
  if (det(cov_matrix) < 1e-10) {
    cov_matrix <- cov_matrix + diag(1e-6, ncol(reference_data))
  }
  
  diff <- as.matrix(new_point - matrix(center, nrow = nrow(new_point), 
                                       ncol = ncol(new_point), byrow = TRUE))
  
  sqrt(rowSums((diff %*% solve(cov_matrix)) * diff))
}

disco_distance <- function(new_point, reference_data, weights = NULL) {
  n_ref <- nrow(reference_data)
  p <- ncol(reference_data)
  
  if (is.null(weights)) weights <- matrix(1, p, p)
  ref_cor <- cor(reference_data)
  
  combined_data <- rbind(reference_data, new_point)
  
  new_cor <- cor(combined_data)
  cor_diff <- ref_cor - new_cor
  
  weighted_diff <- 0.5*cor_diff^2 * weights
  
  sum_diff <- sum(as.numeric(weighted_diff))
  
  log(sum_diff * n_ref^2)
}

entropy_distance <- function(new_point, reference_data) {
  p <- ncol(reference_data)
  cov_ref <- cov(reference_data)
  
  ref_entropy <- 0.5 * (p * log(2 * pi * exp(1)) + log(det(cov_ref)))
  
  sapply(1:nrow(new_point), function(i) {
    combined <- rbind(reference_data, new_point[i, ])
    cov_new <- cov(combined)
    new_entropy <- 0.5 * (p * log(2 * pi * exp(1)) + log(det(cov_new)))
    new_entropy - ref_entropy
  })
}

kl_divergence_cor_specific <- kl_divergence_correct <- function(mu1, sigma1, mu2, sigma2, eps = 1e-10) {
  # 参数验证
  p <- length(mu1)
  if (nrow(sigma1) != p || nrow(sigma2) != p) {
    stop("均值和协方差矩阵维度不匹配")
  }
  
  # 正则化避免奇异矩阵
  sigma2_reg <- sigma2 + diag(eps, p)
  sigma1_reg <- sigma1 + diag(eps, p)
  
  # 计算逆矩阵
  sigma2_inv <- solve(sigma2_reg)
  
  # 正确的迹项计算
  trace_term <- sum(diag(sigma2_inv %*% sigma1_reg))
  
  # 优化行列式计算
  log_det_sigma2 <- determinant(sigma2_reg, logarithm = TRUE)$modulus[1]
  log_det_sigma1 <- determinant(sigma1_reg, logarithm = TRUE)$modulus[1]
  log_det_term <- log_det_sigma2 - log_det_sigma1
  
  # 均值差异项
  mean_diff <- as.vector(mu1 - mu2)
  mean_term <- as.numeric(t(mean_diff) %*% sigma2_inv %*% mean_diff)
  
  # KL散度公式
  0.5 * (trace_term - p + mean_term + log_det_term)
}

kl_distance  <- function(new_point, reference_data) {
  # 从原始数据计算参数（而不是协方差矩阵）
  ref_mean <- colMeans(reference_data)
  ref_cov <- cov(reference_data)
  
  combined_data <- rbind(reference_data, new_point)
  new_mean <- colMeans(combined_data)
  new_cov <- cov(combined_data)
  
  # 计算对称KL距离
  kl_forward <- kl_divergence_correct(new_mean, new_cov, ref_mean, ref_cov)
  #kl_reverse <- kl_divergence_correct(ref_mean, ref_cov, new_mean, new_cov)
  kl_forward
  #0.5 * (kl_forward + kl_reverse)
}
