

library(ggplot2)
library(MASS)
library(ggthemes)
library(patchwork)
library(Matrix)

# 距离计算函数
mahalanobis_distance <- function(new_point, reference_data) {
  center <- colMeans(reference_data)
  cov_matrix <- cov(reference_data)
  
  # 处理协方差矩阵奇异性
  if (det(cov_matrix) < 1e-10) {
    cov_matrix <- cov_matrix + diag(1e-6, ncol(reference_data))
  }
  
  diff <- as.matrix(new_point - center)
  sqrt(t(diff) %*% solve(cov_matrix) %*% diff)
}

# DISCO函数（保持不变）
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
  ref_det <- det(cov_ref)
  
  combined <- rbind(reference_data, new_point)
  cov_new <- cov(combined)
  new_det <- det(cov_new)
  
  # 熵距离与行列式的对数变化相关
  log(new_det) - log(ref_det)
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





























standardize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}











































































































set.seed(123)
n_samples <- 1000  # 增加样本量以更好地表示高维空间
n_dim <- 10       # 10个维度

rho <- 0.1
cor_matrix <- toeplitz(rho^(0:(n_dim-1)))

var_values = runif(n_dim, 0.5, 1.5)
# Convert to covariance matrix
cov_matrix <- diag(sqrt(var_values)) %*% cor_matrix %*% diag(sqrt(var_values))
cov_matrix

# Generate data
mean_vector <- rep(0,n_dim)
data <- mvrnorm(n_samples, mean_vector, cov_matrix)
colnames(data) <- paste0("X", 1:n_dim)
as.data.frame(data)
reference_data=data

# 执行PCA降维
pca <- prcomp(reference_data, scale. = TRUE)

# 计算每个主成分的方差解释比例
variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

# 在PC1-PC2平面上创建网格
grid_size <- 100
pc1_range <- range(pca$x[,1]) * 1.5
pc2_range <- range(pca$x[,2]) * 1.5
pc1_grid <- seq(pc1_range[1], pc1_range[2], length.out = grid_size)
pc2_grid <- seq(pc2_range[1], pc2_range[2], length.out = grid_size)
pc_grid <- expand.grid(PC1 = pc1_grid, PC2 = pc2_grid)

# 固定其他主成分为0（只在PC1-PC2平面上变化）
pc_full_grid <- as.data.frame(matrix(0, nrow = nrow(pc_grid), ncol = n_dim))
colnames(pc_full_grid) <- paste0("PC", 1:n_dim)
pc_full_grid$PC1 <- pc_grid$PC1
pc_full_grid$PC2 <- pc_grid$PC2

# 将PC坐标转回原始特征空间
grid_points <- as.data.frame(
  pca$center + as.matrix(pc_full_grid) %*% t(pca$rotation)
)
dim(grid_points)


# 计算网格中每个点的三种距离
grid_distances <- data.frame(
  PC1 = pc_grid$PC1,
  PC2 = pc_grid$PC2,
  Mahalanobis = apply(grid_points, 1, function(row) mahalanobis_distance(row, reference_data)),
  DISCO = apply(grid_points, 1, function(row) disco_distance(row, reference_data)),
  Entropy = apply(grid_points, 1, function(row) entropy_distance(row, reference_data)),
  KL = apply(grid_points, 1, function(row) kl_distance(row, reference_data))
)


grid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)
grid_distances$DISCO <- standardize(grid_distances$DISCO)
grid_distances$Entropy <- standardize(grid_distances$Entropy)
grid_distances$KL <- standardize(log(grid_distances$KL+1))

# 转换为长格式
grid_distances_long <- reshape2::melt(grid_distances, id.vars = c("PC1", "PC2"))
colnames(grid_distances_long) <- c("PC1", "PC2", "Distance", "Value")

# 计算原始数据在PC1-PC2平面上的投影
pc_projection <- as.data.frame(pca$x[,1:2])
colnames(pc_projection) <- c("PC1", "PC2")

# 绘制三种距离的等高线图（在PC1-PC2平面上）
unique(grid_distances_long$Distance)

grid_distances_long$Distance=factor(grid_distances_long$Distance,
                                    levels = unique(grid_distances_long$Distance)[c(2,1,3,4)])

p_contour_low <- ggplot() +
  geom_raster(data = grid_distances_long, aes(PC1, PC2, fill = Value)) +
  geom_contour(data = grid_distances_long, aes(PC1, PC2, z = Value), 
               color = "black", bins = 10) +
  geom_point(data = pc_projection[sample(nrow(pc_projection),100),], 
             aes(PC1, PC2), 
             size = 0.05, color = "black") + 
  facet_wrap(~ Distance, ncol = 4) +
  scale_fill_distiller(palette = "Spectral") +
  labs(title = "Rho = 0.1",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_few() +
  theme(strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 

p_contour_low

########################
rho <- 0.7
cor_matrix <- toeplitz(rho^(0:(n_dim-1)))

var_values = runif(n_dim, 0.5, 1.5)
# Convert to covariance matrix
cov_matrix <- diag(sqrt(var_values)) %*% cor_matrix %*% diag(sqrt(var_values))
cov_matrix

# Generate data
mean_vector <- rep(0,n_dim)
data <- mvrnorm(n_samples, mean_vector, cov_matrix)
colnames(data) <- paste0("X", 1:n_dim)
as.data.frame(data)
reference_data=data

# 执行PCA降维
pca <- prcomp(reference_data, scale. = TRUE)

# 计算每个主成分的方差解释比例
variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

# 在PC1-PC2平面上创建网格
grid_size <- 100
pc1_range <- range(pca$x[,1]) * 1.5
pc2_range <- range(pca$x[,2]) * 1.5
pc1_grid <- seq(pc1_range[1], pc1_range[2], length.out = grid_size)
pc2_grid <- seq(pc2_range[1], pc2_range[2], length.out = grid_size)
pc_grid <- expand.grid(PC1 = pc1_grid, PC2 = pc2_grid)

# 固定其他主成分为0（只在PC1-PC2平面上变化）
pc_full_grid <- as.data.frame(matrix(0, nrow = nrow(pc_grid), ncol = n_dim))
colnames(pc_full_grid) <- paste0("PC", 1:n_dim)
pc_full_grid$PC1 <- pc_grid$PC1
pc_full_grid$PC2 <- pc_grid$PC2

# 将PC坐标转回原始特征空间
grid_points <- as.data.frame(
  pca$center + as.matrix(pc_full_grid) %*% t(pca$rotation)
)
dim(grid_points)


# 计算网格中每个点的三种距离
grid_distances <- data.frame(
  PC1 = pc_grid$PC1,
  PC2 = pc_grid$PC2,
  Mahalanobis = apply(grid_points, 1, function(row) mahalanobis_distance(row, reference_data)),
  DISCO = apply(grid_points, 1, function(row) disco_distance(row, reference_data)),
  Entropy = apply(grid_points, 1, function(row) entropy_distance(row, reference_data)),
  KL = apply(grid_points, 1, function(row) kl_distance(row, reference_data))
)
head(grid_distances)

grid_distances$Mahalanobis <- standardize(grid_distances$Mahalanobis)
grid_distances$DISCO <- standardize(grid_distances$DISCO)
grid_distances$Entropy <- standardize(grid_distances$Entropy)
grid_distances$KL <- standardize(log(grid_distances$KL+1))

# 转换为长格式
grid_distances_long <- reshape2::melt(grid_distances, id.vars = c("PC1", "PC2"))
colnames(grid_distances_long) <- c("PC1", "PC2", "Distance", "Value")

# 计算原始数据在PC1-PC2平面上的投影
pc_projection <- as.data.frame(pca$x[,1:2])
colnames(pc_projection) <- c("PC1", "PC2")

# 绘制三种距离的等高线图（在PC1-PC2平面上）
unique(grid_distances_long$Distance)

grid_distances_long$Distance=factor(grid_distances_long$Distance,
                                    levels = unique(grid_distances_long$Distance)[c(2,1,3,4)])

p_contour_high <- ggplot() +
  geom_raster(data = grid_distances_long, aes(PC1, PC2, fill = Value)) +
  geom_contour(data = grid_distances_long, aes(PC1, PC2, z = Value), 
               color = "black", bins = 10) +
  geom_point(data = pc_projection[sample(nrow(pc_projection),100),], 
             aes(PC1, PC2), 
             size = 0.05, color = "black") + 
  facet_wrap(~ Distance, ncol = 4) +
  scale_fill_distiller(palette = "Spectral") +
  labs(title = "Rho = 0.7",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_few() +
  theme(strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 

p_contour_high

library(ggpubr)
pdf("Cmp.DISCO_with_KL_high.pdf", height = 6, width = 8)
#print(p_contour_low)
ggarrange(p_contour_low,p_contour_high,nrow = 2,labels = c("a","b"))
dev.off()

#################
# Historical working-directory instruction omitted; the runner supplies output.
library(ggplot2)
library(MASS)
library(ggthemes)
library(patchwork)
plot_KL <- function(rho){
  
  # 生成10维示例数据
  set.seed(123)
  n_samples <- 1000  # 增加样本量以更好地表示高维空间
  n_dim <- 10       # 10个维度
  
  #rho <- 0.7
  cor_matrix <- toeplitz(rho^(0:(n_dim-1)))
  
  var_values = runif(n_dim, 0.5, 1.5)
  # Convert to covariance matrix
  cov_matrix <- diag(sqrt(var_values)) %*% cor_matrix %*% diag(sqrt(var_values))
  cov_matrix
  
  # Generate data
  data <- mvrnorm(n_samples, mean_vector, cov_matrix)
  colnames(data) <- paste0("X", 1:n_dim)
  as.data.frame(data)
  reference_data=data
  
  grid_distances <- data.frame(
    PC1 = pc_grid$PC1,
    PC2 = pc_grid$PC2,
    Mahalanobis = apply(reference_data, 1, function(row) mahalanobis_distance(row, reference_data)),
    DISCO = apply(reference_data, 1, function(row) disco_distance(row, reference_data)),
    Entropy = apply(reference_data, 1, function(row) entropy_distance(row, reference_data)),
    KL=apply(reference_data, 1, function(row) kl_distance(row, reference_data))
  )
  
  p1=ggplot(grid_distances,aes(KL,DISCO)) +
    geom_point(color="#E64B35FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p2=ggplot(grid_distances,aes(KL,Mahalanobis)) +
    geom_point(color="#4DBBD5FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p3=ggplot(grid_distances,aes(KL,Entropy)) +
    geom_point(color="#00A087FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
    
  library(ggpubr)
  
  annotate_figure(ggarrange(p1,NULL,p2,NULL,p3,
                            widths = c(1,0.05,1,0.05,1),nrow = 1),
                  top = text_grob(paste("Rho =",rho),face = "bold",
                                  size = 12, x = 0.02, hjust = 0))
  
}

p1=plot_KL(0.1)
p2=plot_KL(0.3)
p3=plot_KL(0.5)
p4=plot_KL(0.7)

pdf("DM2KL.pdf",height = 8,width = 8)
ggarrange(p1,p2,p3,p4,nrow = 4)
dev.off()

plot_KLL <- function(rho){
  
  # 生成10维示例数据
  set.seed(123)
  n_samples <- 1000  # 增加样本量以更好地表示高维空间
  n_dim <- 10       # 10个维度
  
  #rho <- 0.7
  cor_matrix <- toeplitz(rho^(0:(n_dim-1)))
  
  var_values = runif(n_dim, 0.5, 1.5)
  # Convert to covariance matrix
  cov_matrix <- diag(sqrt(var_values)) %*% cor_matrix %*% diag(sqrt(var_values))
  cov_matrix
  
  # Generate data
  data <- mvrnorm(n_samples, mean_vector, cov_matrix)
  colnames(data) <- paste0("X", 1:n_dim)
  as.data.frame(data)
  reference_data=data
  
  grid_distances <- data.frame(
    PC1 = pc_grid$PC1,
    PC2 = pc_grid$PC2,
    Mahalanobis = apply(reference_data, 1, function(row) mahalanobis_distance(row, reference_data)),
    DISCO = apply(reference_data, 1, function(row) disco_distance(row, reference_data)),
    Entropy = apply(reference_data, 1, function(row) entropy_distance(row, reference_data)),
    KL=apply(reference_data, 1, function(row) kl_distance(row, reference_data))
  )
  
  p1=ggplot(grid_distances,aes(Mahalanobis,DISCO)) +
    geom_point(color="#E64B35FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p2=ggplot(grid_distances,aes(Entropy,DISCO)) +
    geom_point(color="#4DBBD5FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p3=ggplot(grid_distances,aes(Mahalanobis,Entropy)) +
    geom_point(color="#00A087FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  library(ggpubr)
  
  annotate_figure(ggarrange(p1,NULL,p2,NULL,p3,
                            widths = c(1,0.05,1,0.05,1),nrow = 1),
                  top = text_grob(paste("Rho =",rho),face = "bold",
                                  size = 12, x = 0.02, hjust = 0))
  
}
p1=plot_KLL(0.1)
p2=plot_KLL(0.3)
p3=plot_KLL(0.5)
p4=plot_KLL(0.7)

pdf("DM3KL.pdf",height = 8,width = 8)
ggarrange(p1,p2,p3,p4,nrow = 4)
dev.off()
