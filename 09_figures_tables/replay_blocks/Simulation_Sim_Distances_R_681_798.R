plot_KL3 <- function(noise,rho=0.7){
  
  # 生成10维示例数据
  set.seed(123)
  n_samples <- 1000  # 增加样本量以更好地表示高维空间
  n_dim <- 10       # 10个维度
  
  
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
  
  noise=0.1
  grid_distances=data.frame()
  for(noise in seq(0,2,0.2)){
    
    data1 = data + matrix(rnorm(ncol(data)*nrow(data),0,noise),ncol = ncol(data))
    
    distances <- data.frame(
      noise=noise,
      Mahalanobis = apply(data1, 1, function(row) mahalanobis_distance(row, reference_data)),
      DISCO = apply(data1, 1, function(row) disco_distance(row, reference_data)),
      Entropy = apply(data1, 1, function(row) entropy_distance(row, reference_data)),
      KL=apply(data1, 1, function(row) kl_distance(row, reference_data)))
    
    
    grid_distances = rbind(grid_distances,distances)
    
  }
  
  cor.test(grid_distances$noise,grid_distances$DISCO,method = "spearman")
  cor.test(grid_distances$noise,grid_distances$KL,method = "spearman")
  cor.test(grid_distances$noise,grid_distances$Mahalanobis,method = "spearman")
  
  
  #grid_distances$noise=paste("N",grid_distances$noise)
  p1=ggplot(grid_distances,aes(as.factor(noise),DISCO)) +
    geom_boxplot(aes(fill=noise),outlier.colour = NA) +
    theme_classic() +
    stat_cor(aes(noise,DISCO), label.x.npc = 0.1,
             method = "spearman",cor.coef.name = "rho") +
    scale_fill_distiller(palette = "Reds",direction = 1) +
    scale_x_discrete(labels=sprintf("%.1f",seq(0,2,0.2))) +
    xlab("Noise levels") + labs(fill="Var") +theme(legend.position = "None")
  
  p2=ggplot(grid_distances,aes(as.factor(noise),Mahalanobis)) +
    geom_boxplot(aes(fill=noise),outlier.colour = NA) +
    theme_classic() +
    stat_cor(aes(noise,Mahalanobis), label.x.npc = 0.1,
             method = "spearman",cor.coef.name = "rho") +
    scale_x_discrete(labels=sprintf("%.1f",seq(0,2,0.2))) +
    scale_fill_distiller(palette = "Reds",direction = 1) +
    xlab("Noise levels") + labs(fill="Var") +theme(legend.position = "None")
  
  p3=ggplot(grid_distances,aes(as.factor(noise),Entropy)) +
    geom_boxplot(aes(fill=noise),outlier.colour = NA) +
    theme_classic() +
    stat_cor(aes(noise,Entropy), label.x.npc = 0.1,
             method = "spearman",cor.coef.name = "rho") +
    scale_fill_distiller(palette = "Reds",direction = 1) +
    scale_x_discrete(labels=sprintf("%.1f",seq(0,2,0.2))) +
    xlab("Noise levels") + labs(fill="Var") +theme(legend.position = "None")
  
  p11=ggplot(grid_distances,aes(Mahalanobis,DISCO)) +
    geom_point(color="#E64B35FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p21=ggplot(grid_distances,aes(Entropy,DISCO)) +
    geom_point(color="#4DBBD5FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  p31=ggplot(grid_distances,aes(Mahalanobis,Entropy)) +
    geom_point(color="#00A087FF",size=0.5) +
    stat_cor(method = "spearman",cor.coef.name = "rho",
             label.x.npc = 0.7,
             label.y.npc = "bottom",
             aes(label = ..r.label..)) +
    theme_classic()
  
  
  a1=annotate_figure(ggarrange(p1,NULL,p2,NULL,p3,
                            widths = c(1,0.05,1,0.05,1),nrow = 1,ncol = 5),
                  top = text_grob(paste("Rho =",rho),face = "bold",
                                  size = 12, x = 0.05, hjust = 0))
  
  a2=annotate_figure(ggarrange(p11,NULL,p21,NULL,p31,
                            widths = c(1,0.05,1,0.05,1),nrow = 1,ncol = 5),
                  top = text_grob(paste("Rho =",rho),face = "bold",
                                  size = 12, x = 0.05, hjust = 0))
  list(a1,a2)
}

p11=plot_KL3(1,0.1)
p41=plot_KL3(1,0.7)

pdf("Ex.Fig1.pdf",height = 9,width = 9)
ggarrange(p11[[1]],p41[[1]],
          p11[[2]],p41[[2]],
          labels = c("a","b","c","d"),
          nrow = 4,ncol = 1)
dev.off()
