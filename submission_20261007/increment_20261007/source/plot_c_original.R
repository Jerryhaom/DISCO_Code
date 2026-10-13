plot_c <- function(res, filenames = NULL){
  # 如果没有提供文件名，创建占位符
  if (is.null(filenames)) {
    res$filename <- paste0("file_", 1:nrow(res))
  } else {
    res$filename <- filenames
  }
  
  res$type <- res$rho <- res$filename
  res$type <- sapply(res$rho, function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) >= 4) parts[4] else "unknown"
  })
  
  res$rho <- as.numeric(gsub("rho", "", sapply(res$rho, function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) >= 3) parts[3] else "0"
  })))
  
  library(reshape2)
  
  # 提取C-index值（现在有8列）
  c_cols <- c("Mahalanobis_C", "DISCO_C", "Entropy_C", "KL_C")
  se_cols <- c("Mahalanobis_se", "DISCO_se", "Entropy_se", "KL_se")
  
  df_c <- melt(res[, c(c_cols, "rho", "type")], id.vars = c("rho", "type"))
  df_se <- melt(res[, c(se_cols, "rho", "type")], id.vars = c("rho", "type"))
  
  # 合并数据
  df <- data.frame(
    rho = df_c$rho,
    type = df_c$type,
    Distance = df_c$variable,
    value = df_c$value,
    Dse = df_se$value
  )
  
  # 清理距离名称
  df$Distance <- gsub("_C", "", df$Distance)
  
  # 设置因子水平，确保正确的顺序（DISCO, Mahalanobis, Entropy, KL）
  distance_levels <- c("DISCO", "Mahalanobis", "Entropy", "KL")
  df$Distance <- factor(df$Distance, levels = distance_levels)
  
  # 处理类型名称
  df$type <- gsub("exponential", "Exponential", df$type)
  df$type <- gsub("gompertz", "Gompertz", df$type)
  df$type <- gsub("weibull", "Weibull", df$type)
  df$rho <- paste("Corr =", df$rho)
  
  library(ggplot2)
  library(ggsci)
  
  dodge_width <- 0.8
  mx=max(df$value)
  
  ggplot(df, aes(x = type, y = value, fill = Distance)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = dodge_width),
             width = 0.7) +
    geom_errorbar(
      aes(ymin = value - Dse, ymax = value + Dse),
      position = position_dodge(width = dodge_width),
      width = 0.3) +
    facet_wrap(rho ~ .,nrow=1) +
    coord_cartesian(ylim = c(0.5, 1.05*mx)) +
    theme_few() + 
    ylab("C-index") + 
    xlab("") +
    scale_fill_npg(
      labels = c("DISCO", "Mahalanobis", "Entropy", "KL")
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      strip.text = element_text(size = 14)
    )
}
