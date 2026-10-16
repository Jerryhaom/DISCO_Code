  distance_data <- data.frame(
    Mahalanobis = mahalanobis_distance(data[, cl], data[, cl]),
    DISCO = apply(data[, cl], 1, disco_distance, data[, cl]),
    Entropy = entropy_distance(data[, cl], data[, cl]),
    KL = apply(data[, cl], 1, kl_distance, data[, cl]),  # 加入KL距离
    time = data$time,
    status = data$status
  )
