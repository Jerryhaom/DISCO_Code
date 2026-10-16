# Fit the actual original four-distance Cox segments to new synthetic data.
fit_simulation <- function(root, out, manifest) {
  env <- new.env(parent = globalenv())
  source(file.path(root, "source/distance_original.R"), local = env)
  scores <- parse(file.path(root, "source/score_block_original.R"))
  models <- parse(file.path(root, "source/Cox_block_original.R"))
  cache <- records <- long <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i]
    env$data <- as.data.frame(fread(file.path(out, "generated_data", row$file)))
    env$cl <- grep("^X", names(env$data)); x <- as.matrix(env$data[, env$cl])
    key <- paste(row$scenario, row$p, row$rho, sep = "|")
    if (is.null(cache[[key]])) {
      eval(scores, envir = env)
      cache[[key]] <- list(X = x, distances = env$distance_data[, 1:4])
    } else {
      stopifnot(identical(cache[[key]]$X, x))
      env$distance_data <- data.frame(cache[[key]]$distances, time = env$data$time, status = env$data$status)
    }
    result <- eval(models, envir = env)
    wide <- as.data.table(as.list(unname(result)))
    setnames(wide, c("Mahalanobis_C", "Mahalanobis_se", "DISCO_C", "DISCO_se", "Entropy_C", "Entropy_se", "KL_C", "KL_se"))
    wide[, `:=`(file = row$file, scenario = row$scenario, p = row$p, rho = row$rho,
      distribution = row$distribution, N = nrow(env$data), events = sum(env$data$status))]
    records[[row$file]] <- wide
    for (j in seq_along(env$cox_models)) {
      su <- summary(env$cox_models[[j]]); independent <- concordance(env$cox_models[[j]])
      stopifnot(env$cox_models[[j]]$n == 1000L, env$cox_models[[j]]$nevent == row$events,
        abs(independent$concordance - su$concordance[1]) < 1e-12,
        abs(sqrt(independent$var) - su$concordance[2]) < 1e-12)
      metric_name <- c("Mahalanobis", "DISCO", "Entropy", "KL")[j]
      long[[paste(row$file, metric_name)]] <- data.table(file = row$file, metric = metric_name,
        N = 1000L, events = row$events, C = unname(su$concordance[1]), C_se = unname(su$concordance[2]))
    }
    cat("FITTED", row$file, "4 original models\n")
  }
  wide <- rbindlist(records); long <- rbindlist(long)
  fwrite(wide, file.path(out, "models_wide.csv")); fwrite(long, file.path(out, "models_long.csv"))
  reference <- fread(file.path(root, "input/reference_models_384.csv"))
  comparison <- merge(long, reference[, .(file, metric, reference_C = C, reference_SE = C_se)], by = c("file", "metric"))
  stopifnot(nrow(comparison) == nrow(long))
  comparison[, `:=`(C_difference = C - reference_C, SE_difference = C_se - reference_SE)]
  stopifnot(max(abs(comparison$C_difference)) < 1e-10, max(abs(comparison$SE_difference)) < 1e-10)
  fwrite(comparison, file.path(out, "numeric_comparison.csv"))
  write_json(list(status = "PASS", datasets = nrow(wide), newly_fitted_models = nrow(long),
    no_existing_score_input = TRUE, no_existing_time_status_input = TRUE,
    max_C_difference = max(abs(comparison$C_difference)), max_SE_difference = max(abs(comparison$SE_difference)),
    comparison_after_fitting_only = TRUE), file.path(out, "fit_QA.json"), auto_unbox = TRUE, pretty = TRUE)
  wide
}
