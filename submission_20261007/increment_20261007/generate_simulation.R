# Execute recovered original functions. Existing simulated CSVs are never inputs.
generate_simulation <- function(root, out, smoke = FALSE) {
  source(file.path(root, "source/generator_original.R"), local = FALSE)
  log_error <- function(...) stop(paste(..., collapse = " "))
  assign("log_error", log_error, envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  data_out <- file.path(out, "generated_data")
  dir.create(data_out, showWarnings = FALSE)
  grid <- CJ(scenario = c("Standard", "HighNoise"), p = config$samples,
    rho = config$rhos, distribution = c(config$distributions, "Cox"))
  if (smoke) grid <- grid[scenario == "Standard" & p == 10L & rho == 0.1]
  cache <- list(); rows <- list(); warnings <- list(); context <- ""
  logged <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- data.table(file = context, message = conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (i in seq_len(nrow(grid))) {
    g <- grid[i]
    base_seed <- config$seed + g$p + round(g$rho * 10) + ifelse(g$scenario == "Standard", 0L, 1000L)
    transform_seed <- base_seed + 100L
    key <- paste(g$scenario, g$p, g$rho, sep = "|")
    if (is.null(cache[[key]])) cache[[key]] <- generate_base_data(config$n_standard, g$p, g$rho, base_seed)
    noise <- if (g$scenario == "Standard") config$noise_sd else config$high_noise_sd
    context <- sprintf("%s_p%s_rho%.1f_%s_squared_data.csv", g$scenario, g$p, g$rho, g$distribution)
    d <- logged(if (g$distribution == "Cox") {
      generate_cox_data(cache[[key]], function(x) x^2, noise, transform_seed + 500L)
    } else {
      generate_surv_data(cache[[key]], function(x) x^2, g$distribution, 0.1,
        if (g$distribution == "gompertz") 0.1 else 1, noise, transform_seed)
    })
    stopifnot(!is.null(d), nrow(d) == 1000L, all(is.finite(d$time)), all(d$time > 0),
      all(d$status %in% 0:1))
    file <- file.path(data_out, context); fwrite(d, file)
    rows[[context]] <- data.table(file = context, scenario = g$scenario, p = g$p,
      rho = g$rho, distribution = g$distribution, N = nrow(d), events = sum(d$status),
      base_seed = base_seed, predictor_seed = transform_seed + ifelse(g$distribution == "Cox", 800L, 100L),
      survival_seed = transform_seed + ifelse(g$distribution == "Cox", 900L, 200L),
      sha256 = digest(file = file, algo = "sha256"))
    cat("GENERATED", context, "\n")
  }
  manifest <- rbindlist(rows)
  fwrite(manifest, file.path(out, "generation_manifest.csv"))
  if (length(warnings)) fwrite(rbindlist(warnings), file.path(out, "generator_warnings.csv"))
  write_json(list(datasets = nrow(grid), synthetic_rows = nrow(grid) * 1000L,
    mode = if (smoke) "FOUR_DISTRIBUTION_SMOKE_NOT_ALL96" else "FULL96",
    existing_dataset_read = FALSE, rng = RNGkind(), original_algorithm_changes = "None",
    source_sha256 = digest(file = file.path(root, "source/generator_original.R"), algo = "sha256")),
    file.path(out, "generation_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  manifest
}
