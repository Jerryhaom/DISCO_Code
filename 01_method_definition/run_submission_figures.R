#!/usr/bin/env Rscript
# Portable release bootstrap (source/data/output roots are independent).
local({
  frames <- vapply(sys.frames(), function(e) if (is.null(e$ofile)) "" else as.character(e$ofile)[1], character(1))
  frames <- frames[nzchar(frames)]
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  entry <- if (length(frames)) tail(frames, 1) else if (length(arg)) sub("^--file=", "", arg[1]) else ""
  cfg <- file.path(dirname(normalizePath(gsub("~+~", " ", entry, fixed=TRUE), mustWork = TRUE)), "..", "project_config.R")
  source(cfg, local = .GlobalEnv)
})
disco_init_module("01_method_definition")

options(stringsAsFactors = FALSE)
module <- file.path(disco_results_root(), "01_method_definition")
out <- file.path(module, "outputs", "submission_figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressMessages({ library(MASS); library(ggplot2) })

mahalanobis_distance <- function(x, ref) {
  mu <- colMeans(ref); S <- cov(ref)
  if (det(S) < 1e-10) S <- S + diag(1e-6, ncol(ref))
  d <- sweep(as.matrix(x), 2, mu)
  sqrt(rowSums((d %*% solve(S)) * d))
}
disco_distance <- function(x, ref) {
  n <- nrow(ref); cr <- cor(ref); mu <- colMeans(ref); ref_ss <- (n - 1) * cov(ref)
  vapply(seq_len(nrow(x)), function(i) {
    delta <- x[i, ] - mu
    ss <- ref_ss + (n / (n + 1)) * (delta %o% delta)
    cn <- ss / sqrt(diag(ss) %o% diag(ss))
    log(sum(0.5 * (cr - cn)^2) * n^2)
  }, numeric(1))
}
entropy_distance <- function(x, ref) {
  n <- nrow(ref); mu <- colMeans(ref); ref_ss <- (n - 1) * cov(ref)
  lr <- determinant(ref_ss / (n - 1), logarithm=TRUE)$modulus[1]
  vapply(seq_len(nrow(x)), function(i) {
    delta <- x[i, ] - mu
    ss <- ref_ss + (n / (n + 1)) * (delta %o% delta)
    0.5 * (determinant(ss / n, logarithm=TRUE)$modulus[1] - lr)
  }, numeric(1))
}
kl_distance <- function(x, ref, eps = 1e-10) {
  n <- nrow(ref); mu <- colMeans(ref); ref_ss <- (n - 1) * cov(ref); ref_cov <- ref_ss / (n - 1); p <- ncol(ref)
  vapply(seq_len(nrow(x)), function(i) {
    delta <- x[i, ] - mu
    new_cov <- (ref_ss + (n / (n + 1)) * (delta %o% delta)) / n
    s2 <- ref_cov + diag(eps, p); s1 <- new_cov + diag(eps, p); s2i <- solve(s2)
    0.5 * (sum(diag(s2i %*% s1)) - p + as.numeric(t(delta / (n + 1)) %*% s2i %*% (delta / (n + 1))) +
             determinant(s2, logarithm=TRUE)$modulus[1] - determinant(s1, logarithm=TRUE)$modulus[1])
  }, numeric(1))
}

# Figure S5: noise response. Faithful to Sim_Distances.R plot_KL2, with a single
# set.seed(123) base-data generation and one sequential RNG stream for the noise
# draws (identical to the submission script), so the result is bit-reproducible.
n_samples <- 1000; n_dim <- 10; rho <- 0.7
set.seed(123)
cor_matrix <- toeplitz(rho^(0:(n_dim - 1)))
var_values <- runif(n_dim, 0.5, 1.5)
cov_matrix <- diag(sqrt(var_values)) %*% cor_matrix %*% diag(sqrt(var_values))
data <- as.matrix(mvrnorm(n_samples, rep(0, n_dim), cov_matrix))
noise_levels <- seq(0, 2, 0.2)
res <- do.call(rbind, lapply(noise_levels, function(nl) {
  data1 <- data + matrix(rnorm(ncol(data) * nrow(data), 0, nl), ncol = ncol(data))
  data.frame(noise = nl,
             Mahalanobis = mahalanobis_distance(data1, data),
             DISCO = disco_distance(data1, data),
             Entropy = entropy_distance(data1, data),
             KL = kl_distance(data1, data))
}))
write.csv(res, file.path(out, "S5_noise_distances.csv"), row.names = FALSE)

cors <- data.frame(
  metric = c("Mahalanobis", "DISCO", "Entropy", "KL"),
  spearman_with_noise = c(
    cor(res$noise, res$Mahalanobis, method = "spearman"),
    cor(res$noise, res$DISCO, method = "spearman"),
    cor(res$noise, res$Entropy, method = "spearman"),
    cor(res$noise, res$KL, method = "spearman")))
write.csv(cors, file.path(out, "S5_noise_spearman.csv"), row.names = FALSE)

long <- reshape(res, direction = "long", varying = c("Mahalanobis", "DISCO", "Entropy", "KL"),
                v.names = "value", timevar = "metric", times = c("Mahalanobis", "DISCO", "Entropy", "KL"))
long$noise_f <- sprintf("%.1f", long$noise)
g <- ggplot(long, aes(noise_f, value)) +
  geom_boxplot(aes(fill = noise), outlier.colour = NA) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(x = "Noise levels", y = "Distance", fill = "Var") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(out, "S5_noise.pdf"), g, height = 5, width = 8)

cat("S5 Spearman correlations with noise:\n"); print(cors, row.names = FALSE)

# Figures S6/S7: survival benchmark bar plots from the verified C-index table.
suppressMessages({ library(reshape2) })
ci <- read.csv(file.path(module, "outputs", "submission_simulation", "submission_simulation_cindex.csv"))
ci$p <- as.integer(sub(".*_p([0-9]+)_.*", "\\1", ci$file))
ci$rho <- as.numeric(sub(".*_rho([0-9.]+)_.*", "\\1", ci$file))
ci$dist <- sub(".*_(Cox|exponential|gompertz|weibull)_.*", "\\1", ci$file)
ci$scenario <- sub("_p.*", "", ci$file)
ci$dist <- sub("exponential", "Exponential", sub("gompertz", "Gompertz", sub("weibull", "Weibull", ci$dist)))
metrics <- c("DISCO", "Mahalanobis", "Entropy", "KL")
mC <- melt(ci, id.vars = c("scenario", "p", "rho", "dist"), measure.vars = paste0(metrics, "_C"), variable.name = "Distance", value.name = "value")
mse <- melt(ci, id.vars = c("scenario", "p", "rho", "dist"), measure.vars = paste0(metrics, "_se"), variable.name = "Distance", value.name = "se")
mC$Distance <- sub("_C", "", mC$Distance); mse$Distance <- sub("_se", "", mse$Distance)
df <- merge(mC, mse, by = c("scenario", "p", "rho", "dist", "Distance"))
df$Distance <- factor(df$Distance, levels = metrics)
df$rho_f <- paste("Corr =", df$rho)
for (sc in c("Standard", "HighNoise")) {
  for (pv in c(10, 50, 100)) {
    sub <- df[df$scenario == sc & df$p == pv, ]
    g <- ggplot(sub, aes(dist, value, fill = Distance)) +
      geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
      geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge(0.8), width = 0.3) +
      facet_wrap(~ rho_f, nrow = 1) + ylab("C-index") + xlab("") +
      ggtitle(sprintf("%s variables, Noise ~ N(0,%s)", pv, ifelse(sc == "Standard", "0.2", "1"))) +
      theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out, sprintf("S%s_p%s_%s.pdf", ifelse(sc == "Standard", "6", "7"), pv, sc)), g, height = 4, width = 9)
  }
}
cat("S6/S7 survival benchmark figures written to", out, "\n")
