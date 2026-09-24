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


# 01_method simulation lock（可复现版，回应 R2 全面 simulation + R4.3 高维）
# 统一 DGP + 固定 seed + 明确 grid + replicate，输出机器可读 C-index 对比。
# 比较 DISCO / DM(Mahalanobis) / Entropy 在检测 dysregulation 上的 survival 预测能力。
# 不改动历史 Simulation 脚本。

options(stringsAsFactors = FALSE)
suppressMessages({library(MASS); library(survival)})

out_dir <- disco_results_path("01_method_definition/outputs/simulation_lock")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 距离函数 ----
dm_dist <- function(x, ref) {
  mu <- colMeans(ref); S <- cov(ref) + diag(1e-8, ncol(ref))
  sqrt(mahalanobis(x, mu, S))
}
disco_dist <- function(x, ref, weight = NULL) {
  p <- ncol(ref); n0 <- nrow(ref)
  if (is.null(weight)) weight <- matrix(1, p, p)
  ref_cor <- cor(ref)
  sapply(seq_len(nrow(x)), function(i) {
    cc1 <- cor(rbind(ref, x[i, , drop = FALSE]))
    ds <- sum((ref_cor - cc1)^2 * weight)
    log(ds * n0 * n0)
  })
}
entropy_dist <- function(x, ref) {
  mu <- colMeans(ref); S <- cov(ref) + diag(1e-8, ncol(ref))
  ref_det <- det(S)
  sapply(seq_len(nrow(x)), function(i) {
    Sn <- cov(rbind(ref, x[i, , drop = FALSE])) + diag(1e-8, ncol(ref))
    log(det(Sn)) - log(ref_det)
  })
}

# ---- DGP（dysregulation 真相 -> 生存）----
simulate_case <- function(n_ref = 200, n_test = 200, p = 20, rho = 0.5,
                          noise_sd = 0.2, transform = "linear", seed = 1) {
  set.seed(seed)
  sigma <- toeplitz(rho^(0:(p - 1)))
  ref <- mvrnorm(n_ref, rep(0, p), sigma)
  test <- mvrnorm(n_test, rep(0, p), sigma)
  # 潜在 dysregulation 分数（真相）：前 1/3 biomarker 的加权平均 + 噪声
  k <- max(2, round(p / 3))
  latent <- test[, 1:k] %*% rep(1, k) / k
  if (transform == "squared") latent <- latent^2  # 非线性
  # dysregulated 个体：latent 高的，biomarker 偏移（相关结构被扰动）
  dysreg <- latent > quantile(latent, 0.5)
  shift <- 1.2
  test[dysreg, 1:k] <- test[dysreg, 1:k] + shift
  test <- test + matrix(rnorm(n_test * p, 0, noise_sd), n_test, p)
  # 生存：死亡风险与 latent 正相关
  lp <- scale(latent)[, 1]
  haz <- exp(0.6 * lp)
  time <- rexp(n_test, rate = 0.1 * haz)
  status <- ifelse(time < 10, 1, 0)
  time[time >= 10] <- 10
  list(test = test, ref = ref, time = time, status = status, latent = latent)
}

c_index_of <- function(dist_vec, time, status) {
  if (sd(dist_vec) < 1e-12) return(NA)
  d <- data.frame(time = time, status = status, dist = dist_vec)
  fit <- coxph(Surv(time, status) ~ dist, data = d)
  unname(summary(fit)$concordance[1])
}

# ---- grid + replicate ----
grid <- expand.grid(
  p = c(10, 50, 100),
  rho = c(0.1, 0.5),
  noise_sd = c(0.2, 1.0),
  transform = c("linear", "squared"),
  stringsAsFactors = FALSE
)
R <- 20  # replicate per cell

res <- list()
for (g in seq_len(nrow(grid))) {
  p <- grid$p[g]; rho <- grid$rho[g]; noise <- grid$noise_sd[g]; tr <- grid$transform[g]
  for (r in seq_len(R)) {
    seed <- 1000 + g * 1000 + r
    sim <- simulate_case(p = p, rho = rho, noise_sd = noise, transform = tr, seed = seed)
    d_dm <- dm_dist(sim$test, sim$ref)
    d_di <- disco_dist(sim$test, sim$ref)
    d_en <- entropy_dist(sim$test, sim$ref)
    res[[length(res) + 1]] <- data.frame(
      p = p, rho = rho, noise_sd = noise, transform = tr, rep = r, seed = seed,
      C_DISCO = c_index_of(d_di, sim$time, sim$status),
      C_DM = c_index_of(d_dm, sim$time, sim$status),
      C_Entropy = c_index_of(d_en, sim$time, sim$status)
    )
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "simulation_lock_results.csv"), row.names = FALSE)

# 汇总（每格均值）
agg <- aggregate(cbind(C_DISCO, C_DM, C_Entropy) ~ p + rho + noise_sd + transform,
                 data = res, FUN = mean, na.rm = TRUE)
write.csv(agg, file.path(out_dir, "simulation_lock_summary.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "seed_fixed=true",
  "replicates=20_per_cell",
  "grid=p(10,50,100)_rho(0.1,0.5)_noise(0.2,1.0)_transform(linear,squared)",
  "metrics=DISCO_DM_Entropy_Cindex",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))

cat("=== simulation lock 汇总（每格平均 C-index）===\n")
print(agg, digits = 3, row.names = FALSE)
cat("\nDONE\n")
