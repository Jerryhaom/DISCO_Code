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
sim_root <- file.path(disco_project_root(), "Simulation", "simdata")
out <- file.path(module, "outputs", "submission_simulation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(sim_root)) stop("投稿 simulation 附件不存在: ", sim_root)
suppressMessages({ library(survival); library(MASS) })

mahalanobis_distance <- function(x, ref) {
  mu <- colMeans(ref); S <- cov(ref)
  if (det(S) < 1e-10) S <- S + diag(1e-6, ncol(ref))
  d <- sweep(as.matrix(x), 2, mu)
  sqrt(rowSums((d %*% solve(S)) * d))
}
disco_distance <- function(x, ref) {
  n <- nrow(ref); cr <- cor(ref); mu <- colMeans(ref); sdv <- apply(ref, 2, sd)
  ref_ss <- (n - 1) * cov(ref)
  vapply(seq_len(nrow(x)), function(i) {
    # Exact one-row covariance update, algebraically equivalent to cor(rbind(ref, x[i,])).
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
  n <- nrow(ref); mu <- colMeans(ref); ref_ss <- (n - 1) * cov(ref)
  ref_cov <- ref_ss / (n - 1); p <- ncol(ref)
  vapply(seq_len(nrow(x)), function(i) {
    delta <- x[i, ] - mu
    new_cov <- (ref_ss + (n / (n + 1)) * (delta %o% delta)) / n
    s2 <- ref_cov + diag(eps, p); s1 <- new_cov + diag(eps, p)
    s2i <- solve(s2)
    trace_term <- sum(diag(s2i %*% s1))
    log_det_term <- determinant(s2, logarithm=TRUE)$modulus[1] - determinant(s1, logarithm=TRUE)$modulus[1]
    mean_diff <- delta / (n + 1)
    mean_term <- as.numeric(t(mean_diff) %*% s2i %*% mean_diff)
    0.5 * (trace_term - p + mean_term + log_det_term)
  }, numeric(1))
}
cindex <- function(z, time, status) unname(summary(coxph(Surv(time, status) ~ z))$concordance[1])
cindex_se <- function(z, time, status) unname(summary(coxph(Surv(time, status) ~ z))$concordance[2])
files <- list.files(file.path(sim_root, "experiments", "Data"), "_squared_data\\.csv$", recursive=TRUE, full.names=TRUE)
if (!length(files)) stop("未找到投稿 simulation CSV")
rows <- lapply(files, function(f) {
  d <- read.csv(f, check.names=FALSE); x <- as.matrix(d[grepl("^X", names(d))])
  if (!all(c("time", "status") %in% names(d))) stop("缺少 time/status: ", f)
  data.frame(file=basename(f), scenario=basename(dirname(dirname(f))), n=nrow(d), p=ncol(x), events=sum(d$status == 1, na.rm=TRUE),
    Mahalanobis_C=cindex(mahalanobis_distance(x,x), d$time,d$status), DISCO_C=cindex(disco_distance(x,x), d$time,d$status), Entropy_C=cindex(entropy_distance(x,x), d$time,d$status), KL_C=cindex(kl_distance(x,x), d$time,d$status),
    Mahalanobis_se=cindex_se(mahalanobis_distance(x,x), d$time,d$status), DISCO_se=cindex_se(disco_distance(x,x), d$time,d$status), Entropy_se=cindex_se(entropy_distance(x,x), d$time,d$status), KL_se=cindex_se(kl_distance(x,x), d$time,d$status))
})
res <- do.call(rbind, rows)
write.csv(res, file.path(out, "submission_simulation_cindex.csv"), row.names=FALSE)
saved <- readRDS(file.path(sim_root, "Simu.Res.rds")); saved_df <- do.call(rbind, saved); saved_df$file <- saved_df$filename; saved_df$filename <- NULL
key <- match(res$file, saved_df$file); if (anyNA(key)) stop("结果与 Simu.Res.rds 无法完全对应")
cmp <- data.frame(file=res$file, Mahalanobis_abs_diff=abs(res$Mahalanobis_C-saved_df$Mahalanobis_C.C[key]), DISCO_abs_diff=abs(res$DISCO_C-saved_df$DISCO_C.C[key]), Entropy_abs_diff=abs(res$Entropy_C-saved_df$Entropy_C.C[key]))
write.csv(cmp, file.path(out, "submission_simulation_comparison.csv"), row.names=FALSE)
status <- if (all(cmp$Mahalanobis_abs_diff < 1e-10 & cmp$DISCO_abs_diff < 1e-10 & cmp$Entropy_abs_diff < 1e-10)) "EXACT" else "NOT_EXACT_REQUIRES_REVIEW"
writeLines(c(paste0("status=", status), paste0("files=", nrow(res)), "source=Simulation/simdata/Simu.Res.rds", "KL=recomputed_not_in_Simu.Res_rds", paste0("R=", getRversion()), capture.output(sessionInfo())), file.path(out, "audit_status.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt")); cat(status, "\n")
