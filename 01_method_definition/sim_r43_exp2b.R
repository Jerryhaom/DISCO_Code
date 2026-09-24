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

# R4.3 Experiment 2b: DISCO mortality C-index vs p (real UKB proteomics).
# Reference = age<=40 (young); DISCO per individual vs young reference (efficient update).
# Cox Surv(time,status)~DISCO, unadjusted C-index, across p.
suppressMessages({library(ggplot2); library(dplyr); library(survival)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")

d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]; status <- mort$status[mm]; time <- mort$time[mm]
cat("matched:", sum(!is.na(age)), " individuals\n")

# ---- reference = age <= 40 ----
ref_idx <- which(age <= 40 & !is.na(age))
cat("reference (age<=40) N =", length(ref_idx), "\n")

# ---- analysis subset: 10000 non-reference with valid time/status ----
ok <- which(!is.na(time) & !is.na(status) & is.finite(age))
ok <- setdiff(ok, ref_idx)
set.seed(42); analysis <- sample(ok, 10000)

pgrid <- c(50, 100, 200, 500, 1000, 2000, 2920)
res <- list()
for (p in pgrid) {
  set.seed(3000 + p)
  sel <- sample(ncol(X), p)
  Xs <- X[, sel, drop = FALSE]
  Xref <- Xs[ref_idx, , drop = FALSE]; n <- nrow(Xref)
  mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)
  c1 <- (n - 1) / n; c2 <- 1 / (n + 1)
  Ma <- Xs[analysis, , drop = FALSE]
  disco <- vapply(seq_len(nrow(Ma)), function(i) {
    dd <- Ma[i, ] - mu
    Snew <- c1 * S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((Rref - Rnew)^2))
  }, numeric(1))
  dd <- data.frame(time = time[analysis], status = status[analysis], disco = disco)
  fit <- coxph(Surv(time, status) ~ disco, data = dd)
  ci <- summary(fit)$concordance[1]
  res[[length(res)+1]] <- data.frame(p = p, c_index = ci)
  cat(sprintf("p=%d C-index=%.4f\n", p, ci))
}
res <- do.call(rbind, res)
saveRDS(res, file.path(out, "R43_exp2b_cindex_vs_p.rds"))

pfig <- ggplot(res, aes(p, c_index)) +
  geom_point(size = 2.5) + geom_line() +
  scale_x_log10(breaks = pgrid) +
  theme_classic(base_size = 12) +
  labs(x = "Number of features (p)", y = "DISCO mortality C-index (unadjusted)",
       title = "R4.3 (real UKB proteomics): DISCO mortality C-index vs p",
       subtitle = "Reference age<=40; 10000-subject analysis; uniform weight")
ggsave(file.path(out, "R43_exp2b_cindex_vs_p.pdf"), pfig, width = 6.5, height = 5)
cat("wrote R43_exp2b_cindex_vs_p.pdf\n")
