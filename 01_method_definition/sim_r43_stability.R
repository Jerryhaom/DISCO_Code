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

# R4.3 metric 1: two reference lines on one plot.
#   (a) Spearman(DISCO(p), DISCO(50))   -- drift from minimal baseline
#   (b) Spearman(DISCO(p), DISCO(full)) -- convergence to complete score
# Nested selection (50 fixed + additional); 20 random replicates -> SE.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
d <- readRDS(file.path(root, "UKB_Data/UKB.Prot.clean.rds"))
X <- d$X; eid <- d$eid; n_all <- nrow(X); p_full <- ncol(X)
mort <- read.csv(file.path(root, "UKB_Data/UKB.Mortality.csv"))
mm <- match(eid, mort$id)
age <- mort$Age[mm]

disco_vec <- function(Ma, ref) {
  n <- ref$n; c1 <- (n-1)/n; c2 <- 1/(n+1)
  vapply(seq_len(nrow(Ma)), function(i) {
    dd <- Ma[i, ] - ref$mu
    Snew <- c1 * ref$S + c2 * (dd %*% t(dd))
    snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
    log(n^2 * sum((ref$Rref - Rnew)^2))
  }, numeric(1))
}
set_ref <- function(Xr) list(n = nrow(Xr), mu = colMeans(Xr), S = cov(Xr), Rref = cor(Xr))

ref_age <- which(age <= 40 & !is.na(age))
ok <- which(!is.na(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 1000)

# reference scores: DISCO(50) and DISCO(full)
set.seed(1); ref_sel <- sample(p_full, 50)
DISCO_ref50 <- disco_vec(X[analysis, ref_sel, drop = FALSE], set_ref(X[ref_age, ref_sel, drop = FALSE]))
DISCO_full <- disco_vec(X[analysis, , drop = FALSE], set_ref(X[ref_age, , drop = FALSE]))

PGRID <- c(50, 200, 500, 1000, 2000, 2920)
res <- list()
for (p in PGRID) {
  if (p == 50) { sels <- list(ref_sel) }
  else if (p == 2920) { sels <- list(seq_len(p_full)) }
  else {
    sels <- lapply(1:20, function(rep) { set.seed(5000 + p + rep); c(ref_sel, sample(setdiff(seq_len(p_full), ref_sel), p - 50)) })
  }
  r50 <- rfull <- numeric(length(sels))
  for (k in seq_along(sels)) {
    sel <- sels[[k]]
    dp <- disco_vec(X[analysis, sel, drop = FALSE], set_ref(X[ref_age, sel, drop = FALSE]))
    r50[k]   <- cor(dp, DISCO_ref50, method = "spearman")
    rfull[k] <- cor(dp, DISCO_full, method = "spearman")
  }
  se50 <- if (length(sels) > 1) sd(r50) else 0
  sefull <- if (length(sels) > 1) sd(rfull) else 0
  res[[length(res)+1]] <- data.frame(p = p, ref = "DISCO(50)",  rho = mean(r50), se = se50)
  res[[length(res)+1]] <- data.frame(p = p, ref = "DISCO(full)", rho = mean(rfull), se = sefull)
  cat(sprintf("p=%d vs50=%.4f vsfull=%.4f\n", p, mean(r50), mean(rfull)))
}
res <- do.call(rbind, res)
res$p <- factor(res$p, levels = PGRID)
res$ref <- factor(res$ref, levels = c("DISCO(50)", "DISCO(full)"))
saveRDS(res, file.path(out, "R43_score_stability.rds"))

pfig <- ggplot(res, aes(p, rho, colour = ref, group = ref)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se), width = 0.08) +
  coord_cartesian(ylim = c(0.70, 1.02)) +
  scale_colour_manual(values = c("DISCO(50)" = "#E64B35FF", "DISCO(full)" = "#4DBBD5FF")) +
  theme_classic(base_size = 13) +
  labs(x = "Number of features (p)", y = "Spearman correlation", colour = NULL,
       title = "DISCO score stability vs p")
ggsave(file.path(out, "R43_score_stability.pdf"), pfig, width = 6, height = 4.5)
cat("DONE\n")
