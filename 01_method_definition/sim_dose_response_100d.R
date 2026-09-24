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

# 100-dim nonlinearity dose-response (final, FLI1-AKT3 pair within 100 proteins).
# Report: (a) Pearson cor(driver, injected) vs alpha, (b) DISCO deviation vs alpha.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(ggthemes); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
exp  <- readRDS(file.path(out, "proteomics_explore.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X    <- as.matrix(prot[, -1]); colnames(X) <- exp$names
age  <- d$age; ref_idx <- exp$ref_idx

set.seed(7); sel100 <- sample(colnames(X), 100)
if (!"FLI1" %in% sel100) sel100[1] <- "FLI1"
if (!"AKT3" %in% sel100) sel100[2] <- "AKT3"
M <- scale(X[, sel100])
ref <- M[ref_idx, , drop = FALSE]; nref <- nrow(ref); Rc <- cor(ref)

DRIVER <- "FLI1"; TARGET <- "AKT3"
cat("driver-target cor (young ref) =", round(cor(ref[, DRIVER], ref[, TARGET]), 3), "\n")
invnorm <- function(v) qnorm((rank(v) - 0.5) / length(v))

disco <- function(Mx) {
  vapply(seq_len(nrow(Mx)), function(i) {
    ca <- cor(rbind(ref, Mx[i, , drop = FALSE]))
    log(sum((Rc - ca)^2, na.rm = TRUE) * nref^2)
  }, numeric(1))
}

alphas <- seq(0, 1, 0.1); N <- 400
res <- list()
set.seed(20260918)
for (ty in c("monotone","nonmonotone")) {
  for (a in alphas) {
    idx <- sample(ref_idx, N, replace = TRUE)
    Md <- M[idx, , drop = FALSE]
    F <- invnorm(Md[, DRIVER])           # symmetric factor, monotone in FLI1
    if (ty == "monotone") nl <- F + 0.6 * F^3 + 0.3 * rnorm(N)
    else                  nl <- F^2 - 1 + 0.3 * rnorm(N)
    nl <- as.numeric(scale(nl))
    Md[, TARGET] <- (1 - a) * Md[, TARGET] + a * nl
    res[[length(res)+1]] <- data.frame(
      type = ty, alpha = a,
      cor_target = cor(Md[, DRIVER], Md[, TARGET]),
      disco = median(disco(Md))
    )
  }
}
res <- do.call(rbind, res)
base <- res %>% group_by(type) %>% summarise(b = disco[alpha == 0])
res <- merge(res, base, by = "type"); res$dev <- res$disco - res$b
res$type <- factor(res$type, levels = c("monotone","nonmonotone"),
                   labels = c("Monotone nonlinear (cubic)", "Non-monotone nonlinear (U-shape)"))
saveRDS(res, file.path(out, "dose_response_100d.rds"))
print(res[, c("type","alpha","cor_target","dev")])

# ---- figure ----
p1 <- ggplot(res, aes(alpha, cor_target, colour = type)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  scale_colour_manual(values = c("Monotone nonlinear (cubic)" = "#E64B35FF",
                                 "Non-monotone nonlinear (U-shape)" = "#4DBBD5FF")) +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "Nonlinearity strength", y = "Pearson cor (FLI1, AKT3)",
       title = "100-dim: correlation vs nonlinearity strength")

p2 <- ggplot(res, aes(alpha, dev, colour = type)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  scale_colour_manual(values = c("Monotone nonlinear (cubic)" = "#E64B35FF",
                                 "Non-monotone nonlinear (U-shape)" = "#4DBBD5FF")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  labs(x = "Nonlinearity strength", y = "DISCO deviation from linear",
       title = "100-dim: DISCO detection")

p <- p1 / p2 + plot_annotation(
  title = "R2.4 (100 proteins): DISCO misses monotone nonlinear, flags non-monotone",
  subtitle = "Monotone nonlinear keeps Pearson cor high -> DISCO does not flag. Non-monotone collapses cor -> DISCO flags (because correlation vanished).")
ggsave(file.path(out, "R24_dose_response_100d.pdf"), p, width = 7, height = 6.5)
cat("wrote R24_dose_response_100d.pdf\n")
