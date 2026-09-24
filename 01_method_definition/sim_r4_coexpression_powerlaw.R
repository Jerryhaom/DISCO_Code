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

# Co-expression network structure: degree distribution (power-law?), correlation-value density.
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age; p <- ncol(X)

young <- which(age <= 40 & !is.na(age))
Rref <- cor(X[young, , drop = FALSE])

# ---- correlation value density (off-diagonal |r|) ----
off <- abs(Rref[lower.tri(Rref)])
off <- off[off > 0]

# ---- degree distribution at threshold tau ----
degree_dist <- function(tau) {
  A <- abs(Rref) > tau; diag(A) <- 0
  k <- colSums(A)
  tab <- table(k)
  data.frame(k = as.numeric(names(tab)), count = as.numeric(tab), tau = tau)
}

# pick tau = top 10% of |r| as edges (sparse network)
tau10 <- quantile(off, 0.90)
cat("tau (top 10% |r|) =", round(tau10, 3), "\n")
dd <- degree_dist(tau10)
# log-log linear fit (power-law exponent)
pos <- dd$k > 0
fit <- lm(log(dd$count[pos]) ~ log(dd$k[pos]))
gamma <- -coef(fit)[2]
cat("power-law exponent gamma (log-log slope) =", round(gamma, 3), "\n")

saveRDS(list(off = off, dd = dd, tau = tau10, gamma = gamma),
        file.path(out, "R4_coexpression_powerlaw.rds"))

# figures
p1 <- ggplot(dd, aes(k, count)) +
  geom_point(colour = "#4DBBD5FF") +
  scale_x_log10() + scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(x = "Degree k", y = "Count (log)", title = "Co-expression degree distribution") +
  annotate("text", x = max(dd$k)/3, y = max(dd$count)/1.5, size = 4,
           label = sprintf("slope = %.2f", -gamma))
ggsave(file.path(out, "R4_degree_distribution.pdf"), p1, width = 6, height = 4.5)

p2 <- ggplot(data.frame(r = off), aes(r)) +
  geom_histogram(bins = 60, fill = "#E64B35FF", colour = "white") +
  scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(x = "|correlation| (off-diagonal)", y = "Count (log)",
       title = "Co-expression strength distribution")
ggsave(file.path(out, "R4_correlation_density.pdf"), p2, width = 6, height = 4.5)
cat("DONE\n")
