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

# Prove: (1) reference proteome has co-expression structure; (2) dysregulation (aging
# deviation) is structured, NOT random. Therefore the difference DISCO computes
# already carries that structure.
suppressMessages({library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X <- as.matrix(prot[, -1]); X <- scale(X)
age <- d$age
p <- ncol(X)

young <- which(age <= 40 & !is.na(age))
old   <- which(age >= 65 & !is.na(age))
cat("young N =", length(young), " old N =", length(old), "\n")

# ---- Part 1: reference structure (eigen-spectrum of co-expression matrix) ----
Rref <- cor(X[young, , drop = FALSE])
ev <- eigen(Rref, symmetric = TRUE)$values
cat("reference effective rank (eigenvalues > 1):", sum(ev > 1), "/", p, "\n")
cat("top 5 eigenvalues:", round(ev[1:5], 0), "\n")

# ---- Part 2: dysregulation direction is structured, not random ----
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])  # aging direction
delta <- delta / sqrt(sum(delta^2))  # unit direction

# co-expression modes (eigenvectors of Rref)
U <- eigen(Rref, symmetric = TRUE)$vectors

# projection of delta onto each mode (squared, = fraction of variance explained)
c2 <- (t(U) %*% delta)^2          # c2[k] = (u_k' delta)^2
c2 <- c2 / sum(c2)                # normalize: fraction of ||delta||^2 in each mode
cum_real <- cumsum(c2)            # cumulative concentration in top-K modes

# null: shuffle delta (destroy structure), same projection; repeat 200 times
set.seed(1)
cum_null <- matrix(NA, nrow = 200, ncol = p)
for (rep in 1:200) {
  dnull <- sample(delta)          # random permutation
  c2n <- (t(U) %*% dnull)^2
  c2n <- c2n / sum(c2n)
  cum_null[rep, ] <- cumsum(c2n)
}
null_mean <- colMeans(cum_null)

# effective number of modes (participation ratio)
eff <- function(c2v) 1 / sum(c2v^2)
cat(sprintf("effective #modes: real = %.1f, null = %.1f (p=%d)\n",
            eff(c2), eff(rowMeans(t(apply(cum_null, 1, function(r) diff(c(0, r)))))), p))

# how many modes to reach 80% of variance?
k80_real <- which(cum_real >= 0.8)[1]
k80_null <- which(null_mean >= 0.8)[1]
cat(sprintf("modes to reach 80%%: real = %d, null = %d\n", k80_real, k80_null))

saveRDS(list(cum_real = cum_real, null_mean = null_mean, k80_real = k80_real,
             k80_null = k80_null, eff_real = eff(c2)),
        file.path(out, "R4_structure_proof.rds"))

# ---- figure ----
library(ggplot2)
df <- data.frame(k = 1:p, real = cum_real, null = null_mean)
pfig <- ggplot(df, aes(k)) +
  geom_line(aes(y = real, colour = "Real aging deviation")) +
  geom_line(aes(y = null, colour = "Shuffled (random)")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Real aging deviation" = "#E64B35FF",
                                 "Shuffled (random)" = "#4DBBD5FF")) +
  coord_cartesian(xlim = c(0, 60)) +
  theme_classic(base_size = 14) +
  labs(x = "Top-K co-expression modes", y = "Cumulative fraction of deviation",
       colour = NULL, title = "Aging dysregulation is structured, not random") +
  theme(legend.position = c(0.7, 0.3), legend.background = element_blank())
ggsave(file.path(out, "R4_structure_proof.pdf"), pfig, width = 6.5, height = 5)
cat("wrote R4_structure_proof.pdf\n")
