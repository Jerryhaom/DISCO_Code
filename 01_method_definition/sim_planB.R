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

# Plan B: nonlinear age-trajectory injection using real UKB proteins.
# Some proteins change NONLINEARLY with age (quadratic term). Simulate "aged"
# individuals by moving young individuals' nonlinear-age proteins to old-age
# levels along their fitted nonlinear trajectory; check whether DISCO detects them.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(ggthemes)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
exp  <- readRDS(file.path(out, "proteomics_explore.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
X    <- as.matrix(prot[, -1]); colnames(X) <- exp$names
age  <- d$age

# nonlinear-age proteins (top by quadratic significance, with clear age signal)
sel <- c("DPP4","AREG","EGFR","TSPAN1","CCN5","REN","NTproBNP","WFDC2")
M <- X[, sel]
ref_idx <- which(age <= 40)
young_idx <- which(age <= 40)
old_idx   <- which(age >= 65)

# fit nonlinear (quadratic) age trajectory per protein on full data
fit_quad <- function(y) lm(y ~ age + I(age^2))
preds <- vapply(sel, function(p) {
  m <- fit_quad(M[, p]); predict(m, newdata = data.frame(age = 70))
}, numeric(1))
sdres <- vapply(sel, function(p) sd(residuals(fit_quad(M[, p]))), numeric(1))

# simulate "aged" individuals: take young individuals, shift nonlinear proteins to age-70 level + noise
set.seed(20260918)
base <- M[young_idx, , drop = FALSE]
sim_aged <- base
for (p in sel) sim_aged[, p] <- preds[p] + rnorm(nrow(base), 0, sdres[p])

# DISCO per-individual score against young reference (Pearson, uniform weight)
ref <- M[ref_idx, , drop = FALSE]
nref <- nrow(ref); Rc <- cor(ref)
disco <- function(Mx) {
  vapply(seq_len(nrow(Mx)), function(i) {
    ca <- cor(rbind(ref, Mx[i, , drop = FALSE]))
    log(sum((Rc - ca)^2, na.rm = TRUE) * nref^2)
  }, numeric(1))
}

g <- data.frame(
  group = c(rep("Young (real)", 400),
            rep("Simulated aged", 400),
            rep("Old (real)", 400)),
  disco = NA_real_
)
g$disco[g$group == "Young (real)"]     <- sample(disco(M[young_idx, , drop = FALSE]), 400)
g$disco[g$group == "Simulated aged"]    <- disco(sim_aged[sample(nrow(sim_aged), 400), , drop = FALSE])
g$disco[g$group == "Old (real)"]        <- sample(disco(M[old_idx, , drop = FALSE]), 400)
g$group <- factor(g$group, levels = c("Young (real)","Simulated aged","Old (real)"))

print(g %>% group_by(group) %>% summarise(median = median(disco), mean = mean(disco)))
saveRDS(g, file.path(out, "planB_results.rds"))

pB <- ggplot(g, aes(group, disco, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.size = 0.4) +
  scale_fill_manual(values = c("Young (real)" = "#4DBBD5FF",
                               "Simulated aged" = "#E64B35FF",
                               "Old (real)" = "#3C5488FF")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(colour = "black")) +
  labs(x = NULL, y = "DISCO score",
       title = "R2.4: DISCO detects nonlinear aging trajectories",
       subtitle = paste0("Proteins moved along nonlinear age trajectory: ",
                         paste(sel, collapse = ", ")))
ggsave(file.path(out, "R24_planB_nonlinear_age.pdf"), pB, width = 7, height = 5)
cat("wrote planB figure\n")
