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

# Visualize the "continuous co-expression structure captured by DISCO" argument (no clustering).
# A: eigenvalue spectrum (continuous). B: aging deviation mode concentration. C: aging proteins (biology).
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age

tt <- theme(plot.title = element_text(size = 14, face = "bold"))
ax <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11, colour = "black"))

young <- which(age <= 40); old <- which(age >= 65)
Rref <- cor(X[young, , drop = FALSE])
eig <- eigen(Rref, symmetric = TRUE)
lam <- eig$values; U <- eig$vectors

# ---- A: eigenvalue spectrum (continuous) ----
dfA <- data.frame(k = 1:length(lam), lambda = lam)
pA <- ggplot(dfA, aes(k, lambda)) +
  geom_line(colour = "#4DBBD5FF") + geom_point(size = 0.6, colour = "#4DBBD5FF") +
  scale_y_log10() +
  theme_classic() + tt + ax +
  labs(x = "Co-expression mode", y = "Eigenvalue (log)",
       title = "A. Continuous co-expression spectrum") +
  annotate("text", x = 250, y = max(lam), size = 3.5, hjust = 0,
           label = "global factor\nthen continuous tail\n(no discrete modules)")

# ---- B: aging deviation mode concentration ----
delta <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
delta <- delta / sqrt(sum(delta^2))
c2 <- (t(U) %*% delta)^2; c2 <- c2 / sum(c2)
cum_real <- cumsum(c2)
set.seed(1); cum_null <- matrix(NA, 200, length(lam))
for (rep in 1:200) { dn <- sample(delta); cn <- (t(U) %*% dn)^2; cn <- cn/sum(cn); cum_null[rep,] <- cumsum(cn) }
dfB <- data.frame(k = 1:length(lam), real = cum_real, null = colMeans(cum_null))
pB <- ggplot(dfB, aes(k)) +
  geom_line(aes(y = real, colour = "Aging deviation")) +
  geom_line(aes(y = null, colour = "Random")) +
  scale_colour_manual(values = c("Aging deviation" = "#E64B35FF", "Random" = "#4DBBD5FF")) +
  coord_cartesian(xlim = c(0, 80)) +
  theme_classic() + tt + ax +
  theme(legend.position = c(0.75, 0.3), legend.background = element_blank()) +
  labs(x = "Top-K modes", y = "Cumulative deviation explained", colour = NULL,
       title = "B. Aging signal concentrated in few modes") +
  annotate("text", x = 45, y = 0.4, size = 3.5, label = sprintf("~%.0f modes vs %.0f random", 1/sum(c2^2), length(lam)))

# ---- C: aging proteins (biology) ----
bio_cat <- c("GDF15"="Inflammation","IL6"="Inflammation","CXCL9"="Inflammation","MMP12"="Inflammation",
             "CCL19"="Inflammation","TNFRSF10B"="Inflammation","TNFRSF11B"="Inflammation","ADM"="Inflammation",
             "NTproBNP"="Cardiac","NPPB"="Cardiac",
             "NEFL"="Neuro","GFAP"="Neuro",
             "HAVCR1"="Kidney","EDA2R"="Aging")
dlt <- colMeans(X[old, , drop = FALSE]) - colMeans(X[young, , drop = FALSE])
top <- names(sort(abs(dlt), decreasing = TRUE))[1:20]
dfC <- data.frame(protein = factor(top, levels = rev(top)), delta = dlt[top],
                  cat = ifelse(top %in% names(bio_cat), bio_cat[top], "Other"))
colpal <- c("Inflammation"="#E64B35FF","Cardiac"="#00A087FF","Neuro"="#F39B7FFF",
            "Kidney"="#4DBBD5FF","Aging"="#8491B4FF","Other"="grey80")
pC <- ggplot(dfC, aes(protein, delta, fill = cat)) +
  geom_col(width = 0.7, colour = "black") + coord_flip() +
  scale_fill_manual(values = colpal, name = "Biology") +
  theme_classic() + tt + ax +
  labs(x = NULL, y = "Aging deviation (old - young)", title = "C. Aging modes = inflammation / cardiac / neuro")

fig <- (pA | pB) / pC
ggsave(file.path(out, "R4_fig2_modules.pdf"), fig, width = 11, height = 9)
cat("wrote R4_fig2_modules.pdf\n")
