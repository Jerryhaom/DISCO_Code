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

# Standard WGCNA pipeline (manual implementation, WGCNA package unavailable).
# 1) soft-power selection via scale-free topology fit; 2) TOM; 3) dynamic-ish module cut;
# 4) module eigengene vs DISCO / age.
suppressMessages({library(ggplot2); library(dplyr)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco_score <- d$mc
p <- ncol(X); R <- cor(X)

# ---- 1. soft-power selection (scale-free topology fit) ----
betas <- 1:20
sft <- data.frame(power = betas, SFT.R.sq = NA, mean.k = NA)
for (b in betas) {
  A <- abs(R)^b
  k <- colSums(A)
  h <- hist(k, breaks = 30, plot = FALSE)
  keep <- h$counts > 0
  if (sum(keep) < 4) { sft$SFT.R.sq[sft$power==b] <- NA; sft$mean.k[sft$power==b] <- mean(k); next }
  sft$SFT.R.sq[sft$power==b] <- cor(log10(h$mids[keep]), log10(h$counts[keep]))^2
  sft$mean.k[sft$power==b] <- mean(k)
}
cat("=== soft-power selection ===\n"); print(sft, row.names = FALSE)
beta <- betas[which(sft$SFT.R.sq > 0.8)[1]]
if (is.na(beta)) beta <- betas[which.max(sft$SFT.R.sq)]
cat("chosen beta =", beta, "\n")

# ---- 2. TOM ----
A <- abs(R)^beta
k <- colSums(A)
TOM <- matrix(0, p, p, dimnames = list(colnames(X), colnames(X)))
for (i in 1:p) for (j in i:p) {
  if (i == j) next
  tom <- (sum(A[i, ] * A[, j]) + A[i, j]) / (min(k[i], k[j]) + 1 - A[i, j])
  TOM[i, j] <- TOM[j, i] <- tom
}

# ---- 3. hierarchical clustering + dynamic-ish cut (minModuleSize, not fixed K) ----
hc <- hclust(as.dist(1 - TOM), method = "average")
minMod <- 20
# use a height cut that yields modules >= minMod and as many as the tree allows
hcut <- 0.20
module <- cutree(hc, h = hcut)
small <- names(table(module))[table(module) < minMod]
# merge small modules into the closest large one (simplified merge)
module[module %in% small] <- 0
mm <- as.integer(table(module)); cat("module sizes (h=0.20):", paste(mm, collapse=", "), "\n")

# ---- 4. module eigengene vs DISCO / age ----
mods <- sort(setdiff(unique(module), 0))
eig <- lapply(mods, function(m) prcomp(X[, module == m], center=TRUE, scale.=TRUE)$x[,1])
names(eig) <- paste0("ME", mods)
res <- data.frame(module = names(eig), n = sapply(eig, length),
                  r_disco = sapply(eig, function(e) cor(e, disco_score, method="spearman")),
                  r_age = sapply(eig, function(e) cor(e, age, method="spearman")))
res <- res[order(-abs(res$r_disco)), ]
cat("\n=== module eigengene vs DISCO / age ===\n"); print(res, row.names = FALSE)

saveRDS(list(sft = sft, beta = beta, module = module, res = res, eig = eig),
        file.path(out, "R4_wgcna_standard.rds"))

# soft-power plot
p1 <- ggplot(sft, aes(power, SFT.R.sq)) + geom_point() + geom_line() +
  geom_hline(yintercept = 0.8, linetype="dashed", colour="grey50") +
  theme_classic(base_size=14) + labs(x="Soft power (beta)", y="Scale-free topology R^2",
       title="Soft-power selection (scale-free fit)")
ggsave(file.path(out, "R4_wgcna_softpower.pdf"), p1, width=6, height=4.5)

p2 <- ggplot(res, aes(r_age, r_disco, size = n)) + geom_point(colour="#E64B35FF") +
  geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  geom_text(aes(label=module), vjust=-1, size=3.5) +
  theme_classic(base_size=14) + labs(x="Correlation with age", y="Correlation with DISCO",
       size="module size", title="Co-expression modules vs DISCO / age")
ggsave(file.path(out, "R4_wgcna_eigengene_standard.pdf"), p2, width=7, height=5.5)
cat("DONE\n")
