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

# Block-wise: does the cross (clinical x proteomics) block add INDEPENDENT predictive power
# beyond the within blocks? (LRT + C-index increment), crude and age+sex.
suppressMessages({library(survival)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
clin <- readRDS(file.path(root, "Result1/Clean.UKB.Clinical.rds"))
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))

m <- match(prot$eid, clin$sampleID); keep <- which(!is.na(m))
Xc <- scale(as.matrix(clin[m[keep], -1])); Xp <- scale(as.matrix(prot[keep, -1]))
X <- cbind(Xc, Xp); nc <- ncol(Xc); np <- ncol(Xp)
age <- d$age[keep]; sex <- d$sex[keep]; time <- d$time[keep]; status <- d$status[keep]

ref_age <- which(age <= 40 & !is.na(age))
Xref <- X[ref_age, , drop = FALSE]; n <- nrow(Xref)
mu <- colMeans(Xref); S <- cov(Xref); Rref <- cor(Xref)
c1 <- (n-1)/n; c2 <- 1/(n+1)
ic <- 1:nc; ip <- (nc+1):(nc+np)

ok <- which(!is.na(time) & !is.na(status) & is.finite(age)); ok <- setdiff(ok, ref_age)
set.seed(42); analysis <- sample(ok, 5000)

sw <- sc <- numeric(length(analysis))
for (i in seq_along(analysis)) {
  dd <- X[analysis[i], ] - mu
  Snew <- c1 * S + c2 * (dd %*% t(dd))
  snew <- sqrt(diag(Snew)); Rnew <- Snew / (snew %*% t(snew))
  Diff <- (Rref - Rnew)^2
  sw[i] <- sum(Diff[ic, ic]) + sum(Diff[ip, ip])     # within total
  sc[i] <- 2 * sum(Diff[ic, ip])                      # cross
}
disco_within <- log(n^2 * sw)
disco_cross  <- log(n^2 * sc)

dd0 <- data.frame(tm = time[analysis], st = status[analysis], age = age[analysis], sex = sex[analysis])

# crude: does cross add beyond within?
A_c <- coxph(Surv(tm, st) ~ disco_within, data = dd0)
B_c <- coxph(Surv(tm, st) ~ disco_within + disco_cross, data = dd0)
lrt_c <- anova(A_c, B_c)$`Pr(>|Chi|)`[2]
# age+sex: does cross add beyond within?
A_a <- coxph(Surv(tm, st) ~ disco_within + age + sex, data = dd0)
B_a <- coxph(Surv(tm, st) ~ disco_within + disco_cross + age + sex, data = dd0)
lrt_a <- anova(A_a, B_a)$`Pr(>|Chi|)`[2]

cat(sprintf("crude:   within alone C=%.4f -> within+cross C=%.4f (LRT p=%.2e)\n",
            summary(A_c)$concordance[1], summary(B_c)$concordance[1], lrt_c))
cat(sprintf("age+sex: within alone C=%.4f -> within+cross C=%.4f (LRT p=%.2e)\n",
            summary(A_a)$concordance[1], summary(B_a)$concordance[1], lrt_a))
cat(sprintf("\ncross alone C=%.4f (crude), %.4f (age+sex)\n",
            summary(coxph(Surv(tm,st)~disco_cross, dd0))$concordance[1],
            summary(coxph(Surv(tm,st)~disco_cross+age+sex, dd0))$concordance[1]))
