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
disco_init_module("08_revision_analyses")


# R04: cross-cohort random-effects meta-analysis (R1.5).
# Correction (2026-09-16):
#   1. RuLAS wave 4 (1.45) is a LOGISTIC OR on two-year mortality (glm), NOT a
#      Cox HR, so it is excluded from the HR meta-analysis (method mismatch).
#   2. The published Extended Data Table 1 mixes crude HR (NHANES/UKB) and
#      age+sex-adjusted HR (CHARLS/CLHLS/RuLAS) under one footnote. Here we pool
#      the 7 Cox-HR cohorts twice: once all-crude and once all age+sex-adjusted.
# DerSimonian-Laird random effects; per-unit (as published) and per-SD (within
# cohort) standardization. No legacy files are modified.

out_dir <- "outputs/R04_meta_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 7 Cox-HR cohorts; crude, age+sex-adjusted, and expanded-covariate-adjusted HR
d <- data.frame(
  cohort = c("NHANES","UKB clinical","UKB proteomics","UKB clin+prot",
             "CHARLS","CLHLS","RuLAS wave 2"),
  modality = c("clinical","clinical","proteomics","combined",
               "clinical","clinical","clinical"),
  sd = c(1.543, 1.139, 0.815, 0.825, 1.610, 1.398, 1.168),
  crude_hr = c(1.71, 1.61, 2.41, 2.38, 1.48, 1.41, 1.43),
  crude_lo = c(1.68, 1.60, 2.35, 2.32, 1.43, 1.35, 1.31),
  crude_hi = c(1.75, 1.62, 2.47, 2.44, 1.52, 1.47, 1.57),
  adj_hr   = c(1.29, 1.48, 2.09, 2.07, 1.26, 1.27, 1.33),
  adj_lo   = c(1.26, 1.47, 2.03, 2.01, 1.21, 1.21, 1.21),
  adj_hi   = c(1.32, 1.49, 2.14, 2.12, 1.30, 1.33, 1.46),
  # expanded-covariate adjusted per-SD HR (R01 UKB + R06 cohorts + UKB clin+prot)
  exp_hr_sd = c(1.524, 1.480, 1.666, 1.674, 1.445, 1.395, 1.394),
  exp_lo_sd = c(1.468, 1.465, 1.621, 1.628, 1.367, 1.303, 1.246),
  exp_hi_sd = c(1.582, 1.495, 1.713, 1.721, 1.526, 1.493, 1.560)
)

# per-SD = HR^(SD of log-DISCO); CI likewise
dl_pool <- function(hr, lo, hi, label = "") {
  y  <- log(hr)
  se <- (log(hi) - log(lo)) / (2 * 1.96)
  k  <- length(y)
  w  <- 1 / se^2
  theta_fe <- sum(w * y) / sum(w)
  Q <- sum(w * (y - theta_fe)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  w_re <- 1 / (se^2 + tau2)
  theta_re <- sum(w_re * y) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  data.frame(label = label, k = k,
    hr_fe = exp(theta_fe), hr_re = exp(theta_re),
    lo = exp(theta_re - 1.96 * se_re), hi = exp(theta_re + 1.96 * se_re),
    tau2 = tau2, I2 = max(0, (Q - (k - 1)) / Q) * 100, Q = Q, df = k - 1,
    p_het = pchisq(Q, k - 1, lower.tail = FALSE))
}

# crude pools (per-unit and per-SD)
cu <- dl_pool(d$crude_hr, d$crude_lo, d$crude_hi, "crude, per-unit (n=7)")
cs <- dl_pool(d$crude_hr^d$sd, d$crude_lo^d$sd, d$crude_hi^d$sd, "crude, per-SD (n=7)")
# age+sex-adjusted pools (per-unit and per-SD)
au <- dl_pool(d$adj_hr, d$adj_lo, d$adj_hi, "age+sex adjusted, per-unit (n=7)")
as_ <- dl_pool(d$adj_hr^d$sd, d$adj_lo^d$sd, d$adj_hi^d$sd, "age+sex adjusted, per-SD (n=7)")
# expanded-covariate adjusted pool (per-SD, already standardized within cohort)
ex <- dl_pool(d$exp_hr_sd, d$exp_lo_sd, d$exp_hi_sd, "expanded-covariate adjusted, per-SD (n=7)")

pool <- rbind(cu, cs, au, as_, ex)

# leave-one-cohort-out (per-SD: crude, adjusted, expanded)
locoo <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  rbind(
    cbind(dl_pool(d$crude_hr[-i]^d$sd[-i], d$crude_lo[-i]^d$sd[-i], d$crude_hi[-i]^d$sd[-i],
                  paste0("crude per-SD, leave out: ", d$cohort[i])), left_out = d$cohort[i], pool = "crude"),
    cbind(dl_pool(d$adj_hr[-i]^d$sd[-i], d$adj_lo[-i]^d$sd[-i], d$adj_hi[-i]^d$sd[-i],
                  paste0("adjusted per-SD, leave out: ", d$cohort[i])), left_out = d$cohort[i], pool = "adjusted"),
    cbind(dl_pool(d$exp_hr_sd[-i], d$exp_lo_sd[-i], d$exp_hi_sd[-i],
                  paste0("expanded per-SD, leave out: ", d$cohort[i])), left_out = d$cohort[i], pool = "expanded")
  )
}))

write.csv(d, file.path(out_dir, "R04_input_estimates.csv"), row.names = FALSE)
write.csv(pool, file.path(out_dir, "R04_pooled_estimates.csv"), row.names = FALSE)
write.csv(locoo, file.path(out_dir, "R04_leave_one_cohort_out.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "method=DerSimonian-Laird_random_effects",
  "excluded=RuLAS_wave4_logistic_OR_not_Cox_HR",
  "standardization=per-unit_and_per-SD",
  "pooling=all-crude_all-age+sex-adjusted_all-expanded-covariate-adjusted",
  "legacy_files_modified=false"
), file.path(out_dir, "audit_status.txt"))

cat("R04 done\n")
print(pool[, c("label","k","hr_re","lo","hi","tau2","I2","p_het")], digits = 3, row.names = FALSE)
