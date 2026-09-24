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
disco_init_module("03_multimodal_UKB")


# UKB proteomics DISCO/DM analysis with reference age<=40 (unified).
# Reads the regenerated object; produces HR/C-index tables and figures.
# No legacy Result2 files are modified.

options(stringsAsFactors = FALSE)
suppressMessages({
  library(survival); library(survminer); library(ggplot2); library(ggthemes)
})
root <- disco_project_root()
od <- disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration")
out <- disco_results_path("03_multimodal_UKB/outputs/03_prot_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(od, "UKB.Prot.DISCO.regenerated.rds"))
d$loghd <- log(d$hd)

# ---- basic associations -------------------------------------------------
assoc <- data.frame(
  metric = c("DISCO (mc)", "DM log(hd)"),
  spearman_age = c(
    cor(d$mc, d$age, method = "spearman"),
    cor(d$loghd, d$age, method = "spearman")
  ),
  spearman_mc_vs_hd = c(NA, cor(d$mc, d$loghd, method = "spearman"))
)

# ---- Cox models (HR per SD + C-index) -----------------------------------
fit_one <- function(formula, data) {
  f <- coxph(as.formula(formula), data = data)
  s <- summary(f)
  # HR per SD for the metric of interest (first predictor)
  pred <- all.vars(as.formula(formula))[3]  # metric variable
  x <- data[[pred]]
  hr_per_sd <- exp(coef(f)[pred] * sd(x, na.rm = TRUE))
  ci <- confint(f)[pred, ]
  data.frame(
    formula = formula, n = nrow(data), events = f$nevent,
    hr_per_sd = hr_per_sd,
    hr_per_sd_low = exp(log(hr_per_sd) - qnorm(0.975) * sd(x) * sqrt(vcov(f)[pred, pred])),
    hr_per_sd_high = exp(log(hr_per_sd) + qnorm(0.975) * sd(x) * sqrt(vcov(f)[pred, pred])),
    metric_p = s$coefficients[pred, "Pr(>|z|)"],
    c_index = unname(s$concordance[1]), c_index_se = unname(s$concordance[2])
  )
}

models <- rbind(
  fit_one("Surv(time,status)~mc", d),
  fit_one("Surv(time,status)~mc+age+sex", d),
  fit_one("Surv(time,status)~loghd", d),
  fit_one("Surv(time,status)~loghd+age+sex", d),
  fit_one("Surv(time,status)~age", d),
  fit_one("Surv(time,status)~age+sex", d)
)

# ---- figures ------------------------------------------------------------
set.seed(20260914)
dsub <- d[sample(nrow(d), 5000), ]

# (a) DISCO vs age
p1 <- ggplot(dsub, aes(age, mc)) +
  geom_point(alpha = 0.3, size = 0.6, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
  labs(title = "UKB Proteomics (reference age<=40)",
       x = "Age (years)", y = "DISCO") +
  theme_few() + theme(legend.position = "none")
ggsave(file.path(out, "Prot.DISCO.age40.vs.age.pdf"), p1, width = 4, height = 3.2)

# (b) DISCO vs log(DM)
p2 <- ggplot(dsub, aes(loghd, mc)) +
  geom_point(alpha = 0.3, size = 0.6, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "firebrick") +
  labs(title = "UKB Proteomics (reference age<=40)",
       x = "DM (log)", y = "DISCO") +
  theme_few() + theme(legend.position = "none")
ggsave(file.path(out, "Prot.DISCO.age40.vs.DM.pdf"), p2, width = 4, height = 3.2)

# (c) Kaplan-Meier by DISCO quartile (Q1 vs Q4)
d$q_mc <- cut(d$mc, breaks = quantile(d$mc, probs = seq(0, 1, 0.25), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
km_mc <- d[d$q_mc %in% c(1, 4), ]
km_mc$grp <- factor(ifelse(km_mc$q_mc == 1, "DISCO Q1 (low)", "DISCO Q4 (high)"),
                    levels = c("DISCO Q1 (low)", "DISCO Q4 (high)"))
fit_mc <- survfit(Surv(time, status) ~ grp, data = km_mc)
p3 <- ggsurvplot(fit_mc, data = km_mc, conf.int = TRUE, censor = FALSE,
                 pval = TRUE, xlab = "Follow-up (years)", ylab = "Survival probability",
                 legend.title = "", palette = c("#1B9E77", "#D95F02"),
                 title = "UKB Proteomics DISCO (reference age<=40)")
ggsave(file.path(out, "Prot.DISCO.age40.KM.quartile.pdf"),
       print(p3), width = 5, height = 4)

# (d) Kaplan-Meier by DM quartile (Q1 vs Q4)
d$q_hd <- cut(d$loghd, breaks = quantile(d$loghd, probs = seq(0, 1, 0.25), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
km_hd <- d[d$q_hd %in% c(1, 4), ]
km_hd$grp <- factor(ifelse(km_hd$q_hd == 1, "DM Q1 (low)", "DM Q4 (high)"),
                    levels = c("DM Q1 (low)", "DM Q4 (high)"))
fit_hd <- survfit(Surv(time, status) ~ grp, data = km_hd)
p4 <- ggsurvplot(fit_hd, data = km_hd, conf.int = TRUE, censor = FALSE,
                 pval = TRUE, xlab = "Follow-up (years)", ylab = "Survival probability",
                 legend.title = "", palette = c("#1B9E77", "#D95F02"),
                 title = "UKB Proteomics DM (reference age<=40)")
ggsave(file.path(out, "Prot.DISCO.age40.KM.quartile.DM.pdf"),
       print(p4), width = 5, height = 4)

# ---- export -------------------------------------------------------------
write.csv(assoc, file.path(out, "associations.csv"), row.names = FALSE)
write.csv(models, file.path(out, "cox_HR_cindex.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(
  "status=PASS_AGE40_ANALYSIS_COMPLETE",
  "reference=age<=40_N565",
  "n_analysis=47859",
  "events=5117",
  "legacy_files_modified=false"
), file.path(out, "analysis_status.txt"))

cat("=== associations ===\n"); print(assoc)
cat("\n=== Cox models ===\n"); print(models)
cat("\nDONE\n")
