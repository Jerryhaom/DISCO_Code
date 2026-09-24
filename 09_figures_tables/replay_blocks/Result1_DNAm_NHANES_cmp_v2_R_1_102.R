library(ggplot2)
library(ggpubr)
library(survival)

setwd(disco_path("Result1/DNAm/"))

cindex_se <- function(formula, data) {
  fit <- coxph(as.formula(formula), data = data)
  s   <- summary(fit)
  c(cindex = s$concordance[1], se = s$concordance[2])
}

# ── Load data ─────────────────────────────────────────────────────────────────
dat <- read.csv(disco_path("Manuscript/PNAS/cindex_comparison.csv"),
                stringsAsFactors = FALSE)

# Metric selection and display labels
label_map <- c(
  mc            = "DISCO+Age",
  hd            = "DM+Age",
  DunedinPoAm   = "DunedinPACE",
  GrimAgeMort   = "GrimAge",
  PhenoAge      = "PhenoAge",
  HannumAge     = "Hannum",
  HorvathAge    = "Horvath",
  VidalBraloAge = "VidalBralo",
  LinAge        = "Lin",
  ZhangAge      = "Zhang",
  SkinBloodAge  = "SkinBlood",
  WeidnerAge    = "Weidner"
)

dat <- dat[dat$metric %in% names(label_map), ]
dat$label <- label_map[dat$metric]

# ── Age-only C-index ──────────────────────────────────────────────────────────
DNAm  <- readRDS("NHANES.DNAm.res.rds")
DNAm1 <- DNAm[!is.na(DNAm$HorvathAge), ]
ok    <- !is.na(DNAm1$status) & !is.na(DNAm1$time) & !is.na(DNAm1$age)
age_ci <- cindex_se("Surv(time, status) ~ age", DNAm1[ok, ])
cat(sprintf("NHANES Age-only: C=%.4f SE=%.4f\n", age_ci[1], age_ci[2]))

age_row <- data.frame(
  metric   = "age",    label    = "Age",
  adj_c    = age_ci[1], adj_se  = age_ci[2],
  accel_c  = NA,        accel_se = NA
)

# Combine: slim columns, add Age row
d <- rbind(
  data.frame(label    = dat$label,
             adj_c    = dat$adj_c,    adj_se   = dat$adj_se,
             accel_c  = dat$accel_c,  accel_se = dat$accel_se),
  data.frame(label    = age_row$label,
             adj_c    = age_row$adj_c,  adj_se   = age_row$adj_se,
             accel_c  = age_row$accel_c, accel_se = age_row$accel_se)
)

highlight <- c("DISCO+Age", "DM+Age")

# ── Panel 1: Original (adj C-index) ──────────────────────────────────────────
d1 <- d[order(d$adj_c), ]
d1$label <- factor(d1$label, levels = d1$label)
fill1 <- ifelse(d1$label %in% highlight, "#E64B35",
          ifelse(d1$label == "Age",       "#4DBBD5", "grey70"))
ref_adj <- d$adj_c[d$label == "DISCO+Age"]

p1 <- ggplot(d1, aes(adj_c, label)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.9, fill = fill1) +
  geom_errorbar(aes(xmin = adj_c - adj_se, xmax = adj_c + adj_se), width = 0.4) +
  coord_cartesian(xlim = c(0.70, 0.78)) +
  theme_classic() + xlab("C-index") + ylab("Metrics") +
  labs(title = "Original") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        plot.margin = margin(5, 5, 10, 5))

# ── Panel 2: Accel (age-residualised C-index) ─────────────────────────────────
d2 <- d[!is.na(d$accel_c), ]
d2$label <- gsub("\\+Age$", "", d2$label)   # remove "+Age" suffix for Accel panel
d2 <- d2[order(d2$accel_c), ]
d2$label <- factor(d2$label, levels = d2$label)
fill2 <- ifelse(d2$label %in% c("DISCO", "DM"), "#E64B35", "grey70")
ref_accel <- d$accel_c[d$label == "DISCO+Age"]

p2 <- ggplot(d2, aes(accel_c, label)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.9, fill = fill2) +
  geom_errorbar(aes(xmin = accel_c - accel_se, xmax = accel_c + accel_se), width = 0.4) +
  coord_cartesian(xlim = c(0.50, 0.62)) +
  theme_classic() + xlab("C-index") + ylab("Metrics") +
  labs(title = "Accel") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        plot.margin = margin(5, 10, 10, 5))

p_combined <- ggarrange(p1, NULL, p2,
                        widths = c(1, 0.05, 1), ncol = 3)

pdf("NHANES.cmp.v2.pdf", height = 5, width = 7)
print(annotate_figure(p_combined,
                      top = text_grob("NHANES 1999-2002", size = 14,
                                      hjust = 0, x = 0.02)))
dev.off()
cat("Saved NHANES.cmp.v2.pdf\n")
