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

# R05 visualization (2/3): DISCO HR stratified by PRS Q1 vs Q4.

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

wkd <- disco_results_path("08_revision_analyses")
out <- file.path(wkd, "outputs/R05_genetic_interaction/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(wkd, "outputs/R05_genetic_interaction/R05_stratified_disco_hr.csv"),
              stringsAsFactors = FALSE)

lbl <- c(
  amd = "Age-related macular degeneration", ad = "Alzheimer's disease",
  asthma = "Asthma", af = "Atrial fibrillation", bipolar = "Bipolar disorder",
  bowel_ca = "Bowel cancer", breast_ca = "Breast cancer", cvd = "Cardiovascular disease",
  coeliac = "Coeliac disease", cad = "Coronary artery disease", crohn = "Crohn's disease",
  ovarian_ca = "Ovarian cancer", hypertension = "Hypertension",
  ischaemic_stroke = "Ischaemic stroke", melanoma = "Melanoma", ms = "Multiple sclerosis",
  osteoporosis = "Osteoporosis", parkinson = "Parkinson's disease",
  glaucoma = "Glaucoma", prostate_ca = "Prostate cancer", psoriasis = "Psoriasis",
  ra = "Rheumatoid arthritis", schizophrenia = "Schizophrenia", sle = "Systemic lupus",
  t1d = "Type 1 diabetes", t2d = "Type 2 diabetes",
  ulcerative_colitis = "Ulcerative colitis", vte = "Venous thromboembolism")

d$label <- lbl[d$disease]
d$stratum <- factor(d$stratum, levels = c("Q1", "Q4"))
d$version <- ifelse(d$version == "clinical", "Clinical DISCO (N~480k)",
                    "Proteomics DISCO (N~47k)")

# order diseases by Q1 DISCO HR (descending) for readability
ord <- d %>% filter(stratum == "Q1", version == levels(factor(d$version))[1]) %>%
  arrange(desc(disco_hr)) %>% pull(label)
d$label <- factor(d$label, levels = ord)

gg <- ggplot(d, aes(x = disco_hr, y = label, colour = stratum)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = disco_lo, xmax = disco_hi), height = 0.6,
                 position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_point(aes(shape = stratum), size = 1.8, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("Q1" = "#4DBBD5", "Q4" = "#E64B35"),
                      labels = c("Q1" = "Low PRS (Q1)", "Q4" = "High PRS (Q4)")) +
  scale_shape_manual(values = c("Q1" = 16, "Q4" = 17),
                     labels = c("Q1" = "Low PRS (Q1)", "Q4" = "High PRS (Q4)")) +
  scale_x_log10(limits = c(0.5, 6), breaks = c(0.5, 1, 2, 3, 6)) +
  facet_wrap(~ version, ncol = 1, scales = "free_y") +
  labs(x = "DISCO HR per SD (95% CI), within PRS stratum",
       y = NULL, colour = NULL, shape = NULL,
       title = "DISCO-disease association by low vs high genetic risk (PRS Q1 vs Q4)") +
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.position = "top",
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))

ggsave(file.path(out, "R05_stratified_forest.pdf"), gg, width = 6.5, height = 9)
ggsave(file.path(out, "R05_stratified_forest.png"), gg, width = 6.5, height = 9, dpi = 300)
cat("saved stratified forest plot\n")
