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
disco_init_module("08_sensitivity_analyses")

# R05 visualization (1/3): DISCO x PRS interaction-effect forest plot.
# Uses the existing R05_genetic_interaction.csv (no re-fitting).

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

wkd <- disco_results_path("08_sensitivity_analyses")
out <- file.path(wkd, "outputs/R05_genetic_interaction/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(wkd, "outputs/R05_genetic_interaction/R05_genetic_interaction.csv"),
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
d$sig <- d$interaction_p < 0.05
d$version <- ifelse(d$version == "clinical", "Clinical DISCO (N~480k)",
                    "Proteomics DISCO (N~47k)")

# order diseases by interaction HR within each version
d <- d %>% group_by(version) %>%
  mutate(ord = rank(interaction_hr, ties.method = "first")) %>% ungroup()

gg <- ggplot(d, aes(x = interaction_hr, y = ord)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = interaction_lo, xmax = interaction_hi),
                 height = 0.25, linewidth = 0.5,
                 colour = ifelse(d$sig, "#E64B35", "grey55")) +
  geom_point(aes(colour = sig), size = 1.8) +
  scale_colour_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey45"),
                      guide = "none") +
  scale_y_continuous(breaks = d$ord, labels = d$label, expand = c(0.02, 0.02)) +
  scale_x_continuous(limits = c(0.85, 1.10), breaks = seq(0.85, 1.10, 0.05)) +
  facet_wrap(~ version, ncol = 1, scales = "free_y") +
  labs(x = "DISCO \u00d7 PRS interaction HR (95% CI)",
       y = NULL,
       title = "Genetic risk (PRS) modifies the DISCO-disease association") +
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))

ggsave(file.path(out, "R05_interaction_forest.pdf"), gg,
       width = 6.5, height = 9)
ggsave(file.path(out, "R05_interaction_forest.png"), gg,
       width = 6.5, height = 9, dpi = 300)
cat("saved interaction forest plot\n")
