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
disco_init_module("05_organ_network")


# Fig. 4a: expanded-covariate organ/tissue DISCO associations.
# Layout follows the submission figure: organs/tissues on x, outcomes on y,
# with mortality and incident chronic disease arranged as vertical facets.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

root <- disco_results_path("05_organ_network")
input_dir <- file.path(root, "outputs", "organ_expanded_covariates")
output <- file.path(input_dir, "Fig4a_expanded_facet_HR.pdf")

organs <- c(
  "brain", "esophagus", "intestine", "lymphoid tissue", "choroid plexus",
  "liver", "adipose tissue", "heart muscle", "lung", "bone marrow",
  "salivary gland", "kidney", "pancreas", "skin", "stomach",
  "skeletal muscle", "tongue", "retina"
)
mortality_order <- c(
  "All-cause", "Heart disease", "Respiratory disease", "Cancer",
  "Cerebrovascular disease", "Alzheimer's disease", "Others"
)
disease_order <- c(
  "Cancer", "Obesity", "Hypertension", "Diabetes", "MI", "Stroke", "COPD", "Dementia"
)

stars <- function(p) ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))

mortality <- read.csv(file.path(input_dir, "organ_cause_specific_expanded.csv"))
names(mortality)[names(mortality) == "outcome"] <- "outcome"
mortality$panel <- "Mortality"

disease <- read.csv(file.path(input_dir, "organ_chronic_disease_expanded.csv"))
names(disease)[names(disease) == "disease"] <- "outcome"
disease$panel <- "Incident chronic disease"

plot_data <- rbind(
  mortality[, c("organ", "outcome", "Model2_HR", "Model2_P", "panel")],
  disease[, c("organ", "outcome", "Model2_HR", "Model2_P", "panel")]
)
plot_data$organ <- factor(plot_data$organ, levels = organs)
plot_data$panel <- factor(plot_data$panel, levels = c("Mortality", "Incident chronic disease"))
plot_data$outcome_key <- paste(plot_data$panel, plot_data$outcome, sep = "::")
outcome_levels <- c(
  paste("Mortality", rev(mortality_order), sep = "::"),
  paste("Incident chronic disease", rev(disease_order), sep = "::")
)
outcome_labels <- setNames(sub("^.*::", "", outcome_levels), outcome_levels)
plot_data$outcome_key <- factor(plot_data$outcome_key, levels = outcome_levels)
plot_data$label <- paste0(sprintf("%.2f", plot_data$Model2_HR), stars(plot_data$Model2_P))

if (any(is.na(plot_data$organ)) || any(is.na(plot_data$outcome_key))) {
  stop("Unexpected organ or outcome labels in Fig. 4a input data.")
}

pdf(output, width = 10.5, height = 8.5, useDingbats = FALSE)
print(
  ggplot(plot_data, aes(x = organ, y = outcome_key, fill = Model2_HR)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 3.0) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient(
      name = "HR", low = "#fff5f0", high = "#cb181d",
      limits = c(1, 3), breaks = c(1, 1.5, 2, 2.5, 3), oob = scales::squish
    ) +
    scale_y_discrete(labels = outcome_labels) +
    labs(x = "Organ and tissues", y = "Outcomes") +
    theme_classic(base_size = 11) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.spacing = grid::unit(0.7, "lines"),
      legend.position = "right"
    )
)
dev.off()

message("Wrote: ", output)
