plot_ED7a <- function(root, out) {
alld <- fread(file.path(root, "input/ED7a_event_time_curves.csv"))
ann <- fread(file.path(root, "input/ED7a_C_annotations.csv"))
titles <- c("Obesity", "Hypertension", "Diabetes", "Cancer", "Myocardial infarction", "Stroke", "COPD", "Dementia")
glabs <- c("DISCO Q1", "DISCO Q4", "DM Q1", "DM Q4")
gcols <- c("DISCO Q1" = "#E64B35FF", "DISCO Q4" = "#E64B35FF", "DM Q1" = "#4DBBD5FF", "DM Q4" = "#4DBBD5FF")
glty <- c("DISCO Q1" = "solid", "DISCO Q4" = "dashed", "DM Q1" = "solid", "DM Q4" = "dashed")
alld[, group := factor(group, levels = glabs)]
alld[, disease := factor(disease, levels = titles)]
ann[, disease := factor(disease, levels = titles)]
stopifnot(nrow(unique(alld[, .(disease, group)])) == 32L, nrow(ann) == 8L,
  all(is.finite(alld$time)), all(is.finite(alld$surv)))
ymax_df <- alld[, .(ymax = ceiling(max(1 - lower, na.rm = TRUE) / 0.05) * 0.05), by = disease]
p <- ggplot(alld, aes(time, 1 - surv, colour = group, fill = group, linetype = group)) +
  geom_ribbon(aes(ymin = 1 - upper, ymax = 1 - lower), alpha = 0.10, colour = NA, show.legend = FALSE) +
  geom_step(linewidth = 0.45) +
  geom_blank(data = ymax_df, aes(x = 40, y = ymax), inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.2, size = 2.7, colour = "black") +
  facet_wrap(~ disease, ncol = 4, scales = "free") +
  scale_colour_manual(values = gcols) + scale_fill_manual(values = gcols) +
  scale_linetype_manual(values = glty) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Age (Years)", y = "Cumulative incidence") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.margin = margin(t = -6, b = 0, l = 0, r = 0),
    legend.key.size = unit(0.6, "cm"), legend.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"), strip.background = element_blank(),
    panel.spacing = unit(0.9, "lines"), axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 11), axis.title.y = element_text(size = 12))
ggsave(file.path(out, "ED7a_attained_age_1minusKM_NEW_ENDPOINT_CANDIDATE.pdf"), p, width = 8, height = 4.5)
cat("Replayed accepted ED7a drawing from 32 aggregate curves, without private inputs.\n")

}
