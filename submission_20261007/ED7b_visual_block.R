lbl <- c(
  amd = "AMD", ad = "Alzheimer's", asthma = "Asthma", af = "Atrial fibrillation",
  bipolar = "Bipolar", bowel_ca = "Bowel cancer", breast_ca = "Breast cancer", cvd = "CVD",
  coeliac = "Coeliac", cad = "CAD", crohn = "Crohn's", ovarian_ca = "Ovarian cancer",
  hypertension = "Hypertension", ischaemic_stroke = "Isch. stroke", melanoma = "Melanoma",
  ms = "MS", osteoporosis = "Osteoporosis", parkinson = "Parkinson's", glaucoma = "Glaucoma",
  prostate_ca = "Prostate cancer", psoriasis = "Psoriasis", ra = "RA",
  schizophrenia = "Schizophrenia", sle = "SLE", t1d = "T1D", t2d = "T2D",
  ulcerative_colitis = "UC", vte = "VTE")

sys <- c(
  af = "Cardiovascular", cad = "Cardiovascular", cvd = "Cardiovascular",
  hypertension = "Cardiovascular", ischaemic_stroke = "Cardiovascular", vte = "Cardiovascular",
  t1d = "Metabolic", t2d = "Metabolic",
  ad = "Neurological", bipolar = "Neurological",
  ms = "Neurological", parkinson = "Neurological",
  schizophrenia = "Neurological",
  bowel_ca = "Cancer", breast_ca = "Cancer", melanoma = "Cancer",
  ovarian_ca = "Cancer", prostate_ca = "Cancer",
  coeliac = "Autoimmune", crohn = "Autoimmune",
  psoriasis = "Autoimmune", ra = "Autoimmune",
  sle = "Autoimmune", ulcerative_colitis = "Autoimmune",
  asthma = "Respiratory", amd = "Ophthalmic", glaucoma = "Ophthalmic",
  osteoporosis = "Musculoskeletal")

sys_order <- c("Cardiovascular", "Metabolic", "Neurological",
               "Cancer", "Autoimmune", "Respiratory", "Ophthalmic",
               "Musculoskeletal")

d$system <- factor(sys[d$disease], levels = sys_order)
d$label <- lbl[d$disease]

clin <- d[d$version == "clinical", ]
clin <- clin[order(clin$system, clin$interaction_hr), ]
global_order <- clin$disease

build <- function(ver, show_labels) {
  dd <- d[d$version == ver, ]
  Diseases <- c(); HRCI <- c(); Nn <- c(); Ev <- c()
  est <- c(); low <- c(); up <- c(); sysrow <- c()
  for (sy in sys_order) {
    if (show_labels) Diseases <- c(Diseases, sy)
    HRCI <- c(HRCI, ""); Nn <- c(Nn, ""); Ev <- c(Ev, "")
    est <- c(est, NA); low <- c(low, NA); up <- c(up, NA); sysrow <- c(sysrow, TRUE)
    sub <- dd[dd$system == sy, ]
    sub <- sub[match(global_order[global_order %in% sub$disease], sub$disease), ]
    for (i in seq_len(nrow(sub))) {
      if (show_labels) Diseases <- c(Diseases, paste0("   ", sub$label[i]))
      star <- ifelse(sub$interaction_p[i] < 0.001, "***",
              ifelse(sub$interaction_p[i] < 0.01, "**",
              ifelse(sub$interaction_p[i] < 0.05, "*", "")))
      HRCI <- c(HRCI, sprintf("%.2f (%.2f\u2013%.2f)%s", sub$interaction_hr[i],
                              sub$interaction_lo[i], sub$interaction_hi[i], star))
      Nn <- c(Nn, formatC(sub$n[i], big.mark = ",", format = "d"))
      Ev <- c(Ev, formatC(sub$events[i], big.mark = ",", format = "d"))
      est <- c(est, sub$interaction_hr[i]); low <- c(low, sub$interaction_lo[i])
      up <- c(up, sub$interaction_hi[i]); sysrow <- c(sysrow, FALSE)
    }
  }
  if (show_labels) {
    tab <- data.frame(Diseases = Diseases, `Interaction HR (95% CI)` = HRCI, N = Nn, Events = Ev,
                      " " = rep(" ", length(Diseases)),
                      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    tab <- data.frame(`Interaction HR (95% CI)` = HRCI, N = Nn, Events = Ev,
                      " " = rep(" ", length(HRCI)),
                      stringsAsFactors = FALSE, check.names = FALSE)
  }
  list(tab = tab, est = est, lower = low, upper = up, sysrow = sysrow,
       ncol = ifelse(show_labels, 5, 4))
}

make_panel <- function(ver, show_labels) {
  b <- build(ver, show_labels)
  nc <- b$ncol
  tm <- forest_theme(base_size = 15, base_family = "Helvetica",
    ci_pch = 15, ci_fill = "#E64B35", ci_col = "#E64B35",
    ci_lwd = 1.2, ci_lty = 1,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey30"),
    core = list(bg_params = list(fill = "white"), padding = unit(c(2, 3), "mm")),
    colhead = list(fg_params = list(fontface = "bold"), padding = unit(c(2, 3), "mm")))
  p <- forest(b$tab, est = b$est, lower = b$lower, upper = b$upper,
              sizes = .45, ci_column = nc, ref_line = 1,
              xlim = c(0.78, 1.20), ticks_at = c(0.8, 0.9, 1.0, 1.1, 1.2),
              ticks_digits = 1, theme = tm)
  sysidx <- which(b$sysrow)
  if (show_labels) p <- edit_plot(p, row = sysidx, col = 1, gp = gpar(fontface = "bold"))
  p <- add_border(p, part = "header", col = 1:nc, where = "bottom",
                  gp = gpar(lwd = 1.5, col = "black"))
  n <- length(b$sysrow)
  seps <- c(sysidx[-1] - 1, n)
  for (s in seps) p <- add_border(p, row = s, col = 1:nc, where = "bottom",
                                   gp = gpar(lwd = 0.8, col = "grey70"))
  # column widths (shared HR/N/Events/forest; Diseases only in left panel)
  if (show_labels) {
    p$widths[2:(nc + 1)] <- unit(c(1.8, 2.4, 1.0, 1.1, 2.3), "in")
  } else {
    p$widths[2:(nc + 1)] <- unit(c(2.4, 1.0, 1.1, 2.3), "in")
  }
  p$heights[2] <- unit(.38, "in")
  p$heights[3:(nrow(b$tab) + 2)] <- unit(.30, "in")
  list(p = p)
}

left <- make_panel("clinical", show_labels = TRUE)
right <- make_panel("proteomics", show_labels = FALSE)
left$p <- edit_plot(left$p, part='header', col=2, gp=gpar(fontsize=14))
right$p <- edit_plot(right$p, part='header', col=1, gp=gpar(fontsize=14))
pdf(file.path(out, "ED7b_f1_interaction_NEW_ENDPOINT_CANDIDATE.pdf"), width = 20, height = 12, family = "Helvetica",
    title = "ED7b DISCO x PRS interaction forest (raw-date death-first endpoints)")
grid.newpage()
pushViewport(viewport(x = .33, y = .49, width = .43, height = .94)); grid.draw(left$p); popViewport()
pushViewport(viewport(x = .715, y = .49, width = .34, height = .94)); grid.draw(right$p); popViewport()
grid.text("Clinical DISCO", x = .33, y = .985, just = "center",
          gp = gpar(fontsize = 15, fontfamily = "Helvetica", fontface = "bold"))
grid.text("Proteomics DISCO", x = .715, y = .985, just = "center",
          gp = gpar(fontsize = 15, fontfamily = "Helvetica", fontface = "bold"))
dev.off()
cat("saved interaction forest plot (forestploter)\n")
