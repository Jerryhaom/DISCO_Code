plot_simulation <- function(root, out, res, smoke = FALSE, reference_only = FALSE) {
  env <- new.env(parent = globalenv())
  source(file.path(root, "source/plot_c_original.R"), local = env)
  source(file.path(root, "source/plot_ed5_transparent.R"), local = env)
  columns <- c("Mahalanobis_C", "Mahalanobis_se", "DISCO_C", "DISCO_se", "Entropy_C", "Entropy_se", "KL_C", "KL_se")
  if (smoke) {
    stopifnot(nrow(res) == 4L)
    p <- env$plot_c(as.data.frame(res[, ..columns]), res$file) +
      labs(title = "Smoke test: 10 variables, rho=0.1, noise SD=0.2")
    pdf(file.path(out, "SMOKE_four_models_NOT_S6_S7.pdf"), width = 5, height = 4)
    print(p); dev.off()
    return(invisible(NULL))
  }
  stopifnot(nrow(res) == 96L)
  plots <- dat <- list(); k <- 1L
  for (scenario_name in c("Standard", "HighNoise")) for (feature_count in c(10L, 50L, 100L)) {
    z <- res[res$scenario == scenario_name & res$p == feature_count][order(file)]
    stopifnot(nrow(z) == 16L)
    plots[[k]] <- env$plot_c(as.data.frame(z[, ..columns]), z$file) +
      labs(title = sprintf("%s variables, Noise ~ N(0,%s)", feature_count,
        if (scenario_name == "Standard") "0.2" else "1"))
    dat[[k]] <- plots[[k]]$data; k <- k + 1L
  }
  pages <- lapply(list(1:3, 4:6), function(index) ggarrange(plotlist = plots[index],
    common.legend = TRUE, ncol = 1, labels = c("a", "b", "c"), legend = "bottom"))
  for (i in 1:2) {
    pdf(file.path(out, paste0("S", i + 5, "_original_drawing.pdf")), height = 9, width = 9, family = "Helvetica")
    print(pages[[i]]); dev.off()
  }
  ed_dat <- lapply(dat[4:6], function(z) {
    z <- z[as.character(z$Distance) %in% c("DISCO", "Mahalanobis"), , drop = FALSE]
    z$Distance <- factor(as.character(z$Distance), levels = c("DISCO", "Mahalanobis"), labels = c("DISCO", "DM"))
    z$rho <- sub("Corr =", "Rho =", z$rho, fixed = TRUE); z
  })
  ed <- ggarrange(plotlist = Map(env$plot_ed5, ed_dat, c(10, 50, 100)), common.legend = TRUE,
    ncol = 1, labels = c("a", "b", "c"), legend = "bottom")
  pdf(file.path(out, "ED5_transparent_drawing.pdf"), width = 8, height = 10, family = "Helvetica")
  print(ed); dev.off()
  fwrite(rbindlist(Map(function(z, name) as.data.table(z)[, panel := name], dat,
    c("S6a", "S6b", "S6c", "S7a", "S7b", "S7c"))), file.path(out, "S6_S7_plot_cells_384.csv"))
  fwrite(rbindlist(Map(function(z, name) as.data.table(z)[, panel := name], ed_dat,
    c("ED5a", "ED5b", "ED5c"))), file.path(out, "ED5_plot_cells_96.csv"))
  write_json(list(source = if (reference_only) "BUNDLED_MODEL_SUMMARIES_NOT_NEW96_FITS" else "NEW96_GENERATED_DATA384_FITS",
    rows = 96L, S6_S7_cells = 384L, ED5_cells = 96L,
    error_bar = "One model SE, not95%CI", ED5_final_historical_drawing_recovered = FALSE),
    file.path(out, "plot_scope.json"), auto_unbox = TRUE, pretty = TRUE)
}
