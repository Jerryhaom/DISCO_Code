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


# Reproduce organ/network publication figures from LOCKED legacy saved objects.
# Reads objects directly from DISCO_DATA_ROOT/Result3 (do not modify).
# Skips Fig 4a (already done) and Fig 4c (protein-level network, out of scope).

options(stringsAsFactors = FALSE, warn = 1)

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggpubr)
  library(pheatmap)
  library(igraph)
  library(ggnetwork)
  library(survival)
  library(RColorBrewer)
  library(ggthemes)
  library(readxl)
})

root   <- disco_path("Result3")
organ  <- file.path(root, "Organ")
sens   <- file.path(root, "Sensitivity")
outdir <- "outputs/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

panels <- list()
run <- function(name, f) {
  cat("== ", name, " ==\n")
  panels[[name]] <<- tryCatch({ f(); "PASS" },
                              error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(panels[[name]], "\n\n")
}

# ---- shared objects ---------------------------------------------------------
odis    <- readRDS(file.path(organ, "UKB.odis.rds"))
mc2     <- readRDS(file.path(organ, "UKB.tissue.mc2.rds"))
organ18 <- colnames(mc2)                          # 18 organs, canonical (odis) order
cl      <- readRDS(file.path(organ, "Selected.tissue.rds"))   # 18 indices into 36 cols
prot.cla<- readRDS(file.path(organ, "UKB.Prot.class.rds"))
sel     <- readRDS(file.path(organ, "Lasso.select.prot.rds"))
or2     <- readRDS(file.path(organ, "Tissue.order.rds"))
col     <- readRDS(file.path(organ, "Tissue.col.rds"))
col2    <- readRDS(file.path(organ, "Tissue.col2.rds"))
cs1s    <- readRDS(file.path(organ, "Tissue.cs1.rds"))
gmp     <- readRDS(file.path(organ, "GM.res.rds"))
resM    <- readRDS(file.path(organ, "mortality/res.rds"))
plM     <- readRDS(file.path(organ, "mortality/pl.rds"))

age    <- odis$age
sex    <- odis$sex
status <- odis$status
time   <- odis$time

# quantile-based high-risk cut (replicate legacy cutf)
cutf <- function(x, n = 4) {
  yy <- as.numeric(cut(x, breaks = quantile(x, probs = c(0, 1 - 1/n, 1)),
                       include.lowest = TRUE, right = TRUE)) - 1
  yy
}

# significance star helper (breaks auto-sorted by cut())
ft <- function(x) {
  x1 <- cut(as.numeric(x), breaks = c(0, 0.05, 0.01, 1), include.lowest = TRUE)
  levels(x1) <- c("**", "*", "")
  x1
}

# ===========================================================================
# Fig 4b — High-DISCO organs/tissues vs Age (No.at.risk2.pdf)
# ===========================================================================
run("Fig4b", function() {
  odis1 <- odis[, cl]                 # 18 organ scores
  odis2 <- apply(odis1, 2, cutf, n = 10)
  num   <- apply(odis2, 1, sum)
  nd    <- data.frame(age = age, num = num * 10 / 18, time = time, status = status)
  nd1   <- nd %>% group_by(age) %>% summarise(m = mean(num), se = sd(num) / sqrt(length(num)))

  nd11      <- data.frame(nd1[which(nd1$age > 39), ])
  nd11$logm <- log(nd11$m)
  lm_fit    <- lm(logm ~ age, data = nd11)
  start_A   <- exp(coef(lm_fit)[1])
  model     <- nls(m ~ C + A * exp(B * age), data = nd11,
                   start = list(A = start_A, B = 0.1, C = 1))
  s1         <- summary(model)
  fit_values <- predict(model, newdata = data.frame(age = nd11$age))

  pt         <- data.frame(age = 40:70)
  pt$risk1   <- gmp$Estimate[1] + gmp$Estimate[2] * exp(pt$age * gmp$Estimate[4])

  df   <- data.frame(x = 40:70,
                     y2 = pt$risk1[match(40:70, pt$age)],
                     y1 = nd11$m,
                     y3 = fit_values)
  sst  <- (max(df$y3) - min(df$y3)) / (max(df$y2) - min(df$y2))
  ssb  <- -min(df$y2) * sst + min(df$y3)
  df$y4 <- (df$y2) * sst + ssb

  pdf(file.path(outdir, "Fig4b_No_at_risk2.pdf"), height = 3, width = 4)
  p <- ggplot(df, aes(x = x)) +
    geom_point(aes(y = y1), size = 2, shape = 21, color = "black", fill = "black") +
    geom_line(aes(y = y3, color = "Y1")) +
    geom_line(aes(y = y4, color = "Y2")) +
    scale_y_continuous(
      name = "High-DISCO organs/tissues",
      sec.axis = sec_axis(~ (. - ssb) / sst, name = "Mortality hazard",
                          breaks = seq(0.001, 0.007, 0.002))) +
    xlab("Age") +
    scale_color_manual(values = c("Y1" = "blue", "Y2" = "red"),
                       labels = c(paste("Dysregulations, Shape:", sprintf("%0.3f", s1$coefficients[2])),
                                  paste("Mortality, Shape:", sprintf("%0.3f", gmp$Estimate[4])))) +
    theme_few() +
    annotate("text", x = 44, y = 1, color = "blue", label = "rho = 0.98\np < 0.001") +
    theme(legend.position = c(0.4, 0.85),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.title.y.left = element_text(color = "blue"),
          axis.title.y.right = element_text(color = "red"))
  print(p)
  dev.off()
})

# ===========================================================================
# Fig 4d — 18-organ partial-correlation heatmap (overlap; Plot_Fig5C.R)
# ===========================================================================
run("Fig4d", function() {
  prot1 <- readRDS(file.path(organ, "UKB.prot.ogran.rds"))
  res   <- readRDS(file.path(organ, "UKB.tissue.mc.rds"))
  od    <- data.frame(res)
  cc    <- cor(od)
  selp  <- sel
  cl36  <- which(colnames(cc) %in% gsub(" ", ".", organ18))

  pcorr1 <- matrix(0, ncol = length(cl36), nrow = length(cl36))
  k <- 1
  for (i1 in seq_along(cl36)) {
    for (j1 in seq_along(cl36)) {
      i <- cl36[i1]; j <- cl36[j1]
      if (j >= i) { pcorr1[i1, j1] <- 1; next }
      ol <- selp[[i]][selp[[i]] %in% selp[[j]]]
      pd <- data.frame(y = od[, i], x = od[, j], age = age, sex = sex)
      if (length(ol) > 0) {
        ol <- ol[ol %in% colnames(prot1)]
        if (length(ol) > 0) pd <- cbind(pd, prot1[, ol, drop = FALSE])
      }
      l1 <- lm(y ~ ., data = pd[, -2])
      l2 <- lm(x ~ ., data = pd[, -1])
      cval <- cor(l1$residuals, l2$residuals)
      pcorr1[i1, j1] <- cval
      pcorr1[j1, i1] <- cval
      k <- k + 1
    }
    cat("  pair-row", i1, "/", length(cl36), "\n")
  }
  colnames(pcorr1) <- gsub("[.]", " ", colnames(cc)[cl36])
  rownames(pcorr1) <- colnames(pcorr1)

  saveRDS(pcorr1, file.path(outdir, "Fig4d_pcorr_matrix.rds"))

  or21  <- rev(or2)
  pcorr2 <- pcorr1[match(or21, colnames(pcorr1)), match(or21, colnames(pcorr1))]
  pcorr2[upper.tri(pcorr2, diag = TRUE)] <- NA
  pcorr  <- na.omit(melt(pcorr2))
  pcorr  <- data.frame(pcorr)
  colnames(pcorr) <- c("X1", "X2", "X3")
  pcorr$X3 <- as.numeric(pcorr$X3)
  pcorr$X1 <- gsub("[.]", " ", pcorr$X1)
  pcorr$X2 <- gsub("[.]", " ", pcorr$X2)
  pcorr$X1 <- factor(pcorr$X1, levels = rev(or21))
  pcorr$X2 <- factor(pcorr$X2, levels = or21)

  cat("  pcorr range:", range(pcorr$X3), "\n")

  pdf(file.path(outdir, "Fig4d_Organ_pcorr.pdf"), height = 4.5, width = 5.5)
  p <- ggplot(pcorr, aes(X1, X2)) +
    geom_tile(aes(fill = X3), color = "black", alpha = 0.9) +
    geom_text(aes(label = sprintf("%0.2f", X3)), size = 2) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_classic() + labs(fill = "Partial\nCorr") +
    ylab("Organ and tissue DISCO") + xlab("") +
    labs(title = "Age- and sex- adjusted partial correlations\nbetween organ/tissue DISCO scores") +
    theme(legend.position = c(0.8, 0.7),
          legend.key.size = unit(1, "lines"),
          legend.key.height = unit(0.8, "lines"),
          plot.title = element_text(size = 16, lineheight = 0.5, hjust = 1),
          plot.margin = unit(c(0.1, 0.3, 0.1, 0.1) * 2, "lines"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 14, color = "black")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"))
  print(p)
  dev.off()
})

# ===========================================================================
# Fig 4e — K=3 adjacency network (Partial.net2.pdf)
# ===========================================================================
run("Fig4e", function() {
  ig <- readRDS(file.path(organ, "Organ.Net2.rds"))
  coln <- readRDS(file.path(organ, "Tissue.col2.rds"))
  g2 <- fortify(ig)
  pdf(file.path(outdir, "Fig4e_Partial_net2.pdf"), height = 4, width = 6)
  p <- ggplot(data = g2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(colour = "grey80") +
    geom_point(aes(colour = name, size = size)) +
    geom_nodetext(aes(label = name), size = 5) +
    theme_blank() +
    labs(title = "Network of organ/tissue DISCO correlations") +
    scale_color_manual(values = coln) +
    scale_size_continuous(range = c(5, 10)) +
    guides(color = "none") + xlim(-0.1, 1.05) + labs(size = "Degree") +
    theme(plot.title = element_text(size = 16, hjust = 0.5))
  print(p)
  dev.off()
})

# ===========================================================================
# Fig 4f — degree centrality vs age-corr / mortality C-index / disease C-index
# ===========================================================================
run("Fig4f", function() {
  pcorr1 <- readRDS(file.path(organ, "Tissue.pcorr.noverlap2.rds"))
  colnames(pcorr1) <- organ18
  rownames(pcorr1) <- organ18
  cc  <- pcorr1
  cs1 <- cs1s
  deg <- apply(cc, 1, sum) - 1
  cs1$de <- deg[match(cs1$tissue, names(deg))]
  cs2 <- cs1[cs1$variable == "C-index", ]
  or2n <- cs2$tissue[order(cs2$de, cs2$value, decreasing = TRUE)]
  coln <- readRDS(file.path(organ, "Tissue.col.rds"))
  col2n <- coln[match(or2n, names(coln))]

  dci  <- readxl::read_xlsx(file.path(organ, "Disease_cindex.xlsx"))
  dcit <- na.omit(unique(dci[[1]]))
  dcic <- na.omit(dci[[8]])
  dcic <- as.numeric(dcic[2:length(dcic)])
  dcit <- dcit[2:length(dcit)]
  dci1 <- data.frame(tissue = rep(dcit, each = 8),
                     ci = dcic,
                     di = rep(dci[[2]][3:10], 18))
  dci1$di <- factor(dci1$di, levels = dci[[2]][3:10])
  dci1$de <- deg[match(dci1$tissue, names(deg))]
  dci1$tissue <- factor(dci1$tissue, levels = or2n)

  p1 <- ggplot(cs1[cs1$variable == "C-index", ], aes(value, de)) +
    geom_point(aes(fill = tissue), size = 3, shape = 21, color = "black") +
    geom_smooth(method = "lm", color = "navy") +
    stat_cor(cor.coef.name = "r", label.y.npc = 0.95) +
    scale_fill_manual(values = col2n) +
    theme_classic() + ylim(5, 11) + theme(legend.position = "None") +
    xlab("C-index for mortality") + ylab("Degree centrailty")

  p2 <- ggplot(cs1[cs1$variable == "Corr", ], aes(value, de)) +
    geom_point(aes(fill = tissue), size = 3, shape = 21, color = "black") +
    geom_smooth(method = "lm", color = "navy") +
    stat_cor(cor.coef.name = "r", label.y.npc = 0.95) +
    scale_fill_manual(values = col2n) +
    theme_classic() + ylim(5, 11) + theme(legend.position = "None") +
    xlab("Correlation with age") + ylab("Degree centrailty")

  levels(dci1$di) <- c("Obesity", "Hypertension", "Diabetes", "Cancer",
                       "MI", "Stroke", "COPD", "Dementia")

  # ED11a (8-panel facet of degree vs disease C-index)
  pdf(file.path(outdir, "ED11a_Degree_cindex_disease.pdf"), height = 5, width = 8)
  pf <- ggplot(dci1, aes(ci, de)) +
    geom_point(aes(fill = tissue), size = 3, shape = 21, color = "black") +
    geom_smooth(method = "lm", color = "navy") +
    stat_cor(cor.coef.name = "r", label.y.npc = 0.95) +
    facet_wrap(di ~ ., scales = "free_x", nrow = 2) +
    scale_fill_manual(values = coln) +
    theme_few() + labs(fill = "Organ and tissues") +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.key.height = unit(0.5, "lines"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
          legend.box.margin = unit(c(-1, 0.1, 0.1, 0.1), "lines"),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          panel.spacing.x = unit(1, "lines")) +
    ylim(5, 11) +
    xlab("C-index for diseases") + ylab("Degree centrailty")
  print(pf)
  dev.off()

  dci2 <- dci1 %>% group_by(tissue) %>% summarise(cim = mean(ci), de = mean(de))
  p3 <- ggplot(dci2, aes(cim, de)) +
    geom_point(aes(fill = tissue), size = 3, shape = 21, color = "black") +
    geom_smooth(method = "lm", color = "navy") +
    stat_cor(cor.coef.name = "r", label.y.npc = 0.95) +
    scale_fill_manual(values = coln) +
    theme_classic() + ylim(5, 11) +
    xlab("Mean C-index for diseases") + ylab("Degree centrailty") +
    labs(fill = "Organ and tissues") +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.key.height = unit(0.7, "lines"),
          strip.background = element_blank(),
          strip.text = element_text(size = 12))

  pdf(file.path(outdir, "Fig4f_Tissue_imp.pdf"), height = 2.5, width = 8)
  p <- ggarrange(p2, NULL, p1, NULL, p3,
                 ncol = 5, widths = c(1, 0.05, 1, 0.05, 1),
                 common.legend = TRUE, legend = "none")
  print(p)
  dev.off()
})

# ===========================================================================
# ED9a — organ-specific protein count (before / after LASSO)
# ===========================================================================
run("ED9a", function() {
  pcx <- readxl::read_xlsx(file.path(organ, "ukb_simple_specific_1553.xlsx"))
  pcx$specific.tissue <- gsub(" 1", "", pcx$specific.tissue)
  pcx$specific.tissue <- factor(pcx$specific.tissue,
                                levels = names(sort(table(pcx$specific.tissue))))
  pdf(file.path(outdir, "ED9a_Prot_tissue.pdf"), height = 5, width = 4)
  p <- ggplot(pcx, aes(specific.tissue)) +
    geom_bar(stat = "count", aes(fill = sig)) +
    coord_flip() +
    scale_fill_discrete(labels = c("Enhanced", "Enriched")) +
    labs(fill = "") +
    scale_y_continuous(position = "right", expand = c(0.01, 0.01)) +
    theme_classic() + theme(legend.position = c(0.7, 0.1)) +
    xlab("Oragn and tissues") + ylab("Counts")
  print(p)
  dev.off()

  # after LASSO (selected proteins only), Enhanced vs Enriched
  pca <- prot.cla
  pca$sig <- as.factor(pca$sig)
  count <- matrix(0, ncol = 2, nrow = length(cl))
  for (i in seq_along(cl)) {
    i1 <- cl[i]
    l1 <- pca$sig[pca$Gene %in% sel[[i1]]]
    count[i, ] <- c(length(l1[l1 == levels(l1)[1]]),
                    length(l1[l1 == levels(l1)[2]]))
  }
  count <- data.frame(count)
  colnames(count) <- c("Enhanced", "Enriched")
  count$tissue <- organ18
  count$tissue <- factor(count$tissue,
                         levels = count$tissue[order(count$Enhanced + count$Enriched)])
  count <- melt(count, id = "tissue")
  pdf(file.path(outdir, "ED9a_Prot_tissue_sel.pdf"), height = 5, width = 4)
  p <- ggplot(count, aes(value, tissue)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    labs(fill = "") +
    scale_x_continuous(position = "top", expand = c(0.01, 0.01)) +
    theme_classic() + theme(legend.position = c(0.7, 0.1)) +
    ylab("organ and tissue") + xlab("Count") +
    theme(axis.text.y = element_text(size = 14, color = "black"),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14))
  print(p)
  dev.off()
})

# ===========================================================================
# ED9b — shared-protein matrix heatmap (before/after feature selection)
# ===========================================================================
run("ED9b", function() {
  od <- odis
  cc <- cor(od[, 1:36])
  overlapped <- matrix(0, ncol = length(cl), nrow = length(cl))
  for (i in seq_along(cl)) {
    for (j in seq_along(cl)) {
      i1 <- cl[i]; j1 <- cl[j]
      ol <- sel[[i1]][sel[[i1]] %in% sel[[j1]]]
      overlapped[i, j] <- length(ol)
    }
  }
  colnames(overlapped) <- gsub("[.]", " ", colnames(cc)[cl])
  rownames(overlapped) <- colnames(overlapped)

  # before feature selection (all proteins per tissue, ENTREZID)
  pca   <- prot.cla
  tissues <- gsub("[.]", " ", colnames(overlapped))
  pca2  <- pca[pca$specific.tissue %in% tissues, ]
  clt   <- unique(pca2$specific.tissue)
  overlapped2 <- matrix(0, ncol = length(clt), nrow = length(clt))
  for (i in seq_along(clt)) {
    for (j in seq_along(clt)) {
      i1 <- pca2$ENTREZID[pca2$specific.tissue == clt[i]]
      j1 <- pca2$ENTREZID[pca2$specific.tissue == clt[j]]
      ol <- i1[i1 %in% j1]
      overlapped2[i, j] <- length(ol)
    }
  }
  colnames(overlapped2) <- clt
  rownames(overlapped2) <- colnames(overlapped2)

  overlapped3 <- overlapped2[match(or2, rownames(overlapped2)),
                             match(or2, rownames(overlapped2))]
  overlapped31 <- overlapped3
  overlapped31[upper.tri(overlapped3, diag = FALSE)] <- NA
  pd <- na.omit(melt(overlapped31))
  pd$Var1 <- factor(pd$Var1, levels = rev(or2))
  pd$Var2 <- factor(pd$Var2, levels = or2)

  p1 <- ggplot(pd, aes(Var2, Var1)) +
    geom_tile(aes(fill = value), color = "black") +
    geom_text(aes(label = value)) +
    theme_classic() +
    scale_fill_distiller(palette = "Reds", limits = c(0, 256),
                         breaks = c(0, 100, 200), direction = 1) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10)) +
    xlab("") + ylab("Before feature selection") +
    labs(fill = "Counts") +
    theme(legend.direction = "horizontal",
          legend.key.height = unit(0.8, "lines"),
          legend.key.width  = unit(0.8, "lines"),
          legend.position = c(-0.08, 1.1),
          legend.title.position = "top",
          plot.margin = unit(c(7, 7, 0.5, 0.5), "lines"))

  overlapped11 <- overlapped[match(or2, rownames(overlapped)),
                             match(or2, rownames(overlapped))]
  overlapped41 <- overlapped11
  overlapped41[lower.tri(overlapped41, diag = FALSE)] <- NA
  pd <- na.omit(melt(overlapped41))
  pd$Var1 <- factor(pd$Var1, levels = rev(or2))
  pd$Var2 <- factor(pd$Var2, levels = or2)

  p2 <- ggplot(pd, aes(Var2, Var1)) +
    geom_tile(aes(fill = value), color = "black") +
    geom_text(aes(label = value)) +
    theme_classic() +
    scale_fill_distiller(palette = "Reds", limits = c(0, 256),
                         breaks = c(0, 100, 200), direction = 1) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0),
          axis.text.y = element_text(size = 10),
          legend.position = "None") +
    xlab("After feature selection") + ylab("") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"))

  pdf(file.path(outdir, "ED9b_Tissue_overlapped.pdf"), height = 6, width = 7,
      onefile = FALSE)
  text <- "Shared proteins used to calculate organ/tissue DISCO scores"
  g <- ggplotGrob(p2)
  p <- p1 + annotation_custom(g, xmin = 1.5, xmax = 25, ymin = 2, ymax = 26) +
    labs(title = text) +
    theme(plot.title = element_text(hjust = 0, vjust = 32))
  print(p)
  dev.off()
})

# ===========================================================================
# ED9c — 18-organ vs clinical biomarker partial correlation heatmap
# ===========================================================================
run("ED9c", function() {
  r <- readRDS(file.path(organ, "r.rds"))
  p <- readRDS(file.path(organ, "p.rds"))
  pl1 <- p
  p <- apply(p, 2, ft)
  clb <- c(3, 5, 6, 21, 13:20, 26:38)
  r1 <- r[clb, ]
  p1 <- p[clb, ]
  rownames(r1)[17] <- "CHOL"
  rownames(r1)[22] <- "Lp(a)"
  anno1 <- data.frame(Biomarkers = c(rep("Immune", 4), rep("Oxygen", 2),
                                     rep("Liver", 6), rep("Kidney", 4),
                                     rep("Lipid", 7), rep("Glucose", 2)))
  rownames(anno1) <- rownames(r1)
  cc <- brewer.pal(6, "Set2")
  names(cc) <- unique(anno1$Biomarkers)
  cc <- list(Biomarkers = cc)
  r2 <- t(r1)
  r2 <- r2[match(or2, rownames(r2)), ]

  pdf(file.path(outdir, "ED9c_DISCO_biomarker.pdf"), height = 3.5, width = 6)
  pheatmap(r2,
           annotation_col = anno1,
           annotation_colors = cc,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           display_numbers = t(p1),
           number_color = "black",
           fontsize_number = 8,
           color = colorRampPalette(brewer.pal(7, "Reds"))(50),
           treeheight_row = 15)
  dev.off()
})

# ===========================================================================
# ED9d — 3-column bar (corr with age / DISCO C-index / age+DISCO C-index)
# ===========================================================================
run("ED9d", function() {
  od <- odis
  cs <- matrix(0, nrow = 36, ncol = 6)
  for (i in 1:36) {
    fl <- paste("Surv(time,status)~", colnames(od)[i])
    c1 <- coxph(as.formula(fl), data = od); s <- summary(c1)
    cr <- cor.test(od[, i], od$age)
    fl <- paste("Surv(time,status)~", colnames(od)[i], "+age+sex")
    c2 <- coxph(as.formula(fl), data = od); s1 <- summary(c2)
    cs[i, ] <- c(cr$estimate, cr$estimate - cr$conf.int[1],
                 s$concordance[1], s$concordance[2],
                 s1$concordance[1], s1$concordance[2])
  }
  cs <- data.frame(cs)
  colnames(cs) <- c("Corr", "Se", "C-index", "Se1", "C-index1", "Se2")
  cs$tissue <- colnames(od)[1:36]
  cs1 <- data.frame(melt(cs[, c(1, 3, 5, 7)], id = "tissue"),
                    se = melt(cs[, c(2, 4, 6, 7)], id = "tissue")[, 3])
  cs1$tissue <- gsub("[.]", " ", cs1$tissue)
  cs1 <- cs1[cs1$tissue %in% organ18, ]
  cs1$tissue <- factor(cs1$tissue, levels = rev(or2))

  p1 <- ggplot(cs1[cs1$variable == "Corr", ], aes(value, tissue)) +
    geom_bar(stat = "identity", aes(fill = tissue), width = 0.8, color = "black") +
    geom_errorbar(aes(xmin = value - se, xmax = value + se), width = 0.5) +
    scale_fill_manual(values = col2) +
    theme_few() +
    theme(legend.position = "None",
          plot.title = element_text(size = 12, hjust = -0.1),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 0.35)) +
    xlab("Correlation with age") + ylab("Organ and tissue DISCO")

  p2 <- ggplot(cs1[grep("index$", cs1$variable), ], aes(value, tissue)) +
    geom_bar(stat = "identity", aes(fill = tissue), position = "dodge",
             width = 0.8, color = "black") +
    geom_errorbar(aes(xmin = value - se, xmax = value + se), width = 0.5) +
    scale_fill_manual(values = col2) +
    coord_cartesian(xlim = c(0.55, 0.78)) +
    theme_few() +
    theme(legend.position = "None",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 12)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0.55, 0.65, 0.75)) +
    xlab("DISCO C-index")

  c0 <- coxph(Surv(time, status) ~ age + sex, data = od)
  s0 <- summary(c0)

  p3 <- ggplot(cs1[grep("index1$", cs1$variable), ], aes(value, tissue)) +
    geom_bar(stat = "identity", aes(fill = tissue), position = "dodge",
             width = 0.8, color = "black") +
    geom_errorbar(aes(xmin = value - se, xmax = value + se), width = 0.5) +
    geom_vline(xintercept = s0$concordance[1], linetype = "dashed") +
    scale_fill_manual(values = col2) +
    coord_cartesian(xlim = c(0.72, 0.78)) +
    theme_few() +
    theme(legend.position = "None",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 12)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0.72, 0.74, 0.76, 0.78)) +
    xlab("Age + DISCO C-index")

  pdf(file.path(outdir, "ED9d_Organ_corr_cindex2.pdf"), height = 4, width = 7)
  p <- ggarrange(p1, p2, NULL, p3, NULL,
                 nrow = 1, ncol = 5,
                 widths = c(0.8, 0.5, 0.01, 0.5, 0.05))
  print(p)
  dev.off()
})

# ===========================================================================
# ED9e — age-stratified high-risk organ prevalence heatmap
# ===========================================================================
run("ED9e", function() {
  odis1 <- odis[, cl]
  odis2 <- apply(odis1, 2, cutf, n = 10)
  odis21 <- data.frame(odis2)
  odis21$age <- age
  odis21 <- melt(odis21, id = "age")
  stat <- odis21 %>% group_by(age, variable) %>%
    summarise(l = length(value[value > 0]) / length(value))
  stat <- stat[stat$age > 39, ]
  stat$variable <- gsub("[.]", " ", stat$variable)
  stat$variable <- factor(stat$variable, levels = rev(or2))

  pdf(file.path(outdir, "ED9e_Prop_at_risk.pdf"), height = 4, width = 6.5)
  p <- ggplot(stat[stat$age > 39, ], aes(age, variable)) +
    geom_tile(aes(fill = l), color = "black") +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Prop of \nat risk") +
    ylab("Organ and tissues") + xlab("Age") +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 14))
  print(p)
  dev.off()
})

# ===========================================================================
# ED10a — Gompertz M0/Shape/Rate vs high-risk threshold
# ===========================================================================
run("ED10a", function() {
  odis1 <- odis[, cl]
  cal_b <- function(n) {
    odis2 <- apply(odis1, 2, cutf, n)
    num <- apply(odis2, 1, sum)
    nd <- data.frame(age = age, num = num * n / 18, time = time, status = status)
    nd1 <- nd %>% group_by(age) %>% summarise(m = mean(num), se = sd(num) / sqrt(length(num)))
    nd11 <- data.frame(nd1[which(nd1$age > 39), ])
    nd11$logm <- log(nd11$m)
    lm_fit <- lm(logm ~ age, data = nd11)
    start_A <- exp(coef(lm_fit)[1])
    model <- nls(m ~ C + A * exp(B * age), data = nd11,
                 start = list(A = start_A, B = 0.1, C = 1))
    s1 <- summary(model)
    s1$coefficients[, 1:2]
  }
  ns <- c(2:10, 15, 20, 25, 30, 40, 50, 100)
  rb <- sapply(ns, cal_b)
  prb <- data.frame(p = sprintf("%0.2f", 100 / ns), t(rb))
  prb$p <- factor(prb$p, levels = prb$p)
  prb <- cbind(melt(prb[, 1:4], id = "p"),
               se = melt(prb[, c(1, 5:7)], id = "p")[, 3])
  prb$variable <- as.factor(prb$variable)
  levels(prb$variable) <- c("M0", "Shape", "Rate")

  pdf(file.path(outdir, "ED10a_Parameter_est.pdf"), height = 4, width = 7)
  p <- ggplot(prb, aes(value, p)) +
    geom_point() +
    geom_errorbar(aes(xmin = value - se, xmax = value + se), width = 0.3) +
    facet_wrap(variable ~ ., scales = "free_x", ncol = 3) +
    theme_bw() + ylab("Thresholds for high-risks") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12, color = "black"),
          panel.spacing.x = unit(1, "lines")) +
    scale_y_discrete(labels = paste(sprintf("%0.1f",
                                            as.numeric(as.character(unique(prb$p)))), "%", sep = ""))
  print(p)
  dev.off()
})

# ===========================================================================
# ED10 b/c/d + ED11 c/d/e — random-protein sensitivity (plot_sensitivity.R)
# ===========================================================================
run("ED10c", function() {
  cindex <- readRDS(file.path(sens, "UKB.organ.random.cindex.rds"))
  cindex1 <- melt(cindex)
  cindex1$Var1 <- factor(cindex1$Var1, levels = rev(or2))
  coln <- readRDS(file.path(organ, "Tissue.col.rds"))
  pdf(file.path(outdir, "ED10c_Organ_cindex_sensitivity.pdf"), height = 4.5, width = 4)
  p <- ggplot(cindex1, aes(value, Var1)) +
    geom_boxplot(aes(fill = Var1)) +
    scale_fill_manual(values = coln) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.y = element_text(size = 10, color = "black")) +
    xlab("C-index") + ylab("Organ and tissue DISCO") +
    labs(title = "Sensitivity analysis of mortality prediction by\nproteins included in organ/tissue DISCO scores") +
    theme(plot.title = element_text(lineheight = 0.5, hjust = 1, size = 12))
  print(p)
  dev.off()
})

run("ED11c", function() {
  pcor1 <- readRDS(file.path(sens, "Tissue.corr.sensitivity.rds"))
  pcor <- pcor1
  pcor$Var1 <- factor(pcor$Var1, levels = rev(or2))
  pcor$Var2 <- factor(pcor$Var2, levels = rev(or2))
  pcor$value[as.numeric(pcor$Var1) >= as.numeric(pcor$Var2)] <- NA
  pcor$Var2 <- factor(pcor$Var2, levels = or2)

  pdf(file.path(outdir, "ED11c_Organ_pcorr_sensitivity.pdf"), height = 4, width = 5)
  p <- ggplot(pcor[!is.na(pcor$value), ], aes(Var2, Var1)) +
    geom_tile(aes(fill = value), color = "black", alpha = 0.9) +
    geom_text(aes(label = sprintf("%0.2f", value)), size = 2) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_classic() + labs(fill = "Partial\nCorr") +
    ylab("Organ and tissue DISCO") + xlab("") +
    labs(title = "Sensitivity analysis of organ/tissue DISCO correlations") +
    theme(legend.position = c(0.8, 0.7),
          legend.key.size = unit(1, "lines"),
          legend.key.height = unit(0.8, "lines"),
          plot.title = element_text(size = 12, hjust = 1),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, color = "black")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"))
  print(p)
  dev.off()
})

run("ED11d", function() {
  all.pcorr <- readRDS(file.path(organ, "Tissue.pcorr.noverlap.rds"))
  all.pcorr <- data.frame(all.pcorr)
  all.pcorr$X1 <- gsub("[.]", " ", all.pcorr$X1)
  all.pcorr$X2 <- gsub("[.]", " ", all.pcorr$X2)
  pcor1 <- readRDS(file.path(sens, "Tissue.corr.sensitivity.rds"))
  pcor1$raw <- all.pcorr$X3[match(paste(pcor1$Var1, pcor1$Var2),
                                   paste(all.pcorr$X1, all.pcorr$X2))]
  pcor1$raw <- as.numeric(pcor1$raw)

  pdf(file.path(outdir, "ED11d_Organ_pcorr_cmp.pdf"), height = 3, width = 3.5)
  p <- ggplot(pcor1[!is.na(pcor1$raw), ], aes(value, raw)) +
    geom_point(size = 2, color = "grey") +
    geom_smooth(method = "lm", color = "black") +
    theme_classic() +
    ylab("Parital correlation based on\nspecific protein panel") +
    xlab("Parital correlation based on\nrandomly selected protein panel") +
    ylim(0.25, 0.7) +
    labs(title = "Sensitivity of organ/tissue DISCO\ncorrelation to protein composition") +
    stat_cor(cor.coef.name = "r", label.y = 0.7 * 0.95) +
    theme(plot.title = element_text(size = 12, hjust = 1, lineheight = 0.5))
  print(p)
  dev.off()
})

run("ED11e", function() {
  pcor <- readRDS(file.path(sens, "Tissue.corr.sensitivity.rds"))
  pcorr <- dcast(pcor, Var1 ~ Var2, value.var = "value")
  rownames(pcorr) <- pcorr[, 1]
  pcorr <- pcorr[, 2:ncol(pcorr)]
  deg <- apply(pcorr, 1, sum) - 1

  cindex <- readRDS(file.path(sens, "UKB.organ.random.cindex.rds"))
  cindex1 <- melt(cindex)
  cindex2 <- cindex1 %>% group_by(Var1) %>% summarise(m = mean(value))
  cindex2$de <- deg[match(cindex2$Var1, names(deg))]
  coln <- readRDS(file.path(organ, "Tissue.col.rds"))
  cindex2$Var1 <- factor(cindex2$Var1, levels = names(coln))

  pdf(file.path(outdir, "ED11e_Organ_cindex_degree.pdf"), height = 3.5, width = 4)
  p <- ggplot(cindex2, aes(m, de)) +
    geom_point(aes(fill = Var1), size = 3, shape = 21, color = "black") +
    geom_smooth(method = "lm", color = "navy", se = FALSE) +
    stat_cor(cor.coef.name = "r", label.y.npc = 0.95) +
    scale_fill_manual(values = coln) +
    theme_classic() +
    ylim(6, 11) +
    theme(legend.position = "None") +
    xlab("Mean C-index for mortality \n(200 replicates)") +
    ylab("Degree centrailty") +
    labs(title = "Sensitivity analysis of degree centrality-mortality\nassociations of organ/tissue DISCOs") +
    theme(plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.5),
          axis.title = element_text(size = 12))
  print(p)
  dev.off()
})

# ED10b / ED10d need the 200-replicate random DISCO object (1.5 GB)
run("ED10b_d", function() {
  mcl <- readRDS(file.path(sens, "UKB.organ.random.disco.rds"))
  tissue2 <- organ18

  # ED10b — age trend scatter (No.at.risk.sensitivity)
  stm <- function(m) {
    odis <- apply(m, 2, cutf, n = 10)
    apply(odis, 1, sum)
  }
  countm <- lapply(mcl, stm)
  countm1 <- matrix(0, nrow = nrow(mcl[[1]]), ncol = length(mcl))
  for (j in seq_along(mcl)) countm1[, j] <- countm[[j]]
  countm1 <- data.frame(countm1)
  colnames(countm1) <- paste("R", seq(1, 200), sep = "")
  countm1$age <- age
  countm1 <- melt(countm1, id = "age")
  nd1 <- countm1 %>% group_by(age, variable) %>% summarise(m = mean(value * 10 / 18))

  pdf(file.path(outdir, "ED10b_No_at_risk_sensitivity.pdf"), height = 3, width = 4)
  p <- ggplot(data = nd1[nd1$age >= 40, ], aes(age, m)) +
    geom_point(shape = 21, size = 3, color = "grey") +
    geom_smooth(color = "navy", se = FALSE, size = 0.1, aes(group = variable)) +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 0.95) +
    theme_classic() + ylab("High-DISCO organs/tissues") +
    xlab("Age") +
    labs(title = "Sensitivity analysis of age associations by\nproteins included in organ/tissue DISCO scores") +
    theme(plot.title = element_text(lineheight = 0.5, hjust = 1, size = 12))
  print(p)
  dev.off()

  # ED10d — age-stratified prevalence heatmap (Prop.at.risk.sensitivity)
  cal_prop <- function(m) {
    odis21 <- data.frame(m)
    colnames(odis21) <- tissue2
    odis21 <- apply(odis21, 2, cutf, n = 10)
    odis21 <- data.frame(odis21)
    odis21$age <- age
    odis21 <- melt(odis21, id = "age")
    stat <- odis21 %>% group_by(age, variable) %>%
      summarise(l = length(value[value > 0]) / length(value))
    stat[stat$age > 39, ]
  }
  propm <- lapply(mcl, cal_prop)
  propm1 <- matrix(0, nrow = nrow(propm[[1]]), ncol = length(propm))
  for (j in seq_along(mcl)) propm1[, j] <- propm[[j]][[3]]
  propm1 <- data.frame(propm1)
  stat1 <- propm[[1]]
  stat1$l <- apply(propm1, 1, mean)
  stat1$variable <- as.character(stat1$variable)
  stat1$variable <- gsub("[.]", " ", stat1$variable)
  stat1$variable <- factor(stat1$variable, levels = rev(or2))

  pdf(file.path(outdir, "ED10d_Prop_at_risk_sensitivity.pdf"), height = 4, width = 6.5)
  p <- ggplot(stat1[stat1$age > 39, ], aes(age, variable)) +
    geom_tile(aes(fill = l), color = "black") +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Prop of \nat risk") +
    ylab("Organ and tissues") + xlab("Age") +
    labs(title = "Sensitivity analysis of trend of high-DISCO organs/tissues with age") +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 14)) +
    theme(plot.title = element_text(lineheight = 0.5, hjust = 0.5, size = 14))
  print(p)
  dev.off()
})

# ---- summary + MISSING -----------------------------------------------------
cat("\n==================== SUMMARY ====================\n")
for (nm in names(panels)) cat(sprintf("  %-10s %s\n", nm, panels[[nm]]))

miss_lines <- c("# Figure reproduction status",
                "",
                paste0("- ", names(panels), ": ", unname(unlist(panels))),
                "")
fails <- names(panels)[!grepl("^PASS", panels)]
if (length(fails) > 0) {
  writeLines(c(miss_lines,
               "## Failed panels",
               paste0("- ", fails, ": ", unlist(panels[fails]))),
             file.path(outdir, "MISSING.md"))
} else {
  writeLines(c(miss_lines, "## No failed panels"), file.path(outdir, "MISSING.md"))
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nDone. Outputs in", outdir, "\n")
