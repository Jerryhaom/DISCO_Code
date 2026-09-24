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
disco_init_module("01_method_definition")

# Unified R4 figures (fig1, fig2, fig3) with consistent title/font/color.
suppressMessages({library(ggplot2); library(dplyr); library(reshape2); library(patchwork); library(dendextend)})

root <- disco_project_root()
out  <- disco_results_path("01_method_definition/outputs/sensitivity")
prot <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/prot_367_selected.rds"))
d    <- readRDS(disco_results_path("03_multimodal_UKB/outputs/03_prot_regeneration/UKB.Prot.DISCO.regenerated.rds"))
expl <- readRDS(file.path(out, "proteomics_explore.rds"))
X <- as.matrix(prot[, -1]); colnames(X) <- expl$names; X <- scale(X)
age <- d$age; disco <- d$mc

# unified themes
tt <- theme(plot.title = element_text(size = 14, face = "bold"))
ax <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12, colour = "black"))
colpal <- c("Inflammation / TNF-kB"="#E64B35FF","Immune / cytokines"="#4DBBD5FF",
            "Cardiac / natriuretic"="#00A087FF","Senescence / neuro"="#F39B7FFF",
            "Bone / fibrinolysis"="#8491B4FF","Mixed"="grey80")

# ================= FIG 1: structure =================
young <- which(age <= 40); old <- which(age >= 65)
Rref <- cor(X[young, , drop = FALSE])
ev <- eigen(Rref, symmetric = TRUE)$values
p1a <- ggplot(data.frame(k=1:20, lambda=ev[1:20]), aes(k, lambda)) +
  geom_col(fill="#4DBBD5FF", colour="black", width=0.6) +
  geom_text(aes(label=round(lambda)), vjust=-0.3, size=3) +
  theme_classic() + tt + ax +
  labs(x="Co-expression mode", y="Eigenvalue", title="Reference is modular") +
  annotate("text", x=14, y=max(ev[1:20])*0.9, size=3.5, hjust=0, label=sprintf("effective rank %d/%d", sum(ev>1), length(ev)))
delta <- colMeans(X[old,,drop=FALSE]) - colMeans(X[young,,drop=FALSE]); delta <- delta/sqrt(sum(delta^2))
U <- eigen(Rref, symmetric=TRUE)$vectors
c2 <- (t(U)%*%delta)^2; c2 <- c2/sum(c2); cum_real <- cumsum(c2)
set.seed(1); cum_null <- matrix(NA,200,length(ev))
for(rep in 1:200){dn<-sample(delta);c2n<-(t(U)%*%dn)^2;c2n<-c2n/sum(c2n);cum_null[rep,]<-cumsum(c2n)}
dfc <- data.frame(k=1:length(ev), real=cum_real, null=colMeans(cum_null))
p1b <- ggplot(dfc, aes(k)) +
  geom_line(aes(y=real, colour="Aging deviation")) + geom_line(aes(y=null, colour="Random")) +
  scale_colour_manual(values=c("Aging deviation"="#E64B35FF","Random"="#4DBBD5FF")) +
  coord_cartesian(xlim=c(0,60)) + theme_classic() + tt + ax +
  labs(x="Top-K modes", y="Cumulative deviation", colour=NULL, title="Dysregulation is structured") +
  theme(legend.position=c(0.75,0.3), legend.background=element_blank()) +
  annotate("text", x=30, y=0.4, size=3.5, label=sprintf("~%.0f vs %.0f random modes", 1/sum(c2^2), length(ev)))
fig1 <- p1a | p1b
ggsave(file.path(out,"R4_fig1_structure.pdf"), fig1, width=9.5, height=4.3)

# ================= FIG 2: modules =================
r <- readRDS(file.path(out, "R4_modules_eigengene.rds"))
module <- r$module; res <- r$res; eig <- r$eig
bio_map <- c("ME8"="Inflammation / TNF-kB","ME7"="Immune / cytokines","ME5"="Cardiac / natriuretic",
             "ME3"="Senescence / neuro","ME4"="Bone / fibrinolysis")
bio <- function(m) ifelse(m %in% names(bio_map), unname(bio_map[m]), "Mixed")
ord <- res$module[order(-abs(res$r_disco))]

hc <- hclust(as.dist(1 - abs(cor(X))), method="complete")
dend <- as.dendrogram(hc) %>% color_branches(k=12) %>% set("labels_cex", 0.05)
ggd <- as.ggdend(dend)
p2a <- ggplot(ggd$segments) +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), colour=ifelse(is.na(ggd$segments$col),"grey60",ggd$segments$col), linewidth=0.3) +
  coord_flip() + scale_y_reverse(expand=c(0.1,0)) + theme_void() + tt +
  labs(title="A. Module derivation")

mods <- names(eig); E <- do.call(cbind, eig); colnames(E) <- mods
Me <- cor(E)[ord, ord]; Me[upper.tri(Me)] <- NA
hm <- reshape2::melt(Me, na.rm=TRUE); colnames(hm) <- c("m1","m2","cor")
hm$m1 <- factor(hm$m1, levels=ord); hm$m2 <- factor(hm$m2, levels=rev(ord)); lab <- sapply(ord, bio)
p2b <- ggplot(hm, aes(m1,m2,fill=cor)) +
  geom_tile(colour="white") + geom_text(aes(label=sprintf("%.2f",cor)), size=2.5) +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="cor") +
  scale_x_discrete(labels=lab, position="top") + scale_y_discrete(labels=rev(lab)) +
  theme_minimal(base_size=10) + tt +
  theme(axis.text.x=element_text(angle=45,hjust=0,size=8), axis.text.y=element_text(size=8), panel.grid=element_blank()) +
  labs(x=NULL, y=NULL, title="B. Module relationships")

res$bio <- sapply(res$module, bio); res <- res[order(-abs(res$r_disco)),]; res$module <- factor(res$module, levels=res$module)
p2c <- ggplot(res, aes(module, abs(r_disco), fill=bio)) +
  geom_col(width=0.65, colour="black") + geom_text(aes(label=bio), angle=90, hjust=-0.05, size=3) +
  scale_fill_manual(values=colpal, guide="none") + coord_cartesian(ylim=c(0,0.85)) +
  theme_classic() + tt + ax + labs(x="Module", y="|Correlation with DISCO|", title="C. Module vs DISCO")

prot_order <- order(module)
Rsub <- cor(X)[prot_order, prot_order]; Rsub[upper.tri(Rsub)] <- NA
boundaries <- cumsum(as.numeric(table(module)))
hm2 <- reshape2::melt(Rsub, na.rm=TRUE); colnames(hm2) <- c("i","j","cor")
p2d <- ggplot(hm2, aes(i,j,fill=cor)) +
  geom_tile() + scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="cor") +
  geom_vline(xintercept=boundaries+0.5, colour="black", linewidth=0.2) +
  geom_hline(yintercept=boundaries+0.5, colour="black", linewidth=0.2) +
  coord_fixed() + theme_void() + tt + theme(legend.position="right") +
  labs(title="D. Module biology (proteins)")
fig2 <- (p2a | p2b) / (p2c | p2d)
ggsave(file.path(out,"R4_fig2_modules.pdf"), fig2, width=13, height=9.5)

# ================= FIG 3: blockwise =================
cindex <- data.frame(block=c("within-clinical","within-proteomics","cross","full"),
                     ci=c(0.6317,0.6872,0.6701,0.6882))
cindex$block <- factor(cindex$block, levels=c("within-clinical","within-proteomics","cross","full"))
fig3 <- ggplot(cindex, aes(block, ci, fill=block)) +
  geom_col(width=0.6, colour="black") + geom_text(aes(label=sprintf("%.3f",ci)), vjust=-0.4, size=4) +
  scale_fill_manual(values=c("within-clinical"="#E64B35FF","within-proteomics"="#4DBBD5FF","cross"="#00A087FF","full"="grey40"), guide="none") +
  coord_cartesian(ylim=c(0.60,0.71)) + theme_classic() + tt + ax +
  labs(x=NULL, y="Mortality C-index", title="Cross-modality carries signal") +
  annotate("text", x=4.5, y=0.635, size=3.5, hjust=0, label="cross independent: LRT p=0.015")
ggsave(file.path(out,"R4_fig3_blockwise.pdf"), fig3, width=6.5, height=4.8)
cat("DONE\n")
