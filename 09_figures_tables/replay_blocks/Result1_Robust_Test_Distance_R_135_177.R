res$d[4]="Euclidean (unweighted)"

mc=readRDS("NHANES.Clinical.DISCO2.rds")
c=coxph(Surv(time,status)~hd,data = mc)
s1=summary(c)
s1$concordance

c=coxph(Surv(time,status)~mc,data = mc)
s1=summary(c)
res=rbind(res,c(s1$concordance,
                "Euclidean (weighted)"))
res$ci=as.numeric(res$ci)
res$dt=c("Affine Invariant Riemannian Metric",
         "Bhattacharyya distance",
         "Cholesky difference",
         "Euclidean distance",
         "Hellinger distance",
         "Jensen-Bregman Log Determinant Distance",
         "Kullback-Leibler Distance Measure",
         "Log Euclidean Riemannian Metric",
         "Procrustes Size and Shape measure",
         "Procrustes analysis with scale",
         "Weighted Euclidean distance")

res$ci=as.numeric(res$ci)
res$se=as.numeric(res$se)
res$dt=factor(res$dt,levels = res$dt[c(1:3,5:10,4,11)])
p1=ggplot(res,aes(ci,dt)) +
  geom_bar(stat = 'identity',width =0.8,alpha=0.9,
           fill="grey") +
  geom_errorbar(aes(xmin=ci-se,xmax=ci+se),width=0.5) +
  #facet_wrap(Var2~.) +
  coord_cartesian(xlim=c(0.6,0.78)) +
  #geom_vline(xintercept = s$concordance[1],
  #           color="red",linetype="dashed") +
  theme_classic() + xlab("C-index") + ylab("Distance") +
  labs(title = "NHANES") +
  theme(axis.text.y = element_text(size=10,color="black"),
        axis.title = element_text(size=12))

pdf("./Robust/Bench.Distance.pdf",height = 3,width = 5)
p1
dev.off()
