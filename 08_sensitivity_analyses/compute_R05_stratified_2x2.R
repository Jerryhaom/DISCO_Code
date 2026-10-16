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

# R05 visualization data (2/3 & 3/3): stratified DISCO HR and 2x2 group HR.
# Re-fits Cox models; reuses the R05 data-loading logic.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
wkd  <- disco_results_path("08_sensitivity_analyses")
out  <- file.path(wkd, "outputs/R05_genetic_interaction")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressMessages(library(data.table)); suppressMessages(library(survival))

prs <- fread(disco_external_path("DISCO_PRS_FILE", "UKB_Data/UKB_PRS_full.csv"), select = c("eid",
  "26204-0.0","26206-0.0","26210-0.0","26212-0.0","26214-0.0","26218-0.0","26220-0.0",
  "26223-0.0","26225-0.0","26227-0.0","26229-0.0","26232-0.0","26244-0.0","26248-0.0",
  "26252-0.0","26254-0.0","26258-0.0","26260-0.0","26265-0.0","26267-0.0","26269-0.0",
  "26273-0.0","26275-0.0","26278-0.0","26283-0.0","26285-0.0","26287-0.0","26289-0.0"))
setnames(prs, c("eid","AMD","AD","AST","AF","BD","CRC","BC","CVD","CED","CAD","CD","EOC",
  "HT","ISS","MEL","MS","OP","PD","POAG","PC","PSO","RA","SCZ","SLE","T1D","T2D","UC","VTE"))

inc <- readRDS(file.path(root, "UKB_Data/dat_disease_28.rds"))
disc_clin <- readRDS(file.path(root, "Result1/UKB.Clinical.DISCO.rds"))
disc_prot <- readRDS(file.path(root, "Result2/UKB.Prot.DISCO.rds"))

diseases <- list(
  amd=list(prs="AMD",d="amd"), ad=list(prs="AD",d="ad"), asthma=list(prs="AST",d="asthma"),
  af=list(prs="AF",d="af"), bipolar=list(prs="BD",d="bipolar"), bowel_ca=list(prs="CRC",d="bowel_ca"),
  breast_ca=list(prs="BC",d="breast_ca",sex="Female"), cvd=list(prs="CVD",d="cvd"),
  coeliac=list(prs="CED",d="coeliac"), cad=list(prs="CAD",d="cad"), crohn=list(prs="CD",d="crohn"),
  ovarian_ca=list(prs="EOC",d="ovarian_ca",sex="Female"), hypertension=list(prs="HT",d="hypertension"),
  ischaemic_stroke=list(prs="ISS",d="ischaemic_stroke"), melanoma=list(prs="MEL",d="melanoma"),
  ms=list(prs="MS",d="ms"), osteoporosis=list(prs="OP",d="osteoporosis"),
  parkinson=list(prs="PD",d="parkinson"), glaucoma=list(prs="POAG",d="glaucoma"),
  prostate_ca=list(prs="PC",d="prostate_ca",sex="Male"), psoriasis=list(prs="PSO",d="psoriasis"),
  ra=list(prs="RA",d="ra"), schizophrenia=list(prs="SCZ",d="schizophrenia"),
  sle=list(prs="SLE",d="sle"), t1d=list(prs="T1D",d="t1d"), t2d=list(prs="T2D",d="t2d"),
  ulcerative_colitis=list(prs="UC",d="ulcerative_colitis"), vte=list(prs="VTE",d="vte"))

prep <- function(disc, prs) {
  m <- merge(disc[, c("ID","mc","age","sex","time")], prs, by.x="ID", by.y="eid", all.x=TRUE)
  m <- merge(m, inc, by.x="ID", by.y="eid", all.x=TRUE)
  m$mc_sd <- as.numeric(scale(m$mc)); m$sex_num <- as.numeric(m$sex=="Male")
  m
}
clin <- prep(disc_clin, prs); prot <- prep(disc_prot, prs)

# ---- scheme 2: DISCO HR stratified by PRS Q1 vs Q4 ----
strat_hr <- function(m, label) {
  res <- list()
  for (nm in names(diseases)) {
    spec <- diseases[[nm]]
    mm <- m
    if (!is.null(spec$sex)) mm <- mm[mm$sex==spec$sex, ]
    tm <- mm[[paste0("time_",spec$d)]]; st <- mm[[paste0("status_",spec$d)]]
    tm[is.na(tm)] <- mm$time[is.na(tm)]; st[is.na(st)] <- 0
    prev <- tm < 0 & !is.na(tm); tm[prev] <- NA; st[prev] <- NA
    mm$prs_raw <- mm[[spec$prs]]
    cv <- if (is.null(spec$sex)) c("age","sex_num") else "age"
    ok <- complete.cases(mm[, c("mc_sd","prs_raw",cv)]) & !is.na(tm) & !is.na(st) & !is.na(mm$prs_raw)
    dd <- data.frame(tm=tm[ok], st=st[ok], mc=mm$mc_sd[ok], prs=mm$prs_raw[ok], age=mm$age[ok])
    if (!is.null(spec$sex)) {} else dd$sex_num <- mm$sex_num[ok]
    if (nrow(dd) < 500 || sum(dd$st) < 5) next
    q <- quantile(dd$prs, probs = c(0.25, 0.75), na.rm = TRUE)
    for (stratum in c("Q1","Q4")) {
      s <- if (stratum == "Q1") dd[dd$prs <= q[1], ] else dd[dd$prs >= q[2], ]
      rhs <- if (is.null(spec$sex)) "mc+age+sex_num" else "mc+age"
      f <- coxph(as.formula(paste0("Surv(tm,st)~",rhs)), data=s)
      sm <- summary(f)
      res[[length(res)+1L]] <- data.frame(version=label, disease=nm, prs=spec$prs,
        stratum=stratum, n=nrow(s), events=sum(s$st),
        disco_hr=exp(coef(f)["mc"]), disco_lo=exp(confint(f)["mc",1]),
        disco_hi=exp(confint(f)["mc",2]), disco_p=sm$coefficients["mc","Pr(>|z|)"],
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, res)
}
strat <- rbind(strat_hr(clin,"clinical"), strat_hr(prot,"proteomics"))
write.csv(strat, file.path(out,"R05_stratified_disco_hr.csv"), row.names=FALSE)
cat("stratified:", nrow(strat), "rows\n")

# ---- scheme 3: 2x2 group HR (DISCO median x PRS median) ----
twobytwo <- function(m, label, sel_diseases) {
  res <- list()
  for (nm in sel_diseases) {
    spec <- diseases[[nm]]
    mm <- m
    if (!is.null(spec$sex)) mm <- mm[mm$sex==spec$sex, ]
    tm <- mm[[paste0("time_",spec$d)]]; st <- mm[[paste0("status_",spec$d)]]
    tm[is.na(tm)] <- mm$time[is.na(tm)]; st[is.na(st)] <- 0
    prev <- tm < 0 & !is.na(tm); tm[prev] <- NA; st[prev] <- NA
    mm$prs_raw <- mm[[spec$prs]]
    cv <- if (is.null(spec$sex)) c("age","sex_num") else "age"
    ok <- complete.cases(mm[, c("mc_sd","prs_raw",cv)]) & !is.na(tm) & !is.na(st) & !is.na(mm$prs_raw)
    dd <- data.frame(tm=tm[ok], st=st[ok], mc=mm$mc_sd[ok], prs=mm$prs_raw[ok], age=mm$age[ok])
    if (!is.null(spec$sex)) {} else dd$sex_num <- mm$sex_num[ok]
    if (nrow(dd) < 500 || sum(dd$st) < 5) next
    dd$disco_grp <- ifelse(dd$mc <= quantile(dd$mc, 1/3, na.rm = TRUE), "Low",
                    ifelse(dd$mc >= quantile(dd$mc, 2/3, na.rm = TRUE), "High", NA))
    dd$prs_grp  <- ifelse(dd$prs <= quantile(dd$prs, 1/3, na.rm = TRUE), "Low",
                    ifelse(dd$prs >= quantile(dd$prs, 2/3, na.rm = TRUE), "High", NA))
    dd <- dd[!is.na(dd$disco_grp) & !is.na(dd$prs_grp), ]
    dd$grp <- factor(paste(dd$disco_grp, dd$prs_grp, sep="."),
                     levels=c("Low.Low","Low.High","High.Low","High.High"))
    rhs <- if (is.null(spec$sex)) "grp+age+sex_num" else "grp+age"
    f <- coxph(as.formula(paste0("Surv(tm,st)~",rhs)), data=dd)
    ci <- confint(f)
    for (g in c("Low.Low","Low.High","High.Low","High.High")) {
      sub <- dd[dd$grp == g, ]
      if (g == "Low.Low") {
        res[[length(res)+1L]] <- data.frame(version=label, disease=nm, disco=gsub("\\..*","",g),
          prs=gsub(".*\\.","",g), n=nrow(sub), events=sum(sub$st), hr=1, lo=NA, hi=NA, stringsAsFactors=FALSE)
      } else {
        res[[length(res)+1L]] <- data.frame(version=label, disease=nm, disco=gsub("\\..*","",g),
          prs=gsub(".*\\.","",g), n=nrow(sub), events=sum(sub$st), hr=exp(coef(f)[paste0("grp",g)]),
          lo=exp(ci[paste0("grp",g),1]), hi=exp(ci[paste0("grp",g),2]), stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind, res)
}
sel <- c("af","cad","hypertension","t2d","osteoporosis","vte","prostate_ca","amd")
tw <- rbind(twobytwo(clin,"clinical",sel), twobytwo(prot,"proteomics",sel))
write.csv(tw, file.path(out,"R05_2x2_group_hr.csv"), row.names=FALSE)
cat("2x2:", nrow(tw), "rows\n")
