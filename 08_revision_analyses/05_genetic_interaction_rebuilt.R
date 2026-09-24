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

# R05: DISCO x PRS genetic interaction on 28 incident diseases (R1.6).
# Non-case censoring time filled with all-cause follow-up time.

options(stringsAsFactors = FALSE)
root <- disco_project_root()
out_dir <- "outputs/R05_genetic_interaction_rebuilt"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(library(data.table)); suppressMessages(library(survival))

prs <- fread(disco_external_path("DISCO_PRS_FILE", "UKB_Data/UKB_PRS_full.csv"), select = c("eid",
  "26204-0.0","26206-0.0","26210-0.0","26212-0.0","26214-0.0","26218-0.0","26220-0.0",
  "26223-0.0","26225-0.0","26227-0.0","26229-0.0","26232-0.0","26244-0.0","26248-0.0",
  "26252-0.0","26254-0.0","26258-0.0","26260-0.0","26265-0.0","26267-0.0","26269-0.0",
  "26273-0.0","26275-0.0","26278-0.0","26283-0.0","26285-0.0","26287-0.0","26289-0.0"))
setnames(prs, c("eid","AMD","AD","AST","AF","BD","CRC","BC","CVD","CED","CAD","CD","EOC",
  "HT","ISS","MEL","MS","OP","PD","POAG","PC","PSO","RA","SCZ","SLE","T1D","T2D","UC","VTE"))

inc <- readRDS(file.path(root, "UKB_Data/dat_disease_28_rebuilt.rds"))
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

run_version <- function(disc, label) {
  m <- merge(disc[, c("ID","mc","age","sex","time")], prs, by.x="ID", by.y="eid", all.x=TRUE)
  m <- merge(m, inc, by.x="ID", by.y="eid", all.x=TRUE)
  m$mc_sd <- as.numeric(scale(m$mc)); m$sex_num <- as.numeric(m$sex=="Male")
  res <- list()
  for (nm in names(diseases)) {
    spec <- diseases[[nm]]
    mm <- m
    if (!is.null(spec$sex)) mm <- mm[mm$sex==spec$sex, ]
    tm <- mm[[paste0("time_",spec$d)]]; st <- mm[[paste0("status_",spec$d)]]
    tm[is.na(tm)] <- mm$time[is.na(tm)]          # 非病例用全因随访时间作删失
    st[is.na(st)] <- 0
    # 排除患病病例（baseline 前已患病，time<0）
    prev <- tm < 0 & !is.na(tm)
    tm[prev] <- NA; st[prev] <- NA
    mm$prs_sd <- as.numeric(scale(mm[[spec$prs]]))
    cv <- if (is.null(spec$sex)) c("age","sex_num") else "age"
    ok <- complete.cases(mm[, c("mc_sd","prs_sd",cv)]) & !is.na(tm) & !is.na(st)
    dd <- data.frame(tm=tm[ok], st=st[ok], mc=mm$mc_sd[ok], prs=mm$prs_sd[ok], age=mm$age[ok])
    if (!is.null(spec$sex)) {} else dd$sex_num <- mm$sex_num[ok]
    if (nrow(dd) < 500 || sum(dd$st) < 5) next
    rhs0 <- paste(c("mc","prs",cv), collapse="+"); rhs1 <- paste(c("mc","prs","mc:prs",cv), collapse="+")
    f0 <- coxph(as.formula(paste0("Surv(tm,st)~",rhs0)), data=dd)
    f1 <- coxph(as.formula(paste0("Surv(tm,st)~",rhs1)), data=dd)
    s0 <- summary(f0); s1 <- summary(f1); int_row <- which(rownames(s1$coefficients)=="mc:prs")
    res[[length(res)+1L]] <- data.frame(
      version=label, disease=nm, prs=spec$prs, n=nrow(dd), events=sum(dd$st),
      disco_hr=exp(coef(f0)["mc"]), disco_lo=exp(confint(f0)["mc",1]), disco_hi=exp(confint(f0)["mc",2]), disco_p=s0$coefficients["mc","Pr(>|z|)"],
      prs_hr=exp(coef(f0)["prs"]), prs_lo=exp(confint(f0)["prs",1]), prs_hi=exp(confint(f0)["prs",2]), prs_p=s0$coefficients["prs","Pr(>|z|)"],
      interaction_hr=exp(coef(f1)["mc:prs"]), interaction_lo=exp(confint(f1)["mc:prs",1]), interaction_hi=exp(confint(f1)["mc:prs",2]), interaction_p=s1$coefficients[int_row,"Pr(>|z|)"],
      c_index=unname(s0$concordance[1]), c_index_se=unname(s0$concordance[2]),
      stringsAsFactors=FALSE)
  }
  do.call(rbind, res)
}
res <- rbind(run_version(disc_clin,"clinical"), run_version(disc_prot,"proteomics"))
write.csv(res, file.path(out_dir,"R05_genetic_interaction.csv"), row.names=FALSE)
cat("R05 done:", nrow(res), "rows (clinical", sum(res$version=="clinical"), "+ proteomics", sum(res$version=="proteomics"), ")\n")
