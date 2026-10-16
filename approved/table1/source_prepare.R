# Source preparation blocks from run_DISCO_DM_table.R; no legacy model/output loop.
suppressPackageStartupMessages({library(survival);library(openxlsx)})
root <- Sys.getenv("DISCO_DATA_ROOT")
if(!nzchar(root) || !nzchar(Sys.getenv("DISCO_RULAS_QUESTIONNAIRE"))) stop("Set DISCO_DATA_ROOT and DISCO_RULAS_QUESTIONNAIRE")
impute_mode <- function(x) { x[is.na(x)] <- names(sort(table(x),decreasing=TRUE))[1]; x }
# ---------- CHARLS ----------
ch <- readRDS(file.path(root, "Result4/CHARLS.DISCO.rds"))
d <- data.frame(mc=ch$mc, hd=ch$hd, time=ch$time, status=ch$status, age=ch$age, sex=ch$sex)
d$smoke <- impute_mode(ifelse(ch$smoke==1,1,ifelse(ch$smoke==2,0,NA)))
d$drink <- impute_mode(ifelse(ch$drink==1,1,ifelse(ch$drink==2,0,NA)))
edu <- as.numeric(ch$edu); d$education <- impute_mode(ifelse(edu %in% c(1,2,3),0,1))
charls <- d

# ---------- CLHLS ----------
clh <- readRDS(file.path(root, "Result4/CLHLS.DISCO.rds"))
lon <- readRDS(file.path(root, "Result4/CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds"))
d <- data.frame(mc=clh$mc, hd=clh$hd, time=clh$time, status=clh$status, age=clh$age, sex=clh$sex)
d71 <- as.numeric(unclass(lon$d71)); d81 <- as.numeric(unclass(lon$d81))
d$smoke <- impute_mode(ifelse(d71==1,1,ifelse(d71==2,0,NA)))
d$drink <- impute_mode(ifelse(d81==1,1,ifelse(d81==2,0,NA)))
f1 <- as.numeric(unclass(lon$f1)); f1[f1 %in% c(88,99)] <- NA
d$education <- impute_mode(ifelse(f1==0,0,1))
clhls <- d

# ---------- RuLAS wave2 ----------
ru <- readRDS(file.path(root, "Result4/RuLAS.DISCO.rds"))
sh <- read.xlsx(Sys.getenv("DISCO_RULAS_QUESTIONNAIRE"), sheet="Sheet2")
w2 <- match(ru$ID, sh[[1]])
d <- data.frame(mc=ru$mc, hd=ru$hd, time=ru$time, status=ru$status, age=ru$age, sex=ru$sex)
d$smoke <- impute_mode(ifelse(sh[[361]][w2] %in% c("2","3","水烟"),1,ifelse(sh[[361]][w2]=="1",0,NA)))
d$drink <- impute_mode(ifelse(sh[[362]][w2] %in% c("2","3"),1,ifelse(sh[[362]][w2]=="1",0,NA)))
d$education <- impute_mode(ifelse(as.numeric(sh[[101]][w2])==1,0,1))
rulas <- d

# ---------- NHANES ----------
dis2 <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO2.rds"))
dis  <- readRDS(file.path(root, "Result1/NHANES.Clinical.DISCO.rds"))
nhan <- readRDS(file.path(root, "UKB_Data/NHANES_covariates.rds")); nhan <- nhan[match(dis2$id, nhan$sampleID),]
d <- data.frame(mc=dis2$mc, hd=dis2$hd, time=dis2$time, status=dis2$status,
                age=dis$Age[match(dis2$id, dis$id)], sex=nhan$sex)
d$smoke <- impute_mode(ifelse(nhan$smoke==1,1,ifelse(nhan$smoke==0,0,NA)))
d$drink <- impute_mode(nhan$drink)
d$education <- impute_mode(ifelse(nhan$education=="high",1,ifelse(nhan$education=="low",0,NA)))
d$race <- impute_mode(nhan$race)
nhanes <- d

# ---------- UKB（education 修复：-7=无学历→low，-3=拒答→NA，1-6→按阈值）----------
ukb <- readRDS(file.path(root, "UKB_Data/extended_covariates_clinical.rds"))
raw <- readRDS(file.path(root, "UKB_Data/covariates_raw_full.rds"))
edu_raw <- as.numeric(raw$edu_q0[match(ukb$ID, raw$eid)])
race <- read.csv(file.path(root, "UKB_Data/ukb_race_21000.csv"), colClasses="character")
ukb$race <- race$race[match(ukb$ID, race$eid)]; ukb$race[ukb$race==""] <- "White"
d <- data.frame(mc=ukb$mc, hd=ukb$hd, time=ukb$time, status=ukb$status, age=ukb$age,
                sex=ifelse(ukb$sex=="Male",1,0))
d$smoke <- impute_mode(ifelse(ukb$smoke_status %in% c(1,2),1,ifelse(ukb$smoke_status==0,0,NA)))
d$drink <- impute_mode(ifelse(ukb$alcohol_freq==6,0,ifelse(ukb$alcohol_freq %in% 1:5,1,NA)))
d$education <- impute_mode(ifelse(edu_raw %in% c(1,2), 1, ifelse(edu_raw %in% c(-7,3,4,5,6), 0, NA)))
d$race <- impute_mode(ifelse(ukb$race=="White","White","Non-White"))
ukb_clin <- d
cat("UKB education 修复后: low", sum(d$education==0), " high", sum(d$education==1), " NA", sum(is.na(d$education)), "\n")

# ---------- 汇总 ----------
cohorts <- list(NHANES=list(d=nhanes,race=TRUE), "UKB clinical"=list(d=ukb_clin,race=TRUE),
                CHARLS=list(d=charls,race=FALSE), CLHLS=list(d=clhls,race=FALSE),
                "RuLAS wave 2"=list(d=rulas,race=FALSE))

for (nm in names(cohorts)) stop("Preparation-only source; use recompute_table1.R")
