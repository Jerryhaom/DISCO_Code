options(warn=1)
if(nzchar(Sys.getenv('DISCO_R_LIB'))) .libPaths(c(Sys.getenv('DISCO_R_LIB'),.libPaths()))
library(survival);library(survminer);library(ggplot2)
base <- file.path(Sys.getenv('DISCO_DATA_ROOT'),'Result6')
task <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
arg <- commandArgs(trailingOnly=TRUE);pass<-if(length(arg))arg[1] else 'pass_1'
out <- file.path(task,'nmr_corrected',pass);dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs<-file.path(base,c('Rugao.Met.DISCO.rds','Rugao.Met.DISCO2.rds','Rugao.2019.mortality.rds','Plot_surv.R'))
before<-tools::md5sum(inputs)
candidate<-readRDS(file.path(task,'NMR_selected23_dimension_matched_candidate.rds'))
legacy<-readRDS(inputs[2]);corrected<-legacy
stopifnot(!anyDuplicated(candidate$ID),setequal(legacy$ID,candidate$ID))
corrected$mc<-candidate$score[match(corrected$ID,candidate$ID)]
stopifnot(all(is.finite(corrected$mc)),identical(legacy[,setdiff(names(legacy),'mc')],corrected[,setdiff(names(corrected),'mc')]))
saveRDS(corrected,file.path(out,'Rugao.Met.DISCO2_CORRECTED.rds'))
write.csv(data.frame(ID=corrected$ID,mc_old=legacy$mc,mc_corrected=corrected$mc),file.path(out,'corrected_NMR_scores.csv'),row.names=FALSE)

prepare<-function(pd){
 mortality<-readRDS(inputs[3]);mortality<-mortality[match(pd$ID,mortality$x2019ID),]
 mortality$time[is.na(mortality$status)]<-24;mortality$status[is.na(mortality$status)]<-0
 pd$time<-mortality$time;pd$status<-mortality$status;pd
}
legacy<-prepare(legacy);corrected<-prepare(corrected)
stopifnot(identical(legacy[,setdiff(names(legacy),'mc')],corrected[,setdiff(names(corrected),'mc')]))
results<-list()
for(version in c('legacy','corrected'))for(model in c('crude','age_sex'))for(units in c('per_raw_score_unit','per_input_population_SD')){
 pd<-get(version);pd$score<-if(units=='per_input_population_SD')as.numeric(scale(pd$mc)) else pd$mc
 form<-if(model=='crude')Surv(time,status)~score else Surv(time,status)~age+sex+score
 fit<-coxph(form,data=pd);s<-summary(fit);ci<-s$conf.int['score',]
 results[[length(results)+1]]<-data.frame(version,model,units,N=fit$n,deaths=fit$nevent,
   HR=ci[1],lower=ci[3],upper=ci[4],P=s$coefficients['score','Pr(>|z|)'],C=s$concordance[1],C_SE=s$concordance[2],
   score_SD=sd(pd$mc),row.names=NULL)
}
result<-do.call(rbind,results);write.csv(result,file.path(out,'mortality_model_comparison.csv'),row.names=FALSE)
crude_meta<-result[result$version=='corrected'&result$model=='crude'&result$units=='per_raw_score_unit',]
adjusted_meta<-result[result$version=='corrected'&result$model=='age_sex'&result$units=='per_raw_score_unit',]
stopifnot(nrow(crude_meta)==1L,nrow(adjusted_meta)==1L,
  all(is.finite(as.matrix(result[c('N','deaths','HR','C','C_SE')]))),all(result$N>=result$deaths),all(result$deaths>0))

# Evaluate only the exact original NMR plotting segment, redirecting I/O.
# The only scientific alteration is readRDS(Met.DISCO2)$mc after the proven dimension repair.
env<-new.env(parent=.GlobalEnv)
env$readRDS<-function(file,...){if(basename(file)=='Rugao.Met.DISCO2.rds')return(corrected);base::readRDS(file.path(base,file),...)}
env$saveRDS<-function(object,file,...)base::saveRDS(object,file.path(out,basename(file)),...)
env$pdf<-function(file,...)grDevices::pdf(file.path(out,basename(file)),...)
lines<-readLines(inputs[4]);block<-lines[125:235]
writeLines(block,file.path(out,'executed_original_Plot_surv_lines125_235.R'))
for(ex in parse(text=block)){
 # The unassigned interactive preview at original line204 is outside a PDF;
 # it does not contribute to either saved panel and looks in .GlobalEnv for data.
 if(is.call(ex) && identical(ex[[1]],as.name('ggsurvplot')) && grDevices::dev.cur()==1)next
 value<-withVisible(eval(ex,env))
 if(grDevices::dev.cur()>1 && value$visible && (inherits(value$value,'ggplot')||inherits(value$value,'ggsurvplot')))print(value$value)
}
saveRDS(env$fit1,file.path(out,'age_timescale_survfit.rds'))
saveRDS(env$p1,file.path(out,'Fig2d_original_ggplot_corrected.rds'))
write.csv(as.data.frame(table(env$pd.vd1$group1)),file.path(out,'extreme_quartile_group_sizes.csv'),row.names=FALSE)
write.csv(data.frame(source_inputs_unchanged=identical(before,tools::md5sum(inputs)),
                     score_dimension=length(candidate$selected_PCs),reference_N=length(candidate$reference_ID),
                     score_population_N=nrow(corrected),score_population_deaths=sum(corrected$status),
                     crude_Cox_N=crude_meta$N,crude_Cox_deaths=crude_meta$deaths,
                     age_sex_Cox_N=adjusted_meta$N,age_sex_Cox_deaths=adjusted_meta$deaths,
                     correction='dmc(met1[,bb],ref[,bb],weight); original pp<0.05 selection and reference/weights retained'),
          file.path(out,'verification.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(result)
