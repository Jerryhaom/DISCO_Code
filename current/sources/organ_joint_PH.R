# Final common-test PH diagnostics; no participant rows exported.
suppressPackageStartupMessages({library(survival);library(data.table)})
folder<-file.path(Sys.getenv('DISCO_DATA_ROOT'),'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild')
m<-readRDS(file.path(folder,'PRIVATE/common_prediction_models.rds'))
rows<-list()
for(o in names(m))for(method in c('Oh','Goeminne')){
 d<-m[[o]]$test;k<-if(method=='Oh')'O' else 'G'
 d$Dz<-as.numeric(scale(resid(lm(D~age,d))))
 d$Pz<-as.numeric(scale(resid(lm(reformulate('age',k),d))))
 f<-coxph(Surv(time,status)~age+sex+Dz+Pz,d,x=TRUE,y=TRUE,model=TRUE,
     control=coxph.control(eps=1e-10,iter.max=50))
 z<-cox.zph(f,transform='km')$table
 for(term in rownames(z))rows[[length(rows)+1]]<-data.frame(organ=o,method=method,term=term,
   N=f$n,deaths=f$nevent,chisq=z[term,'chisq'],df=z[term,'df'],P=z[term,'p'])
}
r<-rbindlist(rows);r[,q_BH_18:=p.adjust(P,'BH'),by=term]
fwrite(r,file.path(Sys.getenv('DISCO_RUN_OUTPUT'),'qa/final_joint_PH.csv'))
print(r[term%in%c('Dz','Pz','GLOBAL'),.(nominal_under05=sum(P<.05),BH_under05=sum(q_BH_18<.05),models=.N),by=term])
