#!/usr/bin/env Rscript
# Extract numeric objects, not PDF pixels. No model/source object is overwritten.
options(stringsAsFactors=FALSE,device=function(...)grDevices::pdf(file=NULL,...))
.libPaths(c(Filter(nzchar,Sys.getenv('DISCO_R_LIBRARY')),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(jsonlite);library(digest);library(survival);library(reshape2);library(igraph)})
root<-Sys.getenv('DISCO_RESULTS_INPUT');disco<-Sys.getenv('DISCO_INPUT_ROOT')
out<-Sys.getenv('DISCO_RUN_MODULE')
entries<-list();errors<-list()
sha<-function(p)digest(file=p,algo='sha256')
logical_path<-function(p){
 p<-normalizePath(p,mustWork=TRUE)
 if(nzchar(Sys.getenv('DISCO_SIMULATION_SOURCE'))&&p==normalizePath(Sys.getenv('DISCO_SIMULATION_SOURCE')))return('caller:simulation_source')
 for(k in c('DISCO_RUN_MODULE','DISCO_INPUT_ROOT','DISCO_RESULTS_INPUT')){
  b<-Sys.getenv(k);if(nzchar(b)){b<-normalizePath(b,mustWork=FALSE);if(startsWith(p,paste0(b,'/')))return(paste0('caller:',k,'/',substring(p,nchar(b)+2L)))}
 }
 basename(p)
}
emit<-function(panel,component,d,src,def,privacy='aggregate_public_candidate',complete=TRUE,selector='explicit object member'){
 if(is.null(d)||!length(d))stop('Empty actual source: ',panel,'/',component)
 if(is.matrix(d))d<-as.data.frame(d)
 if(!is.data.frame(d))d<-data.frame(name=if(is.null(names(d)))seq_along(d)else names(d),value=as.vector(d))
 d<-as.data.frame(d)
 for(n in names(d))if(inherits(d[[n]],'haven_labelled'))d[[n]]<-as.numeric(unclass(d[[n]]))
 for(n in names(d))if(is.list(d[[n]]))d[[n]]<-vapply(d[[n]],function(x)toJSON(x,auto_unbox=TRUE,na='null'),character(1))
 if(privacy=='aggregate_public_candidate'&&any(tolower(names(d))%in%c('eid','seqn','id','id1','model_ids')))stop('Participant identifier in public data')
 folder<-file.path(out,'SourceData',if(privacy=='private_participant')'private'else'public',panel)
 dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 if(privacy=='private_participant')Sys.chmod(folder,mode='0700')
 dest<-file.path(folder,paste0(component,'.csv'));write.csv(d,dest,row.names=FALSE,na='')
 if(privacy=='private_participant')Sys.chmod(dest,mode='0600')
 entries[[length(entries)+1L]]<<-list(panel=panel,component=component,source=logical_path(src),source_sha256=sha(src),output=logical_path(dest),output_sha256=sha(dest),rows=nrow(d),columns=names(d),selection='current',privacy=privacy,permission=if(privacy=='private_participant')'Restricted individual data; NOT for public ZIP/XLSX or upload without cohort/data-owner authorization'else'Aggregate-only publication candidate; data-owner publication permission still required',definition=def,filter=selector,coverage=if(complete)'complete_component'else'partial_component')
}
run<-function(panel,f){scope<-Sys.getenv('DISCO_EXPORT_SCOPE');if(!nzchar(scope))stop('Set explicit DISCO_EXPORT_SCOPE');if(scope=='ED1_original'&&panel!='ED1_original')return(invisible(NULL));tryCatch({f();cat('EXPORTED',panel,'\n')},error=function(e){errors[[panel]]<<-conditionMessage(e);cat('UNRESOLVED',panel,conditionMessage(e),'\n')})}
read<-function(p)readRDS(file.path(root,p))
km<-function(s){
 fields<-c('time','n.risk','n.event','n.censor','surv','std.err','lower','upper')
 d<-as.data.frame(s[intersect(fields,names(s))]);d$stratum<-rep(names(s$strata),s$strata)
 d
}
matlong<-function(m){d<-as.data.frame(as.table(m));names(d)<-c('row','column','value');d}
private_points<-function(panel,component,d,src,def,complete=TRUE){
 d<-d[,!tolower(names(d))%in%c('id','id1','seqn','eid'),drop=FALSE]
 emit(panel,component,d,src,def,'private_participant',complete=complete)
}
raw_def<-'Actual retained individual plotting coordinates/score fields; restricted, not public SourceData. No IDs exported; medical measurements still participant-level.'

run('Fig1_b',function(){src<-file.path(root,'09_figures_tables/outputs/evidence_reconstruction_sim/20260917T070825Z_verified/pass_2/numeric_objects.rds');z<-readRDS(src)
 emit('Fig1_b','noise_simulation_points',z$Fig1b_drawing,src,'Synthetic noise grid actual reconstructed drawing; not historical final source recovered.',selector='Fig1b_drawing')
 emit('Fig1_b','boxplot_statistics',z$Fig1b_built$box[[1]],src,'Actual ggplot_build boxplot statistics; synthetic values; whiskers/outliers retained.',selector='Fig1b_built$box[[1]]')
 emit('Fig1_b','correlation_annotations',z$Fig1b_rho,src,'Actual synthetic Spearman annotations.',selector='Fig1b_rho')})
run('Fig1_c',function(){src<-file.path(root,'02_primary_clinical/reproduction/Fig1c/numeric_objects.rds');z<-readRDS(src)
 private_points('Fig1_c','individual_age_DM_DISCO_points',z$plot_data,src,raw_def)
 emit('Fig1_c','Spearman_annotations',c(age_DISCO=z$age_correlation,logDM_DISCO=z$hd_correlation),src,'Actual saved Spearman plot annotations; full scatter coordinates separately restricted.',complete=FALSE)})
run('Fig1_d',function(){src<-file.path(root,'02_primary_clinical/reproduction/Fig1d/numeric_objects.rds');z<-readRDS(src)
 emit('Fig1_d','attained_age_KM_curves',km(z$survival),src,'Actual original4marginalQ1/Q4 curves; attained-age delayedentry;95%logKMCI. Quartile cohorts can overlap between metrics.',selector='survival')
 emit('Fig1_d','quartile_thresholds',data.frame(metric=rep(c('mc','hd'),each=5),quantile=rep(seq(0,1,.25),2),threshold=unlist(z$quartiles,use.names=FALSE)),src,'Original marginal quartile cutoffs; not standardization SD.')})
run('ED2_a',function(){src<-file.path(root,'09_figures_tables/outputs/strict_extensions_clinical/20260917T104101_clinical_original_blocks_verified/pass_2/ED2a/numeric_objects.rds');z<-readRDS(src)
 private_points('ED2_a','individual_DM_entropy_points',z$plot_data[,c('hd','en')],src,raw_def)
 emit('ED2_a','Spearman_annotation',z$hd_entropy_correlation,src,'Actual Spearman(logDM,Entropy) correlation annotation; pointcloud private.',complete=FALSE)})
for(letter in c('d','e','f'))local({letter<-letter;run(paste0('ED2_',letter),function(){src<-file.path(root,paste0('00_manifest/figures_two_pass/20260928/numeric/pass_2/original_ED2',letter,'/ED2',letter,'_plot_data.rds'));z<-readRDS(src)
 private_points(paste0('ED2_',letter),'individual_DM_entropy_points',z$data,src,raw_def)
 emit(paste0('ED2_',letter),'Spearman_annotation',cor(log(z$data$hd),z$data$en,method='spearman',use='pairwise.complete.obs'),src,'Source-input Spearman(logDM,Entropy) annotation recomputed from actualsource, not extracted PDF; actualcoordinates restricted.',complete=FALSE)})})
for(letter in c('b','c','d'))local({letter<-letter;run(paste0('ED3_',letter),function(){src<-file.path(root,paste0('04_external_validation/reproduction/ED3',letter,'/numeric_objects.rds'));z<-readRDS(src)
 private_points(paste0('ED3_',letter),'age_DISCO_individual_points',z$age_plot_data,src,raw_def)
 private_points(paste0('ED3_',letter),'DM_DISCO_individual_points',z$hd_plot_data,src,raw_def)
 emit(paste0('ED3_',letter),'attained_age_KM_curves',km(z$survival),src,if(letter=='b')'Current historicalCHARLSwave-assigned2/4/6death,8censor; precise followupNOTverified. Attained-age4marginalcurves,95%KMCI.20261008 interval endpoints separately notselected.'else'Actual original attained-age delayed-entry marginal4KMcurves;95%logKMCI.',selector='survival')})})
run('Fig3_b',function(){src<-file.path(root,'02_primary_clinical/reproduction/Fig3b/Fig3b_KM/numeric_objects.rds');z<-readRDS(src)
 emit('Fig3_b','six_KM_curves',km(z$survival),src,'Actual six riskgroups, NHANES mortality, age-time delayedentry;95%logKMCI.',complete=FALSE,selector='survival')
 private_points('Fig3_b','individual_model_plot_coordinates',z$plot_data,src,raw_def)})
run('ED4_b',function(){src<-file.path(root,'01_method_definition/ED4_sensitivity/20260928/ED4b_smaller_v4/assets/plot_source.rds');z<-readRDS(src)
 for(n in names(z$plots)){
  p<-z$plots[[n]];emit('ED4_b',paste0(n,'_coordinates'),p$data,src,'Actual four-column distribution/correlation plot coordinates; m/se summary quantities; dualC-axis supplementarylinebelow.',selector=paste0('plots$',n,'$data'))
  for(i in seq_along(p$layers))if(is.data.frame(p$layers[[i]]$data))emit('ED4_b',paste0(n,'_layer',i),p$layers[[i]]$data,src,'Explicit layer actualsource data; dual-axis C uses producer affine display transformation, nativeC summary separately.',selector=paste0('plots$',n,'$layers[[',i,']]$data'))
 }
 emit('ED4_b','correlation_Cindex_native_summary',z$summary,src,'Mean acrossdraws; Corr_SE/C_SE repeated-estimate SE; mean_model_SE explicitlydistinct; dualaxis display not new scientific variable.')})
run('ED4_c',function(){src<-file.path(root,'01_method_definition/ED4_sensitivity/20260928/ED4abcd_typography_v3/assets/ED4c_plot_source.rds');z<-readRDS(src)
 emit('ED4_c','reference_age_Cindex',z$data,src,'SameNHANES evaluation; ageweighted/unweighted crudeC,95%normalCI C+/-1.96SE; thresholds equallyspaced display; actualmodel fieldretained.')})
run('ED6_a',function(){src<-file.path(root,'04_external_validation/HRS_imputed_reference/20261001/median/plot_source.rds');z<-readRDS(src)
 for(n in c('original','accel'))emit('ED6_a',paste0('HRS35_',n,'_Cindex'),z$data$HRS35[[n]],src,'Actual selected HRS35 internal<60median-imputed reference; Originalmetric+age vs clocks; Accel age residuals. Errorbar +/-modelSE, not95%CI.',selector=paste0('data$HRS35$',n))})
run('ED6_b',function(){src<-file.path(root,'04_external_validation/HRS_imputed_reference/20261001/median/ED6_update/plot_source.rds');z<-readRDS(src)
 emit('ED6_b','HRS35_six_KM_curves',z$b$data,src,'Actual selected HRS35 Trackerverified endpoint attained-age delayedentry6groups;95%logKMCI.',selector='b$data')})
run('Fig4_b',function(){src<-file.path(root,'05_organ_network/reproduction/Fig4b/Fig4b_plot_data.rds');z<-readRDS(src)
 emit('Fig4_b','dual_axis_model_coordinates',z$df,src,'Actual plotted Gompertz/age summary dualaxiscoordinates, NOTindividualdata. y-field semantics fromPlot_Fig4C.R327–425.',selector='df')
 emit('Fig4_b','age_summary',z$nd11,src,'Actual age-summary mean,se,logmean fromoriginal producer.',selector='nd11')
 emit('Fig4_b','Gompertz_parameters',z$gmp,src,'Actual original curve parameters.',selector='gmp')})
run('Fig4_c',function(){src<-file.path(root,'09_figures_tables/outputs/original_block_replay/Fig4c_CANDIDATE_plot_data.rds');z<-readRDS(src)
 emit('Fig4_c','protein_graph_scatter',z,src,'Actual saved-network per-protein mean/age-correlation/degree scatterdata. DoesNOTcertify currentrawproteinupstreamnetwork chain;20261008newgraph separatelycandidate.',selector='whole protein-level table')})
run('Fig4_d',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/shared_selected_protein_partial_correlations.csv');z<-read.csv(src)
 emit('Fig4_d','organ_pair_partial_correlations',z,src,'153uniqueorganpairs; age+sex+sharedselectedproteins adjusted partialPearson, exactactualsharedmatrix branch; trianglefilters byoriginalordering; noCI plotted.')})
run('Fig4_e',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/Fig4e__THRESHOLD_0p5_RECALCULATED.rds');z<-readRDS(src)
 emit('Fig4_e','network_edges',as_data_frame(z,what='edges'),src,'Actual undirected71edge network fromsharedproteinadjusted r>.5; not K3alternative; graphlayoutnot reestimated.')
 nodes<-as_data_frame(z,what='vertices');nodes$degree<-degree(z);nodes$layout_x<-z$layout[,1];nodes$layout_y<-z$layout[,2]
 emit('Fig4_e','network_nodes',nodes,src,'Actual18organ nodes degree and savedlayout; aggregateonly.')})
run('Fig4_f',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/organ36_provenance_and_models.csv');z<-read.csv(src);org<-readRDS(file.path(root,'05_organ_network/corrected_organ/20260926/locked_design.rds'))$organs
 z<-z[match(org,z$organ),];m<-read('05_organ_network/corrected_organ/20260926/pass_2/specific_panel_nooverlap_partial_exact.rds');z$strength<-rowSums(m)-1
 dci<-readxl::read_xlsx(file.path(disco,'Result3/Organ/Disease_cindex.xlsx'));dn<-na.omit(unique(dci[[1]]))[-1];cv<-as.numeric(na.omit(dci[[8]])[-1]);da<-aggregate(ci~tissue,data.frame(tissue=rep(dn,each=8),ci=cv),mean);z$disease_C<-da$ci[match(z$organ,da$tissue)]
 emit('Fig4_f','three_organ_association_coordinates',z[,c('organ','agecor','C_crude','disease_C','strength')],src,'Actual three organ-scatter coordinates fromcorrectedscore/nooverlapmatrix; diseaseC remains historicalpreservedworkbook, not endpointrecomputed. strength=sumr minusdiag. OriginalpointwiseLMbanddefault95%; notacrossparticipantcurve.')})
run('ED9_a',function(){src<-file.path(root,'00_manifest/figures_two_pass/20260928/numeric/pass_2/original_ED9a_selected/ED9a_selected_plot_data.rds');emit('ED9_a','after_selection_counts',readRDS(src),src,'Actual18organ before/after selectedcount subregion; beforepanel rawtissuecount separatelyrequired.',complete=FALSE)
 src<-file.path(disco,'Result3/Organ/UKB.Prot.class.rds');z<-readRDS(src);emit('ED9_a','before_selection_tissue_counts',as.data.frame(table(z$specific.tissue)),src,'Tissuecount source; finaloriginalPlot_Fig4A ordering/filter mustbecrosschecked.',complete=FALSE)})
run('ED9_c',function(){src<-file.path(root,'05_organ_network/corrected_annotations/20260926/pass_2/ED9c_results.rds');z<-readRDS(src);d<-matlong(z$r);d$P<-as.vector(z$p);d$stars<-as.vector(z$stars)
 emit('ED9_c','clinical_organ_heatmap_cells',d,src,'Current saved18x25 r/P/stars withcorrectedannotation alignment; actualcorrelationadjustment fromproducer remainsnotfull20261008rawrebuild.',selector='r,p,stars matrices sameposition')})
run('ED9_d',function(){src<-file.path(root,'05_organ_network/corrected_annotations/20260926/pass_2/ED9d_results.rds');z<-readRDS(src)
 emit('ED9_d','three_forest_series',z$results,src,'Actualcorrected54bars; CorrFisherz95%CI; C-indexbars +/-oneSE (not95%CI); score-onlyvsage+sex+score exactvariables retained.',selector='results')
 emit('ED9_d','age_sex_baseline_Cindex',z$baseline,src,'Actualage+sex baseline dashedreference.')})
run('ED9_e',function(){src<-file.path(root,'05_organ_network/reproduction/ED9e/ED9e_plot_data.rds');z<-readRDS(src)
 emit('ED9_e','age_organ_fraction_cells',z$stat,src,'558actualorganage-group fractioncells; thresholdhighrisk rule fromproducer, not perindividualdata.')})
run('S10_a',function(){src<-file.path(root,'05_organ_network/corrected_annotations/20260926/pass_2/ED10a_results.rds');z<-readRDS(src)
 emit('S10_a','threshold_Gompertz_parameters',z$parameters,src,'Corrected48threshold/3parameter mappings; parameters actualvalue/se; errorbar +/-SE; original16thresholds retained.')})
run('S10_b',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/ED10bd_recomputed.rds');z<-readRDS(src)
 emit('S10_b','random_panel_age_replicate_mean',z$age_replicate,src,'Actual6400 aggregateagebin/randompanelmeans, NOTparticipantpoints; N peragebin, 200replicate scheme; exactproducer correlationannotation must distinguish means.',selector='age_replicate')})
run('S10_c',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/random_summary.rds');z<-readRDS(src)
 emit('S10_c','random_panel_Cindex_boxplot_values',z$stats[,c('replicate','organ','N','deaths','C_crude','SE_crude')],src,'3600originalcrudeC modelreplicates=18organsx200; boxplot empiricaldistribution not95%CI; noindividualscore rows.')})
run('S10_d',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/ED10bd_recomputed.rds');z<-readRDS(src)
 emit('S10_d','age_organ_fraction_cells',z$age_organ_mean,src,'Actual agebin/organ averaged highriskfraction across200replicates; definition retained inobject; not individualpoints.')})
run('S11_c',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/random_summary.rds');z<-readRDS(src)
 emit('S11_c','random_mean_score_partial_matrix',matlong(z$partial_of_mean),src,'Actual residualage+sex partialmatrix of mean randomscores; triangleoriginalfilter. Notmeanof200partialmatrices.',selector='partial_of_mean')})
run('S11_d',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/random_summary.rds');z<-readRDS(src);m<-read('05_organ_network/corrected_organ/20260926/pass_2/specific_panel_nooverlap_partial_exact.rds');ix<-upper.tri(m)
 d<-data.frame(organ1=rownames(m)[row(m)[ix]],organ2=colnames(m)[col(m)[ix]],random_partial=z$partial_of_mean[ix],specific_partial=m[ix])
 emit('S11_d','organ_pair_scatter153',d,src,'Actual153dyads; plotQAP9999organlabelpermutations preservesdyadicdependence; notnominalindependentpairP.',selector='upper.tri(partial_of_mean), upper.tri(specific nooverlap matrix)')})
run('S11_e',function(){src<-file.path(root,'05_organ_network/corrected_organ/20260926/pass_2/organ_meanC_strength.csv');emit('S11_e','organ_connectivity_mortality_scatter18',read.csv(src),src,'Actual18aggregate organ meansfrom200randommodels; currentplot uses C_crude andsignedstrength; adjustedC notsubstituted.')})
run('ED4_g',function(){for(branch in c('clinical','clinicalprot')){
 src<-if(branch=='clinical')file.path(disco,'Result1/UKB.clinical.sim2.rds')else file.path(disco,'Result1/Num/UKB.clinicalprot.sim2.rds');z<-readRDS(src)
 emit('ED4_g',paste0(branch,'_dimension_Cindex_replicates'),z,src,'ActualoriginalTest2.R retainedreplicate C-index vsfeaturecount; scoremodelreplicates aggregateonly; loess defaultCI95%notmodelC CI.',selector='entire panel-specific input table')
}})
for(num in 1:5)local({num<-num;run(paste0('S',num,'_all'),function(){src<-file.path(root,'01_method_definition/corrected_simulation/20261007_signed_half/pass_2/numeric_objects.rds');z<-readRDS(src);keys<-grep(paste0('^S',num,'(_|$)'),names(z),value=TRUE)
 for(n in keys)emit(paste0('S',num,'_all'),n,z[[n]],src,'Currentcorrectedsyntheticgeometry values; signedEntropyhalf convention, correctedinversePCA/noisewise rowmap/SDlabel. Syntheticcoordinates public; actualSnumber-specificmembers only.',selector=n)
})})
write_json(entries,file.path(out,'results/r_exports.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
write_json(errors,file.path(out,'results/r_export_errors.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'logs/R_sessionInfo.txt'))
cat('TOTAL components',length(entries),'ERRORS',length(errors),'\n')
