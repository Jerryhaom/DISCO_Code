#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(jsonlite);library(digest);library(data.table)})
args<-commandArgs(TRUE);m<-if(length(args))normalizePath(args[1])else normalizePath('05_organ_network/S11b_full_pair_recompute_20261008');root<-Sys.getenv('DISCO_INPUT_ROOT');g<-readRDS(file.path(root,'Result3/Organ/UKB.prot.imputed.rds'));od<-readRDS(file.path(root,'Result3/Organ/UKB.odis.rds'));raw<-read.csv(file.path(m,'input/private/organ_full_scores.csv'),check.names=FALSE);bridge<-read_json(file.path(m,'results/server_core_fingerprint.json'),simplifyVector=TRUE);ix<-match(as.character(raw$eid),as.character(g[,'eid']));stopifnot(!anyNA(ix),!anyDuplicated(raw$eid));org<-gsub('.',' ',bridge$score_fields,fixed=TRUE);delta<-as.matrix(raw[,org])-as.matrix(od[ix,bridge$score_fields]);stopifnot(max(abs(delta))<1e-10)
d<-data.frame(eid=as.character(g[,'eid']),od[,bridge$score_fields],check.names=FALSE)
# Reproduce the actual original export serialization, not an invented higher-
# precision representation. write.table() rounds doubles when forming the CSV;
# this affects a few exact rounding-boundary fingerprints at decimal place10.
# All CSV text stays in local memory; no records or scores are uploaded.
tc<-textConnection('csv_text','w',local=TRUE);write.table(d,tc,sep=',',row.names=FALSE,na='NA');close(tc)
rt<-as.data.frame(fread(paste(csv_text,collapse='\n'),colClasses=c(eid='character')))
stopifnot(identical(as.character(d$eid),as.character(rt$eid)),identical(names(d),names(rt)))
serialization_error<-max(abs(as.matrix(d[,-1])-as.matrix(rt[,-1])));stopifnot(serialization_error<1e-12)
rt<-rt[order(rt$eid),];z<-as.matrix(rt[,-1]);txt<-vapply(seq_len(nrow(rt)),function(i)paste(c(rt$eid[i],sprintf('%.10f',z[i,])),collapse='|'),'');h<-digest(paste(txt,collapse='\n'),algo='sha256',serialize=FALSE);stopifnot(h==bridge$core_sha256,nrow(d)==bridge$N)
write_json(list(status='PASS',N=nrow(d),organs=18,individual_score_pairs_checked=nrow(d)*18,local_raw_to_original_ODIS_id_matched_max_error=max(abs(delta)),original_export_write_table_csv_roundtrip_max_error=serialization_error,server_to_reproduced_original_export_all_id_score_tuples_round10_equal=TRUE,canonical_core_sha256=h,exact_bitwise_identity_NOT_claimed=TRUE,server_comparison_bound='Every same eid-score tuple agrees after10decimal rounding with the actual original write.table CSV export; absolute tolerance<1e-10. Local raw versus saved original source directly measured separately.',server_fingerprint_file=file.path(m,'results/server_core_fingerprint.json'),private_score_upload_performed=FALSE),file.path(m,'results/id_score_bridge_acceptance.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
