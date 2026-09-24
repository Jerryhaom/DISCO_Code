bb=colnames(biomarker1)[2:ncol(biomarker1)]
en = entropy_distance_parallel_fast(biomarker1[,bb],
                                    biomarker1[which(phe$age2016<=75),bb])
pn=data.frame(hd=hd,en=en)
