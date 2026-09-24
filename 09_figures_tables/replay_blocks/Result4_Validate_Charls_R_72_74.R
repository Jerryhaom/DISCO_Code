en = entropy_distance_parallel_fast(biomarker1[,bb],
             biomarker1[which(info$age<=50),bb])
pn=data.frame(hd=hd,en=en)
