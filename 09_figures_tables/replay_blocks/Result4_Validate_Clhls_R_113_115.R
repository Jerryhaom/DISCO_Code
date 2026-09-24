en = entropy_distance_parallel_fast(biomarker1[,bb],
                                    biomarker1[which(biomarker$trueage<=70),bb])
pn=data.frame(hd=hd,en=en)
