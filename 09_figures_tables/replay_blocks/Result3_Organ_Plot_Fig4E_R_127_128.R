odis$cause[odis$cause=="Diabetes mellitus"]="Others"
cause=na.omit(unique(unique(odis$cause)))[c(2,1,6,5,3,4)]
