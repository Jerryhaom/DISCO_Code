or2=readRDS("Tissue.order.rds")

overlapped3=overlapped2[match(or2,rownames(overlapped2)),
                        match(or2,rownames(overlapped2))]
