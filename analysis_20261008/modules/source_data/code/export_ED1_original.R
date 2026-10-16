#!/usr/bin/env Rscript
# Replay the actual original producer, exposing local numeric objects only.
source(file.path(Sys.getenv('DISCO_RUN_MODULE'),'code/export_more_sources.R'))
run('ED1_original', function() {
  src <- Sys.getenv('DISCO_SIMULATION_SOURCE')
  lines <- readLines(src, warn = FALSE)
  e <- ctx()
  for (ex in parse(text = lines[3:101])) eval(ex, e)
  # The original global mean is established at lines 346–347.
  e$mean_vector <- rep(0, 10)
  for (ex in parse(text = lines[681:788])) eval(ex, e)
  stopifnot(is.function(e$plot_KL3))
  original_body <- body(e$plot_KL3)
  statements <- as.list(original_body)
  e$.capture <- function(d, p) { e$.captured <- list(data = d, plots = p) }
  # Insert a read-only observer immediately before the unchanged return value.
  body(e$plot_KL3) <- as.call(append(statements,
    list(quote(.capture(grid_distances, list(p1,p2,p3,p11,p21,p31)))),
    after = length(statements)-1L))
  for (sp in list(list(rho=.1, box='a', scatter='c'),
                  list(rho=.7, box='b', scatter='d'))) {
    invisible(e$plot_KL3(1, sp$rho))
    z <- e$.captured
    stopifnot(nrow(z$data)==11000, all(is.finite(as.matrix(z$data))))
    def <- paste('Exact original set.seed(123), n=1000,p=10, toeplitz rho=',sp$rho,
      '; 11 repeated noise levels on the same synthetic subjects. Entropy is the original full logdet difference (2*Gaussian signed DeltaH), not absolute entropy. Observer changes no calculation or RNG. No biological records.',sep='')
    emit(paste0('ED1_',sp$box),'original_noise_score_rows',z$data,src,def,
      selector='local grid_distances from actual plot_KL3; captured before unchanged return')
    emit(paste0('ED1_',sp$scatter),'original_metric_scatter_rows',z$data,src,def,
      selector='same local grid_distances used by p11,p21,p31')
    for(j in 1:3) plotlayers(paste0('ED1_',sp$box),paste0('original_box_',j),z$plots[[j]],src,def)
    emit(paste0('ED1_',sp$scatter),'Spearman_pair_annotations',data.frame(
      pair=c('Mahalanobis_DISCO','Entropy_DISCO','Mahalanobis_Entropy'),
      rho=c(cor(z$data$Mahalanobis,z$data$DISCO,method='spearman'),
            cor(z$data$Entropy,z$data$DISCO,method='spearman'),
            cor(z$data$Mahalanobis,z$data$Entropy,method='spearman'))),src,def)
  }
})
write_json(entries,file.path(out,'results/r_exports.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
write_json(errors,file.path(out,'results/r_export_errors.json'),pretty=TRUE,auto_unbox=TRUE)
cat('FINAL components',length(entries),'ERRORS',length(errors),'\n')

stopifnot(length(errors)==0L,length(entries)>0L,all(vapply(entries,function(x)x$privacy=='aggregate_public_candidate',TRUE)))
