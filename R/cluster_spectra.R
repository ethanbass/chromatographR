setClass("cluster", representation(peaks = "character", pval = "numeric"))

cluster_spectra <- function(pkTab, chrom_list, peak_no = c(5,100),
                            alpha=0.95, nboot=1000, plot_dend=T, plot_spectra=T,
                            verbose=T, save=T, parallel=T){
  if (verbose==T) print('...collecting representative spectra')
  rep <- sapply(1:ncol(pkTab), function(j){
    sp <- plot_spectrum(peak=j, peak_table=pkTab, chrom_list = dat.pr,
                        scale_spectrum=T, plot_trace=F, export_spectrum = T, plot_spectrum=F,verbose=F)
  })
  rep <- data.frame(do.call(cbind,rep))
  names(rep)<-paste0('V',1:ncol(rep))
  d<-1-abs(cor(rep,method="pearson"))
  
  if (verbose==T) print('...clustering spectra')
  result <- pvclust::pvclust(rep, method.dist="cor", method.hclust="average", nboot=nboot, parallel=parallel)
  
  if (plot_dend==T){
  plot(result,labels=F, cex.pv=0.5, print.pv='au',print.num = F)
  pvclust::pvrect(result, alpha=alpha, max.only = F)
  }
  if (save==T) saveRDS(result, 'pvclust.RDS')
  p <- pvclust::pvpick(result, alpha=alpha, max.only=F)
  l <- sapply(p$clusters, length)
  sub <- p$clusters[which(l > mn & l < mx)]
  pval<-1-result$edges[p$edges[which(l > mn & l < mx)],'au']
  sub <- lapply(1:length(sub), function(i) new("cluster", peaks=sub[[i]], pval=pval[i]))
  pval=format(round(result$edges[p$edges[which(l > mn & l < mx)],'au'],2), nsmall=2)
  names(sub)<-paste0('c',1:length(sub))
  
  if (plot_spectra==T){
    if (verbose==T) print('...plotting clustered spectra')
    sapply(1:length(sub), function(i){ 
      matplot(new.lambdas,rep[,as.numeric(gsub('V','',sub[[i]]@peaks))],
              type='l', ylab='', yaxt='n', xlab=expression(lambda),
              main=paste0('cluster ', i, '; p = ', format(round(sub[[i]]@pval,2),nsmall=2))
                          )})
  }
  return(sub)
}

dend<-as.dendrogram(result)
result %>%
  as.dendrogram() %>%
  hang.dendrogram() %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
result %>% text()
dend %>% set("labels_cex",0.5) %>% plot()
result %>% as.dendrogram %>% plot(labels=F)
result %>% pvrect(alpha = 0.95)

# dend %>%
#   pvclust_show_signif(result, signif_type = 'au', max.only=F, hang=-1) %>%
#   plot()
# dend %>%
#   pvclust_show_signif(result, show_type = "lwd") %>%
#   plot()
# result %>% text()
# dend %>% plot()
# result %>% pvrect(alpha = 0.95, max.only=F)
# 
# dendextend::pvclust_show_signif(dend, result)