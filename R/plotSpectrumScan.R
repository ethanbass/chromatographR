## Plot spectrum at a particular time point.

plot_spectrum_scan <- function(rt, chrom_list, chr = 'max', lambda = 'max',
                               ts = new.ts, lambdas = new.lambdas, plot_spectrum = T, 
                               plot_chrom = T,export_spectrum=FALSE,
                               spectrum_labels=T, ...){
  RT <- round(as.numeric(rt),2)
  t <- which(elementwise.all.equal(RT,new.ts))
  if (chr == 'max'){
    chr <- which.max(sapply(chrom_list,function(x) max(x[RT,])))
  }
  y=chrom_list[[chr]][t,]
  if (lambda == 'max'){
    lambda = names(which.max(y))
  }
  if (plot_spectrum == T){
    matplot(x=lambdas, y=y, type='l',
            #main=paste(peak, '\n', names(chrom_list)[chr], ';','\n', 'RT = ', RT, 'mins','; ',chr),
            #ylab = 'Intensity', xlab = 'Wavelength (nm)',
            ylim=c(0,max(y)*1.2), ...)
    if (spectrum_labels == T){
      pks <- alsace::findpeaks(y,span=3)
      pks <- data.frame(names(y)[pks], y[pks],stringsAsFactors = F)
      text(pks[,1],pks[,2],pks[,1],pos=3,offset=.3,cex = .8)
    }
  }
  if (plot_chrom == T){
    matplot(x=ts, y=chrom_list[[chr]][,lambda],type='l',
            #main=paste(names(chrom_list)[chr], ';','\n', 'RT = ', RT,
            #          '; Wavelength = ', lambda, 'nm'))
            ylab='', xlab='')
    abline(v=RT,col='red',lty=3)
  } 
  if (export_spectrum==TRUE){
    return(data.frame(lambdas,y))}
}