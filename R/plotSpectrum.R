elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

pkTab = pkTab.PAs
CList = rts.warp6
new.ts <- ts
peak='X27'
plotSpectrum <- function(peak, pkTab, CList, chr = 'max', lambda = 'max',
                          ts = new.ts, lambdas = new.lambdas, plot_spectrum = T, 
                          plot_chrom = T, export_spectrum=FALSE,
                          spectrum_labels=T, ...){
  RT <- round(pkTab['RT',peak],2)
  pkTab<-pkTab[4:(nrow(pkTab)),]
  t <- which(elementwise.all.equal(RT,ts))
  if (chr == 'max'){
    chr <- which.max(pkTab[,peak])
  }
  y=CList[[chr]][t,]
  if (lambda == 'max'){
    lambda = names(which.max(y))
  }
  if (plot_spectrum == T){
    matplot(x=lambdas, y=y, type='l',
            #main=paste(names(CList)[chr], ';    ', 'RT = ', RT),
            #ylab = 'Intensity', xlab = 'Wavelength (nm)',
            ylim=c(0,max(y)*1.2), ...)
    if (spectrum_labels == T){
      pks <- alsace::findpeaks(y,span=3)
      pks <- data.frame(names(y)[pks], y[pks],stringsAsFactors = F)
      text(pks[,1],pks[,2],pks[,1],pos=3,offset=.3,cex = .8)
    }
  }
  if (plot_chrom == T){
    matplot(x=ts, y=CList[[chr]][,lambda],type='l',
            main=bquote(.(names(CList)[chr])*';' ~~~~ lambda ~ '=' ~ .(lambda) ~ 'nm'),
            ylab='', xlab='')
    abline(v=RT,col='red',lty=3)
  } 
  if (export_spectrum==TRUE){
    return(data.frame(lambdas,y))}
}
