## Function to plot spectra of peaks in peak table. Returns chromatogram
## with "highlighting" of selected peak and/or spectrum.

## Elementwise all equal function from Brian Diggs
## (https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal)

elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

plot_spectrum <- function(peak, peak_table, chrom_list, chr = 'max', lambda = 'max',
                          plot_spectrum = T, plot_trace = T, export_spectrum=FALSE,
                          spectrum_labels=T, verbose=T, scale_spectrum=F, ...){
  RT <- round(peak_table['RT',peak],2)
  peak_table<-peak_table[4:(nrow(peak_table)-3),]
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  t <- which(elementwise.all.equal(RT,new.ts))
  if (chr == 'max'){
    chr <- which.max(peak_table[,peak])
  }
  y=chrom_list[[chr]][t,]
  if (lambda == 'max'){
    lambda = names(which.max(y))
  } else lambda <- as.character(lambda)
  if (verbose==T){
    print(paste0("chrome no. ", chr, "; RT: ", RT, "; lambda = ", lambda, " nm"))
    }
  if (scale_spectrum == T){
    y<-scales::rescale(y)
  }
  if (plot_spectrum == T){
    matplot(x=new.lambdas, y=y, type='l',
            #main=paste(peak, '\n', names(chrom_list)[chr], ';','\n', 'RT = ', RT, 'mins','; ',chr),
            #ylab = 'Intensity', xlab = 'Wavelength (nm)',
            ylim=c(0,max(y)*1.2), ...)
    if (spectrum_labels == T){
      pks <- alsace::findpeaks(y)
      pks <- data.frame(names(y)[pks], y[pks],stringsAsFactors = F)
      text(pks[,1],pks[,2],pks[,1],pos=3,offset=.3,cex = .8)
    }
  }
  if (plot_trace == T){
    matplot(x=new.ts, y=chrom_list[[chr]][,lambda],type='l',
            #main=paste(names(chrom_list)[chr], ';','\n', 'RT = ', RT,
            #          '; Wavelength = ', lambda, 'nm'))
            ylab='', xlab='')
    abline(v=RT,col='red',lty=3)
  } 
  if (export_spectrum==TRUE){
    return(data.frame(new.lambdas,y))}
}
