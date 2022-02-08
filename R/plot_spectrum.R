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
    y<-rescale(y)
  }
  if (plot_spectrum == T){
    matplot(x=new.lambdas, y=y, type='l',
            #main=paste(peak, '\n', names(chrom_list)[chr], ';','\n', 'RT = ', RT, 'mins','; ',chr),
            ylab = 'Intensity', xlab = 'Wavelength (nm)',
            ylim=c(0,max(y)*1.2), ...)
    if (spectrum_labels == T){
      pks <- find_peaks(y,slope_thresh=.00001, bounds=F)
      if (length(pks)>0){
      pks <- data.frame(names(y)[pks], y[pks],stringsAsFactors = F)
      text(pks[,1],pks[,2],pks[,1],pos=3,offset=.3,cex = .8)
      }
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
    return(data.frame(y))}
}

## Function to plot all spectra of chosen peaks in peak table.

plot_all_spectra <- function(peak, peak_table, chrom_list, plot_spectrum = T,
                             export_spectrum=T, scale_spectrum=T, overlapping=T, verbose=F, ...){
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  sp <- sapply(1:length(chrom_list), function(chr){
    plot_spectrum(peak=peak, peak_table=peak_table, chrom_list=chrom_list, chr=chr,
                  plot_spectrum=F, plot_trace=F, export_spectrum = T,
                  scale_spectrum=scale_spectrum, verbose=verbose)
  })
  sp<-do.call(cbind, sp)
  if(overlapping==T){
    matplot(new.lambdas, sp, type='l', xlab='wavelength', ylab='intensity',las=2)
  } else{
    apply(sp, 2,function(spp){
      plot(new.lambdas, spp, type='l', xlab='', ylab='',las=2)
    })
  }
  if(export_spectrum==T){
    sp
  }
}