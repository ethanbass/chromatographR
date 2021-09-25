## Function to plot all spectra of chosen peaks in peak table.

plot_all_spectra <- function(peak, peak_table, chrom_list, plot_spectrum = T, plot_chrom = F,
                             export_spectrum=T, add=F, scale_spectrum=T, ...){
  #par(mfrow=c(3,2))
  if (add == T){
    #par(mfrow=c(1,1))
    scale_spectrum <- T
  }
  sp <- sapply(1:length(chrom_list), function(chr){
    plot_spectrum(peak=peak, peak_table=peak_table, chrom_list=chrom_list, chr=chr,
                  plot_spectrum=plot_spectrum, plot_trace=plot_chrom, export_spectrum = export_spectrum,
                  scale_spectrum=scale_spectrum, ...)
  })
  if(plot_spectrum==T){
    do.call(cbind, sp)
  }
}
