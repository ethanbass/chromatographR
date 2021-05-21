## Function to plot all spectra of chosen peaks in peak table.

plot_all_spectra <- function(peak=i, peak_table = pkTab.lvs2.s, chrom_list = mz.lvs.warp,
                             lambdas=new.lambdas,plot_spectrum=T,plot_chrom=F,export_spectrum=T, add=F, scale_spectrum=T, ...){
  #par(mfrow=c(3,2))
  if (add == T){
    par(mfrow=c(1,1))
    scale_spectrum <- T
  }
  for (i in 1:length(chrom_list)){
    plot_spectrum(peak=peak, peak_table=peak_table, chrom_list=chrom_list, lambdas=lambdas, chr=i,
                  plot_spectrum=plot_spectrum, plot_chrom=plot_chrom, export_spectrum = export_spectrum, ...)
  }
}
