#' Plot spectrum from peak table
#' 
#' A function to plot the trace and/or spectrum for a given peak in peak table.
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text
#' @param loc The name of the peak or retention time for which you wish to
#' extract spectral data.
#' @param peak_table The peak table (output from \code{\link{getPeakTable}}
#' function).
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param chr Numerical index of chromatogram you wish to plot, or "max"" if
#' you want to automatically plot the chromatogram with the largest signal.
#' @param lambda The wavelength you wish to plot the trace at if plot_chrom ==
#' T and/or the wavelength to be used for the determination of signal
#' abundance.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param plot_trace Logical. If TRUE, plots the trace of the chosen peak at
#' lambda. Defaults to TRUE.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param verbose Logical. If TRUE, prints verbose output to console. Defaults
#' to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param what What to look for. Either "peak" to extract spectral information
#' for a certain peak, "rt" to scan by retention time, or "click" to manually
#' select retention time by clicking on the chromatogram. Defaults to "peak"
#' mode.
#' @param ... Additional arguments.
#' @author Ethan Bass
#' @export plot_spectrum
plot_spectrum <- function(loc, peak_table=NULL, chrom_list, chr = 'max', lambda = 'max',
                          plot_spectrum = T, plot_trace = T, export_spectrum=FALSE,
                          spectrum_labels=T, verbose=T, scale_spectrum=F, what=c("peak", "rt", "click"), ...){
  what <- match.arg(what, c("peak", "rt", "click"))
  if (what=="peak"){
    if(is.null(peak_table)){
      stop("Peak table is required to locate peak spectrum.")}
    if (!(loc %in% colnames(peak_table))){
      stop(paste0("No match found for peak \'", loc, "\' in peak table."))}
  }
  if (what=="rt" & chr=="max"){
    stop("Must specify chromatogram for scan function.")}
  if (what=="click" & (!is.numeric(chr) | !(lambda %in% new.lambdas))){
      stop("Chromatogram ('chr') and wavelength ('lambda') must be specified for manual selection.")}
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  sig <- max(sapply(strsplit(rownames(chrom_list[[1]]),".",fixed=T),function(x) nchar(x[2])),na.rm=T)
  if (what=="click"){
    y<-chrom_list[[chr]][,lambda]
    matplot(x=new.ts, y=y, type='l', ylab='', xlab='')
    print("click trace to select timepoint")
    t <- identify(new.ts, y,n=1, plot=F)
    RT <- new.ts[t]
    abline(v=RT,col='red',lty=3)
    title(paste0("\n\n chr ", chr,  " ;   rt: ", RT, " ;  abs: ", round(y[t],1)))
    y=chrom_list[[chr]][t,]
    # closest match
    if (verbose){
      print(paste0("chrome no. ", chr, "; RT: ", RT, "; lambda = ", lambda, " nm"))
    if (!is.null(peak_table)){
    pk <- names(which.min(abs(peak_table["RT",]-RT)))
    print(paste("closest match in peak table is",pk))
    }}
  } else {
    if (what=="peak"){
    RT <- round(peak_table['RT',loc], sig)
    peak_table<-peak_table[4:(nrow(peak_table)-3),]
    } else if (what=="rt"){
      RT <- round(as.numeric(loc), sig)
      }
    t <- which(elementwise.all.equal(RT,new.ts))
    if (chr == 'max'){
      chr <- which.max(peak_table[,loc])
    }
    y=chrom_list[[chr]][t,]
    if (lambda == 'max'){
      lambda = names(which.max(y))
    } else lambda <- as.character(lambda)
    if (plot_trace){
      matplot(x=new.ts, y=chrom_list[[chr]][,lambda],type='l',
              #main=paste(names(chrom_list)[chr], ';','\n', 'RT = ', RT,
              #          '; Wavelength = ', lambda, 'nm'))
              ylab='', xlab='')
      abline(v=RT,col='red',lty=3)
      if (verbose==T){
        print(paste0("chrome no. ", chr, "; RT: ", RT, "; lambda = ", lambda, " nm"))
      }
    }
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
  if (export_spectrum){
    return(data.frame(y))}
}


## Elementwise all equal function from Brian Diggs
## (https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal)
#' @noRd
elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

#' Plot all spectra for chosen peak
#' 
#' Plot multiple for a given peak in peak table. Wrapper for
#' \code{\link{plot_spectrum}}.
#' 
#' @param peak The name of the peak you wish to plot (must be in character
#' format)
#' @param peak_table The peak table (output from \code{\link{getPeakTable}}
#' function)
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints x components).
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' @param overlapping Logical. If TRUE, plot spectra in single plot.
#' @param verbose Logical. If TRUE, print verbose output to console.
#' @param \dots Additional arguments to plot_spectrum.
#' @author Ethan Bass
#' @seealso \code{\link{plot_spectrum}}
#' @export plot_all_spectra
plot_all_spectra <- function(peak, peak_table, chrom_list, plot_spectrum = T,
                             export_spectrum=T, scale_spectrum=T, overlapping=T, verbose=F, ...){
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  sp <- sapply(1:length(chrom_list), function(chr){
    plot_spectrum(loc=peak, peak_table=peak_table, chrom_list=chrom_list, chr=chr,
                  plot_spectrum=F, plot_trace=F, export_spectrum = T,
                  scale_spectrum=scale_spectrum, verbose=verbose, what="peak")
  })
  sp<-do.call(cbind, sp)
  colnames(sp) <- names(chrom_list)
  if (plot_spectrum==T){
    if(overlapping==T){
      matplot(new.lambdas, sp, type='l', xlab='wavelength', ylab='intensity',las=2)
    } else{
      apply(sp, 2,function(spp){
        plot(new.lambdas, spp, type='l', xlab='', ylab='',las=2)
      })
  }
  }
  if(export_spectrum==T){
    sp
  }
}
