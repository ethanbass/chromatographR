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
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param chr Numerical index of chromatogram you wish to plot, or "max" to
#' automatically plot the chromatogram with the largest signal.
#' @param lambda The wavelength you wish to plot the trace at if plot_trace ==
#' TRUE and/or the wavelength to be used for the determination of signal
#' abundance.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param plot_trace Logical. If TRUE, plots the trace of the chosen peak at
#' lambda. Defaults to TRUE.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param verbose Logical. If TRUE, prints verbose output to console. Defaults
#' to TRUE.
#' @param what What to look for. Either "peak" to extract spectral information
#' for a certain peak, "rt" to scan by retention time, or "click" to manually
#' select retention time by clicking on the chromatogram. Defaults to "peak"
#' mode.
#' @param ... Additional arguments.
#' @author Ethan Bass
#' @export plot_spectrum
plot_spectrum <- function(loc, peak_table, chrom_list=NULL,
                          chr = 'max', lambda = 'max',
                          plot_spectrum = TRUE, plot_trace = TRUE,
                          spectrum_labels=TRUE, scale_spectrum=FALSE,
                          export_spectrum = FALSE, verbose=TRUE, 
                          what=c("peak", "rt", "click"), ...){
  what <- match.arg(what, c("peak", "rt", "click"))
  if (what == "peak"){
    if(is.null(peak_table)){
      stop("Peak table is required to locate peak spectrum.")}
    if (!(loc %in% colnames(peak_table$tab))){
      stop(paste0("No match found for peak \'", loc, "\' in peak table."))}
    if (is.null(chrom_list)){
      chrom_list <- try(get(peak_table$args["chrom_list"]))
      if (class(chrom_list)=="try-error") stop("Chromatograms not found!")
    }
  }
  if (what == "rt" | what == "click"){
    if (chr == "max")
      stop("Chromatogram must be specified for scan function.")
    if (!(lambda %in% new.lambdas))
      stop("Lambda must be specified for scan function")
    if (is.null(chrom_list)){
      stop("List of chromatograms must be provided for scan function.")
    }
  }
  tab <- peak_table$tab
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  sig <- max(sapply(strsplit(rownames(chrom_list[[1]]),".",fixed=T),function(x) nchar(x[2])),na.rm=T)
  if (what == "click"){
    y<-chrom_list[[chr]][,lambda]
    matplot(x=new.ts, y=y, type='l', ylab='', xlab='')
    print("Click trace to select timepoint")
    time <- identify(new.ts, y,n=1, plot=F)
    RT <- new.ts[time]
    abline(v=RT,col='red',lty=3)
    title(paste0("\n\n chr ", chr,  " ;   rt: ", RT, " ;  abs: ", round(y[t],1)))
    y <- chrom_list[[chr]][time,]
    # closest match
    if (verbose){
      print(paste0("chrome no. ", chr, "; RT: ", RT, "; lambda = ", lambda, " nm"))
    if (!is.null(tab)){
    pk <- names(which.min(abs(tab["RT",] - RT)))
    print(paste("closest match in peak table is", pk))
    }}
  } else {
    if (what == "peak"){
    RT <- round(peak_table$pk_meta['RT',loc], sig)
    } else if (what == "rt"){
      RT <- round(as.numeric(loc), sig)
      }
    time <- which(elementwise.all.equal(RT, new.ts))
    if (chr == 'max'){
      chr <- which.max(tab[,loc])
    }
    y <- chrom_list[[chr]][time,]
    if (lambda == 'max'){
      lambda <- names(which.max(y))
    } else lambda <- as.character(lambda)
    if (plot_trace){
      matplot(x=new.ts, y=chrom_list[[chr]][,lambda],type='l',
              #main=paste(names(chrom_list)[chr], ';','\n', 'RT = ', RT,
              #          '; Wavelength = ', lambda, 'nm'))
              ylab='', xlab='')
      abline(v=RT,col='red',lty=3)
      if (verbose){
        print(paste0("chrome no. ", chr, "; RT: ", RT, "; lambda = ", lambda, " nm"))
      }
    }
  }
  if (scale_spectrum){
    y<-rescale(y)
  }
  if (plot_spectrum){
    matplot(x=new.lambdas, y=y, type='l',
            #main=paste(peak, '\n', names(chrom_list)[chr], ';','\n', 'RT = ', RT, 'mins','; ',chr),
            ylab = 'Intensity', xlab = 'Wavelength (nm)',
            ylim=c(0,max(y)*1.2), ...)
    if (spectrum_labels){
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
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function)
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints x components).
#' @param chrs Vector of chromatograms to plot.
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
plot_all_spectra <- function(peak, peak_table, chrom_list=NULL, chrs="all", 
                             plot_spectrum = T, export_spectrum=TRUE,
                             scale_spectrum=TRUE, overlapping=TRUE,
                             verbose=FALSE, ...){
  if (is.null(chrom_list)){
    chrom_list <- get(peak_table$args["chrom_list"])
  }
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  if ("all" %in% chrs) chrs <- seq_along(chrom_list)
  sp <- sapply(chrs, function(chr){
    plot_spectrum(loc = peak, peak_table = peak_table, chrom_list=chrom_list,
                  chr=chr, plot_spectrum=FALSE, plot_trace=FALSE, export_spectrum = TRUE,
                  scale_spectrum=scale_spectrum, verbose=verbose, what="peak")
  })
  sp<-do.call(cbind, sp)
  colnames(sp) <- names(chrom_list)[chrs]
  if (plot_spectrum){
    if(overlapping){
      matplot(new.lambdas, sp, type='l', xlab='wavelength', ylab='intensity',las=2)
    } else{
      apply(sp, 2,function(spp){
        plot(new.lambdas, spp, type='l', xlab='', ylab='',las=2)
      })
  }
  }
  if(export_spectrum){
    sp
  }
}

#' Plot spectrum from peak table
#' 
#' A function to plot the trace and/or spectrum for a given peak in peak table.
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text boxplot
#' @importFrom stats as.formula
#' @param x The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param loc The name of the peak or retention time for which you wish to
#' extract spectral data.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param chr Numerical index of chromatogram you wish to plot; "max" to
#' plot the chromatogram with the largest signal; or "all" to plot spectra
#' for all chromatograms.
#' @param lambda The wavelength you wish to plot the trace at if
#' plot_chrom ==T and/or the wavelength to be used for the determination
#' of signal abundance.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param plot_trace Logical. If TRUE, plots the trace of the chosen peak at
#' lambda. Defaults to TRUE.
#' @param box_plot Logical. If TRUE, plots box plot using categories
#' defined by \code{vars}.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param verbose Logical. If TRUE, prints verbose output to console. Defaults
#' to TRUE.
#' @param vars Independent variables for boxplot.
#' @param ... Additional arguments.
#' @author Ethan Bass
#' @rdname plot.peak_table
#' @export

plot.peak_table <- function(x, ..., loc=NULL, chrom_list=NULL,
                          chr = 'max', lambda = 'max',
                          plot_spectrum = TRUE, plot_trace = TRUE,
                          box_plot = FALSE,
                          spectrum_labels=TRUE, scale_spectrum=FALSE,
                          export_spectrum=FALSE, verbose=TRUE, vars=NULL){
  if (is.null(loc)){
    loc <- readline(prompt="Which peak would you like to plot? \n")
    loc <- gsub('\\"','',loc)
    loc <- gsub("\\'","",loc)
    if (!(loc %in% colnames(x$tab)))
      stop("Peak not found.")
  }
  if (plot_spectrum == TRUE | plot_trace == TRUE){
    if (chr == "all"){
      plot_all_spectra(loc, x, chrom_list=NULL,
                       plot_spectrum = plot_spectrum,
                       export_spectrum=export_spectrum,
                       verbose=verbose, ...)
    } else{
    plot_spectrum(loc, x, chrom_list, chr=chr,
                  lambda=lambda, plot_spectrum=plot_spectrum,
                  plot_trace=plot_trace, spectrum_labels=spectrum_labels,
                  scale_spectrum=scale_spectrum, export_spectrum=export_spectrum,
                  verbose=verbose, what="peak")
    }
  }
  if (box_plot == T){
    if (is.null(vars))
      stop("Must provide independent variable or variables for boxplot")
    boxplot(as.formula(paste("x$tab[,loc]",vars,sep="~")), data = x$sample_meta,
            main = paste(loc, '\n', 'RT = ', round(x$pk_meta['RT', loc],2)),
            ylab="abs", xlab="", ...)
  }
}