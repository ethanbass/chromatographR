#' Extract all peaks from the chromatographic profiles of a list of
#' chromatograms.
#' 
#' Function to extract all peaks in a list of chromatograms. Peaks are located
#' as local maxima within the given span (function \code{\link{find_peaks}})
#' and at the given positions a gaussian curve is fit (function
#' \code{\link{fit_peaks}}).
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints times wavelengths).
#' @param lambdas Character vector of wavelengths to find peaks at.
#' @param fit What type of fit to use. Current options are gaussian and
#' exponential gaussian hybrid.
#' @param sd.max Maximum width (standard deviation) for peaks. Defaults to 50.
#' @param max.iter Maximum number of interations for non-linear least squares
#' in \code{\link{fit_peaks}}.
#' @param \dots Additional arguments to \code{\link{find_peaks}}.
#' @return The result is a list, with each element corresponding to one data
#' file, and containing data for the fitted peaks for each of the ALS
#' components. Note that this function presents the "rt", "sd" and "FWHM"
#' fields in real time units.
#' @note Function is adapted from the
#' \url{https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R}{getAllPeaks}
#' function authored by Ron Wehrens.
#' @author Ethan Bass & Ron Wehrens
#' @export get_peaks

get_peaks <- function (chrom_list, lambdas, fit = c("egh", "gaussian"),
                       sd.max=50, max.iter=100, ...){
  fit <- match.arg(fit, c("egh", "gaussian"))
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  chrom_list_str <- deparse(substitute(chrom_list))
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=F])
  peak_positions <- lapply(chrom_list, function(c_mat){
    apply(c_mat, 2, function(x) find_peaks(x, ...))})
  result <- lapply(seq_along(chrom_list), function(smpl){
    ptable <- lapply(seq_along(peak_positions[[smpl]]), function(cmpnd){
      fit_peaks(chrom_list[[smpl]][,cmpnd], peak_positions[[smpl]][[cmpnd]],
                fit = fit, max.iter = max.iter, sd.max = sd.max)
    })
    names(ptable) <- names(peak_positions[[smpl]])
    ptable
  })
  names(result) <- names(peak_positions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks){
    pks[apply(pks, 1, function(x) !any(is.na(x))), , drop = FALSE]
    }))
  timepoints <- as.numeric(rownames(chrom_list[[1]]))
  tdiff <- median(diff(timepoints))
  result <- lapply(result, function(smpl) lapply(smpl, function(cmpnd){
    x <- cmpnd
    x[, c('rt', 'start', 'end')] <- sapply(c('rt', 'start', 'end'),
                                           function(j) timepoints[x[,j]])
    x[, c('sd', 'FWHM')] <- x[, c('sd', 'FWHM')] * tdiff
    if (!is.null(x$tau)){x[, c('tau')] <- x[, c('tau')] * tdiff} 
    x
  }))
  structure(result,
            chrom_list = chrom_list_str,
            lambdas = deparse(substitute(lambdas)), fit=fit, sd.max=sd.max,
            max.iter=max.iter,
            class="peak_list")
}

#' Function to visually assess accuracy of integration
#' 
#' Visually assess integration accuracy by fitted peaks onto chromatogram.
#'
#' @importFrom stats median
#' @importFrom graphics polygon arrows
#' @param x Peak_list object. Output from the \code{get_peaks} function.
#' @param chrom_list List of chromatograms (retention time x wavelength
#' matrices)
#' @param index Index or name of chromatogram to be plotted.
#' @param lambda Wavelength for plotting.
#' @param points Logical. If TRUE, plot peak maxima. Defaults to FALSE.
#' @param ticks Logical. If TRUE, mark beginning and end of each peak. Defaults
#' to FALSE.
#' @param a Select "alpha"" parameter controlling transparency of shapes.
#' @param cex.points Size of points. Defaults to 0.5
#' @param \dots Additional arguments to plot function.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @rdname plot.peak_list
#' @export
#' 
plot.peak_list <- function(x, ..., chrom_list=NULL, index=1, lambda=NULL,
                       points=FALSE, ticks=FALSE, a=0.5, cex.points=0.5){
  if (is.null(chrom_list)){
    chrom_list <- get(attr(x, "chrom_list"))
  }
  if (is.null(lambda)){
    lambda <- names(x[[1]])[1]
  }
  if (!(lambda %in% names(x[[1]]))){
    stop('Error: lambda must match one of the wavelengths in your peak list')
  }
  if (is.numeric(lambda)){lambda <- as.character(lambda)}
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  y <- chrom_list[[index]][,lambda]
  pks <- data.frame(x[[index]][[lambda]])
  fit <- ifelse("tau" %in% colnames(pks), "egh", "gaussian")
  plot(new.ts, y, type='l', xlab='', ylab='', xaxt='n', yaxt='n', ...)
  if (points){
    points(pks$rt, pks$height, pch=20, cex=cex.points, col='red')
  }
  if (ticks){
    arrows(pks$start, y[which(new.ts %in% pks$start)]-5, pks$start, y[which(new.ts %in% pks$start)]+5, col="blue", length=0)
    arrows(pks$end, y[which(new.ts %in% pks$end)]-5, pks$end, y[which(new.ts %in% pks$end)]+5, col="blue", length=0)
  }
  for (i in seq_len(nrow(pks))){
    peak.loc<-seq.int((pks$start[i]),(pks$end[i]), by = .01)
      if (fit == "gaussian"){
        yvals <- gaussian(peak.loc, center=pks$rt[i], width=pks$sd[i], height = pks$height[i])
        color <- "red"
      }
      else if (fit == "egh"){
        yvals <- egh(x=peak.loc, center=pks$rt[i], width=pks$sd[i], height = pks$height[i], tau=pks$tau[i])
        color <- "purple"
      }
      sapply(1:(length(peak.loc)-1), function(i){
        polygon(peak.loc[c(i,i,(i+1),(i+1))], c(0,yvals[i:(i+1)],0),col=scales::alpha(color,a), lty=3,border=NA)
      })
  }
}
