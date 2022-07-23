#' Get peak list.
#' 
#' Finds and fits peaks and extracts peak parameters from a list of chromatograms
#' at the specified wavelengths.
#' 
#' Peaks are located by finding zero-crossings in the smoothed first derivative
#' of the specified chromatographic traces (function \code{\link{find_peaks}}).
#' At the given positions, an exponential-gaussian hybrid (or regular gaussian)
#' function is fit to the signal using \code{\link{fit_peaks}}). The area is then
#' calculated using a trapezoidal approximation.
#' 
#' The \code{sd}, \code{FWHM}, \code{tau}, and \code{area} are returned in units
#' determined by \code{time.units}. By defaults the units are in minutes.
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints x wavelengths).
#' @param lambdas Character vector of wavelengths to find peaks at.
#' @param fit What type of fit to use. Current options are exponential-gaussian
#' hybrid (\code{egh}), gaussian or raw. The \code{raw} setting performs
#' trapezoidal integration directly on the raw data without fitting a peak shape.
#' @param sd.max Maximum width (standard deviation) for peaks. Defaults to 50.
#' @param max.iter Maximum number of iterations for non-linear least squares
#' in \code{\link{fit_peaks}}.
#' @param time.units Units of \code{sd}, \code{FWHM}, \code{area}, and \code{tau}
#' (if applicable). Options are minutes \code{"min"}, seconds (\code{"s"}, or 
#' milliseconds \code{"ms"}.
#' @param \dots Additional arguments to \code{\link{find_peaks}}.
#' @return The result is an S3 object of class \code{peak_list}, containing a nested
#' list of data.frames containing information about the peaks fitted for each
#' chromatogram at each specified wavelength. The data.frame includes information
#' about the retention time (\code{rt}), \code{start} and \code{end} of each peak,
#' as well as the standard deviation (\code{sd}), \code{tau} (if \code{egh} is 
#' selected), full width at half maximum (\code{FWHM}), \code{height}, \code{area},
#' and \code{r.squared} (coefficient of determination). (*Note:* This last
#' parameter is determined from a linear model of the fitted peak values to the
#' raw data. This approach is not really statistically valid but it can be useful
#' as a rough metric for "goodness-of-fit").
#' @author Ethan Bass
#' @note The function is adapted from the
#' \href{https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R}{getAllPeaks}
#' function authored by Ron Wehrens (though the underlying algorithms for peak
#' identification and peak-fitting are not the same).
#' @references 
#' Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}
#' 
#' #' Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. \emph{Journal of
#' Chromatography A} \bold{915}:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions - A
#' good model for chromatographic peaks in isocratic HPLC? \emph{Chromatographia},
#' /bold{26}: 285-296. \doi{10.1007/BF02268168}.
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'), sd.max=50, fit="egh")
#' @seealso \code{\link{find_peaks}}, \code{\link{fit_peaks}}
#' @export get_peaks

get_peaks <- function(chrom_list, lambdas, fit = c("egh", "gaussian", "raw"),
                       sd.max=50, max.iter=100, time.units = c("min", "s", "ms"), ...){
  time.units <- match.arg(time.units, c("min", "s", "ms"))
  tfac <- switch(time.units, "min" = 1, "s" = 60, "ms" = 60*1000)
  fit <- match.arg(fit, c("egh", "gaussian", "raw"))
  chrom_list_string <- deparse(substitute(chrom_list))
  if (class(chrom_list)[1] == "matrix")
    chrom_list <- list(chrom_list)
  if (missing(lambdas)){
    if (ncol(chrom_list[[1]]) == 1){
      lambdas <- colnames(chrom_list[[1]])
    } else stop("Wavelengths (`lambdas`) must be provided.")
  }
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  if (is.null(names(chrom_list))){
    warning("Sample names not found. It is recommended to include names for your samples.", immediate. = TRUE)
    names(chrom_list) <- seq_along(chrom_list)
  }
  peaks<-list()
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=FALSE])
  result <- lapply(seq_along(chrom_list), function(sample){
    suppressWarnings(ptable <- lapply(lambdas, function(lambda){
      cbind(sample = names(chrom_list)[sample], lambda,
            fit_peaks(chrom_list[[sample]][,lambda], fit = fit,
                      max.iter = max.iter, sd.max = sd.max, ...))
    }))
    names(ptable) <- lambdas
    ptable
  })
  names(result) <- names(chrom_list)
  result <- lapply(result, function(sample) lapply(sample, function(pks){
    pks[apply(pks, 1, function(x){
      !any(is.na(x)) & (x["rt"] > x["start"]) & x["rt"] < x["end"]}), , drop=FALSE]
  }))
  timepoints <- as.numeric(rownames(chrom_list[[1]]))
  tdiff <- median(diff(timepoints))
  result <- lapply(result, function(smpl) lapply(smpl, function(lambda){
    x <- lambda
    x[, c('rt', 'start', 'end')] <- sapply(c('rt', 'start', 'end'),
                                           function(j) timepoints[x[,j]])
    x[, c('sd', 'FWHM', 'area')] <- x[, c('sd', 'FWHM', 'area')] * tdiff * tfac
    if (!is.null(x$tau)){x[, c('tau')] <- x[, c('tau')] * tdiff * tfac} 
    x
  }))
  structure(result,
            chrom_list = chrom_list_string,
            lambdas = lambdas, fit=fit, sd.max=sd.max,
            max.iter = max.iter,
            time.units = time.units,
            class = "peak_list")
}

#' Plot fitted peak shapes.
#' 
#' Visually assess integration accuracy by plotting fitted peaks over trace.
#'
#' @importFrom stats median
#' @importFrom graphics polygon arrows
#' @importFrom scales alpha
#' @param x Peak_list object. Output from the \code{get_peaks} function.
#' @param chrom_list List of chromatograms (retention time x wavelength
#' matrices)
#' @param index Index or name of chromatogram to be plotted.
#' @param lambda Wavelength for plotting.
#' @param points Logical. If TRUE, plot peak maxima. Defaults to FALSE.
#' @param ticks Logical. If TRUE, mark beginning and end of each peak. Defaults
#' to FALSE.
#' @param a Alpha parameter controlling the transparency of fitted shapes.
#' @param color The color of the fitted shapes.
#' @param cex.points Size of points. Defaults to 0.5
#' @param \dots Additional arguments to plot function.
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots a chromatographic trace from the specified chromatogram (\code{chr})
#' at the specified wavelength (\code{lambda}) with fitted peak shapes from the
#' provided \code{peak_list} drawn underneath the curve. 
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @rdname plot.peak_list
#' @export
#' 
plot.peak_list <- function(x, ..., chrom_list=NULL, index=1, lambda=NULL,
                       points=FALSE, ticks=FALSE, a=0.5, color=NULL,
                       cex.points=0.5){
  time.units <- attributes(x)$time.units
  tfac <- switch(time.units, "min" = 1, "s" = 1/60, "ms" = 1/60000)
  if (is.null(chrom_list)){
    chrom_list <- get_chrom_list(x)
  }
  if (is.null(lambda)){
    lambda <- names(x[[1]])[1]
  }
  if (!(lambda %in% names(x[[1]]))){
    stop('Error: lambda must match one of the wavelengths in your peak list')
  }
  if (is.numeric(lambda)){
    lambda <- as.character(lambda)
  }
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  y <- chrom_list[[index]][,lambda]
  pks <- data.frame(x[[index]][[lambda]])
  if ("r.squared" %in% colnames(pks)){
    fit <- ifelse("tau" %in% colnames(pks), "egh", "gaussian")
  } else{
    fit <- "raw"
  }
  plot(new.ts, y, type='l', xlab='', ylab='', xaxt='n', yaxt='n', ...)
  if (points){
    points(pks$rt, pks$height, pch=20, cex=cex.points, col='red')
  }
  if (ticks){
    arrows(pks$start, y[which(new.ts %in% pks$start)]-5,
           pks$start, y[which(new.ts %in% pks$start)]+5,
           col="blue", length=0)
    arrows(pks$end, y[which(new.ts %in% pks$end)]-5,
           pks$end,y[which(new.ts %in% pks$end)]+5,
           col="blue", length=0)
  }
  res <- median(diff(as.numeric(rownames(chrom_list[[1]]))))
  for (i in seq_len(nrow(pks))){
    peak.loc<-seq.int((pks$start[i]),(pks$end[i]), by = res)
      if (fit == "gaussian"){
        yvals <- gaussian(peak.loc, center=pks$rt[i],
                          width=pks$sd[i]*tfac, height = pks$height[i])
        if (is.null(color))
          color <- "red"
      }
      else if (fit == "egh"){
        yvals <- egh(x = peak.loc, center = pks$rt[i],
                     width=pks$sd[i]*tfac, height = pks$height[i],
                     tau=pks$tau[i]*tfac)
        if (is.null(color))
          color <- "purple"
      }
      else if (fit == "raw"){
        yvals <- chrom_list[[index]][as.character(peak.loc), lambda]
        if (is.null(color))
          color <- "hotpink"
      }
      sapply(1:(length(peak.loc) - 1), function(i){
        polygon(peak.loc[c(i, i, (i+1), (i+1))], c(0, yvals[i:(i+1)], 0),
                col=alpha(color, a), lty=3, border = NA)
      })
  }
}
