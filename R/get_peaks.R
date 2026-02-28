#' Get peak list.
#' 
#' Finds and fits peaks and extracts peak parameters from a list of chromatograms
#' at the specified wavelengths.
#' 
#' Peaks are located by finding zero-crossings in the smoothed first derivative
#' of the specified chromatographic traces (function \code{\link{find_peaks}}).
#' At the given positions, a peak model is fit to the signal using 
#' \code{\link{fit_peaks}} according to the value of the \code{fit} argument.
#' Available models include the bidirectional exponentially modified gaussian
#' (\code{bemg}),exponential-gaussian hybrid (\code{egh} or regular 
#' \code{gaussian}.  Finally, the area is calculated using a trapezoidal 
#' approximation.
#'
#' Additional arguments can be provided to \code{\link{find_peaks}} to fine-tune
#' the peak-finding algorithm. For example, the \code{smooth_window} can be
#' increased to prevent peaks from being split into multiple features. Overly
#' aggressive smoothing may cause small peaks to be overlooked. 
#' 
#' The standard deviation (\code{sd}), full-width at half maximum (\code{FWHM}),
#' \eqn{\tau} (\code{tau}), \eqn{\tau_{R}} (\code{tau_right}), 
#' \eqn{\tau_{L}} (\code{tau_left}), and \code{area} are returned in units determined by 
#' \code{time_unit}. By default, the units are in minutes. To compare directly
#' with 'ChemStation' integration results, the time units should be changed to
#' seconds.
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints × wavelengths).
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' find peaks at. For one-dimensional chromatograms, this argument can be ignored.
#' @param fit What type of fit to use. Current options are bidirectional 
#' exponentially modified gaussian (\code{bemg}), exponential-gaussian
#' hybrid (\code{egh}), \code{gaussian}, or \code{raw}. The \code{raw} setting 
#' performs trapezoidal integration directly on the raw data without fitting a 
#' model to the data. Defaults to \code{bemg}.
#' @param sd_max Maximum width (standard deviation) for peaks. Defaults to 
#' \code{50}.
#' @param sd.max The \code{sd.max} argument is deprecated. Please use 
#' \code{sd_max} instead.
#' @param max_iter Maximum number of iterations for non-linear least squares
#' in \code{\link{fit_peaks}}.
#' @param max.iter The \code{max.iter} argument is deprecated. Please use 
#' \code{max_iter} instead.
#' @param time_unit Units of \code{sd}, \code{FWHM}, \code{area}, and \code{tau}
#' (if applicable). Options are minutes (\code{"min"}), seconds (\code{"s"}), or 
#' milliseconds (\code{"ms"}).
#' @param time.units The \code{time.units} argument is deprecated. Please use 
#' \code{time_unit} instead.
#' @param estimate_purity Logical. Whether to estimate purity or not. Defaults
#' to \code{FALSE}. (If \code{TRUE}, this will slow down the function significantly).
#' @param noise_threshold Noise threshold. Argument to \code{get_purity}.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' \code{TRUE} if \code{\link[pbapply]{pbapply}} is installed.
#' @param cl Argument to \code{\link[pbapply]{pblapply}} or \code{\link[parallel]{mclapply}}.
#' Either an integer specifying the number of clusters to use for parallel
#' processing or a cluster object created by \code{\link[parallel]{makeCluster}}.
#' Defaults to 2. On Windows integer values will be ignored.
#' @param collapse Logical. Whether to collapse multiple peak lists per sample
#' into a single list when multiple wavelengths (\code{lambdas}) are provided.
#' @param \dots Additional arguments to \code{\link{find_peaks}}. Arguments
#' provided to \code{\link{find_peaks}} can be used to fine-tune the peak-finding
#' algorithm. Most importantly, the \code{smooth_window} should be increased if
#' features are being split into multiple bins. Other arguments that can be used
#' here include \code{smooth_type}, \code{slope_thresh}, and \code{amp_thresh}.
#' @return The result is an S3 object of class \code{peak_list}, containing a 
#' nested list of data.frames containing information about the peaks fitted for 
#' each chromatogram at each of wavelengths specified by the \code{lamdas}
#' argument. Each row in these data.frames is a peak and the columns contain 
#' information about various peak parameters:
#' * `rt`: The retention time of the peak maximum.
#' * \code{start}: The retention time where the peak is estimated to begin.
#' * \code{end}: The retention time where the peak is estimated to end.
#' * \code{sd}: The standard deviation (\eqn{\sigma}) of the fitted peak shape.
#' * \code{tau} Empirical skewness parameter (in units of time) determining peak
#' asymmetry when \code{fit = "egh"}.
#' * \code{tau_right} Exponential rate constant controlling right-sided tailing
#' when \code{fit = "bemg"}. Note that unlike \code{tau} in the EGH model,
#' this parameter has units of 1/time.
#' * \code{tau_left} Exponential rate constant controlling left-sided tailing
#' when \code{fit = "bemg"}. Note that unlike \code{tau} in the EGH model,
#' this parameter has units of 1/time.
#' * \code{FWHM}: The full-width at half maximum calculated as \eqn{2.355 \sigma}.
#' * \code{height}: The height of the peak.
#' * \code{area}: The area of the peak as determined by trapezoidal approximation.
#' * \code{r.squared} The coefficient of determination (\eqn{R^2}) of the fitted
#' model to the raw data. (\strong{Note}: this value is calculated by fitting a
#' linear model of the fitted peak values to the raw data. This approach is
#' statistically questionable, since the models are fit using non-linear least
#' squares. Nevertheless, it can still be useful as a rough metric for 
#' "goodness-of-fit").
#' * \code{purity} The peak purity as estimated by \code{\link{get_purity}}.
#' @author Ethan Bass
#' @note The bones of this function are adapted from the
#' \href{https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R}{getAllPeaks}
#' function authored by Ron Wehrens (though the underlying algorithms for peak
#' identification and peak-fitting are not the same).
#' @references 
#' * Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. \emph{Journal of
#' Chromatography A} \strong{915}:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' * Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions 
#' — A good model for chromatographic peaks in isocratic HPLC? \emph{Chromatographia},
#' \strong{26}: 285-296. \doi{10.1007/BF02268168}.
#' 
#' * Arase, Shuntaro, Kanta Horie, Takashi Kato, et al. 2016. Intelligent Peak 
#' Deconvolution through In-Depth Study of the Data Matrix from Liquid 
#' Chromatography Coupled with a Photo-Diode Array Detector Applied to 
#' Pharmaceutical Analysis. *Journal of Chromatography. A*, **1469**: 35–47. 
#' \doi{doi.org/10.1016/j.chroma.2016.09.037}.
#'
#' * O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement.
#' \url{https://terpconnect.umd.edu/~toh/spectrum/} (Accessed January, 2022).
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'), sd_max = 50, fit = "egh")
#' @seealso \code{\link{find_peaks}}, \code{\link{fit_peaks}}
#' @export get_peaks
#' @md

get_peaks <- function(chrom_list, lambdas, 
                      fit = c("bemg", "egh", "gaussian", "raw"),
                      sd_max = 50, max_iter = 100,
                      time_unit = c("min", "s", "ms"),
                      estimate_purity = FALSE,  noise_threshold = .001,
                      show_progress = NULL, cl = 2, collapse = FALSE, 
                      time.units = NULL, sd.max = NULL, max.iter = NULL, ...){
  max_iter <- resolve_deprecated(max.iter, max_iter)
  sd_max <- resolve_deprecated(sd.max, sd_max)
  time_unit <- resolve_deprecated(time.units, time_unit)
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
  fit <- match.arg(fit, c("bemg", "egh", "gaussian", "raw"))
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
    warning("Sample names not found. It is recommended to include names for your samples.",
            immediate. = TRUE)
    names(chrom_list) <- seq_along(chrom_list)
  }
  peaks <- list()
  laplee <- choose_apply_fnc(show_progress, cl = cl)
  
  result <- laplee(seq_along(chrom_list), function(sample){
    suppressWarnings(ptable <- lapply(lambdas, function(lambda){
      pks <- fit_peaks(chrom_list[[sample]], lambda = lambda, fit = fit,
                       max_iter = max_iter, sd_max = sd_max,
                       estimate_purity = estimate_purity,
                       noise_threshold = noise_threshold, ...)
      pks <- cbind(sample = names(chrom_list)[sample], lambda, pks)
      pks <- remove_bad_peaks(pks)
      pks <- convert_indices_to_times(pks, chrom_list = chrom_list, 
                                      idx = sample, tfac = tfac)
      pks
    }))
    names(ptable) <- lambdas
    if (collapse){
      ptable <- do.call(rbind, ptable)
    }
    ptable
  })
  names(result) <- names(chrom_list)
  
  structure(result,
            chrom_list = chrom_list_string,
            lambdas = lambdas, fit = fit, sd_max = sd_max,
            max_iter = max_iter,
            time_unit = time_unit,
            intensity_unit = get_metadata_attribute(chrom_list[[1]], "detector_y_unit"),
            class = c("peak_list", "list"))
}

#' Remove bad peaks
#' This function is called internally by \code{get_peaks}.
#' @author Ethan Bass
#' @noRd
remove_bad_peaks <- function(pks){
  pks[which(apply(pks, 1, function(x)!all(is.na(x))) & 
              apply(pks[,c("rt","start","end")], 1, function(x)all(!is.na(x))) &
              pks[,"rt"] >= 1 
              # pks[,"rt"] > pks[,"start"] &
              # pks[,"rt"] < pks[,"end"]
            ), , drop = FALSE]
}

#' Convert indices to times
#' This function is called internally by \code{get_peaks}.
#' @author Ethan Bass
#' @noRd
convert_indices_to_times <- function(x, chrom_list, idx, tfac){
  timepoints <- get_times(chrom_list, idx = idx)
  tdiff <- get_time_resolution(chrom_list, idx = idx)
  x[, c('rt', 'start', 'end')] <- sapply(c('rt', 'start', 'end'),
                                         function(j) timepoints[x[, j]])
  x[, c('sd', 'FWHM', 'area')] <- x[, c('sd', 'FWHM', 'area')] * tdiff * tfac
  if (!is.null(x$tau)){
    x[, c('tau')] <- x[, c('tau')] * tdiff * tfac
  } 
  if (!is.null(x$tau_right)){
    x[, c('tau_right', 'tau_left')] <- x[, c('tau_right', 'tau_left')] / (tdiff * tfac)
  } 
  x
}
