#' Get peak list.
#' 
#' Detects chromatographic peaks, fits peak models, and extracts peak parameters
#' at one or more wavelengths.
#' 
#' Peaks are detected by finding zero-crossings in the smoothed first derivative
#' of the specified chromatographic traces (using [`find_peaks`]).
#' Peak models are then fit to the detected features using [`fit_peaks`] 
#' according to the value of `fit`. Available models include the bidirectional
#' exponentially modified gaussian (`"bemg"`), exponential-gaussian hybrid 
#' (`"egh"`) or `"gaussian"`. Peak areas are then calculated using a 
#' trapezoidal approximation.
#' 
#' Additional arguments can be provided to `find_peaks` to fine-tune
#' the peak-finding algorithm. For example, the `smooth_window` can be
#' increased to reduce splitting of broad peaks into multiple features. Overly
#' aggressive smoothing may cause small peaks to be overlooked. 
#' 
#' The parameters `sd`, `FWHM`, `tau`, `tau_right`, `tau_left`, and
#' `area` are returned in units determined by `time_unit`, which defaults to
#' `minutes`. To compare directly with 'ChemStation' integration results, 
#' the time units should be changed to seconds.
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths).
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' find peaks at. For one-dimensional chromatograms, this argument can be ignored.
#' @param fit What type of fit to use. Current options are bidirectional 
#' exponentially modified gaussian (`"bemg"`), exponential-gaussian
#' hybrid (`"egh"`), `"gaussian"`, or `"raw"`. The `raw` setting 
#' performs trapezoidal integration directly on the raw data without fitting a 
#' model to the data. Defaults to `egh`.
#' @param sd_max Maximum width (standard deviation) for peaks. Defaults to `50`.
#' @param sd.max The `sd.max` argument is deprecated. Please use `sd_max` instead.
#' @param max_iter Maximum number of iterations for non-linear least squares
#' in [`fit_peaks`].
#' @param max.iter The `max.iter` argument is deprecated. Please use `max_iter` 
#' instead.
#' @param time_unit Specifies the units for `sd`, `FWHM`, `area`, and `tau` (if
#' applicable). Options are minutes (`"min"`), seconds (`"s"`), or 
#' milliseconds (`"ms"`). Defaults to `"min"`.
#' @param time.units The `time.units` argument is deprecated. Please use 
#' `time_unit` instead.
#' @param estimate_purity Logical. Whether to estimate purity or not. Defaults
#' to `FALSE`. (If `TRUE`, this will slow down the function significantly).
#' @param noise_threshold Noise threshold. Argument to `get_purity`.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed.
#' @param cl Either an integer specifying the number of cores to use for
#' parallel processing or a cluster object created by
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows integer
#' values will be ignored.
#' @param baseline Whether to fit a baseline offset beneath each peak model.
#' Options are `"none"` (no baseline, default), `"flat"` (constant offset per
#' peak), or `"sloped"` (linearly varying offset per peak, useful for gradient
#' elution). Not available when `fit = "raw"`.
#' @param collapse Logical. Whether to collapse multiple peak lists per sample
#' into a single list when multiple wavelengths (`lambdas`) are provided.
#' @param \dots Additional arguments to [`find_peaks`]. Arguments provided to
#' `find_peaks` can be used to fine-tune the peak-finding algorithm. Most
#' importantly, the `smooth_window` should be increased if features are being
#' split into multiple bins. Other arguments that can be used
#' here include `smooth_type`, `slope_thresh`, and `amp_thresh`.
#' @return The result is an S3 object of class `peak_list`, containing a
#' nested list of data.frames containing information about the peaks fitted for
#' each chromatogram at each of wavelengths specified by the `lambdas`
#' argument. Each row in these data.frames is a peak and the columns contain
#' information about various peak parameters:
#' * `rt`: The retention time of the peak maximum.
#' * `start`: The retention time where the peak is estimated to begin.
#' * `end`: The retention time where the peak is estimated to end.
#' * `sd`: The standard deviation (\eqn{\sigma}) of the fitted peak shape.
#' * `tau`: Empirical skewness parameter (in units of time) determining peak
#' asymmetry when `fit == "egh"`.
#' * `tau_right`: Exponential rate constant controlling right-sided tailing
#' when `fit = "bemg"`. Note that unlike `tau` in the EGH model,
#' this parameter has units of 1/time.
#' * `tau_left`: Exponential rate constant controlling left-sided tailing
#' when `fit = "bemg"`. Note that unlike `tau` in the EGH model,
#' this parameter has units of 1/time.
#' * `h`: The `height` parameter of the fitted BEMG shape function (only when
#' `fit = "bemg"`). This is the value of the pure shape at `center`, before
#' any baseline is added, and is distinct from `height`.
#' * `center`: The center parameter of the fitted BEMG shape function (only
#' when `fit = "bemg"`). This is the location around which the exponential
#' tails are applied. Because BEMG peaks can be asymmetric
#' (`tau_right ≠ tau_left`) and a sloped baseline shifts the fitted curve,
#' the actual peak maximum (`rt`) may differ from `center`.
#' * `FWHM`: The full-width at half maximum calculated as \eqn{2.355 \sigma}.
#' * `height`: The height of the peak above the fitted baseline.
#' * `area`: The area of the peak above the fitted baseline, estimated by
#' trapezoidal approximation.
#' * `r.squared`: The coefficient of determination (\eqn{R^2}) of the fitted
#' model to the raw data. (**Note**: this value is calculated by fitting a
#' linear model of the fitted peak values to the raw data. This approach is
#' statistically dubious, since the models are fit using non-linear least
#' squares. Nevertheless, it can still be useful as a rough metric for
#' "goodness-of-fit").
#' * `purity`: The peak purity as estimated by [`get_purity`].
#' * `floor`: Constant baseline level (only when `baseline = "flat"`).
#' * `floor_start`: Baseline level at the left edge of the peak window (only
#' when `baseline = "sloped"`).
#' * `floor_end`: Baseline level at the right edge of the peak window (only
#' when `baseline = "sloped"`). Both `floor_start` and `floor_end` are
#' constrained to be non-negative.
#' @author Ethan Bass
#' @note The bones of this function are adapted from the
#' `getAllPeaks` function (<https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R>)
#' by Ron Wehrens (though the underlying algorithms for finding and integrating
#' peaks are not the same.
#' @references 
#' * Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. *Journal of
#' Chromatography A* 915:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' * Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions 
#' — A good model for chromatographic peaks in isocratic HPLC? *Chromatographia*,
#' 26: 285-296. \doi{10.1007/BF02268168}.
#' 
#' * Arase, Shuntaro, Kanta Horie, Takashi Kato, et al. 2016. Intelligent Peak 
#' Deconvolution through In-Depth Study of the Data Matrix from Liquid 
#' Chromatography Coupled with a Photo-Diode Array Detector Applied to 
#' Pharmaceutical Analysis. *Journal of Chromatography. A*, 1469: 35–47. 
#' \doi{doi.org/10.1016/j.chroma.2016.09.037}.
#'
#' * O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement.
#' <https://terpconnect.umd.edu/~toh/spectrum/> (Accessed January, 2022).
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. 
#' *Metabolomics* 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'), sd_max = 50, fit = "egh")
#' @seealso [`find_peaks`], [`fit_peaks`], [`print.peak_list`], 
#' \code{\link{[.peak_list}}
#' @export get_peaks
#' @md

get_peaks <- function(chrom_list, lambdas,
                      fit = c("egh", "bemg", "gaussian", "raw"),
                      sd_max = 50, max_iter = 100,
                      time_unit = c("min", "s", "ms"),
                      baseline = c("none", "flat", "sloped"),
                      estimate_purity = FALSE,  noise_threshold = .001,
                      show_progress = NULL, cl = 2, collapse = FALSE,
                      time.units = NULL, sd.max = NULL, max.iter = NULL, ...){
  check_duplicated_names(names(chrom_list))
  max_iter <- resolve_deprecated(max.iter, max_iter)
  sd_max <- resolve_deprecated(sd.max, sd_max)
  time_unit <- resolve_deprecated(time.units, time_unit)
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  baseline <- match.arg(baseline, c("none", "flat", "sloped"))
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
                       noise_threshold = noise_threshold,
                       baseline = baseline, ...)
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
    ptable <- transfer_metadata(ptable, chrom_list[[sample]], 
                                exclude = c('names','row.names','dim',
                                              'dimnames', "format_out", 
                                              "data_format", "parser"))
    ptable
  })
  names(result) <- names(chrom_list)
  
  structure(result,
            chrom_list = chrom_list_string,
            lambdas = lambdas, fit = fit, sd_max = sd_max,
            max_iter = max_iter,
            time_unit = time_unit,
            baseline = baseline,
            intensity_unit = get_metadata_attribute(chrom_list[[1]], "detector_y_unit"),
            class = c("peak_list", "list"))
}

#' Remove bad peaks
#' This function is called internally by `get_peaks`.
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
#' This function is called internally by `get_peaks`.
#' @author Ethan Bass
#' @noRd
convert_indices_to_times <- function(x, chrom_list, idx, tfac){
  timepoints <- get_times(chrom_list, idx = idx)
  tdiff <- get_time_resolution(chrom_list, idx = idx)
  idx_cols <- intersect(c('rt', 'start', 'end', 'center'), colnames(x))
  x[, idx_cols] <- sapply(idx_cols, function(j) timepoints[x[, j]])
  x[, c('sd', 'FWHM', 'area')] <- x[, c('sd', 'FWHM', 'area')] * tdiff * tfac
  if (!is.null(x$tau)){
    x[, c('tau')] <- x[, c('tau')] * tdiff * tfac
  } 
  if (!is.null(x$tau_right)){
    x[, c('tau_right', 'tau_left')] <- x[, c('tau_right', 'tau_left')] / (tdiff * tfac)
  } 
  x
}

#' Subset peak list
#'
#' Subsets a `peak_list` while preserving metadata attributes.
#'
#' @param x A `peak_list` object.
#' @param i Indices specifying elements to extract.
#' @param j Not used.
#' @param ... Additional arguments (ignored).
#' @param drop Not used.
#' @return A `peak_list` object.
#' @author Ethan Bass
#' @method [ peak_list
#' @export
`[.peak_list` <- function(x, i, j, ..., drop = FALSE) {
  meta <- attributes(x)
  result <- .subset(x, i)
  exclude = c('names','row.names','dim','dimnames')
  meta[exclude] <- NULL
  attributes(result) <- c(attributes(result), meta)
  result
}

#' Print peak list
#'
#' Prints a summary of a `peak_list` object, including the number of
#' samples, wavelengths, fit method, and total number of peaks.
#'
#' @param x A `peak_list` object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @author Ethan Bass
#' @method print peak_list
#' @export
print.peak_list <- function(x, ...) {
  meta <- attributes(x)
  n_samples <- length(x)
  lambdas <- paste(meta$lambdas, collapse = ", ")
  n_peaks <- sum(sapply(x, function(s) sum(sapply(s, nrow))))
  mean_peaks <- round(n_peaks / n_samples)
  
  cat(sprintf("A peak_list with %d samples and %d wavelength(s) (%s nm)\n",
              n_samples, length(meta$lambdas), lambdas))
  cat(sprintf("Fit method: %s | Time unit: %s | sd_max: %s\n",
              meta$fit, meta$time_unit, meta$sd_max))
  cat(sprintf("Total peaks: %d (mean %s per sample)\n",
              n_peaks, mean_peaks))
  if (!is.na(meta$chrom_list))
    cat(sprintf("Source: %s\n", meta$chrom_list))
  invisible(x)
}
