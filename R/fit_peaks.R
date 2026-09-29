#' Find peaks
#' 
#' Detects peaks in a chromatographic profile by finding zero-crossings in the
#' smoothed first derivative of the signal (`y`). Smoothing suppresses spurious
#' local extrema. Peaks can be filtered by minimum curvature (`min_curvature`),
#' minimum amplitude (`min_amplitude`), or minimum prominence relative to the
#' baseline noise (`min_prominence_snr`).
#' 
#' @details
#' Preprocess with [`preprocess`] before peak detection: high chromatographic
#' resolution can cause peak splitting. Increase `smooth_window` or reduce
#' time-axis resolution via the `dim1` argument to [`preprocess`].
#'
#' Three filters remove small or spurious peaks after detection:
#' `min_curvature`, `min_amplitude`, and `min_prominence_snr`. All three are
#' off by default (`0`) and can be combined.
#'
#' `min_curvature` and `min_amplitude` are in absolute instrument units, so a
#' usable value must be chosen for each dataset. `min_prominence_snr` measures
#' each peak against the noise instead. A peak's prominence is its height above
#' the higher of the two valleys bounding it — equal to its height when the peak
#' stands alone, but less when a neighboring peak elevates one valley. The noise
#' is a single number for the whole chromatogram: `find_peaks` masks every
#' detected peak, then takes the standard deviation of what remains by iterative
#' sigma clipping. `min_prominence_snr = 4` keeps a peak that rises four noise
#' standard deviations above its own valleys, regardless of signal scale. The
#' prominence is local; the noise is not — every peak uses the same estimate.
#'
#' Because the unmasked timepoints include baseline drift, the estimated
#' standard deviation can reflect drift rather than noise, so
#' `min_prominence_snr` is only meaningful on baseline-corrected data.
#' Estimation also requires some timepoints to survive masking: too narrow a
#' `smooth_window` leaves hundreds of noise-driven apices whose boundaries tile
#' the whole trace, in which case the noise cannot be estimated and the filter
#' is skipped with a warning. See *Choosing the smoothing window* for guidance
#' on when to use it.
#' 
#' @section Choosing the smoothing window:
#' Smoothing trades noise reduction against resolution. Too narrow a window
#' fills the peak list with spurious apices from noise-driven sign changes in
#' the derivative. Too wide a window merges neighboring peaks into one
#' detection. Neither raises an error, and a merge is easier to overlook:
#' the merged peak is still fitted and the fit can look reasonable. A merged
#' peak shows up as a fitted width well above the rest and a lower `r.squared`,
#' so plot the peak list over the trace rather than reading the table alone.
#'
#' Widening the window cuts spurious detections by orders of magnitude while
#' costing at most the faintest peak or two, so resolution — losing a peak to
#' a merge — is the failure mode to design around, not recall. Set
#' `smooth_window` to roughly the full width at half maximum of a typical peak,
#' in timepoints, then narrow it if closely-eluting peaks disappear. Give
#' the width in timepoints rather than as a proportion, which follows run length
#' rather than peak width: lengthening a gradient widens the window even though
#' the peaks have not changed.
#'
#' `min_prominence_snr` measures the noise for you rather than asking for a
#' threshold in instrument units, so its cutoff does not depend on signal scale.
#' Set `smooth_window` first, though. Under-smoothing lets
#' noise notch the summit of a small peak, which then registers as two apices
#' with a shallow valley between them; `find_peaks` measures prominence from
#' that interior notch rather than from the foot of the peak, so the peak
#' reads as noise-sized and the filter discards it. A tall peak's summit is
#' too steep to notch, so the loss falls on the small peaks the filter exists
#' to judge. Once the window is wide enough to keep each peak a single apex,
#' the filter clears spurious apices at no cost to the real ones.
#'
#' **Method-specific caveats.** `"savgol"` requires a window around the peak
#' width, and errors below 7 timepoints rather than silently widening:
#' Savitzky-Golay fits a 4th-order polynomial, so a 5-point window is an exact
#' fit that returns the signal unchanged, quietly reducing the method to
#' `"none"`. `"mva"` and `"tmva"` round the window to the nearest odd integer,
#' floored at three (an even-width running mean is not centered on a sample and
#' reports the apex one timepoint off).
#'
#' @importFrom caTools runmean
#' @importFrom minpack.lm nlsLM
#' @importFrom stats deriv lm ksmooth uniroot optimize pnorm
#' @importFrom utils tail
#' @param y A chromatographic signal (as a numeric vector).
#' @param smooth_type Type of smoothing. Either gaussian kernel (`"gaussian"`),
#' box kernel (`"box"`), Savitzky-Golay smoothing (`"savgol"`),
#' moving average (`"mva"`), triangular moving average (`"tmva"`), or 
#' no smoothing (`"none"`).
#' @param smooth_window Width of the smoothing window, used to suppress the
#' spurious local extrema that noise produces in the first derivative. A value
#' between 0 and 1 is interpreted as a proportion of the length of `y`; a value
#' of 1 or greater is interpreted as an absolute number of timepoints. Defaults
#' to `0.001`, a proportion, so the window you get follows the length of the
#' run rather than the width of the peaks: 3 timepoints for a 3,000-point trace
#' but 100 for a 100,000-point one. Set it deliberately, and see *Choosing the
#' smoothing window* for how the moving-average and Savitzky-Golay methods
#' round a non-integer or even-width value to an odd window.
#' @param min_curvature Minimum curvature at the apex, as the magnitude of the
#' second derivative of the smoothed signal. Curvature scales as height over
#' width squared, so larger values exclude broader peaks. Defaults to `0`.
#' @param min_amplitude Minimum peak amplitude, as the intensity of the
#' unsmoothed signal at the apex. Larger values exclude smaller peaks. Defaults
#' to `0`.
#' @param min_prominence_snr Minimum peak prominence, expressed as a multiple
#' of the baseline noise standard deviation. Unlike `min_curvature` and
#' `min_amplitude`, this threshold is scale-free and transfers between datasets.
#' Requires `bounds = TRUE`, since the noise estimate masks peak regions using
#' the peak boundaries. Only meaningful on baseline-corrected data (see
#' *Details*). Defaults to `0`, which applies no filtering.
#' @param bounds Logical. If `TRUE` (default), includes peak boundaries in a
#' `data.frame`.
#' @param slope_thresh The `slope_thresh` argument is deprecated. Please use
#' `min_curvature` instead. (The name is misleading: the filter applies to the
#' second derivative, not the slope.)
#' @param amp_thresh The `amp_thresh` argument is deprecated. Please use
#' `min_amplitude` instead.
#' @return If `bounds == TRUE`, a `data.frame` with the center, start, end, and
#' prominence of each peak. Otherwise, a numeric vector of peak centers. All
#' locations are expressed as indices.
#' @note The `find_peaks` function is adapted from MATLAB code included in
#' Prof. Tom O'Haver's [Pragmatic Introduction to Signal Processing](
#' https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm).
#' @keywords internal
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' find_peaks(Sa_pr[[1]][,"220"])
#' @seealso [`fit_peaks`], [`get_peaks`]
#' @references Oboril, J., Haas, C.P., Lubbesmeyer, M., Nicholls, R.,
#' Gressling, T., Jensen, K.F., Volpin, G., and Hillenbrand, J. 2024. Automated
#' processing of chromatograms: a comprehensive python package with a GUI for
#' intelligent peak identification and deconvolution in chemical reaction
#' analysis. *Digital Discovery* \doi{10.1039/D4DD00214H}. (The noise
#' estimate behind `min_prominence_snr` is adapted from 'MOCCA2'.)
#'
#' O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement. <https://terpconnect.umd.edu/~toh/spectrum/>
#' (Accessed January, 2022).
#' @export find_peaks

find_peaks <- function(y, smooth_type = c("gaussian", "box", "savgol", "mva", 
                                          "tmva", "none"), smooth_window = .001,
                       min_curvature = 0, min_amplitude = 0,
                       min_prominence_snr = 0, bounds = TRUE,
                       slope_thresh = NULL, amp_thresh = NULL){
  min_curvature <- resolve_deprecated(slope_thresh, min_curvature)
  min_amplitude <- resolve_deprecated(amp_thresh, min_amplitude)
  if (!is.vector(y)){
    stop("Please provide a vector to argument `y` to proceed.")
  }
  if (min_prominence_snr > 0 && !bounds){
    stop("`min_prominence_snr` requires `bounds = TRUE`, since the noise ",
         "estimate masks peak regions using the peak boundaries.")
  }
  smooth_type <- match.arg(smooth_type, c("gaussian", "box", "savgol", "mva",
                                          "tmva", "none"))
  if (smooth_window < 1){
    smooth_window <- max(length(y) * smooth_window, 1)
  }
  # compute derivative (with or without smoothing)
  if (smooth_type == "savgol"){
    if (smooth_window < 7){
      stop("`smooth_window` resolves to ", signif(smooth_window, 3),
           " timepoints, which is too narrow for `smooth_type = \"savgol\"`. ",
           "Savitzky-Golay fits a 4th-order polynomial, so a 5-point window ",
           "is an exact fit that returns the signal unchanged and a narrower ",
           "one is underdetermined. Supply `smooth_window` on the order of ",
           "the peak width in timepoints (at least 7).")
    }
    smooth_window <- odd_window(smooth_window)
    d <- savgol(diff(y), fl = smooth_window)
  } else if (smooth_type == "mva"){
    smooth_window <- odd_window(smooth_window)
    d <- caTools::runmean(diff(y), k = smooth_window)
  } else if (smooth_type == 'gaussian'){
    d <- diff(ksmooth(seq_along(y), y, kernel = "normal",
                      bandwidth = smooth_window)$y)
  }  else if (smooth_type == "box"){
    d <- diff(ksmooth(seq_along(y), y, kernel = "box",
                      bandwidth = smooth_window)$y)
  } else if (smooth_type == "tmva"){
    smooth_window <- odd_window(smooth_window)
    d <- runmean(runmean(diff(y), k = smooth_window), k = smooth_window)
  } else{
    d <- diff(y)
  }
  
  p1 <- which(sign(d[1:(length(d)-1)]) > sign(d[2:length(d)])) + 1L
  
  p2 <- which(abs(diff(d)) > min_curvature) + 1L
  
  # detect y-vals exceeding amplitude threshold
  p3 <- which(y > min_amplitude) 
  p <- intersect(intersect(p1,p2), p3)
  if (bounds){
    p4 <- which(sign(d[1:(length(d)-1)]) < sign(d[2:length(d)])) + 1L
    
    # `p4` is ascending, so the flanking valleys are a binary search rather
    # than a scan of every valley per apex
    lo <- findInterval(p - 1L, p4)
    bl <- ifelse(lo == 0L, 1, p4[pmax(lo, 1L)])

    hi <- findInterval(p, p4) + 1L
    bu <- ifelse(hi > length(p4), length(y), p4[pmin(hi, length(p4))])
    na <- which(is.na(y))
    if (length(na)){
      k <- findInterval(p, na)
      bl <- pmax(bl, ifelse(k == 0L, 1L, na[pmax(k, 1L)] + 1L))
      bu <- pmin(bu, ifelse(k >= length(na), length(y),
                            na[pmin(k + 1L, length(na))] - 1L))
    }
    p <- data.frame(pos = p, lower = bl, upper = bu)
    p$prominence <- get_prominence(y, p)
    if (nrow(p) > 0 && min_prominence_snr > 0){
      noise <- estimate_noise(y, p)
      if (noise > 0){
        p <- p[p$prominence > min_prominence_snr * noise, , drop = FALSE]
        rownames(p) <- NULL
      }
    }
  }
  p
}

#' Fit chromatographic peaks
#' 
#' Peak parameters are estimated by fitting one of several models to the data
#' using nonlinear least squares ([`nlsLM`][minpack.lm::nlsLM]). Supported 
#' models include Gaussian, bidirectional exponentially modified Gaussian,
#' or exponential-Gaussian hybrid functions. The area under each fitted curve is 
#' estimated using trapezoidal approximation.
#' 
#' @param x A single chromatogram in matrix format.
#' @param lambda Wavelength at which to fit peaks. Must correspond to a 
#' wavelength present in `x`.
#' @param pos A `data.frame` of peak positions and boundaries, as returned by
#' [`find_peaks`], with columns `pos`, `lower` and `upper` holding index
#' positions into the trace at `lambda`. A bare vector of apex positions is not
#' enough, since every peak model needs a window to fit over. If `NULL` (the
#' default), [`find_peaks`] is run on that trace. Cannot be combined with
#' `windows`.
#' @param windows Optional two-column matrix or `data.frame` of start and end
#' retention times (or a single pair of times), giving manual integration
#' windows. If provided, peak detection is skipped entirely and one peak is
#' returned for every window, whether or not a peak is detectable there. Two
#' adjacent windows sharing a boundary time act as a vertical drop line: their
#' areas sum exactly to the area of their union, which is how a peak is split at
#' a chosen time. `max_fwhm` is ignored, since the width of each window is
#' chosen by the user. Cannot be combined with `pos`. Defaults to
#' `NULL`.
#' @param max_fwhm Maximum peak width, as a full width at half maximum, in
#' index units. Peaks whose `FWHM` could not be measured (`NA`) are retained.
#' Defaults to `NULL`, which applies no width filter. Ignored if `windows` is
#' provided.
#' @param sd_max The `sd_max` argument is deprecated, because a standard
#' deviation is not a peak width for any model except `"gaussian"`. Please use
#' `max_fwhm` instead. If supplied, it is applied as a cutoff on `FWHM` of
#' `2.355 * sd_max`.
#' @param fit Peak model to use. Options include `"bemg"` (bidirectional 
#' exponentially modified Gaussian), `"egh"` (exponential-Gaussian hybrid),
#' `"gaussian"`, and `"raw"`. If `raw` is selected, peaks are integrated using
#' trapezoidal integration without model fitting. Defaults to `"egh"`.
#' @param max_iter Maximum number of iterations to use in nonlinear least
#' squares peak-fitting. Defaults to `1000`.
#' @param estimate_purity Logical. Whether to estimate purity. Defaults
#' to `TRUE`, but is set to `FALSE` when `x` holds a single wavelength.
#' @param purity_args A named list of additional arguments passed to
#' [`get_purity`], e.g. `list(background = "none", metric = "max")`.
#' @param noise_variance Variance of noise, passed to [`get_purity`]. If not
#' provided, it is estimated once for the whole chromatogram rather than once
#' per peak.
#' @param baseline Whether to subtract a local baseline beneath each peak.
#' Options are `"none"` (no baseline, default), `"flat"` (constant offset), or
#' `"sloped"` (linearly varying offset). For the peak models, the offset is
#' estimated by non-linear least squares along with the shape parameters. When
#' `fit = "raw"`, it is instead read directly off the boundaries of the
#' integration window: `"flat"` draws a horizontal line at the lower of the two
#' boundary intensities, while `"sloped"` draws a line between them
#' (valley-to-valley integration). A baseline lying above the signal produces a
#' negative area, which is reported as such, with a warning.
#' @param ... Additional arguments to `find_peaks`.
#' @return A `data.frame` with one row per peak and the following columns:
#' * `rt`: Peak maximum location.
#' * `start`: Peak start.
#' * `end`: Peak end.
#' * `sd`: Standard deviation of the Gaussian component of the fitted model,
#' which is not the width of the peak unless `fit = "gaussian"`. `NA` when
#' `fit = "raw"`.
#' * `tau`: Exponential decay constant (only when `fit = "egh"`).
#' * `tau_right`: Exponential decay constant for right side of peak (only when
#' `fit = "bemg"`).
#' * `tau_left`: Exponential decay constant for left side of peak (only when
#' `fit = "bemg"`).
#' * `h`: The `height` parameter of the fitted BEMG shape function: the value
#' of the pure shape at `center`, before any baseline is added, and distinct
#' from `height` (only when `fit = "bemg"`).
#' * `center`: The center parameter of the fitted BEMG shape function, the
#' location around which the exponential tails are applied. May differ from
#' `rt` (only when `fit = "bemg"`).
#' * `FWHM`: The full width at half maximum of the peak, comparable across
#' values of `fit`. `NA` for a `"raw"` peak whose integration window does not
#' fall below half its height on both sides.
#' * `height`: Peak height above the fitted baseline.
#' * `area`: Peak area above the fitted baseline.
#' * `r.squared`: Coefficient of determination (\eqn{R^2}), from a linear
#' regression of fitted values on the raw data. Statistically informal, since
#' the model was fit by nonlinear least squares rather than OLS, but useful as
#' a fit diagnostic.
#' * `purity`: The spectral purity of peak as assessed by [`get_purity`].
#' * `window_limited`: `1` if a fitted shape parameter came to rest on a bound
#' derived from the integration window, meaning the window rather than the data
#' determined the reported shape, and `0` otherwise. Not returned when
#' `fit = "raw"`.
#' * `converged`: `1` if the optimizer reported convergence for the final stage
#' of the fit, and `0` otherwise. A `0` also covers a stage skipped for holding
#' fewer timepoints than the model has parameters, in which case the reported
#' parameters are starting estimates rather than a measurement. Not returned
#' when `fit = "raw"`.
#' * `floor`: Constant baseline offset (only when `baseline = "flat"`).
#' * `floor_start`: Baseline level at the left edge of the peak window (only
#' when `baseline = "sloped"`).
#' * `floor_end`: Baseline level at the right edge of the peak window (only
#' when `baseline = "sloped"`).
#' The table is in index units rather than retention times: `rt`, `start`,
#' `end` and `center` are index positions; `sd`, `FWHM`, `area` and `tau` are
#' widths on the index scale; and `tau_right` and `tau_left` are in inverse
#' index units. [`get_peaks`] converts them to the chromatogram's time axis.
#' @note The `fit_peaks` function is adapted from Dr. Robert Morrison's
#' [DuffyTools](https://github.com/robertdouglasmorrison/DuffyTools) package
#' as well as code published in Ron Wehrens' [alsace](
#' https://github.com/rwehrens/alsace) package.
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' fit_peaks(Sa_pr[[1]], lambda = 220)
#' @seealso [`find_peaks`], [`get_peaks`]
#' @references
#' * Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. *Journal of
#' Chromatography A* 915:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' * Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions - A
#' good model for chromatographic peaks in isocratic HPLC? *Chromatographia*,
#' 26: 285-296. \doi{10.1007/BF02268168}.
#' @export fit_peaks
#' @md
#' @keywords internal

fit_peaks <- function (x, lambda, pos = NULL, windows = NULL,
                       max_fwhm = NULL,
                       sd_max = NULL,
                       fit = c("egh", "bemg", "gaussian", "raw"),
                       max_iter = 1000, estimate_purity = TRUE,
                       purity_args = list(), noise_variance = NULL,
                       baseline = c("none", "flat", "sloped"), ...){
  lambda.idx <- get_lambda_idx(lambda, get_lambdas(x), allow_max = FALSE)
  y <- x[,lambda.idx]
  # `find_peaks` returns a bare vector of apex positions when `bounds = FALSE`,
  # which cannot be integrated, since every peak model needs a window.
  dots <- list(...)
  if ("bounds" %in% names(dots) && !isTRUE(dots$bounds)){
    stop("`bounds = FALSE` cannot be used with `fit_peaks`, since peak boundaries are required in order to integrate. It is only meaningful when calling `find_peaks` directly.")
  }
  if (!is.null(sd_max)){
    if (!missing(max_fwhm))
      stop("Supply either `max_fwhm` or `sd_max`, not both. `max_fwhm` is a ",
           "full width at half maximum in timepoints, whereas `sd_max` is a ",
           "standard deviation in timepoints.")
    warning("The `sd_max` argument is deprecated. It is expressed as a ",
            "standard deviation, which is not a peak width for any model ",
            "except `\"gaussian\"`. Please use `max_fwhm` instead, which is ",
            "expressed as a full width at half maximum.", immediate. = TRUE)
    max_fwhm <- FWHM_CONST * sd_max
  }
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  baseline <- match.arg(baseline, c("none", "flat", "sloped"))
  if (!is.null(pos) && !(is.data.frame(pos) &&
                         all(c("pos", "lower", "upper") %in% names(pos)))){
    stop("`pos` must be a data.frame with columns `pos`, `lower` and `upper`, ",
         "as returned by `find_peaks`.")
  }
  if (!is.null(windows)){
    if (!is.null(pos)){
      stop("Please supply either `windows` or `pos`, not both.")
    }
    pos <- windows_to_pos(windows, y = y, times = get_times(x))
    # the width of a manual window is chosen by the user, so it must not be
    # filtered.
    max_fwhm <- NULL
  }
  if (is.null(pos)){
    pos <- find_peaks(y = y, ...)
  }
  if (ncol(x) == 1){
    estimate_purity <- FALSE
  }
  tabnames <- switch(fit,
                     "gaussian" = c("rt", "start", "end", "sd", "FWHM",
                                    "height", "area", "r.squared", "purity",
                                    "window_limited", "converged"),
                     "egh" = c("rt", "start", "end", "sd", "tau", "FWHM",
                               "height", "area", "r.squared", "purity",
                               "window_limited", "converged"),
                     "bemg" = c("rt", "start", "end", "sd", "tau_right",
                                "tau_left", "h", "center", "FWHM", "height",
                                "area", "r.squared", "purity",
                                "window_limited", "converged"),
                     "raw" = c("rt", "start", "end", "sd", "FWHM", "height",
                               "area", "purity")
  )
  tabnames <- c(tabnames, floor_cols(baseline))
  noPeaksMat <- as.data.frame(matrix(numeric(0), nrow = 0,
                                     ncol = length(tabnames),
                                     dimnames = list(NULL, tabnames)),
                              check.names = FALSE)
  if (nrow(pos) == 0)
    return(noPeaksMat)
  on.edge <- sapply(pos$pos, function(x){
    x <= 1 || is.na(y[x + 1]) || is.na(y[x - 1])
  })
  if (!is.null(windows) && any(on.edge)){
    warning(sprintf(
      "The apex of the following window%s falls on the edge of the chromatogram, so %s dropped: %s.",
      plural(sum(on.edge)),
      ifelse(sum(on.edge) > 1, "they were", "it was"),
      paste(sprintf("%s-%s", signif(get_times(x)[pos$lower[on.edge]], 4),
                    signif(get_times(x)[pos$upper[on.edge]], 4)), collapse = ", ")),
      immediate. = TRUE)
  }
  pos <- pos[!on.edge, , drop = FALSE]

  if (nrow(pos) == 0)
    return(noPeaksMat)
  fitpk <- switch(fit,
                  "gaussian" = fitpk_gaussian,
                  "egh" = fitpk_egh,
                  "raw" = fitpk_raw,
                  "bemg" = fitpk_bemg)
  if (estimate_purity && is.null(noise_variance)){
    noise_variance <- hoist_noise_variance(x, purity_args)
  }
  huhn <- data.frame(t(apply(pos, 1, fitpk, x = x,
                             lambda = lambda.idx, max_iter = max_iter,
                             estimate_purity = estimate_purity,
                             noise_variance = noise_variance,
                             purity_args = purity_args,
                             baseline = baseline)))
  colnames(huhn) <- tabnames
  huhn <- data.frame(sapply(huhn, as.numeric, simplify = FALSE))
  huhn <- drop_unfittable_peaks(huhn)
  huhn <- huhn[keep_by_width(huhn$FWHM, max_fwhm), , drop = FALSE]
  df <- try(huhn[which(huhn$rt > 0), ], silent = TRUE)
  if(inherits(df, "try-error")) NA else df
}

#' Drop peaks whose model did not yield a usable intensity
#'
#' A model that diverges can leave finite coefficients alongside a non-finite
#' fitted curve, so `height` and `area` are the columns that reveal it. Such a
#' row is a failed fit rather than a measurement, and a non-finite `response`
#' propagates into [`get_peaktable`], so it is removed here with a warning
#' rather than passed on.
#'
#' @param huhn A `data.frame` of fitted peak parameters.
#' @return `huhn`, without rows whose `height` or `area` is non-finite.
#' @author Ethan Bass
#' @noRd
drop_unfittable_peaks <- function(huhn){
  cols <- intersect(c("height", "area"), colnames(huhn))
  if (length(cols) == 0 || nrow(huhn) == 0) return(huhn)
  bad <- !Reduce(`&`, lapply(huhn[cols], is.finite))
  if (any(bad)){
    warning(sprintf(
      "Discarded %d peak%s whose fitted %s not finite, indicating the peak model did not converge.",
      sum(bad), plural(sum(bad)),
      if (length(cols) > 1) "height or area was" else paste(cols, "was")),
      call. = FALSE)
  }
  huhn[!bad, , drop = FALSE]
}

#' Gaussian function
#' @note: Adapted from <https://github.com/robertdouglasmorrison/DuffyTools/blob/master/R/gaussian.R>
#' @noRd
gaussian <- function(x, center = 0, width = 1, height = NULL, floor = 0, slope = 0){

  # adapted from Earl F. Glynn;  Stowers Institute for Medical Research, 2007
  twoVar <- 2 * width * width
  sqrt2piVar <- sqrt( pi * twoVar)
  y <- exp( -( x - center)^2 / twoVar) / sqrt2piVar

  # by default, the height is such that the curve has unit volume
  if ( ! is.null (height)) {
    scalefactor <- sqrt2piVar
    y <- y * scalefactor * height
  }
  y + floor + slope * (x - x[1])
}

#' Fit gaussian peak
#' @importFrom stats coef fitted lm quantile residuals
#' @noRd
fit_gaussian <- function(x, y, start.center = NULL,
                         start.width = NULL, start.height = NULL,
                         start.floor = NULL, baseline = "none",
                         max_iter = 1000){
  st <- fit_starts(x, y, start.center, start.width, start.height)
  start.center <- st$center
  start.width <- st$width
  start.height <- st$height
  peak_width <- st$span
  shape_lower <- st$lower
  shape_upper <- st$upper

  controlList <- minpack.lm::nls.lm.control(maxiter = max_iter,
                                           maxfev = 10 * max_iter)
  starts <- list("center" = start.center, "width" = start.width, "height" = start.height)

  # Stage 1: shape only
  nlsAns <- skip_underdetermined(length(x), 3) %||%
    try(minpack.lm::nlsLM(y ~ gaussian(x, center, width, height),
                       start = starts, control = controlList,
                       lower = shape_lower, upper = shape_upper), silent = TRUE)
  nls_stage1 <- nlsAns

  upper_floor <- max(min(y), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- skip_underdetermined(length(x), 4) %||%
      try(minpack.lm::nlsLM(
      y ~ gaussian(x, center, width, height, floor),
      start = starts, control = controlList,
      lower = c(shape_lower, floor = 0),
      upper = c(shape_upper, floor = upper_floor)
    ), silent = TRUE)
  }
  nls_stage2 <- nlsAns

  # Stage 3: slope via (floor_s, floor_e) so both endpoints are bounded >= 0
  # and m$y stays consistent with baseline_vec in fitpk_*
  if (baseline == "sloped") {
    f0 <- if (!inherits(nlsAns, "try-error")) unname(coef(nlsAns)["floor"])
          else if (!is.null(start.floor)) start.floor else 0
    s3 <- c(as.list(if (!inherits(nlsAns, "try-error"))
                      coef(nlsAns)[c("center", "width", "height")]
                    else c(center = start.center, width = start.width, height = start.height)),
            floor_s = f0, floor_e = f0)
    nlsAns <- skip_underdetermined(length(x), 5) %||%
      try(minpack.lm::nlsLM(
      y ~ gaussian(x, center, width, height, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(shape_lower, floor_s = 0, floor_e = 0),
      upper = c(shape_upper, floor_s = Inf, floor_e = Inf)
    ), silent = TRUE)
  }

  if (inherits(nlsAns, "try-error")) {
    if (baseline == "sloped" && !inherits(nls_stage2, "try-error")) {
      fb <- coef(nls_stage2)
      floorAns <- unname(fb["floor"]); slopeAns <- 0
      yAns <- gaussian(x, fb[1], fb[2], fb[3], floorAns, slopeAns)
      out <- list("center" = fb[1], "width" = fb[2], "height" = fb[3],
                  "y" = yAns, "residual" = y - yAns)
    } else if (!inherits(nls_stage1, "try-error")) {
      fb <- coef(nls_stage1)
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- gaussian(x, fb[1], fb[2], fb[3], floorAns, slopeAns)
      out <- list("center" = fb[1], "width" = fb[2], "height" = fb[3],
                  "y" = yAns, "residual" = y - yAns)
    } else {
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- gaussian(x, start.center, start.width, start.height, floorAns, slopeAns)
      out <- list("center" = start.center, "width" = start.width,
                  "height" = start.height, "y" = yAns, "residual" = y - yAns)
    }
  } else {
    coefs <- coef(nlsAns)
    out <- list("center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                "y" = fitted(nlsAns), "residual" = residuals(nlsAns))
    if (baseline == "sloped") {
      floorAns <- coefs["floor_s"]
      slopeAns  <- (coefs["floor_e"] - coefs["floor_s"]) / peak_width
    } else {
      floorAns <- if (baseline != "none") coefs["floor"] else 0; slopeAns <- 0
    }
  }
  out$window_limited <- is_window_limited(out, shape_lower, shape_upper)
  out$converged <- !inherits(nlsAns, "try-error") &&
    isTRUE(nlsAns$convInfo$isConv)
  if (baseline != "none")    out <- c(out, "floor" = unname(floorAns))
  if (baseline == "sloped")  out <- c(out, "slope" = unname(slopeAns))
  return(out)
}

###############################################################################

#' Expontential-gaussian hybrid
#' @noRd
egh <- function(x, center, width, height, tau, floor = 0, slope = 0){
  result <- rep(0, length(x))
  index <- which(2*width^2 + tau*(x-center)>0)
  result[index] <- height*exp(-(x[index] - center)^2/(2*width^2 + tau*(x[index] - center)))
  return(result + floor + slope * (x - x[1]))
}


#' Fit exponential-gaussian hybrid peak
#' @importFrom stats coef fitted lm quantile residuals
#' @noRd
fit_egh <- function(x1, y1, start.center = NULL, start.width = NULL,
                    start.tau = NULL, start.height = NULL, start.floor = NULL,
                    baseline = "none", max_iter = 1000) {

  if (is.null(start.tau))    start.tau    <- 0

  st <- fit_starts(x1, y1, start.center, start.width, start.height)
  start.center <- st$center
  start.width <- st$width
  start.height <- st$height
  peak_width <- st$span
  shape_lower <- c(st$lower, tau = -Inf)
  shape_upper <- c(st$upper, tau = Inf)

  controlList <- minpack.lm::nls.lm.control(maxiter = max_iter,
                                           maxfev = 10 * max_iter)
  starts <- list("center" = start.center, "width" = start.width,
                 "height" = start.height, "tau" = start.tau)

  # Stage 1: shape only
  nlsAns <- skip_underdetermined(length(x1), 4) %||%
    try(minpack.lm::nlsLM(y1 ~ egh(x1, center, width, height, tau),
                      start = starts, control = controlList,
                      lower = shape_lower, upper = shape_upper), silent = TRUE)
  nls_stage1 <- nlsAns

  upper_floor <- max(min(y1), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- skip_underdetermined(length(x1), 5) %||%
      try(minpack.lm::nlsLM(
      y1 ~ egh(x1, center, width, height, tau, floor),
      start = starts, control = controlList,
      lower = c(shape_lower, floor = 0),
      upper = c(shape_upper, floor = upper_floor)
    ), silent = TRUE)
  }
  nls_stage2 <- nlsAns

  # Stage 3: slope via (floor_s, floor_e) so both endpoints are bounded >= 0
  if (baseline == "sloped") {
    f0 <- if (!inherits(nlsAns, "try-error")) unname(coef(nlsAns)["floor"])
          else if (!is.null(start.floor)) start.floor else 0
    s3 <- c(as.list(if (!inherits(nlsAns, "try-error"))
                      coef(nlsAns)[c("center", "width", "height", "tau")]
                    else c(center = start.center, width = start.width,
                           height = start.height, tau = start.tau)),
            floor_s = f0, floor_e = f0)
    nlsAns <- skip_underdetermined(length(x1), 6) %||%
      try(minpack.lm::nlsLM(
      y1 ~ egh(x1, center, width, height, tau, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(shape_lower, floor_s = 0, floor_e = 0),
      upper = c(shape_upper, floor_s = Inf, floor_e = Inf)
    ), silent = TRUE)
  }

  if (inherits(nlsAns, "try-error")) {
    if (baseline == "sloped" && !inherits(nls_stage2, "try-error")) {
      fb <- coef(nls_stage2)
      floorAns <- unname(fb["floor"]); slopeAns <- 0
      yAns <- egh(x1, fb[1], fb[2], fb[3], fb[4], floorAns, slopeAns)
      out <- list("center" = fb[1], "width" = fb[2], "height" = fb[3],
                  "tau" = fb[4], "y" = yAns, "residual" = y1 - yAns)
    } else if (!inherits(nls_stage1, "try-error")) {
      fb <- coef(nls_stage1)
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- egh(x1, fb[1], fb[2], fb[3], fb[4], floorAns, slopeAns)
      out <- list("center" = fb[1], "width" = fb[2], "height" = fb[3],
                  "tau" = fb[4], "y" = yAns, "residual" = y1 - yAns)
    } else {
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- egh(x1, start.center, start.width, start.height, start.tau, floorAns, slopeAns)
      out <- list("center" = start.center, "width" = start.width,
                  "height" = start.height, "tau" = start.tau,
                  "y" = yAns, "residual" = y1 - yAns)
    }
  } else {
    coefs <- coef(nlsAns)
    out <- list("center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                "tau" = coefs[4], "y" = fitted(nlsAns), "residual" = residuals(nlsAns))
    if (baseline == "sloped") {
      floorAns <- coefs["floor_s"]
      slopeAns  <- (coefs["floor_e"] - coefs["floor_s"]) / peak_width
    } else {
      floorAns <- if (baseline != "none") coefs["floor"] else 0; slopeAns <- 0
    }
  }

  out$window_limited <- is_window_limited(out, shape_lower, shape_upper)
  out$converged <- !inherits(nlsAns, "try-error") &&
    isTRUE(nlsAns$convInfo$isConv)
  if (baseline != "none")    out <- c(out, "floor" = unname(floorAns))
  if (baseline == "sloped")  out <- c(out, "slope" = unname(slopeAns))
  return(out)
}

## Conversion from Gaussian standard deviation to full width at half maximum:
## 2 * sqrt(2 * log(2)). Previously hard-coded as the rounded 2.35 while the
## documentation stated 2.355.
FWHM_CONST <- 2 * sqrt(2 * log(2))

#' Full width at half maximum of an exponential-gaussian hybrid
#' @noRd
egh_fwhm <- function(width, tau){
  unname(sqrt((log(2) * tau)^2 + 8 * log(2) * width^2))
}

#' Half-height width of a unimodal shape function
#' @noRd
shape_fwhm <- function(f, apex, scale, max_mult = 1024){
  half <- f(apex) / 2
  if (!is.finite(half) || half <= 0) return(NA_real_)
  if (!is.finite(scale) || scale <= 0) return(NA_real_)
  side <- function(direction){
    d <- scale
    while (d <= max_mult * scale && isTRUE(f(apex + direction * d) > half)){
      d <- 2 * d
    }
    if (!isTRUE(f(apex + direction * d) <= half)) return(NA_real_)
    uniroot(function(t) f(apex + direction * t) - half,
                   c(0, d), tol = 1e-8)$root
  }
  left <- side(-1)
  right <- side(1)
  if (is.na(left) || is.na(right)) NA_real_ else unname(left + right)
}

#' Half-height width of a fitted bi-exponentially modified gaussian
#' @noRd
bemg_fwhm <- function(m){
  scale <- max(unname(m$width), 1 / abs(unname(m$a)), 1 / abs(unname(m$b)))
  span <- 20 * scale
  if (!is.finite(span) || span <= 0) return(NA_real_)
  shape <- function(t) bemg(t, m$center, m$width, m$height, m$a, m$b)
  center <- unname(m$center)
  apex <- optimize(shape, c(center - span, center + span),
                          maximum = TRUE)$maximum
  shape_fwhm(shape, apex, scale)
}

#' Half-height width of a sampled peak, by linear interpolation
#' @noRd
measured_fwhm <- function(y){
  n <- length(y)
  if (n < 3) return(NA_real_)
  apex <- which.max(y)
  half <- y[[apex]] / 2
  if (!is.finite(half) || half <= 0) return(NA_real_)
  below <- which(y[seq_len(apex)] < half)
  if (length(below) == 0) return(NA_real_)
  i <- max(below)
  left <- i + (half - y[[i]]) / (y[[i + 1]] - y[[i]])
  below <- which(y[apex:n] < half)
  if (length(below) == 0) return(NA_real_)
  j <- apex + min(below) - 1
  right <- (j - 1) + (y[[j - 1]] - half) / (y[[j - 1]] - y[[j]])
  unname(right - left)
}

#' Whether a fitted shape came to rest on a window-derived limit
#'
#' Most bounds passed to `nlsLM` are derived from the integration window, so a
#' coefficient sitting on one means the window, rather than the data, decided
#' the reported shape. Infinite bounds are skipped, since they constrain
#' nothing, and so are the bounds that do not come from the window: the lower
#' bounds on `width` and on the `bemg` rates `a` and `b`, which are positivity
#' constraints, and the bounds on `height`.
#'
#' @param out A fitted model, as returned by `fit_gaussian`/`fit_egh`/`fit_bemg`.
#' @param lower,upper The named bound vectors given to `nlsLM`.
#' @param tol Relative tolerance for treating a value as resting on its bound.
#' @return `TRUE` if any bounded coefficient lies on its bound.
#' @author Ethan Bass
#' @noRd
is_window_limited <- function(out, lower, upper, tol = 1e-6){
  on_bound <- function(b){
    p <- intersect(names(b), names(out))
    if (length(p) == 0) return(FALSE)
    v <- vapply(out[p], function(z) unname(as.numeric(z))[1], numeric(1))
    b <- b[p]
    f <- is.finite(b)
    any(abs(v[f] - b[f]) <= tol * pmax(1, abs(b[f])), na.rm = TRUE)
  }
  on_bound(lower[!names(lower) %in% c("width", "height", "a", "b")]) ||
    on_bound(upper[names(upper) != "height"])
}

#' A skipped fit, in the shape of a failed one
#'
#' `minpack.lm::nls.lm` writes out of bounds when given fewer observations than
#' parameters, which corrupts the heap and eventually segfaults. Such a window
#' cannot determine the model anyway, so the stage is skipped and reported as a
#' failure, letting the caller fall back as it would for any fit that did not
#' converge.
#'
#' @param n Number of observations.
#' @param p Number of parameters.
#' @return A `try-error` if the fit must be skipped, otherwise `NULL`.
#' @author Ethan Bass
#' @noRd
skip_underdetermined <- function(n, p){
  if (n >= p) return(NULL)
  structure(sprintf("Error : window of %d points cannot determine %d parameters\n", n, p),
            class = "try-error")
}

#' Reject a bi-exponential fit whose tail grows away from the apex
#'
#' `a` and `b` are the rates of the two exponentials convolved into the
#' gaussian, so a negative one describes a tail that grows without bound
#' instead of decaying. Over a finite integration window that can imitate a
#' rising baseline or an unresolved shoulder and so fit well, but it is not a
#' peak: it has no half-height crossing on that side, which leaves `FWHM`
#' unmeasurable and lets the peak slip past a width filter. The bounds cannot
#' prevent it, because a finite bound changes the search for every peak and not
#' only for those near it, so the fit is rejected afterwards and handled like
#' any other that did not converge.
#'
#' @param ans The result of an `nlsLM` call, or a `try-error`.
#' @return `ans`, or a `try-error` if either rate is negative.
#' @author Ethan Bass
#' @noRd
reject_growing_tail <- function(ans){
  if (inherits(ans, "try-error")) return(ans)
  cf <- try(coef(ans), silent = TRUE)
  if (inherits(cf, "try-error")) return(ans)
  if (any(cf[c("a", "b")] < 0, na.rm = TRUE)){
    return(structure("Error : fitted tail grows away from the apex\n",
                     class = "try-error"))
  }
  ans
}

#' Peaks retained by a width filter
#' @noRd
keep_by_width <- function(fwhm, cutoff){
  if (is.null(cutoff) || !is.finite(cutoff)) return(rep(TRUE, length(fwhm)))
  is.na(fwhm) | fwhm <= cutoff
}

#' Spacing of a fitting abscissa
#'
#' The heuristic start widths count timepoints above half height, which is a
#' width only when the abscissa is in index units. `fitpk_*` passes indices, but
#' the fitters are also called directly with a time axis, so the count is scaled
#' by the spacing to give a width in the units of `x`.
#'
#' @param x A numeric vector of timepoints or indices.
#' @return The median spacing of `x`, or `1` if it cannot be determined.
#' @author Ethan Bass
#' @noRd
x_spacing <- function(x){
  d <- axis_dt(x)
  if (!is.finite(d) || d <= 0) 1 else d
}

## Largest fitted peak width, as a multiple of the width of the integration
## window, admitted by the peak models. `find_peaks` bounds a peak at its
## nearest valleys, which for a resolved peak sit near +/- 2 sigma, so a window
## is typically about 1.7 FWHM and a limit of exactly one window leaves almost
## no headroom. Two windows admits a peak whose bounds were clipped by a
## shoulder or by the end of the trace while still excluding a fit many times
## wider than the data supporting it.
##
## The same multiple also bounds the `egh` tail, which is what actually sets the
## width of a skewed fit: `log(2) * tau` is tau's share of the `egh` FWHM (see
## `egh_fwhm`), so holding it to two windows caps the whole peak at 2.83 windows
## (`sqrt(2) * mult`, with the gaussian and tail terms both at their bound).
##
## It deliberately does NOT bound the `bemg` rates, for two measured reasons.
## `bemg` trades its gaussian core against its two exponential tails, so a rate
## floor does not bound the peak: the optimizer routes the width through the
## core instead and reaches a worse fit. On the `Sa` test data a floor of
## `1/(2 * span)` inflated `sd` by twentyfold or more on 37 of 128 peaks and
## doubled the number of fits below an r-squared of 0.8. Even a floor of zero,
## which is a physical constraint rather than a width one, moved 9 of 128 peak
## areas by more than 10% -- none of which had a negative rate to begin with --
## because a finite bound switches `nls.lm` into its bounded code path and
## changes the search for every peak, not only for those near the bound. A
## negative rate does describe an exponentially growing tail rather than a peak,
## but that has to be caught after the fit, not by constraining the optimizer.
FIT_WIDTH_MAX_MULT <- 2

## Fraction of detected peaks whose integration window cannot determine the
## peak model, above which `get_peaks` reports that the peak list is dominated
## by windows too small to fit. On a well sampled run essentially no peak should
## be underdetermined: across two datasets the healthy case ran from 1.4% to
## 10.1% while a run in which noise was being detected as peaks ran from 28.6%
## to 43.4%, depending on the model. A fifth sits between the two with room on
## either side, but it is a judgement calibrated on those two datasets rather
## than a property of the method.
UNDERDETERMINED_WARN_FRAC <- 0.2

## Fraction of detected peaks discarded by `max_fwhm`/`sd_max`, above which
## `get_peaks` warns that the cutoff may be too tight. Not independently
## calibrated; reuses `UNDERDETERMINED_WARN_FRAC`'s threshold since both flag
## the same failure shape (peaks silently missing from the output).
WIDTH_FILTER_WARN_FRAC <- 0.2

#' Parameter bounds derived from the integration window
#'
#' @param x1 The timepoints (or indices) of the integration window.
#' @param mult Largest admissible fitted `FWHM`, as a multiple of the span of
#' `x1`.
#' @return A list with the window `span`, the upper bound on `width`, and the
#' bounds on `center`.
#' @author Ethan Bass
#' @noRd
fit_window_bounds <- function(x1, mult = FIT_WIDTH_MAX_MULT){
  span <- x1[length(x1)] - x1[1]
  if (!is.finite(span) || span <= 0) span <- 1
  margin <- (mult - 1) / 2 * span
  list(span = span, width_max = mult * span / FWHM_CONST,
       center_min = x1[1] - margin, center_max = x1[length(x1)] + margin)
}

#' Clamp a starting value into its bounds
#'
#' `nlsLM` projects a start value that lies outside `lower` or `upper` onto the
#' bound, reaching a worse optimum than a feasible start, and the heuristic
#' start widths are counts of timepoints above half height, which can exceed
#' the window-derived bound.
#'
#' @param v A starting value.
#' @param lo,hi Bounds.
#' @return `v`, moved inside `[lo, hi]`.
#' @author Ethan Bass
#' @noRd
clamp_start <- function(v, lo, hi){
  if (!is.finite(v)) v <- (lo + hi) / 2
  min(max(v, lo), hi)
}

#' Starting values and shape bounds shared by the peak models
#'
#' The heuristic start width counts timepoints above half height, so it is
#' scaled by the abscissa spacing, then clamped into the window bounds along
#' with the center. Models with a tail parameter append its bounds to `lower`
#' and `upper` themselves, in the position the parameter occupies in `start`,
#' since `nlsLM` matches bounds by position rather than by name.
#' @noRd
fit_starts <- function(x, y, start.center, start.width, start.height,
                       height_mult = NULL){
  who.max <- which.max(y)
  if (is.null(start.center)) start.center <- x[who.max]
  if (is.null(start.height)) start.height <- y[who.max]
  if (is.null(start.width))
    start.width <- sum(y > (start.height / 2)) / 2 * x_spacing(x)
  bnd <- fit_window_bounds(x)
  list(center = clamp_start(start.center, bnd$center_min, bnd$center_max),
       width = clamp_start(abs(start.width), 0, bnd$width_max),
       height = start.height,
       span = bnd$span,
       lower = c(center = bnd$center_min, width = 0, height = -Inf),
       upper = c(center = bnd$center_max, width = bnd$width_max,
                 height = if (is.null(height_mult)) Inf
                          else height_mult * start.height))
}

#' Number of free parameters in a peak model
#'
#' The shape parameters of the model, plus the baseline terms implied by
#' `baseline`. A window holding fewer timepoints than this cannot determine the
#' model, whatever the data look like.
#'
#' @param fit One of `"gaussian"`, `"egh"` or `"bemg"`.
#' @param baseline One of `"none"`, `"flat"` or `"sloped"`.
#' @return The number of free parameters, or `NA` for a model without a fit.
#' @author Ethan Bass
#' @noRd
n_fit_params <- function(fit, baseline){
  shape <- switch(fit, "gaussian" = 3L, "egh" = 4L, "bemg" = 5L, NA_integer_)
  if (is.na(shape)) return(NA_integer_)
  shape + length(floor_cols(baseline))
}

floor_cols <- function(baseline){
  switch(baseline, "none" = character(0), "flat" = "floor",
         "sloped" = c("floor_start", "floor_end"))
}

#' Purity of one peak
#' @noRd
peak_purity <- function(x, pos, estimate_purity, noise_variance, purity_args){
  purity_args[c("x", "pos", "try", "noise_variance")] <- NULL
  do.call(get_purity, c(list(x = x, pos = pos, try = estimate_purity,
                             noise_variance = noise_variance), purity_args))
}

#' Trapezoidal area above a local baseline
#' @noRd
integrate_above_baseline <- function(y_fit, peak.loc, floor_val, slope_val){
  y_above <- y_fit - (floor_val + slope_val * (peak.loc - peak.loc[1]))
  list(y_above = y_above,
       area = sum(diff(peak.loc) * (y_above[-length(y_above)] + y_above[-1]) / 2))
}

#' Local baseline implied by a fitted peak model
#' @noRd
get_local_baseline <- function(m, baseline){
  list(floor = if (baseline != "none") unname(m$floor) else 0,
       slope = if (baseline == "sloped") unname(m$slope) else 0)
}

#' Append the local baseline columns to a peak's result vector
#' @noRd
append_floor_cols <- function(result, baseline, floor_val, slope_val, peak.loc){
  switch(baseline,
         "none" = result,
         "flat" = c(result, floor = floor_val),
         "sloped" = c(result, floor_start = floor_val,
                      floor_end = floor_val +
                        slope_val * (peak.loc[length(peak.loc)] - peak.loc[1])))
}

fitpk_bemg <- function(x, pos, lambda, max_iter,
                       estimate_purity = TRUE, noise_variance = NULL,
                       purity_args = list(), baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_bemg(peak.loc, y[peak.loc], start.center = xloc,
                                 start.height = y[[xloc]], max_iter = max_iter,
                                 baseline = baseline)
  )
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- peak_purity(x, pos, estimate_purity, noise_variance, purity_args)
  bl <- get_local_baseline(m, baseline)
  int <- integrate_above_baseline(m$y, peak.loc, bl$floor, bl$slope)
  y_above <- int$y_above
  area <- int$area
  apex_idx <- which.max(y_above)
  rt_val <- if (length(apex_idx) == 1L) peak.loc[apex_idx] else unname(round(m$center))
  result <- c("rt" = rt_val, "start" = pos[[2]], "end" = pos[[3]],
    "sd" = unname(m$width), "tau_right" = unname(m$b), "tau_left" = unname(m$a),
    "h" = unname(m$height), "center" = unname(round(m$center)),
    "FWHM" = bemg_fwhm(m),
    "height" = unname(max(y_above)),
    "area" = area, "r.squared" = unname(r.squared),
    purity = purity, "window_limited" = as.numeric(m$window_limited),
    "converged" = as.numeric(m$converged))
  append_floor_cols(result, baseline, bl$floor, bl$slope, peak.loc)
}

#' Fit peak (gaussian)
#' @noRd
fitpk_gaussian <- function(x, pos, lambda, max_iter,
                           estimate_purity = TRUE, noise_variance = NULL,
                           purity_args = list(), baseline = "none", ...){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_gaussian(peak.loc, y[peak.loc],
                                     start.center = xloc,
                                     start.height = y[[xloc]],
                                     max_iter = max_iter,
                                     baseline = baseline)
  )
  bl <- get_local_baseline(m, baseline)
  area <- integrate_above_baseline(m$y, peak.loc, bl$floor, bl$slope)$area
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- peak_purity(x, pos, estimate_purity, noise_variance, purity_args)
  result <- c("rt" = unname(m$center), "start" = pos[2], "end" = pos[3],
    "sd" = unname(m$width), "FWHM" = FWHM_CONST * unname(m$width),
    "height" = unname(m$height), "area" = area,
    "r.squared" = unname(r.squared), purity = purity,
    "window_limited" = as.numeric(m$window_limited),
    "converged" = as.numeric(m$converged))
  append_floor_cols(result, baseline, bl$floor, bl$slope, peak.loc)
}

#' Fit peak (exponential-gaussian hybrid)
#' @noRd
fitpk_egh <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_variance = NULL,
                      purity_args = list(), baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc,
                                start.height = y[[xloc]], max_iter = max_iter,
                                baseline = baseline)
  )
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- peak_purity(x, pos, estimate_purity, noise_variance, purity_args)
  bl <- get_local_baseline(m, baseline)
  area <- integrate_above_baseline(m$y, peak.loc, bl$floor, bl$slope)$area
  result <- c("rt" = unname(m$center), "start" = pos[2], "end" = pos[3],
    "sd" = unname(m$width), "tau" = unname(m$tau),
    "FWHM" = egh_fwhm(m$width, m$tau),
    "height" = unname(m$height), "area" = area,
    "r.squared" = unname(r.squared), purity = purity,
    "window_limited" = as.numeric(m$window_limited),
    "converged" = as.numeric(m$converged))
  append_floor_cols(result, baseline, bl$floor, bl$slope, peak.loc)
}

#' Fit peak (raw)
#' @noRd
fitpk_raw <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_variance = NULL,
                      purity_args = list(), baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  yp <- y[peak.loc]

  # unlike the model fits, which estimate `floor` by non-linear least squares,
  # the local baseline is read directly off the boundaries of the window, since
  # those were chosen by the user. `min(yp)` is not used for the flat baseline
  # because `preprocess` clips negative intensities to zero, so a corrected
  # trace typically touches zero somewhere in a wide window, which would make
  # `baseline = "flat"` indistinguishable from `baseline = "none"`.
  n <- length(yp)
  floor_val <- switch(baseline,
                      "none" = 0,
                      "flat" = min(yp[1], yp[n]),
                      "sloped" = yp[1])
  slope_val <- if (baseline == "sloped" && n > 1){
    (yp[n] - yp[1]) / (peak.loc[n] - peak.loc[1])
  } else 0

  # perform trapezoidal integration on raw signal
  int <- integrate_above_baseline(yp, peak.loc, floor_val, slope_val)
  area <- int$area
  if (!is.na(area) && area < 0){
    warning(sprintf(
      "The local baseline lies above the signal in the window from %s to %s, giving a negative area.",
      signif(get_times(x)[pos[2]], 4), signif(get_times(x)[pos[3]], 4)),
      immediate. = TRUE)
  }
  purity <- peak_purity(x, pos, estimate_purity, noise_variance, purity_args)
  result <- c("rt" = pos[1], "start" = pos[2], "end" = pos[3],
    "sd" = NA_real_, "FWHM" = measured_fwhm(int$y_above),
    "height" = y[xloc] - (floor_val + slope_val * (xloc - peak.loc[1])),
    "area" = area, purity = purity)
  append_floor_cols(result, baseline, floor_val, slope_val, peak.loc)
}


#' Savitsky Golay Smoothing ported from pracma
#' @author Hans W. Borchers
#' @param T Vector of signals to be filtered
#' @param fl Filter length (for instance fl = 51..151), has to be odd.
#' @param forder filter order Filter order (2 = quadratic filter, 4 = quartic).
#' @param dorder Derivative order (0 = smoothing, 1 = first derivative, etc.).
#' @note This function is bundled from [pracma](
#' https://cran.r-project.org/web/packages/pracma/index.html), where it is 
#' licensed under GPL (>= 3).
#' @importFrom stats convolve
#' @noRd
savgol <- function(T, fl, forder = 4, dorder = 0) {
  stopifnot(is.numeric(T), is.numeric(fl))
  if (fl <= 1 || fl %% 2 == 0)
    stop("Argument 'fl' must be an odd integer greater than 1.")
  n <- length(T)
  
  # -- calculate filter coefficients --
  fc <- (fl-1)/2                          # index: window left and right
  X <- outer(-fc:fc, 0:forder, FUN = "^")   # polynomial terms and coeffs
  Y <- pinv(X);                           # pseudoinverse
  
  # -- filter via convolution and take care of the end points --
  T2 <- convolve(T, rev(Y[(dorder + 1),]), type = "o")   # convolve(...)
  T2 <- T2[(fc+1):(length(T2)-fc)]
  
  Tsg <- (-1)^dorder * T2
  return( Tsg )
}

#' 'pinv' function from pracma
#' @author Hans W. Borchers
#' @note This function is ported from [pracma](
#' https://cran.r-project.org/web/packages/pracma/index.html), where it is 
#' licensed under GPL (>= 3).
#' @noRd
pinv <- function (A, tol = .Machine$double.eps^(2/3)) {
  stopifnot(is.numeric(A) || is.complex(A), is.matrix(A))
  
  s <- svd(A)

  if (is.complex(A)) s$u <- Conj(s$u)
  
  p <- ( s$d > max(tol * s$d[1], 0) )
  if (all(p)) {
    mp <- s$v %*% (1/s$d * t(s$u))
  } else if (any(p)) {
    mp <- s$v[, p, drop = FALSE] %*% (1/s$d[p] * t(s$u[, p, drop = FALSE]))
  } else {
    mp <- matrix(0, nrow = ncol(A), ncol = nrow(A))
  }
  
  return(mp)
}

#' Bi-Exponentially Modified Gaussian
#' @noRd
bemg <- function(x, center, width, height, a, b, floor = 0, slope = 0) {
  tm <- x - center
  as2 <- a * width^2
  bs2 <- b * width^2
  sq2s2 <- sqrt(2 * width^2)

  tau_left  <- exp(a * as2 / 2 + a * tm + log_erfc((as2 + tm) / sq2s2))
  tau_right <- exp(b * bs2 / 2 - b * tm + log_erfc((bs2 - tm) / sq2s2))

  return(floor + slope * (x - x[1]) + height * 0.5 * (tau_left + tau_right))
}

#' Fit Bi-Exponentially Modified Gaussian peak
#' @importFrom stats coef fitted quantile residuals
#' @noRd
fit_bemg <- function(x1, y1, start.center = NULL, start.width = NULL,
                     start.a = NULL, start.b = NULL, start.height = NULL,
                     start.floor = NULL, baseline = "none", max_iter = 1000) {

  if (is.null(start.a))      start.a      <- 1.0
  if (is.null(start.b))      start.b      <- 1.0

  st <- fit_starts(x1, y1, start.center, start.width, start.height,
                   height_mult = 5)
  start.center <- st$center
  start.width <- st$width
  start.height <- st$height
  peak_width <- st$span
  shape_lower <- c(st$lower, a = -Inf, b = -Inf)
  shape_upper <- c(st$upper, a = Inf, b = Inf)

  controlList <- minpack.lm::nls.lm.control(maxiter = max_iter,
                                           maxfev = 10 * max_iter)
  starts <- list(center = start.center, width = start.width,
                 height = start.height, a = start.a, b = start.b)

  # Stage 1: shape only
  nlsAns <- skip_underdetermined(length(x1), 5) %||%
    try(minpack.lm::nlsLM(
    y1 ~ bemg(x1, center, width, height, a, b),
    start = starts, control = controlList,
    lower = shape_lower, upper = shape_upper), silent = TRUE)
  nlsAns <- reject_growing_tail(nlsAns)
  nls_stage1 <- nlsAns

  upper_floor <- max(min(y1), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, floor = start.floor)
    nlsAns <- skip_underdetermined(length(x1), 6) %||%
      try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b, floor),
      start = starts, control = controlList,
      lower = c(shape_lower, floor = 0),
      upper = c(shape_upper, floor = upper_floor)
    ), silent = TRUE)
  }
  nlsAns <- reject_growing_tail(nlsAns)
  nls_stage2 <- nlsAns

  # Stage 3: slope via (floor_s, floor_e) so both endpoints are bounded >= 0
  if (baseline == "sloped") {
    f0 <- if (!inherits(nlsAns, "try-error")) unname(coef(nlsAns)["floor"])
          else if (!is.null(start.floor)) start.floor else 0
    s3 <- c(as.list(if (!inherits(nlsAns, "try-error"))
                      coef(nlsAns)[c("center", "width", "height", "a", "b")]
                    else c(center = start.center, width = start.width,
                           height = start.height, a = start.a, b = start.b)),
            floor_s = f0, floor_e = f0)
    nlsAns <- skip_underdetermined(length(x1), 7) %||%
      try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(shape_lower, floor_s = 0, floor_e = 0),
      upper = c(shape_upper, floor_s = Inf, floor_e = Inf)
    ), silent = TRUE)
  }

  nlsAns <- reject_growing_tail(nlsAns)
  if (inherits(nlsAns, "try-error")) {
    if (baseline == "sloped" && !inherits(nls_stage2, "try-error")) {
      fb <- coef(nls_stage2)
      floorAns <- unname(fb["floor"]); slopeAns <- 0
      yAns <- bemg(x1, fb["center"], fb["width"], fb["height"],
                   fb["a"], fb["b"], floorAns, slopeAns)
      out <- list(center = fb["center"], width = fb["width"],
                  height = fb["height"], a = fb["a"], b = fb["b"],
                  y = yAns, residual = y1 - yAns)
    } else if (!inherits(nls_stage1, "try-error")) {
      fb <- coef(nls_stage1)
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- bemg(x1, fb["center"], fb["width"], fb["height"],
                   fb["a"], fb["b"], floorAns, slopeAns)
      out <- list(center = fb["center"], width = fb["width"],
                  height = fb["height"], a = fb["a"], b = fb["b"],
                  y = yAns, residual = y1 - yAns)
    } else {
      floorAns <- if (baseline != "none") start.floor else 0; slopeAns <- 0
      yAns <- bemg(x1, start.center, start.width, start.height,
                   start.a, start.b, floorAns, slopeAns)
      out <- list(center = start.center, width = start.width,
                  height = start.height, a = start.a, b = start.b,
                  y = yAns, residual = y1 - yAns)
    }
  } else {
    coefs <- coef(nlsAns)
    out <- list(center = coefs["center"], width = coefs["width"],
                height = coefs["height"], a = coefs["a"], b = coefs["b"],
                y = fitted(nlsAns), residual = residuals(nlsAns))
    if (baseline == "sloped") {
      floorAns <- coefs["floor_s"]
      slopeAns  <- (coefs["floor_e"] - coefs["floor_s"]) / peak_width
    } else {
      floorAns <- if (baseline != "none") coefs["floor"] else 0; slopeAns <- 0
    }
  }

  out$window_limited <- is_window_limited(out, shape_lower, shape_upper)
  out$converged <- !inherits(nlsAns, "try-error") &&
    isTRUE(nlsAns$convInfo$isConv)
  if (baseline != "none")    out <- c(out, floor = unname(floorAns))
  if (baseline == "sloped")  out <- c(out, slope = unname(slopeAns))
  return(out)
}

#' @noRd
erfc <- function(x) 2 * pnorm(-x * sqrt(2))

#' Log of the complementary error function
#'
#' Returns `log(erfc(x))` without forming `erfc(x)` first, so it stays accurate
#' in the far upper tail where `erfc` underflows to zero. Callers that need
#' `exp(A) * erfc(x)` should compute `exp(A + log_erfc(x))`: the product is
#' bounded even when `exp(A)` alone overflows to `Inf`, which is what makes the
#' exponentially modified Gaussians numerically stable for wide peaks and large
#' tailing parameters.
#'
#' @param x A numeric vector.
#' @return A numeric vector of the same length as `x`.
#' @author Ethan Bass
#' @noRd
log_erfc <- function(x) log(2) + pnorm(-x * sqrt(2), log.p = TRUE)
