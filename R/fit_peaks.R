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
#' @param lambda Wavelength at which to fit peaks. Must correspond to a wavelength present in `x`.
#' @param pos Locations of peaks in vector `y`. If `NULL`, `find_peaks` will run
#' automatically to find peak positions.
#' @param sd_max Maximum width (standard deviation) for peaks. Defaults to `50`.
#' @param fit Peak model to use. Options include `"bemg"` (bidirectional 
#' exponentially modified Gaussian), `"egh"` (exponential-Gaussian hybrid),
#' `"gaussian"`, and `"raw"`. If `raw` is selected, peaks are integrated using
#' trapezoidal integration without model fitting. Defaults to `"egh"`.
#' @param max_iter Maximum number of iterations to use in nonlinear least
#' squares peak-fitting. Defaults to `1000`.
#' @param estimate_purity Logical. Whether to estimate purity or not. Defaults
#' to `TRUE`.
#' @param noise_threshold Noise threshold. Input to `get_purity`.
#' @param baseline Whether to fit a baseline offset beneath each peak model.
#' Options are `"none"` (no baseline, default), `"flat"` (constant offset), or
#' `"sloped"` (linearly varying offset). Not available when `fit = "raw"`.
#' @param ... Additional arguments to `find_peaks`.
#' @return A matrix with one row per peak and the following columns:
#' * `rt`: Peak maximum location.
#' * `start`: Peak start (only included in table if `bounds = TRUE`).
#' * `end`: Peak end (only included in table if `bounds = TRUE`).
#' * `sd`: Peak standard deviation.
#' * `tau`: Exponential decay constant (only included in table if
#' `fit = "egh"`).
#' * `tau_right`: Exponential decay constant for right side of peak (when
#' `bemg` is selected).
#' * `tau_left`: Exponential decay constant for left side of peak (when
#' `bemg` is selected).
#' * `FWHM`: The full width at half maximum.
#' * `height`: Peak height above the fitted baseline.
#' * `area`: Peak area above the fitted baseline.
#' * `r.squared`: The R<sup>2</sup> value for linear fit of the model to the data.
#' * `purity`: The spectral purity of peak as assessed by [`get_purity`].
#' * `floor`: Constant baseline offset (only when `baseline = "flat"`).
#' * `floor_start`: Baseline level at the left edge of the peak window (only
#' when `baseline = "sloped"`).
#' * `floor_end`: Baseline level at the right edge of the peak window (only
#' when `baseline = "sloped"`).
#' The first five elements (`rt`, `start`, `end`, `sd` and `FWHM`) are expressed
#' as index positions rather than absolute retention times. The transformation
#' to real time is done in `get_peaks`.
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

fit_peaks <- function (x, lambda, pos = NULL, sd_max = 50,
                       fit = c("egh", "bemg", "gaussian", "raw"),
                       max_iter = 1000, estimate_purity = TRUE,
                       noise_threshold = .001,
                       baseline = c("none", "flat", "sloped"), ...){
  lambda.idx <- get_lambda_idx(lambda, get_lambdas(x))
  y <- x[,lambda.idx]
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  baseline <- match.arg(baseline, c("none", "flat", "sloped"))
  if (is.null(pos)){
    pos <- find_peaks(y = y, ...)
  }
  if (ncol(x) == 1){
    estimate_purity <- FALSE
  }
  tabnames <- switch(fit,
                     "gaussian" = c("rt", "start", "end", "sd", "FWHM",
                                    "height", "area", "r-squared", "purity"),
                     "egh" = c("rt", "start", "end", "sd", "tau", "FWHM",
                               "height", "area", "r.squared", "purity"),
                     "bemg" = c("rt", "start", "end", "sd", "tau_right",
                                "tau_left", "h", "center", "FWHM", "height",
                                "area", "r.squared", "purity"),
                     "raw" = c("rt", "start", "end", "sd", "FWHM", "height",
                               "area", "purity")
  )
  if (baseline == "flat" && fit != "raw")    tabnames <- c(tabnames, "floor")
  if (baseline == "sloped" && fit != "raw") tabnames <- c(tabnames, "floor_start", "floor_end")
  noPeaksMat <- matrix(rep(NA, length(tabnames)),
                       nrow = 1, dimnames = list(NULL, tabnames))
  on.edge <- sapply(pos$pos, function(x){
    x <= 1 || is.na(y[x + 1]) || is.na(y[x - 1])
  })
  pos <- pos[!on.edge,]

  if (nrow(pos) == 0)
    return(noPeaksMat)
  fitpk <- switch(fit,
                  "gaussian" = fitpk_gaussian,
                  "egh" = fitpk_egh,
                  "raw" = fitpk_raw,
                  "bemg" = fitpk_bemg)
  huhn <- data.frame(t(apply(pos, 1, fitpk, x = x,
                             lambda = lambda.idx, max_iter = max_iter,
                             estimate_purity = estimate_purity,
                             noise_threshold = noise_threshold,
                             baseline = baseline)))
  colnames(huhn) <- tabnames
  huhn <- data.frame(sapply(huhn, as.numeric, simplify = FALSE))
  if (!is.null(sd_max)) {
    huhn <- huhn[huhn$sd < sd_max, ]
  }
  df <- try(huhn[huhn$rt > 0,], silent = TRUE)
  if(inherits(df, "try-error")) NA else df
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
#' @importFrom stats coef fitted lm nls.control quantile residuals
#' @noRd
fit_gaussian <- function(x, y, start.center = NULL,
                         start.width = NULL, start.height = NULL,
                         start.floor = NULL, baseline = "none",
                         max_iter = 1000){
  who.max <- which.max(y)
  if (is.null(start.center)) start.center <- x[who.max]
  if (is.null(start.height)) start.height <- y[who.max]
  if (is.null(start.width))  start.width  <- sum(y > (start.height / 2)) / 2

  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512, warnOnly = TRUE)
  starts <- list("center" = start.center, "width" = start.width, "height" = start.height)

  # Stage 1: shape only
  nlsAns <- try(minpack.lm::nlsLM(y ~ gaussian(x, center, width, height),
                       start = starts, control = controlList), silent = TRUE)
  nls_stage1 <- nlsAns

  peak_width <- x[length(x)] - x[1]
  if (peak_width == 0) peak_width <- 1
  upper_floor <- max(min(y), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- try(minpack.lm::nlsLM(
      y ~ gaussian(x, center, width, height, floor),
      start = starts, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, floor = 0),
      upper = c(center = Inf, width = Inf, height = Inf, floor = upper_floor)
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
    nlsAns <- try(minpack.lm::nlsLM(
      y ~ gaussian(x, center, width, height, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, floor_s = 0, floor_e = 0)
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
#' @importFrom stats coef fitted lm nls.control quantile residuals
#' @noRd
fit_egh <- function(x1, y1, start.center = NULL, start.width = NULL,
                    start.tau = NULL, start.height = NULL, start.floor = NULL,
                    baseline = "none", max_iter = 1000) {

  who.max <- which.max(y1)
  if (is.null(start.center)) start.center <- x1[who.max]
  if (is.null(start.height)) start.height <- y1[who.max]
  if (is.null(start.width))  start.width  <- sum(y1 > (start.height / 2)) / 2
  if (is.null(start.tau))    start.tau    <- 0

  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512, warnOnly = TRUE)
  starts <- list("center" = start.center, "width" = start.width,
                 "height" = start.height, "tau" = start.tau)

  # Stage 1: shape only
  nlsAns <- try(minpack.lm::nlsLM(y1 ~ egh(x1, center, width, height, tau),
                      start = starts, control = controlList), silent = TRUE)
  nls_stage1 <- nlsAns

  peak_width <- x1[length(x1)] - x1[1]
  if (peak_width == 0) peak_width <- 1
  upper_floor <- max(min(y1), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ egh(x1, center, width, height, tau, floor),
      start = starts, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, tau = -Inf, floor = 0),
      upper = c(center = Inf, width = Inf, height = Inf, tau = Inf, floor = upper_floor)
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
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ egh(x1, center, width, height, tau, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, tau = -Inf, floor_s = 0, floor_e = 0)
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

  if (baseline != "none")    out <- c(out, "floor" = unname(floorAns))
  if (baseline == "sloped")  out <- c(out, "slope" = unname(slopeAns))
  return(out)
}

fitpk_bemg <- function(x, pos, lambda, max_iter,
                       estimate_purity = TRUE, noise_threshold = .001,
                       baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_bemg(peak.loc, y[peak.loc], start.center = xloc,
                                 start.height = y[[xloc]], max_iter = max_iter,
                                 baseline = baseline)
  )
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  floor_val <- if (baseline != "none") unname(m$floor) else 0
  slope_val <- if (baseline == "sloped") unname(m$slope) else 0
  baseline_vec <- floor_val + slope_val * (peak.loc - peak.loc[1])
  y_above <- m$y - baseline_vec
  area <- sum(diff(peak.loc) * (y_above[-length(y_above)] + y_above[-1]) / 2)
  apex_idx <- which.max(y_above)
  rt_val <- if (length(apex_idx) == 1L) peak.loc[apex_idx] else unname(round(m$center))
  result <- c("rt" = rt_val, "start" = pos[[2]], "end" = pos[[3]],
    "sd" = unname(m$width), "tau_right" = unname(m$a), "tau_left" = unname(m$b),
    "h" = unname(m$height), "center" = unname(round(m$center)),
    "FWHM" = 2.35 * unname(m$width),
    "height" = unname(max(y_above)),
    "area" = area, "r.squared" = unname(r.squared),
    purity = purity)
  if (baseline == "flat")    result <- c(result, floor = floor_val)
  if (baseline == "sloped") result <- c(result,
    floor_start = floor_val,
    floor_end   = floor_val + slope_val * (tail(peak.loc, 1) - peak.loc[1]))
  result
}

#' Fit peak (gaussian)
#' @noRd
fitpk_gaussian <- function(x, pos, lambda, max_iter,
                           estimate_purity = TRUE, noise_threshold = .001,
                           baseline = "none", ...){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_gaussian(peak.loc, y[peak.loc],
                                     start.center = xloc,
                                     start.height = y[[xloc]],
                                     max_iter = max_iter,
                                     baseline = baseline)
  )
  floor_val <- if (baseline != "none") unname(m$floor) else 0
  slope_val <- if (baseline == "sloped") unname(m$slope) else 0
  baseline_vec <- floor_val + slope_val * (peak.loc - peak.loc[1])
  y_above <- m$y - baseline_vec
  area <- sum(diff(peak.loc) * (y_above[-length(y_above)] + y_above[-1]) / 2)
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  result <- c("rt" = unname(m$center), "start" = pos[2], "end" = pos[3],
    "sd" = unname(m$width), "FWHM" = 2.35 * unname(m$width),
    "height" = unname(m$height), "area" = area,
    "r.squared" = unname(r.squared), purity = purity)
  if (baseline == "flat")    result <- c(result, floor = floor_val)
  if (baseline == "sloped") result <- c(result,
    floor_start = floor_val,
    floor_end   = floor_val + slope_val * (tail(peak.loc, 1) - peak.loc[1]))
  result
}

#' Fit peak (exponential-gaussian hybrid)
#' @noRd
fitpk_egh <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_threshold = .001,
                      baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc,
                                start.height = y[[xloc]], max_iter = max_iter,
                                baseline = baseline)
  )
  r.squared <- suppressWarnings(try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE))
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  floor_val <- if (baseline != "none") unname(m$floor) else 0
  slope_val <- if (baseline == "sloped") unname(m$slope) else 0
  baseline_vec <- floor_val + slope_val * (peak.loc - peak.loc[1])
  y_above <- m$y - baseline_vec
  area <- sum(diff(peak.loc) * (y_above[-length(y_above)] + y_above[-1]) / 2)
  result <- c("rt" = unname(m$center), "start" = pos[2], "end" = pos[3],
    "sd" = unname(m$width), "tau" = unname(m$tau), "FWHM" = 2.35 * unname(m$width),
    "height" = unname(m$height), "area" = area,
    "r.squared" = unname(r.squared), purity = purity)
  if (baseline == "flat")    result <- c(result, floor = floor_val)
  if (baseline == "sloped") result <- c(result,
    floor_start = floor_val,
    floor_end   = floor_val + slope_val * (tail(peak.loc, 1) - peak.loc[1]))
  result
}

#' Fit peak (raw)
#' @noRd
fitpk_raw <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_threshold = .001,
                      baseline = "none"){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  
  # perform trapezoidal integration on raw signal
  yp <- y[peak.loc]
  area <- sum(diff(peak.loc) * (yp[-length(yp)] + yp[-1]) / 2)
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  c("rt" = pos[1], "start" = pos[2], "end" = pos[3], 
    "sd" = pos[3] - pos[2], "FWHM" = 2.35 * (pos[3] - pos[2]),
    "height" = y[xloc], "area" = area, purity = purity)
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

  tau_right <- exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / sq2s2)
  tau_left  <- exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / sq2s2)

  return(floor + slope * (x - x[1]) + height * 0.5 * (tau_right + tau_left))
}

#' Fit Bi-Exponentially Modified Gaussian peak
#' @importFrom stats coef fitted nls.control quantile residuals
#' @noRd
fit_bemg <- function(x1, y1, start.center = NULL, start.width = NULL,
                     start.a = NULL, start.b = NULL, start.height = NULL,
                     start.floor = NULL, baseline = "none", max_iter = 1000) {

  who.max <- which.max(y1)
  if (is.null(start.center)) start.center <- x1[who.max]
  if (is.null(start.height)) start.height <- y1[who.max]
  if (is.null(start.width))  start.width  <- sum(y1 > (start.height / 2)) / 2
  if (is.null(start.a))      start.a      <- 1.0
  if (is.null(start.b))      start.b      <- 1.0

  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512, warnOnly = TRUE)
  starts <- list(center = start.center, width = start.width,
                 height = start.height, a = start.a, b = start.b)

  # Stage 1: shape only
  upper <- c(center = Inf, width = Inf, height = 5 * start.height, a = Inf, b = Inf)
  nlsAns <- try(minpack.lm::nlsLM(
    y1 ~ bemg(x1, center, width, height, a, b),
    start = starts, control = controlList, upper = upper), silent = TRUE)
  nls_stage1 <- nlsAns

  peak_width <- x1[length(x1)] - x1[1]
  if (peak_width == 0) peak_width <- 1
  upper_floor <- max(min(y1), 0)

  # Stage 2: add constant floor
  if (baseline != "none") {
    if (!inherits(nlsAns, "try-error")) starts <- as.list(coef(nlsAns))
    if (is.null(start.floor)) start.floor <- 0
    starts <- c(starts, floor = start.floor)
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b, floor),
      start = starts, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, a = -Inf, b = -Inf, floor = 0),
      upper = c(center = Inf, width = Inf, height = 5 * start.height,
                a = Inf, b = Inf, floor = upper_floor)
    ), silent = TRUE)
  }
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
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b, floor_s, (floor_e - floor_s) / peak_width),
      start = s3, control = controlList,
      lower = c(center = -Inf, width = -Inf, height = -Inf, a = -Inf, b = -Inf,
                floor_s = 0, floor_e = 0),
      upper = c(center = Inf, width = Inf, height = 5 * start.height,
                a = Inf, b = Inf, floor_s = Inf, floor_e = Inf)
    ), silent = TRUE)
  }

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

  if (baseline != "none")    out <- c(out, floor = unname(floorAns))
  if (baseline == "sloped")  out <- c(out, slope = unname(slopeAns))
  return(out)
}

#' @noRd
erfc <- function(x) 2 * stats::pnorm(-x * sqrt(2))
