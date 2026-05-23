#' Find peaks
#' 
#' Detects peaks in chromatographic profile.
#' 
#' Detects peaks by looking for zero-crossings in the smoothed first derivative of
#' the signal (`y`) that exceed the specified slope threshold (`slope_thresh`).
#' Additionally, peaks can be filtered by supplying a minimal
#' amplitude threshold (`amp_thresh`), filtering out peaks below the
#' specified height. Smoothing is intended to prevent the algorithm from
#' getting caught up on local minima and maxima that do not represent true
#' features. Several smoothing options are available, including `"gaussian"`,
#' box kernel (`"box"`), savitzky-golay smoothing (`"savgol"`),
#' moving average (`"mva"`), triangular moving average (`"tmva"`), or 
#' no smoothing (`"none"`).
#' 
#' It is recommended to do pre-processing using the [`preprocess`]
#' function before peak detection. Overly high chromatographic resolution can 
#' sometimes cause peaks to be split into multiple segments. In this case,
#' it is recommended to increase the `smooth_window` or reduce the
#' resolution along the time axis by adjusting the `dim1` argument during
#' preprocessing.
#' 
#' @importFrom caTools runmean
#' @importFrom minpack.lm nlsLM
#' @importFrom stats deriv lm ksmooth
#' @importFrom utils tail
#' @param y A chromatographic signal (as a numeric vector).
#' @param smooth_type Type of smoothing. Either gaussian kernel (`"gaussian"`),
#' box kernel (`"box"`), Savitzky-Golay smoothing (`"savgol"`),
#' moving average (`"mva"`), triangular moving average (`"tmva"`), or 
#' no smoothing (`"none"`).
#' @param smooth_window Smoothing window. Larger values of this parameter will 
#' exclude sharp, narrow features. If the supplied value is between 0 and
#' 1, the window will be interpreted as a proportion of points to include. 
#' Otherwise, the window will be interpreted as the absolute number of points to
#' include in the window. Defaults to `0.001`.
#' @param slope_thresh Minimum threshold for slope of the smoothed first
#' derivative. This parameter filters peaks on the basis of their width, such
#' that larger values will exclude broad peaks from the peak list. Defaults to
#' `0`.
#' @param amp_thresh Minimum threshold for peak amplitude. This parameter
#' filters on the basis of peak height, such that larger values will exclude 
#' small peaks from the peak list. Defaults to `0`.
#' @param bounds Logical. If `TRUE` (default), includes peak boundaries in 
#' data.frame.
#' @return If `bounds == TRUE`, returns a `data.frame` containing the center,
#' start, and end of each identified peak. Otherwise, returns a numeric vector
#' of peak centers. All locations are expressed as indices.
#' @note The `find_peaks` function is adapted from MATLAB code included in
#' Prof. Tom O'Haver's [Pragmatic Introduction to Signal Processing](
#' http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm).
#' @keywords internal
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' find_peaks(Sa_pr[[1]][,"220"])
#' @seealso [`fit_peaks`], [`get_peaks`]
#' @references O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement. <https://terpconnect.umd.edu/~toh/spectrum/>
#' (Accessed January, 2022).
#' @export find_peaks

find_peaks <- function(y, smooth_type = c("gaussian", "box", "savgol", "mva", 
                                          "tmva", "none"), smooth_window = .001,
                       slope_thresh = 0, amp_thresh = 0,
                       bounds = TRUE){
  if (!is.vector(y)){
    stop("Please provide a vector to argument `y` to proceed.")
  }
  smooth_type <- match.arg(smooth_type, c("gaussian", "box", "savgol", "mva",
                                          "tmva", "none"))
  if (smooth_window < 1){
    smooth_window <- max(length(y) * smooth_window, 1)
  }
  # compute derivative (with or without smoothing)
  if (smooth_type == "savgol"){
    if ((smooth_window %% 2) == 0){
      smooth_window <- smooth_window + 1
    }
    d <- savgol(diff(y), fl = smooth_window)
  } else if (smooth_type == "mva"){
    d <- caTools::runmean(diff(y), k = smooth_window)
  } else if (smooth_type == 'gaussian'){
    d <- diff(ksmooth(seq_along(y), y, kernel = "normal",
                      bandwidth = smooth_window)$y)
  }  else if (smooth_type == "box"){
    d <- diff(ksmooth(seq_along(y), y, kernel = "box",
                      bandwidth = smooth_window)$y)
  } else if (smooth_type == "tmva"){
    d <- runmean(runmean(diff(y), k = smooth_window), k = smooth_window)
  } else{
    d <- diff(y)
  }
  
  # detect zero-crossing of first derivative (peak apex)
  p1 <- which(sign(d[1:(length(d)-1)]) > sign(d[2:length(d)]))
  
  # detect second derivative exceeding slope threshold
  p2 <- which(abs(diff(d)) > slope_thresh)
  
  # detect y-vals exceeding amplitude threshold
  p3 <- which(y > amp_thresh) 
  p <- intersect(intersect(p1,p2), p3)
  if (bounds){
    p4 <- which(sign(d[1:(length(d)-1)]) < sign(d[2:length(d)]))
    
    # find lower bound
    suppressWarnings(bl <- sapply(p, function(v) max(p4[p4 < v])))
    bl[which(bl == -Inf)] <- 1
    
    # find upper bound
    suppressWarnings(bu <- sapply(p, function(v) min(p4[p4 > v])))
    bu[which(bu == Inf)] <- length(y)
    p <- data.frame(pos = p, lower = bl, upper = bu)
  }
  p
}

#' Fit chromatographic peaks
#' 
#' Fit peak parameters using the specified peak model
#' 
#' Peak parameters are calculated by fitting the data to a gaussian, 
#' bidirectional exponentially modified gaussian, or exponential-gaussian hybrid
#' curve using non-linear least squares estimation as implemented in 
#' [`nlsLM`][minpack.lm::nlsLM]. The area under the fitted curve is then 
#' estimated using trapezoidal approximation.
#' 
#' @param x A single chromatogram in matrix format.
#' @param lambda Wavelength to fit peaks at.
#' @param pos Locations of peaks in vector `y`. If `NULL`, `find_peaks` will run
#' automatically to find peak positions.
#' @param sd_max Maximum width (standard deviation) for peaks. Defaults to `50`.
#' @param fit Function for peak fitting. Currently, bidirectional exponentially 
#' modified Gaussian (`"bemg"`), exponential-gaussian hybrid (`"egh"`),
#' `"gaussian"` and `"raw"` settings are supported. If `raw` mode is
#' selected, trapezoidal integration will be performed on raw data
#' without fitting a model to the peak. Defaults to `"bemg"`.)
#' @param max_iter Maximum number of iterations to use in nonlinear least
#' squares peak-fitting. Defaults to `1000`.
#' @param estimate_purity Logical. Whether to estimate purity or not. Defaults
#' to `TRUE`.
#' @param noise_threshold Noise threshold. Input to `get_purity`.
#' @param ... Additional arguments to `find_peaks`.
#' @return Returns a matrix with columns containing the following information 
#' about each peak:
#' * `rt`: Location of the peak maximum.
#' * `start`: Start of peak (only included in table if `bounds = TRUE`).
#' * `end`: End of peak (only included in table if `bounds = TRUE`).
#' * `sd`: The standard deviation of the peak.
#' * `tau`: Exponential decay constant (only included in table if 
#' `fit = "egh"`).
#' * `tau_right`: Exponential decay constant for right side of peak (when 
#' `bemg` is selected).
#' * `tau_left`: Exponential decay constant for left side of peak (when 
#' `bemg` is selected).
#' * `FWHM`: The full width at half maximum.
#' * `height`: Peak height.
#' * `area`: Peak area.
#' * `r.squared`: The R<sup>2</sup> value for linear fit of the model to the data.
#' * `purity`: The spectral purity of peak as assessed by [`get_purity`].
#' Again, the first five elements (rt, start, end, sd and FWHM) are expressed
#' as indices, so not in terms of the real retention times. The transformation
#' to "real" time is done in function the `get_peaks` function.
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
#' @keywords internal
#' @md

fit_peaks <- function (x, lambda, pos = NULL, sd_max = 50,
                       fit = c("bemg", "egh", "gaussian", "raw"), 
                       max_iter = 1000, estimate_purity = TRUE, 
                       noise_threshold = .001, ...){
  lambda.idx <- get_lambda_idx(lambda, get_lambdas(x))
  y <- x[,lambda.idx]
  fit <- match.arg(fit, c("bemg", "egh", "gaussian", "raw"))
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
                                "tau_left", "h", "FWHM", "height", "area", 
                                "r.squared", "purity"),
                     "raw" = c("rt", "start", "end", "sd", "FWHM", "height",
                                                      "area", "purity")
  )
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
                             noise_threshold = noise_threshold)))
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
gaussian <- function(x, center = 0, width = 1, height = NULL, floor = 0){
  
  # adapted from Earl F. Glynn;  Stowers Institute for Medical Research, 2007
  twoVar <- 2 * width * width
  sqrt2piVar <- sqrt( pi * twoVar)
  y <- exp( -( x - center)^2 / twoVar) / sqrt2piVar
  
  # by default, the height is such that the curve has unit volume
  if ( ! is.null (height)) {
    scalefactor <- sqrt2piVar
    y <- y * scalefactor * height
  }
  y + floor
}

#' Fit gaussian peak
#' @importFrom stats coef fitted lm nls.control quantile residuals
#' @noRd
fit_gaussian <- function(x, y, start.center = NULL,
                         start.width = NULL, start.height = NULL,
                         start.floor = NULL, fit.floor = FALSE,
                         max_iter = 1000){
  # estimate starting values
  who.max <- which.max(y)
  if (is.null(start.center)) start.center <- x[who.max]
  if (is.null(start.height)) start.height <- y[who.max]
  if (is.null(start.width)) start.width <- sum( y > (start.height/2)) / 2
  
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512,
                              warnOnly = TRUE)
  starts <- list("center" = start.center, "width" = start.width,
                  "height" = start.height)
  if (!fit.floor) {
    nlsAns <- try(nlsLM(y ~ gaussian(x, center, width, height),
                         start = starts, control = controlList), silent = TRUE)
  } else{
    if (is.null(start.floor)) start.floor <- quantile(y, seq(0, 1, 0.1))[2]
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- try(nlsLM( y ~ gaussian(x, center, width, height, floor),
                         start = starts, control = controlList), silent = TRUE)
  }
  # package up the results to pass back
    if (inherits(nlsAns, "try-error")){
      yAns <- gaussian(x, start.center, start.width, start.height, start.floor)
      out <- list("center" = start.center, "width" = start.width,
                  "height" = start.height,
                  "y" = yAns, "residual" = y - yAns)
      floorAns <- if (fit.floor) start.floor else 0
    } else {
      coefs <-coef(nlsAns)
      out <- list( "center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                   "y" = fitted(nlsAns), "residual" = residuals(nlsAns))
      floorAns <- if (fit.floor) coefs[4] else 0
    }
    if (fit.floor) {
      out <- c( out, "floor" = floorAns)
    }
  return( out)
}

###############################################################################

#' Expontential-gaussian hybrid
#' @noRd
egh <- function(x, center, width,  height, tau, floor = 0){
  result <- rep(0, length(x))
  index <- which(2*width^2 + tau*(x-center)>0)
  result[index] <- height*exp(-(x[index] - center)^2/(2*width^2 + tau*(x[index] - center)))
  return(result)
}


#' Fit exponential-gaussian hybrid peak
#' @importFrom stats coef fitted lm nls.control quantile residuals
#' @noRd
fit_egh <- function(x1, y1, start.center = NULL, start.width = NULL,
                    start.tau = NULL, start.height = NULL, start.floor = NULL,
                    fit.floor = FALSE, max_iter = 1000) {
  
  # try to find the best egh to fit the given data
  
  # make some rough estimates from the values of Y
  who.max <- which.max(y1)
  if (is.null(start.center)){
    start.center <- x1[who.max]
  }
  if (is.null(start.height)){
    start.height <- y1[who.max]
  }
  if (is.null(start.width)){
    start.width <- sum(y1 > (start.height/2)) / 2
  }
  if (is.null(start.tau)){
    start.tau <- 0
  }
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512,
                             warnOnly = TRUE)
  starts <- list("center" = start.center, "width" = start.width, 
                 "height" = start.height, "tau" = start.tau)
  if (!fit.floor){
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau),
                        start = starts, control = controlList), silent = TRUE)
  } else{
    if (is.null( start.floor)) start.floor <- quantile( y1, seq(0, 1, 0.1))[2]
    starts <- c(starts, "floor" = start.floor)
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau, floor),
                        start = starts, control = controlList), silent = TRUE)
  }
  
  # package up the results to pass back
  if (inherits(nlsAns, "try-error")) {
    yAns <- egh(x1, start.center, start.width, start.height,
                start.tau, start.floor)
    out <- list("center" = start.center, "width" = start.width,
                "height" = start.height, "tau" = start.tau,
                "y" = yAns, "residual" = y1 - yAns)
    floorAns <- if (fit.floor) start.floor else 0
  } else {
    coefs <-coef(nlsAns)
    out <- list( "center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                 "tau" = coefs[4], "y" = fitted(nlsAns),
                 "residual" = residuals(nlsAns))
    floorAns <- if (fit.floor) coefs[5] else 0
  }
  
  if (fit.floor) {
    out <- c( out, "floor" = floorAns)
  }
  return(out)
}

fitpk_bemg <- function(x, pos, lambda, max_iter,
                       estimate_purity = TRUE, noise_threshold = .001){
  y <- x[,lambda]
  xloc <- pos[[1]]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_bemg(peak.loc, y[peak.loc], start.center = xloc,
                                 start.height = y[[xloc]], max_iter = max_iter)
  )
  r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE)
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  # trapezoidal integration
  area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y, -1))))
  c("rt" = unname(round(m$center)), "start" = pos[[2]], "end" = pos[[3]], 
    "sd" = unname(m$width), "tau_right" = unname(m$a), "tau_left" = unname(m$b), 
    "h" = unname(m$height), "FWHM" = 2.35 * unname(m$width),
    "height" = unname(max(m$y)), "area" = area, "r.squared" = unname(r.squared), 
    purity = purity)
}

#' Fit peak (gaussian)
#' @noRd
fitpk_gaussian <- function(x, pos, lambda, max_iter,
                           estimate_purity = TRUE, noise_threshold = .001, ...){
  
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_gaussian(peak.loc, y[peak.loc],
                                     start.center = xloc,
                                     start.height = y[xloc],
                                     max_iter = max_iter)
  )
  area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1))))
  r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE)
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
    "sd" = m$width, "FWHM" = 2.35 * m$width,
    "height" = m$height, "area" = area, "r.squared" = r.squared, purity = purity)
}

#' Fit peak (exponential-gaussian hybrid)
#' @noRd
fitpk_egh <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_threshold = .001){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc,
                                start.height = y[xloc], max_iter = max_iter)
  )
  r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE)
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  # trapezoidal integration
  area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y, -1))))
  c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
    "sd" = m$width, "tau" = m$tau, "FWHM" = 2.35 * m$width,
    "height" = m$height, "area" = area, "r.squared" = r.squared, 
    purity = purity)
}

#' Fit peak (raw)
#' @noRd
fitpk_raw <- function(x, pos, lambda, max_iter,
                      estimate_purity = TRUE, noise_threshold = .001){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  
  # perform trapezoidal integration on raw signal
  area <- sum(diff(peak.loc) * mean(c(y[peak.loc][-1], tail(y[peak.loc],-1))))
  purity <- get_purity(x = x, pos = pos, try = estimate_purity,
                       noise_threshold = noise_threshold)
  c("rt" = pos[1], "start" = pos[2], "end" = pos[3], 
    "sd" = pos[3] - pos[2], "FWHM" = 2.35 * pos[3] - pos[2],
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
bemg <- function(x, center, width, height, a, b, floor = 0) {
  tm <- x - center
  as2 <- a * width^2
  bs2 <- b * width^2
  sq2s2 <- sqrt(2 * width^2)
  
  tau_right <- exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / sq2s2)
  tau_left <- exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / sq2s2)
  
  return(floor + height * 0.5 * (tau_right + tau_left))
}

#' Fit Bi-Exponentially Modified Gaussian peak
#' @importFrom stats coef fitted nls.control quantile residuals
#' @noRd
fit_bemg <- function(x1, y1, start.center = NULL, start.width = NULL,
                     start.a = NULL, start.b = NULL, start.height = NULL,
                     start.floor = NULL, fit.floor = FALSE, max_iter = 1000) {
  
  who.max <- which.max(y1)
  if (is.null(start.center)) start.center <- x1[who.max]
  if (is.null(start.height)) start.height <- y1[who.max]
  if (is.null(start.width))  start.width  <- sum(y1 > (start.height / 2)) / 2
  if (is.null(start.a))      start.a      <- 1.0
  if (is.null(start.b))      start.b      <- 1.0
  
  controlList <- nls.control(maxiter = max_iter, minFactor = 1/512,
                             warnOnly = TRUE)
  starts <- list(center = start.center, width = start.width,
                 height = start.height, a = start.a, b = start.b)
  
  if (!fit.floor) {
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b),
      start = starts, control = controlList), silent = TRUE)
  } else {
    if (is.null(start.floor)) start.floor <- quantile(y1, seq(0, 1, 0.1))[2]
    starts <- c(starts, floor = start.floor)
    nlsAns <- try(minpack.lm::nlsLM(
      y1 ~ bemg(x1, center, width, height, a, b, floor),
      start = starts, control = controlList), silent = TRUE)
  }
  
  if (inherits(nlsAns, "try-error")) {
    yAns <- bemg(x1, start.center, start.width, start.height,
                 start.a, start.b, if (fit.floor) start.floor else 0)
    out <- list(center = start.center, width = start.width,
                height = start.height, a = start.a, b = start.b,
                y = yAns, residual = y1 - yAns)
    floorAns <- if (fit.floor) start.floor else 0
  } else {
    coefs <- coef(nlsAns)
    out <- list(center = coefs["center"], width = coefs["width"],
                height = coefs["height"], a = coefs["a"], b = coefs["b"],
                y = fitted(nlsAns), residual = residuals(nlsAns))
    floorAns <- if (fit.floor) coefs["floor"] else 0
  }
  
  if (fit.floor) out <- c(out, floor = floorAns)
  return(out)
}

#' @noRd
erfc <- function(x) 2 * stats::pnorm(-x * sqrt(2))
