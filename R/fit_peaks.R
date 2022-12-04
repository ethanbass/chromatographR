#' Find peaks in chromatographic profile
#' 
#' Find peaks in chromatographic profile.
#' 
#' Find peaks by looking for zero-crossings in the smoothed first derivative of
#' the signal (\code{y}) that exceed the specified slope threshold
#' (\code{slope_thresh}). Peaks can also be filtered by supplying a minimal
#' amplitude threshold (\code{amp_thresh}), to filter out peaks below the
#' specified height. 
#' 
#' It is recommended to do pre-processing using the \code{\link{preprocess}}
#' function before peak detection. Overly high chromatographic resolution can 
#' sometimes lead to peaks being split into multiple segments. In this case,
#' it is recommended to reduce the chromatogram resolution using the \code{dim2}
#' argument before proceeding with further analyses.
#' 
#' @importFrom caTools runmean
#' @importFrom minpack.lm nlsLM
#' @importFrom stats deriv lm ksmooth
#' @importFrom utils tail
#' @param y response (numerical vector)
#' @param smooth_type Type of smoothing. Either gaussian kernel (\code{gaussian}),
#' box kernel (\code{box}), savitzky-golay smoothing (\code{savgol}),
#' moving average (\code{mva}), triangular moving average (\code{tmva}), or 
#' no smoothing (\code{none}).
#' @param smooth_window Smoothing window. If the supplied value is between 0 and
#' 1, window will be interpreted as a proportion of points to include. Otherwise,
#' the window will be the absolute number of points to include in the window.
#' (Defaults to .001).
#' @param slope_thresh Minimum threshold for peak slope. (Defaults to 0).
#' @param amp_thresh Minimum threshold for peak amplitude. (Defaults to 0).
#' @param bounds Logical. If TRUE, includes peak boundaries in data.frame.
#' (Defaults to TRUE).
#' @return If bounds == TRUE, returns a data.frame containing the center, start,
#' and end of each identified peak. Otherwise, returns a numeric vector of peak
#' centers. All locations are expressed as indices.
#' @note The \code{find_peaks} function is adapted from matlab code in Prof.
#' Tom O'Haver's
#' \href{http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm}{
#' Pragmatic Introduction to Signal Processing}.
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' find_peaks(Sa_pr[[1]][,"220"])
#' @seealso \code{\link{fit_peaks}}, \code{\link{get_peaks}}
#' @references O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement.
#' /href{https://terpconnect.umd.edu/~toh/spectrum/} (Accessed January, 2022).
#' @export find_peaks
find_peaks <- function(y, smooth_type=c("gaussian", "box", "savgol", "mva","tmva","none"),
                       smooth_window = .001, slope_thresh = 0, amp_thresh = 0,
                       bounds = TRUE){
  smooth_type <- match.arg(smooth_type, c("gaussian", "box","savgol", "mva","tmva","none"))
  if (smooth_window < 1){
    smooth_window <- length(y)*smooth_window
  }
  #compute derivative (with or without smoothing)
  if (smooth_type == "savgol"){
    if ((smooth_window %% 2) == 0){smooth_window <- smooth_window + 1}
    d <- savgol(diff(y), fl = smooth_window)
  } else if (smooth_type == "mva"){
    d <- caTools::runmean(diff(y), k = smooth_window)
  } else if (smooth_type == 'gaussian'){
    d <- diff(ksmooth(seq_along(y), y, kernel = "normal", bandwidth = smooth_window)$y)
  }  else if (smooth_type == "box"){
    d <- diff(ksmooth(seq_along(y), y, kernel = "box", bandwidth = smooth_window)$y)
  } else if (smooth_type == "tmva"){
    d <- runmean(runmean(diff(y), k=smooth_window), k=smooth_window)
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

#' Fit chromatographic peaks to an exponential-gaussian hybrid or gaussian
#' profile
#' 
#' Fit peak parameters using exponential-gaussian hybrid or gaussian function.
#' 
#' Peak parameters are calculated using \code{fit_peaks}, which fits the data
#' to a gaussian or exponential-gaussian hybrid curve using non-linear least
#' squares estimation as implemented in \code{\link[minpack.lm:nlsLM]{nlsLM}}.
#' Area under the fitted curve is estimated using trapezoidal estimation.
#' 
#' @param x A chromatogram in matrix format.
#' @param lambda Wavelength to fit peaks at.
#' @param pos Locations of peaks in vector y. If NULL, \code{find_peaks} will
#' run automatically to find peak positions.
#' @param sd.max Maximum width (standard deviation) for peaks. Defaults to 50.
#' @param fit Function for peak fitting. (Currently exponential-gaussian hybrid
#' \code{egh}, \code{gaussian} and \code{raw} settings are supported). If \code{
#' raw} is selected, trapezoidal integration will be performed on raw data
#' without fitting a peak shape. Defaults to \code{egh}.)
#' @param max.iter Maximum number of iterations to use in nonlinear least
#' squares peak-fitting. (Defaults to 1000).
#' @param noise_threshold Noise threshold. Input to \code{get_purity}.
#' @param ... Additional arguments to \code{find_peaks}.
#' @return Function \code{fit_peaks} returns a matrix, whose columns contain
#' the following information:
#' \item{rt}{location of the maximum of the peak (x)}
#' \item{start}{start of peak (only included in table if `bounds==TRUE`)}
#' \item{end}{end of peak (only included in table if `bounds==TRUE`)}
#' \item{sd}{width of the peak (x)} \item{tau}{tau parameter (only included in
#' table if `fit=="egh"`)}
#' \item{FWHM}{full width at half maximum (x)}
#' \item{height}{height of the peak (y)}
#' \item{area}{peak area}
#' \item{r.squared}{r-squared value for linear fit of model to data.}
#' #' \item{purity}{spectral purity of peak as assessed by \code{\link{get_purity}}.}
#' Again, the first five elements (rt, start, end, sd and FWHM) are expressed
#' as indices, so not in terms of the real retention times. The transformation
#' to "real" time is done in function \code{get_peaks}.
#' @note The \code{\link{fit_peaks}} function is adapted from Dr. Robert
#' Morrison's
#' \href{https://github.com/robertdouglasmorrison/DuffyTools}{DuffyTools package}
#' as well as code published in Ron Wehrens'
#' \href{https://github.com/rwehrens/alsace}{alsace} package.
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' fit_peaks(Sa_pr[[1]], lambda="220")
#' @seealso \code{\link{find_peaks}}, \code{\link{get_peaks}}
#' @references
#' Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. \emph{Journal of
#' Chromatography A} \bold{915}:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions - A
#' good model for chromatographic peaks in isocratic HPLC? \emph{Chromatographia},
#' /bold{26}: 285-296. \doi{10.1007/BF02268168}.
#' @export fit_peaks
fit_peaks <- function (x, lambda, pos = NULL, sd.max = 50,
                       fit = c("egh", "gaussian", "raw"), 
                       max.iter = 1000, noise_threshold=.001, ...){
  y <- x[,lambda]
  fit <- match.arg(fit, c("egh", "gaussian", "raw"))
  if (is.null(pos)){
    pos <- find_peaks(y, ...)
  }
  tabnames <- switch(fit,
                     "gaussian" = c("rt", "start", "end", "sd", "FWHM", 
                                    "height", "area", "r-squared", "purity"),
                     "egh" = c("rt", "start", "end", "sd", "tau", "FWHM", 
                               "height", "area", "r.squared", "purity"),
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
                  "raw" = fitpk_raw)
  
  huhn <- data.frame(t(apply(pos, 1, fitpk, x = x,
                             lambda = lambda, max.iter = max.iter,
                             noise_threshold = noise_threshold)))
  colnames(huhn) <- tabnames
  huhn <- data.frame(sapply(huhn, as.numeric, simplify = FALSE))
  if (!is.null(sd.max)) {
    huhn <- huhn[huhn$sd < sd.max, ]
  }
  x <- try(huhn[huhn$rt > 0,], silent=TRUE)
  if(inherits(x,  "try-error")){NA} else {x}
}

#' Gaussian function
#' @note: Adapted from \href{https://github.com/robertdouglasmorrison/DuffyTools/blob/master/R/gaussian.R}
#' @noRd
gaussian <- function(x, center=0, width=1, height=NULL, floor=0) {
  
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

fit_gaussian <- function(x, y, start.center=NULL, start.width=NULL, start.height=NULL,
                         start.floor=NULL, fit.floor=FALSE, max.iter=1000) {
  # estimate starting values
  who.max <- which.max(y)
  if ( is.null( start.center)) start.center <- x[ who.max]
  if ( is.null( start.height)) start.height <- y[ who.max]
  if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2
  
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter = max.iter, minFactor=1/512, warnOnly=TRUE)
  starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
  if ( ! fit.floor) {
    nlsAns <- try(nlsLM( y ~ gaussian(x, center, width, height),
                         start=starts, control=controlList), silent=TRUE)
  } else{
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- c(starts,"floor"=start.floor)
    nlsAns <- try(nlsLM( y ~ gaussian( x, center, width, height, floor),
                         start=starts, control=controlList), silent=TRUE)
  }
  
  # package up the results to pass back
  
    if (inherits(nlsAns, "try-error")){
      yAns <- gaussian(x, start.center, start.width, start.height, start.floor)
      out <- list("center"=start.center, "width"=start.width, "height"=start.height,
                  "y"=yAns, "residual"= y - yAns)
      floorAns <- if ( fit.floor) start.floor else 0
    } else {
      coefs <-coef(nlsAns)
      out <- list( "center"=coefs[1], "width"=coefs[2], "height"=coefs[3],
                   "y"=fitted( nlsAns), "residual"=residuals(nlsAns))
      floorAns <- if ( fit.floor) coefs[4] else 0
    }
    if (fit.floor) {
      out <- c( out, "floor"=floorAns)
    }
  return( out)
}

###########################################################################################

#' Expontential-gaussian hybrid
#' @noRd
egh <- function(x, center, width,  height, tau, floor=0){
    result <- rep(0, length(x))
    index <- which(2*width^2 + tau*(x-center)>0)
    result[index] <- height*exp(-(x[index]-center)^2/(2*width^2 + tau*(x[index]-center)))
    return(result)
  }


#' Fit exponential-gaussian hybrid peak
#' @importFrom stats coef fitted lm nls.control quantile residuals
#' @noRd
fit_egh <- function(x1, y1, start.center=NULL, start.width=NULL, start.tau=NULL,
                    start.height=NULL, start.floor=NULL, fit.floor=FALSE,
                    max.iter=1000) {
  
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
  controlList <- nls.control(maxiter=max.iter, minFactor=1/512, warnOnly=TRUE)
  starts <- list("center"=start.center, "width"=start.width, 
                 "height"=start.height, "tau"=start.tau)
  if (!fit.floor){
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau),
                        start=starts, control=controlList), silent=TRUE)
  } else{
    if (is.null( start.floor)) start.floor <- quantile( y1, seq(0,1,0.1))[2]
    starts <- c(starts, "floor"=start.floor)
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau, floor),
                        start=starts, control=controlList), silent=TRUE)
  }
  
  # package up the results to pass back
  if (inherits(nlsAns, "try-error")) {
    yAns <- egh(x1, start.center, start.width, start.height, start.tau, start.floor)
    out <- list("center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau,
                "y"=yAns, "residual"= y1 - yAns)
    floorAns <- if ( fit.floor) start.floor else 0
  } else {
    coefs <-coef(nlsAns)
    out <- list( "center"=coefs[1], "width"=coefs[2], "height"=coefs[3], "tau"=coefs[4],
                 "y"=fitted( nlsAns), "residual"=residuals(nlsAns))
    floorAns <- if ( fit.floor) coefs[5] else 0
  }
  
  if (fit.floor) {
    out <- c( out, "floor"=floorAns)
  }
  return(out)
}

#' @noRd
fitpk_gaussian <- function(x, pos, lambda, max.iter, noise_threshold = .001, ...){
  
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_gaussian(peak.loc, y[peak.loc],
                                     start.center = xloc,
                                     start.height = y[xloc],
                                     max.iter = max.iter)
  )
  area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
  r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=TRUE)
  purity <- ifelse(dim(x)[2] > 1, try(get_purity(x, pos, noise_threshold = noise_threshold)), NA)
  c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
    "sd" = m$width, "FWHM" = 2.35 * m$width,
    "height" = y[xloc], "area" = area, "r.squared" = r.squared, purity=purity)
}

#' @noRd
fitpk_egh <- function(x, pos, lambda, max.iter, noise_threshold = .001){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  suppressWarnings(m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc,
                                start.height = y[xloc], max.iter = max.iter)
  )
  r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=TRUE)
  purity <- ifelse(dim(x)[2] > 1, try(get_purity(x, pos, noise_threshold)), NA)
  area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
  c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
    "sd" = m$width, "tau" = m$tau, "FWHM" = 2.35 * m$width,
    "height" = y[xloc], "area" = area, "r.squared" = r.squared, purity=purity)
}

#' @noRd
fitpk_raw <- function(x, pos, lambda, max.iter, noise_threshold = .001){
  y <- x[,lambda]
  xloc <- pos[1]
  peak.loc <- seq.int(pos[2], pos[3])
  area <- sum(diff(peak.loc) * mean(c(y[peak.loc][-1], tail(y[peak.loc],-1)))) # trapezoidal integration
  purity <- ifelse(dim(x)[2] > 1, try(get_purity(x, pos, noise_threshold)), NA)
  c("rt" = pos[1], "start" = pos[2], "end" = pos[3], 
    "sd" = pos[3]-pos[2], "FWHM" = 2.35 * pos[3]-pos[2],
    "height" = y[xloc], "area" = area, purity = purity)
}


#' Savitsky Golay Smoothing from pracma
#' @author Hans W. Borchers
#' @param T Vector of signals to be filtered
#' @param fl Filter length (for instance fl = 51..151), has to be odd.
#' @param forder filter order Filter order (2 = quadratic filter, 4 = quartic).
#' @param dorder Derivative order (0 = smoothing, 1 = first derivative, etc.).
#' @importFrom stats convolve
savgol <- function(T, fl, forder = 4, dorder = 0) {
  stopifnot(is.numeric(T), is.numeric(fl))
  if (fl <= 1 || fl %% 2 == 0)
    stop("Argument 'fl' must be an odd integer greater than 1.")
  n <- length(T)
  
  # -- calculate filter coefficients --
  fc <- (fl-1)/2                          # index: window left and right
  X <- outer(-fc:fc, 0:forder, FUN="^")   # polynomial terms and coeffs
  Y <- pinv(X);                           # pseudoinverse
  
  # -- filter via convolution and take care of the end points --
  T2 <- convolve(T, rev(Y[(dorder+1),]), type="o")   # convolve(...)
  T2 <- T2[(fc+1):(length(T2)-fc)]
  
  Tsg <- (-1)^dorder * T2
  return( Tsg )
}
#' pinv port from pracma
#' @noRd
pinv <- function (A, tol = .Machine$double.eps^(2/3)) {
  stopifnot(is.numeric(A) || is.complex(A), is.matrix(A))
  
  s <- svd(A)
  # D <- diag(s$d); Dinv <- diag(1/s$d)
  # U <- s$u; V <- s$v
  # A = U D V'
  # X = V Dinv U'
  if (is.complex(A)) s$u <- Conj(s$u)
  
  p <- ( s$d > max(tol * s$d[1], 0) )
  if (all(p)) {
    mp <- s$v %*% (1/s$d * t(s$u))
  } else if (any(p)) {
    mp <- s$v[, p, drop=FALSE] %*% (1/s$d[p] * t(s$u[, p, drop=FALSE]))
  } else {
    mp <- matrix(0, nrow=ncol(A), ncol=nrow(A))
  }
  
  return(mp)
}