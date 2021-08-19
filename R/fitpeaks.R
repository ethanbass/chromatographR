

## legacy peak-finding function. detects local maximum within windows of width determined by "span"
# 
# findpeaks <- function(y, span = NULL)
# {
#   if (is.null(span)) span <- round(.2 * length(y))
#   
#   z <- embed(y, span)
#   s <- span %/% 2
#   v <- max.col(z, ties.method = "first") == 1 + s
# 
#   which(c(rep(FALSE, s), v, rep(FALSE, s)))
# }

## new peak finding function ported from matlab
## (see http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm)

findpeaks <- function(y, smooth_type='gaussian', smooth_window = 1, smooth_width = 0.1,
                               slope_thresh=.05, amp_thresh=0){
  if (smooth_type=='gaussian'){
    d <- smoother::smth.gaussian(diff(y),window = smooth_window, alpha=smooth_width)
  } else{
    d=deriv(y)
  }
  p1 <- which(sign(d[1:(length(d)-1)])>sign(d[2:length(d)])) #detects zero-crossing
  p2 <- which(abs(diff(d)) > slope_thresh) # detects second derivative exceeding slope threshold
  p3 <- which(y > amp_thresh) # detects y vals exceeding amplitude threshold
  p <- intersect(intersect(p1,p2), p3)
  p
}

# fit peaks using gaussian distribution (egh setting (exponential gaussian hybrid) doesn't work yet).

fitpeaks <- function (y, pos, w=1, sd.max=50, fit=c("gaussian","egh","emg")){
  #names(y) <- NULL
  fit <- match.arg(fit,c("gaussian","egh","emg"))
  if (fit=="gaussian"){
    tabnames <- c("rt", "sd", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 5), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      m <- fit.gaussian(peak.loc, y[peak.loc], start.center = xloc, 
                        start.height = y[which(peak.loc==xloc)])
      c(m$center, m$width, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                m$width))
    }
  } else if(fit == "egh"){
    tabnames <- c("rt", "sd","tau", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 6), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      m <- chromatographR:::fit.egh(peak.loc, y[peak.loc], start.center = xloc,
                                    start.height = y[which(peak.loc==xloc)])
      c(m$center, m$width, m$tau, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                       m$width))
    }
  } else if(fit == "emg"){
    tabnames <- c("rt", "sd","lambda", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 6), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      #m <- emg::emg.mle(y[peak.loc])
      m <- fit.EMG(peak.loc,y[peak.loc], start.center = xloc,
                   start.height = y[which(peak.loc==xloc)])
      c(m$center, m$width, m$tau, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                       m$width))
    }
  }
  huhn <- data.frame(t(sapply(pos, fitpk)))
  colnames(huhn) <- tabnames
  huhn[huhn$sd<sd.max,]
}
####################################
# from https://github.com/robertdouglasmorrison/DuffyTools/blob/master/R/gaussian.R

gaussian <- function( x, center=0, width=1, height=NULL, floor=0) {
  
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


fit.gaussian <- function(x, y, start.center=NULL, start.width=NULL, start.height=NULL,
                         start.floor=NULL, fit.floor=FALSE) {
  
  # try to find the best gaussian to fit the given data
  
  # make some rough estimates from the values of Y
  who.max <- which.max(y)
  if ( is.null( start.center)) start.center <- x[ who.max]
  if ( is.null( start.height)) start.height <- y[ who.max]
  if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2
  
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
  if ( ! fit.floor) {
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
    nlsAns <- try( nls( y ~ gaussian( x, center, width, height), start=starts, control=controlList))
  } else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
                    "floor"=start.floor)
    nlsAns <- try(nls( y ~ gaussian( x, center, width, height, floor), start=starts, control=controlList))
  }
  
  # package up the results to pass back
  if ( class( nlsAns) == "try-error") {
    centerAns <- start.center
    widthAns <- start.width
    heightAns <- start.height
    floorAns <- if ( fit.floor) start.floor else 0
    yAns <- gaussian( x, centerAns, widthAns, heightAns, floorAns)
    residualAns <- y - yAns
  } else {
    coefs <-coef(nlsAns)
    centerAns <- coefs[1]
    widthAns <- coefs[2]
    heightAns <- coefs[3]
    floorAns <- if ( fit.floor) coefs[4] else 0
    yAns <- fitted( nlsAns)
    residualAns <- residuals( nlsAns)
  }
  
  # always report the SD as a possitive value
  widthAns <- abs( widthAns) #width = SD
  
  out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "y"=yAns,
               "residual"=residualAns)
  if ( fit.floor) {
    out <- c( out, "floor"=floorAns)
  }
  
  return( out)
}


####################################

egh <- 
  function(x, center, width,  height, tau, floor=0){
    result <- rep(0, length(x))
    index <- which(2*width^2 + tau*(x-center)>0)
    result[index] <- height*exp(-(x[index]-center)^2/(2*width^2 + tau*(x[index]-center)))
    return(result)
  }

fit.egh <- function(x1, y1, start.center=NULL, start.width=NULL, start.tau=NULL, start.height=NULL,
                    start.floor=NULL, fit.floor=FALSE) {
  
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
  #A<-which.max(which(y[1:start.center] < (start.height/2)))
  #B<-which(y[start.center:length(y)] < (start.height/2))[1]+start.center
  #if (is.null(start.tau)) start.tau <- y<(start.height/2)
  #if (is.null(start.tau)) start.tau <- (B-A)/(-log(.5))
  if (is.null(start.tau)){
    start.tau <- 0
  }
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
  if (!fit.floor){
    starts <- list("center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau)
    nlsAns <- try(nls(y1 ~ egh(x1, center, width, height, tau), start=starts, control=controlList))
  }else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau, 
                    "floor"=start.floor)
    nlsAns <- try(nls( y ~ egh(x, center, width, height, tau, floor), start=starts, control=controlList))
  }
  
  # package up the results to pass back
  if ( class( nlsAns) == "try-error") {
    centerAns <- start.center
    widthAns <- start.width
    heightAns <- start.height
    tauAns <- start.tau
    floorAns <- if ( fit.floor) start.floor else 0
    yAns <- egh( x, centerAns, widthAns, heightAns, tauAns, floorAns)
    residualAns <- y - yAns
  } else {
    coefs <-coef(nlsAns)
    centerAns <- coefs[1]
    widthAns <- coefs[2]
    heightAns <- coefs[3]
    tauAns <- coefs[4]
    floorAns <- if ( fit.floor) coefs[5] else 0
    yAns <- fitted( nlsAns)
    residualAns <- residuals( nlsAns)
  }
  
  # always report the SD as a possitive value
  widthAns <- abs( widthAns) #width = SD
  
  out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "tau"=tauAns, "y"=yAns,
               "residual"=residualAns)
  if ( fit.floor) {
    out <- c( out, "floor"=floorAns)
  }
  
  return(out)
}

####################################

EMG <- function (x, mu = 0, sigma = 1, height = 1, alpha = 1){
  f1 <- function (x, alpha, sigma, mu) 
    exp(sigma^2/(2 * alpha^2) + (mu - x)/alpha)
  f2 <- function (alpha, sigma) 
    sqrt(2 * pi) * sigma/alpha
  f3 <- function (x, alpha, sigma, mu) 
  {
    z <- (sigma/alpha + (mu - x)/sigma)/sqrt(2)
    erfc(z)
  }
  erfc <- function (z) 
    pnorm(z * sqrt(2), lower.tail = FALSE)
  height * f1(x, alpha, sigma, mu) * f2(alpha, sigma) * f3(x, alpha, sigma, 
                                                           mu)
}

fit.EMG <- function(x1, y1, start.center=NULL, start.width=NULL, start.alpha=NULL, start.height=NULL,
                    start.floor=NULL, fit.floor=FALSE) {
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
  if (is.null(start.alpha)){
    start.alpha <- 1
  }
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
  if (!fit.floor){
    starts <- list("mu"=start.center, "sigma"=start.width, "height"=start.height, "alpha"=start.alpha)
    #nlsAns <- try(nls(y1 ~ EMG(x1, mu, sigma, height, alpha), start=starts, control=controlList))
    nlsAns <- try(nls(y1 ~ EMG(x1, mu, sigma, height, alpha), start=starts, control=controlList,
                      lower = lowers, upper=c(Inf,Inf,Inf,Inf), algorithm="port"))
  }else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "mu"=start.center, "sigma"=start.width, "height"=start.height, "alpha"=start.alpha, 
                    "floor"=start.floor)
    #nlsAns <- try(nls(y ~ EMG(x, mu, sigma, height, alpha, floor), start=starts, control=controlList))
    nlsAns <- try(nls(y ~ EMG(x, mu, sigma, height, alpha), start=starts, control=controlList,
               lower = lowers, upper=c(Inf,Inf,Inf,Inf), algorithm="port"))
  }
  
  # package up the results to pass back
  if ( class( nlsAns) == "try-error") {
    centerAns <- start.center
    widthAns <- start.width
    heightAns <- start.height
    alphaAns <- start.alpha
    floorAns <- if ( fit.floor) start.floor else 0
    yAns <- EMG( x, centerAns, widthAns, heightAns, alphaAns, floorAns)
    residualAns <- y - yAns
  } else {
    coefs <-coef(nlsAns)
    centerAns <- coefs[1]
    widthAns <- coefs[2]
    heightAns <- coefs[3]
    tauAns <- coefs[4]
    floorAns <- if ( fit.floor) coefs[5] else 0
    yAns <- fitted( nlsAns)
    residualAns <- residuals( nlsAns)
  }
  
  # always report the SD as a possitive value
  widthAns <- abs( widthAns) #width = SD
  
  out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "tau"=tauAns, "y"=yAns,
               "residual"=residualAns)
  if ( fit.floor) {
    out <- c( out, "floor"=floorAns)
  }
  
  return(out)
}

fitpeaks_at_max <- function (mat, pos, w=5, sd.max=50, fit=c("gaussian","egh")){
  #names(y) <- NULL
  if(fit=="gaussian"){
    tabnames <- c("rt", "lambda", "sd", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 6), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    y<-mat[,1]
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      lambda <- which.max(mat[xloc,])
      y <- mat[,lambda]
      peak.loc <- seq.int(xloc-w, xloc+w)
      m <- fit.gaussian(peak.loc, y[peak.loc])
      #fit <- fit.egh(peak.loc, y[peak.loc])
      #fit <- egh_optim(peak.loc, y[peak.loc], par=c(1,1,1,1), upper=c(Inf,Inf,Inf,Inf), lower=c(0,0,0,0))
      c(m$center, colnames(mat)[lambda], m$width, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                m$width))
    }
  }
  if(fit=="egh"){
    tabnames <- c("rt", "lambda","sd","tau", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 7), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      m <- fit.egh(peak.loc, y[peak.loc])
      c(m$center, m$width, m$tau, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                       m$width))
    }
  }
  huhn <- data.frame(t(sapply(pos, fitpk)))
  colnames(huhn) <- tabnames
  huhn[huhn$sd<sd.max,]
}

