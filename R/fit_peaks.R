# peak finding function adapted from matlab function by Prof. Tom O'Haver
## (see http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm)

find_peaks <- function(y, smooth_type="gaussian", smooth_window = 1, smooth_width = 0.1,
                               slope_thresh=.05, amp_thresh=0, bounds=F){
  #compute derivative (with or without smoothing)
  if (smooth_type=='gaussian'){
    d <- smoother::smth.gaussian(diff(y),window = smooth_window, alpha=smooth_width)
  } else{
    d=deriv(y)
  }
  p1 <- which(sign(d[1:(length(d)-1)])>sign(d[2:length(d)])) # detects zero-crossing of first derivative (peak apex)
  p2 <- which(abs(diff(d)) > slope_thresh) # detects second derivative exceeding slope threshold
  p3 <- which(y > amp_thresh) # detects y-vals exceeding amplitude threshold
  p <- intersect(intersect(p1,p2), p3)
  if (bounds==T){
    p4 <- which(sign(d[1:(length(d)-1)]) < sign(d[2:length(d)]))
    bl <- sapply(p, function(v) max(p4[p4 < v])) # lower bound
    bl[which(bl == -Inf)]<-0
    bu <- sapply(p,function(v) min(p4[p4 > v])) # upper bound
    bu[which(bu==Inf)] <- length(y)
    data.frame(pos=p, lower=bl, upper=bu)
  } else 
  p
}

# fit peaks to gaussian or exponential-gaussian hybrid ('egh') function using nlsLM

fit_peaks <- function (y, pos, sd.max = 50, fit = c("egh", "gaussian", "emg"), 
                      max.iter = 100) 
{
  fit <- match.arg(fit, c("egh", "gaussian", "emg"))
  if (fit == "gaussian") {
    tabnames <- c("rt","start","end", "sd", "FWHM", "height", "area", "r-squared")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    # on.edge <- sapply(pos$pos, function(x) y[x + 1] == 0 | 
    #                     y[x - 1] == 0)
    # pos$pos <- pos$pos[!on.edge]
    if (nrow(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(pos) {
      xloc <- pos[1]
      peak.loc <- seq.int(pos[2], pos[3])
      m <- fit_gaussian(peak.loc, y[peak.loc], start.center = xloc, 
                        start.height = y[xloc], max.iter = max.iter)
      area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
      # area <- y[xloc]/dnorm(m$center, m$center, m$width)
      r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=T)
      c("rt" = m$center, "start" = pos[1], "end" = pos[2], "sd" = m$width, "FWHM" = 2.35 * m$width,
        "height" = y[xloc], "area" = area, "r.squared" = r.squared)
    }
  }
  else if (fit == "egh") {
    tabnames <- c("rt","start","end", "sd", "tau", "FWHM", "height", "area", 
                  "r.squared")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    # on.edge <- sapply(pos$pos, function(x) y[x + 1] == 0 | 
    #                     y[x - 1] == 0)
    # pos$pos <- pos$pos[!on.edge]
    if (nrow(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(pos){
      xloc <- pos[1]
      peak.loc <- seq.int(pos[2], pos[3])
      m <- fit_egh(peak.loc, y[peak.loc], 
                                    start.center = xloc, start.height = y[xloc])
      r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=T)
      area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
      c("rt" = m$center, "start" = pos[1], "end" = pos[2], "sd" = m$width, "tau" = m$tau, "FWHM" = 2.35 * m$width,
        "height" = y[xloc], "area" = area, "r.squared" = r.squared)
    }
  }
  else if (fit == "emg") {
    tabnames <- c("rt","start","end", "sd", "tau", "FWHM", "height", "area", 
                  "r.squared")
    #tabnames <- c("rt", "sd", "lambda", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    # on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
    #                                                        1] == 0)
    # pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc) {
      xloc <- pos[1]
      peak.loc <- seq.int(pos[2], pos[3])
      m <- fit_EMG(peak.loc, y[peak.loc], start.center = xloc, 
                   start.height = y[xloc])
      r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=T)
      c(m$center, pos[2], pos[3], m$width, m$tau, 2.35 * m$width, y[xloc], 
        y[xloc]/dnorm(m$center, m$center, m$width), r.squared)
    }
  }
  huhn <- data.frame(t(apply(pos, 1, fitpk)))
  colnames(huhn) <- tabnames
  huhn <- data.frame(sapply(huhn, as.numeric))
  if (!is.null(sd.max)) {
    huhn <- huhn[huhn$sd < sd.max, ]
  }
  huhn[huhn$rt>0,]
}
#################################################################################################
### gaussian
## from https://github.com/robertdouglasmorrison/DuffyTools/blob/master/R/gaussian.R

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


fit_gaussian <- function(x, y, start.center=NULL, start.width=NULL, start.height=NULL,
                         start.floor=NULL, fit.floor=FALSE, max.iter) {
  # estimate starting values
  who.max <- which.max(y)
  if ( is.null( start.center)) start.center <- x[ who.max]
  if ( is.null( start.height)) start.height <- y[ who.max]
  if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2
  
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter = max.iter, minFactor=1/512, warnOnly=TRUE)
  if ( ! fit.floor) {
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
    nlsAns <- try(nlsLM( y ~ gaussian( x, center, width, height), start=starts, control=controlList), silent=T)
  } else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
                    "floor"=start.floor)
    nlsAns <- try(nlsLM( y ~ gaussian( x, center, width, height, floor), start=starts, control=controlList), silent=T)
  }
  
  # package up the results to pass back
  if (class( nlsAns) == "try-error") {
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


###########################################################################################
### expontential-gaussian hybrid
egh <- function(x, center, width,  height, tau, floor=0){
    result <- rep(0, length(x))
    index <- which(2*width^2 + tau*(x-center)>0)
    result[index] <- height*exp(-(x[index]-center)^2/(2*width^2 + tau*(x[index]-center)))
    return(result)
  }

fit_egh <- function(x1, y1, start.center=NULL, start.width=NULL, start.tau=NULL, start.height=NULL,
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
  if (is.null(start.tau)){
    start.tau <- 0
  }
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
  if (!fit.floor){
    starts <- list("center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau)
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau), start=starts, control=controlList), silent=T)
  }else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau, 
                    "floor"=start.floor)
    nlsAns <- try(nlsLM( y ~ egh(x, center, width, height, tau, floor), start=starts, control=controlList), silent=T)
  }
  
  # package up the results to pass back
  if ( class( nlsAns) == "try-error") {
    centerAns <- start.center
    widthAns <- start.width
    heightAns <- start.height
    tauAns <- start.tau
    floorAns <- if ( fit.floor) start.floor else 0
    yAns <- egh(x, centerAns, widthAns, heightAns, tauAns, floorAns)
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
  
  # always report the SD as a positive value
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

fit_EMG <- function(x1, y1, start.center=NULL, start.width=NULL, start.alpha=NULL, start.height=NULL,
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
    nlsAns <- try(nlsLM(y1 ~ EMG(x1, mu, sigma, height, alpha), start=starts, control=controlList,
                      lower = c(0,0,0,0), upper=c(Inf,Inf,Inf,Inf), algorithm="port"), silent=T)
  }else {
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- list( "mu"=start.center, "sigma"=start.width, "height"=start.height, "alpha"=start.alpha, 
                    "floor"=start.floor)
    nlsAns <- try(nlsLM(y ~ EMG(x, mu, sigma, height, alpha), start=starts, control=controlList,
               lower = c(0,0,0,0), upper=c(Inf,Inf,Inf,Inf), algorithm="port"), silent=T)
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
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
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
      m <- fit_gaussian(peak.loc, y[peak.loc])
      #fit <- fit_egh(peak.loc, y[peak.loc])
      #fit <- egh_optim(peak.loc, y[peak.loc], par=c(1,1,1,1), upper=c(Inf,Inf,Inf,Inf), lower=c(0,0,0,0))
      c(m$center, colnames(mat)[lambda], m$width, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                m$width))
    }
  }
  if(fit=="egh"){
    tabnames <- c("rt", "lambda","sd","tau", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      m <- fit_egh(peak.loc, y[peak.loc])
      c(m$center, m$width, m$tau, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                       m$width))
    }
  }
  huhn <- data.frame(t(sapply(pos, fitpk)))
  colnames(huhn) <- tabnames
  huhn[huhn$sd<sd.max,]
}


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


#emg::emg.mle()