# peak finding function adapted from matlab function by Prof. Tom O'Haver
## (see http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm)

find_peaks <- function(y, smooth_type="gaussian", smooth_window = 1, smooth_width = 0.1,
                               slope_thresh=.05, amp_thresh=0, bounds=T){
  #compute derivative (with or without smoothing)
  if (smooth_type=='gaussian'){
    d <- smth.gaussian(diff(y),window = smooth_window, alpha=smooth_width)
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

fit_peaks <- function (y, pos, sd.max = 50, fit = c("egh", "gaussian"), max.iter = 1000){
  fit <- match.arg(fit, c("egh", "gaussian"))
  if (fit == "gaussian") {
    tabnames <- c("rt", "start", "end", "sd", "FWHM", "height", "area", "r-squared")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos$pos, function(x) y[x + 1] == 0 |
                        y[x - 1] == 0)
    pos <- pos[!on.edge,]
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
      c("rt" = m$center, "start" = pos[2], "end" = pos[3], "sd" = m$width, "FWHM" = 2.35 * m$width,
        "height" = y[xloc], "area" = area, "r.squared" = r.squared)
    }
  }
  else if (fit == "egh") {
    tabnames <- c("rt", "start", "end", "sd", "tau", "FWHM", "height", "area", 
                  "r.squared")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos$pos, function(x) y[x + 1] == 0 |
                        y[x - 1] == 0)
    pos <- pos[!on.edge,]
    if (nrow(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(pos){
      xloc <- pos[1]
      peak.loc <- seq.int(pos[2], pos[3])
      m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc, start.height = y[xloc], max.iter = max.iter)
      r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent=T)
      area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
      c("rt" = m$center, "start" = pos[2], "end" = pos[3], "sd" = m$width, "tau" = m$tau, "FWHM" = 2.35 * m$width,
        "height" = y[xloc], "area" = area, "r.squared" = r.squared)
    }
  }
  huhn <- data.frame(t(apply(pos, 1, fitpk)))
  colnames(huhn) <- tabnames
  huhn <- data.frame(sapply(huhn, as.numeric))
  if (!is.null(sd.max)) {
    huhn <- huhn[huhn$sd < sd.max, ]
  }
  x <- try(huhn[huhn$rt>0,],silent=T)
  if(inherits(x,  "try-error")){NA} else {x}
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
    nlsAns <- try(nlsLM( y ~ gaussian( x, center, width, height), start=starts, control=controlList), silent=T)
  } else{
    if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
    starts <- c(starts,"floor"=start.floor)
    nlsAns <- try(nlsLM( y ~ gaussian( x, center, width, height, floor), start=starts, control=controlList), silent=T)
  }
  
  # package up the results to pass back
  
    if (class( nlsAns) == "try-error") {
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
### expontential-gaussian hybrid
egh <- function(x, center, width,  height, tau, floor=0){
    result <- rep(0, length(x))
    index <- which(2*width^2 + tau*(x-center)>0)
    result[index] <- height*exp(-(x[index]-center)^2/(2*width^2 + tau*(x[index]-center)))
    return(result)
  }

fit_egh <- function(x1, y1, start.center=NULL, start.width=NULL, start.tau=NULL, start.height=NULL,
                    start.floor=NULL, fit.floor=FALSE, max.iter=1000) {
  
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
  starts <- list("center"=start.center, "width"=start.width, "height"=start.height, "tau"=start.tau)
  if (!fit.floor){
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau), start=starts, control=controlList), silent=T)
  } else{
    if (is.null( start.floor)) start.floor <- quantile( y1, seq(0,1,0.1))[2]
    starts <- c(starts, "floor"=start.floor)
    nlsAns <- try(nlsLM(y1 ~ egh(x1, center, width, height, tau, floor), start=starts, control=controlList), silent=T)
  }
  
  # package up the results to pass back
  if ( class( nlsAns) == "try-error") {
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

#########################

fitpeaks_at_max <- function (mat, pos, w=5, sd.max=50, fit=c("gaussian","egh")){
  #names(y) <- NULL
  if(fit=="gaussian"){
    tabnames <- c("rt", "lambda", "sd", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, length(tabnames)), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    y<-mat[,1]
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge,]
    if (length(pos) == 0) 
      return(noPeaksMat)
    fitpk <- function(xloc){
      lambda <- which.max(mat[xloc,])
      y <- mat[,lambda]
      peak.loc <- seq.int(xloc-w, xloc+w)
      m <- fit_gaussian(peak.loc, y[peak.loc])
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
    pos <- pos[!on.edge,]
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
