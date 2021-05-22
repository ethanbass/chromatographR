

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
  p <- intersect(p1,p2) %>% intersect(p3)
  p
}

# fit peaks using gaussian (egh setting doesn't work yet)

fitpeaks <- function (y, pos, w=5, sd.max=50, fit=c("gaussian","egh")){
  #names(y) <- NULL
  if(fit=="gaussian"){
    tabnames <- c("rt", "sd", "FWHM", "height", "area")
    noPeaksMat <- matrix(rep(NA, 5), nrow = 1, dimnames = list(NULL, 
                                                               tabnames))
    on.edge <- sapply(pos, function(x) y[x + 1] == 0 | y[x - 
                                                           1] == 0)
    pos <- pos[!on.edge]
    if (length(pos) == 0) 
      return(noPeaksMat)
    #xloc<-pks[1]
    fitpk <- function(xloc){
      peak.loc<-seq.int(xloc-w, xloc+w)
      m <- fit.gaussian(peak.loc, y[peak.loc])
      #fit <- fit.egh(peak.loc, y[peak.loc])
      #fit <- egh_optim(peak.loc, y[peak.loc], par=c(1,1,1,1), upper=c(Inf,Inf,Inf,Inf), lower=c(0,0,0,0))
      c(m$center, m$width, 2.35*m$width, y[xloc], y[xloc]/dnorm(m$center, m$center, 
                                                                m$width))
    }
  }
  if(fit=="egh"){
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
      m <- fit.egh(peak.loc, y[peak.loc])
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

## function to visually check integration accuracy
## fit is output of getallpeaks for chrome

plot_peaks <- function(chrome_list, pks, index=1, lambda='230', w=100, slope=.01, points=F, a=0.5, time=c('rt','raw'),
                       fit=c("gaussian","egh")){
  y<-chrome_list[[index]][,lambda]
  pks<-data.frame(pks[[index]][[lambda]])
  plot(new.ts,y,type='l',xlab='',ylab='',xaxt='n',yaxt='n')
  if (time=='rt'){
    pks$rt <- (pks$rt - new.ts[1])*100
    pks$sd <- pks$sd*100
  }
  if (points==T){
    points(new.ts[pks$rt],pks$height,pch=20,cex=0.5,col='red')
  }
  for (i in 1:nrow(pks)){
    #print(i/nrow(pks))
    xs<-seq.int((pks$rt[i]-w),(pks$rt[i]+w))
    mi <- xs[min(which(abs(diff(gaussian(xs, center=pks$rt[i], width=pks$sd[i],height = pks$height[i])))>h*slope))]
    ma <- xs[max(which(abs(diff(gaussian(xs, center=pks$rt[i], width=pks$sd[i],height = pks$height[i])))>h*slope))]
    if(is.na(mi)|is.na(ma)){
      next #skip bad peaks
    } else{
      xs2<-seq.int(mi,ma)
      if(fit=="gaussian"){
        polygon(new.ts[xs2], gaussian(xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i]),
                col=alpha('red',a), border=NA)
      }
      if(fit=="egh"){
        
        polygon(new.ts[xs2], egh(x=xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i], tau=pks$tau[i]),
                col=alpha('red',a), border=NA)
      }
    }
  }
}