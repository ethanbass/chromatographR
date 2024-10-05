fit_gaussian <- function(x, y, start.center = NULL,
                         start.width = NULL, start.height = NULL,
                         start.floor = NULL, fit.floor = FALSE,
                         max.iter = 1000){
  # estimate starting values
  who.max <- which.max(y)
  if (is.null(start.center)) start.center <- x[who.max]
  if (is.null(start.height)) start.height <- y[who.max]
  if (is.null(start.width)) start.width <- sum( y > (start.height/2)) / 2
  
  # call the Nonlinear Least Squares, either fitting the floor too or not
  controlList <- nls.control(maxiter = max.iter, minFactor = 1/512,
                             warnOnly = TRUE)
  starts <- list( "center" = start.center, "width" = start.width,
                  "height" = start.height)
  if (!fit.floor) {
    nlsAns <- try(nlsLM( y ~ gaussian(x = x, center = center,
                                      width = width, height = height),
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