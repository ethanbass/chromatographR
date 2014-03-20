## wrapper for the als function from the ALS package, setting some
## always-used default values, and also calculating some summary
## statistics that for larger datasets take too long to be done in the
## summary.ALS function

doALS <- function(Xl, PureS, maxiter = 100) {
  Cini <- lapply(Xl, function(xl) xl[,1:ncol(PureS)])
  
  capture.output(result <- als(PsiList = Xl, CList = Cini, S = PureS,
                               maxiter = maxiter, normS = .5,
                               nonnegS = TRUE, optS1st = FALSE,
                               nonnegC = TRUE, uniC = FALSE,
                               baseline = FALSE) )

  colnames(result$S) <- paste("Component", 1:ncol(PureS))
  for (i in 1:length(result$CList))
      colnames(result$CList[[i]]) <- colnames(result$S)

  predicted.values <- lapply(1:length(Xl),
                             function(ii)
                             Xl[[ii]] - result$resid[[ii]])
  predvals2 <- sum(unlist(predicted.values)^2)
  npoints <- prod(c(length(result$CList),
                    nrow(result$CList[[1]]),
                    nrow(result$S)))
  result$summ.stats <-
      list(lof = 100 * sqrt(result$rss / predvals2),
           rms = sqrt(result$rss / npoints),
           r2 = 1 - result$rss / predvals2)
  
  class(result) <- "ALS"
  result
}

summary.ALS <- function(object, ...) {
  npoints <- prod(c(length(object$CList),
                    nrow(object$CList[[1]]),
                    nrow(object$S)))
  cat("ALS object fitting", length(object$CList),
      "samples with", ncol(object$S), "components.",
      "\nEach data matrix contains", nrow(object$S),
      "wavelengths and", nrow(object$CList[[1]]), "timepoints\n")

  cat("\tRMS fit error:", round(object$summ.stats$rms, 5), "\n")
  cat("\tLOF: ", round(object$summ.stats$lof, 2), "%\n", sep = "")
  cat("\tR2: ", round(object$summ.stats$r2, 5), "\n", sep = "")
  
  invisible()
}

print.ALS <- function(x, ...) {
  cat("ALS object fitting", length(x$CList), "samples; RMS fit error",
      round(x$summ.stats$rms, 5), "\n")
}

plot.ALS <- function(x,
                     what = c("spectra", "profiles"),
                     showWindows = TRUE,
                     mat.idx,
                     comp.idx,
                     xlab, ylab, main, ...) {
  what <- match.arg(what)

  if (missing(ylab)) ylab = "Response"
  if (missing(comp.idx)) comp.idx <- 1:ncol(x$S)
  if (missing(mat.idx)) mat.idx <- 1:length(x$CList)
  
  if (what == "spectra") {
    lambdas <- getWavelength(x)
    if (missing(xlab)) xlab <- "Wavelength"
    if (missing(main)) main <- NULL

    matplot(lambdas, x$S[,comp.idx], type = "l", lty = 1, 
            xlab = xlab, ylab = ylab, main = main, ...)
  }

  if (what == "profiles") {
    tpoints <- getTime(x)
    CL <- x$CList[mat.idx]
      
    if (missing(xlab)) xlab <- "Time"
    if (missing(main)) main <- names(CL)

    for (i in 1:length(CL)) {
      matplot(tpoints, CL[[i]][,comp.idx], xlab = xlab, ylab = ylab,
              type = "l", lty = 1, main = main[i], ...)

      if (showWindows) {
        if (!is.null(splitpoints <- attr(x, "splitpoints"))) {
          abline (v = splitpoints, col = "gray")
          
          if (!is.null(overlap <- attr(x, "overlap"))) {
            splitpoints.idx <- which(tpoints %in% splitpoints)
            overlap.boundaries <- c(tpoints[splitpoints.idx - overlap],
                                    tpoints[splitpoints.idx + overlap])
            abline(v = overlap.boundaries, col = "gray", lty = 2)
          }
        }
      }
    }
  }
  
  invisible()
}

getTime <- function(x) {
  if (is.null(rownames(x[[1]][[1]]))) {
    1:nrow(x[[1]][[1]])
  } else {
    as.numeric(rownames(x[[1]][[1]]))
  }
}

getWavelength <- function(x) {
  if (is.null(rownames(x$S))) {
    1:nrow(x$S)
  } else {
    as.numeric(rownames(x$S))
  }
}
