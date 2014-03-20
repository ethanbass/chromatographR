findpeaks <- function(y, span = NULL)
{
  if (is.null(span)) span <- round(.2 * length(y))
  
  z <- embed(y, span)
  s <- span %/% 2
  v <- max.col(z, ties.method = "first") == 1 + s

  which(c(rep(FALSE, s), v, rep(FALSE, s)))
}

fitpeaks <- function(y, pos) {
  names(y) <- NULL
  tabnames <- c("rt", "sd", "FWHM", "height", "area")
  noPeaksMat <- matrix(rep(NA, 5), nrow = 1,
                       dimnames = list(NULL, tabnames))

  ## when combining windows, sometimes profiles show abrupt changes to
  ## pure zero values. Such points should not be seen as individual
  ## peaks. Given the algorithm for finding peak positions, pos can
  ## never be the extreme point of a profile, so x+1 and x-1 should
  ## always exist. Also spikes are taken out like this.
  on.edge <- sapply(pos,
                    function(x)
                    y[x+1] == 0 | y[x-1] == 0)
  pos <- pos[!on.edge]
  if (length(pos) == 0)
      return(noPeaksMat)

  fitpk <- function(xloc) {
    ## find all areas higher than half the current max
    peak.loc <- which(y > 0.5*y[xloc])
    peak.loc.diff <- diff(peak.loc)
    boundaries <- c(0, which(diff(peak.loc) != 1), length(peak.loc)+1)
    peaknrs <- rep(1:length(boundaries),
                   c(boundaries[1], diff(c(boundaries))))
    peaknrs[boundaries[-1]] <- NA
    current.peak <- peaknrs[peak.loc == xloc]
    current.peak <- current.peak[!is.na(current.peak)]
    if (length(current.peak) == 0)
        return(rep(NA, 5))
    
    ## only retain those points adjacent to the current max
    FWHM <- diff(range(peak.loc[peaknrs == current.peak],
                       na.rm = TRUE))
    pksd <- FWHM / (2*sqrt(2*log(2)))
    
    c(xloc, pksd, FWHM, y[xloc], y[xloc] / dnorm(xloc, xloc, pksd))
  }

  huhn <- t(sapply(pos, fitpk))
  colnames(huhn) <- tabnames

  huhn
}
