filterPeaks <- function(peakList, minHeight, minArea, minWHM, maxWHM) {
  if (missing(minHeight) & missing(minArea) &
      missing(maxWHM) & missing(minWHM)) {
    warning("Nothing to filter...")
    return(peakList)
  }

  if (missing(maxWHM)) ## find max value in peakList and add 1 to be sure...
      maxWHM <- 1 + max(unlist(sapply(peakList,
                                      function(samp)
                                      sapply(samp,
                                             function(comp) comp[,"FWHM"]))))
  if (missing(minWHM)) minWHM <- 0
  if (missing(minHeight)) minHeight <- 0
  if (missing(minArea)) minArea <- 0

  lapply(peakList,
         function(smpl)
         lapply(smpl,
                function(comp)
                comp[comp[,"FWHM"] < maxWHM &
                     comp[,"FWHM"] > minWHM &
                     comp[,"height"] > minHeight &
                     comp[,"area"] > minArea,,drop = FALSE]))
}
