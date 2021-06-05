## Function to extract all peaks from the concentration profiles of an
## ALS object. For each sample, a list of matrices is returned,
## corresponding to the peaks found in each of the components.

getAllPeaks <- function (CList, wavelengths, ...){
  peaks<-list()
  CList2 <- lapply(CList, function(Cmat) Cmat[,wavelengths])
  peakPositions <- lapply(CList2, function(Cmat){
    apply(Cmat, 2, function(x) findpeaks(x))})
  Cmat <- CList2[5]
  result <- lapply(1:length(CList2), function(smpl) {
    ptable <- lapply(1:length(peakPositions[[smpl]]), function(cmpnd) fitpeaks(CList2[[smpl]][,cmpnd], peakPositions[[smpl]][[cmpnd]],...))
    names(ptable) <- names(peakPositions[[smpl]])
    ptable
  })
  names(result) <- names(peakPositions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
                                                                               1, function(x) !any(is.na(x))), , drop = FALSE]))
  timepoints <- as.numeric(rownames(CList[[1]]))
  tdiff <- median(diff(timepoints))
  lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
    x <- cmpnd
    x[, 1] <- timepoints[x[, 1]]
    x[, 2:3] <- x[, 2:3] * tdiff
    x
  }))
}

getAllPeaks2 <- function (CList, wavelengths, span = NULL, findpeaks=findpeaks_by_slope,
                         fitpeaks=fitpeaks, ...){
  peaks<-list()
  CList2 <- lapply(CList, function(Cmat) Cmat[,wavelengths])
  peakPositions <- lapply(CList2, function(Cmat){
    apply(Cmat, 2, function(x) findpeaks(x))})
  Cmat <- CList2[5]
  result <- lapply(1:length(CList2), function(smpl) {
    ptable <- lapply(1:length(peakPositions[[smpl]]), function(cmpnd) fitpeaks(CList2[[smpl]], peakPositions[[smpl]][[cmpnd]],...))
    names(ptable) <- names(peakPositions[[smpl]])
    ptable
  })
  names(result) <- names(peakPositions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
                                                                               1, function(x) !any(is.na(x))), , drop = FALSE]))
  timepoints <- as.numeric(rownames(CList[[1]]))
  tdiff <- median(diff(timepoints))
  lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
    x <- cmpnd
    x[, 1] <- timepoints[x[, 1]]
    x[, 2:3] <- x[, 2:3] * tdiff
    x
  }))
}
