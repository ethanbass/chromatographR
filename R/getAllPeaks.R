## Function to extract all peaks from the concentration profiles of an
## ALS object. For each sample, a list of matrices is returned,
## corresponding to the peaks found in each of the components.

getAllPeaks <- function(CList, span = NULL, eps = 1e-1) {
  peakPositions <- lapply(CList,
                          function(Cmat)
                            lapply(1:ncol(Cmat),
                                   function(ii)
                                     findpeaks(Cmat[,ii], span = span)))
  
  result <-
    lapply(1:length(CList),
           function(smpl) {
             ptable <- lapply(1:length(peakPositions[[smpl]]),
                              function(cmpnd)
                                fitpeaks(CList[[smpl]][,cmpnd],
                                         peakPositions[[smpl]][[cmpnd]]))
             names(ptable) <- paste("Component",
                                    1:length(peakPositions[[smpl]]))
             ptable
           })
  names(result) <- names(peakPositions)
  ## sometimes NA values are returned, e.g. when a local max does not
  ## really correspond to a peak. Remove them!
  result <- lapply(result,
                   function(smpl)
                     lapply(smpl,
                            function(pks) 
                              pks[apply(pks, 1,
                                        function(x)
                                          !any(is.na(x))),,drop = FALSE]))
  
  timepoints <- as.numeric(rownames(CList[[1]]))
  tdiff <- median(diff(timepoints))
  
  ## convert results in indices to results in real time: the column
  ## mean must be replaced by the corresponding timepoints, and
  ## columns sd and FWHM must be multiplied by the difference between
  ## consecutive time points
  output <- lapply(result,
                   function(smpl)
                     lapply(smpl,
                            function(cmpnd) {
                              x <- cmpnd
                              x[,1] <- timepoints[x[,1]]
                              x[,2:3] <- x[,2:3] * tdiff
                              
                              x
                            }))
  
  ## only return those lines that have a width larger
  ## than zero. Sometimes we see zeros at sd, in those cases also area
  ## is zero. Perhaps we should check for empty matrices... let's see.
  lapply(output,
         function(sample)
           lapply(sample,
                  function(peakmat)
                    peakmat[peakmat[,"sd"] > eps,, drop = FALSE]))
                
}
