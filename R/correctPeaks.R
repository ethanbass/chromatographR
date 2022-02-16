#' Correct peak positions according to a ptw warping model
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in a separate column
#' in the list of peak tables.
#' 
#' 
#' @param peakList A nested list of peak tables: the first level is the sample,
#' and the second level is the component. Every component is described by a
#' matrix where every row is one peak, and the columns contain information on
#' retention time, full width at half maximum (FWHM), peak width, height, and
#' area.
#' @param modList A list of ptw models.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention time.
#' @author Ron Wehrens
#' @seealso \code{\link{correctRT}}
#' @keywords manip
correctPeaks <- function(peakList, modList) {
  mapply(function(samp, mod) {
    lapply(samp,
           function(profile) {
             if (nrow(profile) > 0) {
               cbind(profile,
                     rt.cor = c(predict(mod, profile[,1], what = "time")))
             } else {
               cbind(profile, rt.cor = rep(0, 0))
             }
           })},
         peakList, modList, SIMPLIFY = FALSE)
}
