#' Filter peak lists
#' 
#' Utility function to remove peaks from a peak list, e.g. because their
#' intensity is too low. Currently one can filter on peak height, peak area,
#' standard deviation, and/or retention time.
#' 
#' @param peak_list A peak_list object, consisting of a nested list of peak
#' tables, where the first level is the sample, and the second level is the 
#' spectral component. Every component is described by a matrix where every row 
#' is one peak, and the columns contain information on retention time, 
#' full width at half maximum (FWHM), peak width, height, and area.
#' @param min_height Minimum peak height.
#' @param min_area Minimum peak area.
#' @param min_sd Minimal standard deviation.
#' @param max_sd Maximum standard deviation.
#' @param min_rt Minimum retention time.
#' @param max_rt Maxmimum retention time.
#' @return A peak list similar, with all rows removed
#' from the peak tables that are not satisfying the criteria.
#' @author Ron Wehrens, Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @export filter_peaks
filter_peaks <- function(peak_list, min_height, min_area,
                         min_sd, max_sd,min_rt, max_rt) {
  if (missing(min_height) & missing(min_area) &
      missing(min_sd) & missing(max_sd) &
      missing(min_rt) & missing(max_rt)) {
    warning("Nothing to filter...")
    return(peak_list)
  }
  x <- do.call(rbind, do.call(rbind, peak_list))
  if (!missing(min_height)){
    if (min_height < min(x$height))
      warning("'min_height' is less than minimum peak height.")}
  if (!missing(min_area)){
    if (min_area < min(x$area))
      warning("'min_area' is less than minimum peak area.")}
  if (!missing(min_sd)){
    if (min_sd < min(x$sd))
      warning("'min_sd' is less than minimum peak standard deviation.")}
  if (!missing(max_sd)){
    if (max_sd > max(x$sd))
      warning("'max_sd' is greater than maximum peak standard deviation,")}
  if (missing(max_sd)) ## find max value in peak_list and add 1 to be sure...
    max_sd <- 1 + max(unlist(sapply(peak_list,
                                    function(samp)
                                      sapply(samp,
                                             function(comp) comp[,"sd"]))))
  if (missing(max_rt)) ## find max value in peak_list and add 1 to be sure...
    max_rt <- 1 + max(unlist(sapply(peak_list,
                                    function(samp)
                                      sapply(samp,
                                             function(comp) comp[,"rt"]))))
  if (missing(min_sd)) min_sd <- 0
  if (missing(min_height)) min_height <- 0
  if (missing(min_area)) min_area <- 0
  if (missing(min_rt)) min_rt <- 0
  
  result <- lapply(peak_list,
         function(smpl)
           lapply(smpl,
                  function(comp)
                    comp[which(comp[,"sd"] < max_sd &
                           comp[,"sd"] > min_sd &
                           comp[,"height"] > min_height &
                           comp[,"area"] > min_area &
                           comp[,"rt"] > as.numeric(min_rt) &
                           comp[,"rt"] < as.numeric(max_rt)), , drop = FALSE]))
  att <- attributes(peak_list)
  structure(result,
            chrom_list = att$chrom_list,
            lambdas = att$lambdas, fit=att$fit, sd.max=att$sd.max,
            max.iter=att$max.iter,
            class="peak_list")
}
