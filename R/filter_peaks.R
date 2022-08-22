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
#' @param max_rt Maximum retention time.
#' @return A peak list similar to the input, with all rows removed
#' from that do not satisfy the specified criteria.
#' @author Ron Wehrens, Ethan Bass
#' @seealso \code{\link{get_peaks}}, \code{\link{filter_peaktable}}
#' @export filter_peaks
filter_peaks <- function(peak_list, min_height, min_area,
                         min_sd, max_sd, min_rt, max_rt){
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

#' Filter peak table
#' 
#' Utility function to remove peaks from peak table, e.g. because their
#' intensity is too low. Currently one can filter on mean or median peak intensity,
#' or retention time.
#' 
#' @param peak_table A peak_table object from \code{\link{get_peaktable}}.
#' @param rts Vector of retention times to include in the peak table.
#' @param min_rt Minimum retention time to include in the peak table.
#' @param max_rt Maximum retention time to include in the peak table.
#' @param min_value Minimal cutoff for average peak intensity.
#' @param what Whether to average intensities using \code{mean} or \code{median}.
#' @param comp Component(s) to include in peak table (e.g. wavelengths if you
#' are using HPLC-DAD/UV).
#' @param tol Tolerance for matching of retention times to \code{rts}.
#' @return A peak table similar to the input, with all columns removed
#' from the peak table that do not satisfy the specified criteria.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}}, \code{\link{filter_peaks}}
#' @examples
#' data(pk_tab)
#' pk_tab <- filter_peaktable(pk_tab, min_rt = 10, max_rt = 16)
#' @export filter_peaktable
filter_peaktable <- function(peak_table, rts, min_rt, max_rt, min_value, comp,
                              what = c("median","mean"), tol = 0){
  if (missing(rts) & missing(min_rt) &
      missing(max_rt) & missing(min_value) & missing(comp)) {
    warning("Nothing to filter...")
    return(peak_table)
  }
  what <- match.arg(what, c("median","mean"))
  if (!missing(rts)){
    rts <- as.numeric(rts)
    if (!inherits(rts, c("numeric"))){
      stop("`rts` should be a vector of retention times.")
    }
    idx.rt <- as.numeric(sapply(rts, function(x){
      which(elementwise.all.equal(x, peak_table$pk_meta["rt",],
                                  tolerance = tol, scale = 1))
    }))
    nas <- is.na(idx.rt)
    if (any(nas)){
      warning(paste0("The following retention times were not identified in the peak table: ",
                     paste(rts[nas], collapse = ', '),
              ". \n", "     You can try increasing the tolerance (`tol`) to permit fuzzier matching."),
              immediate. = TRUE)
      idx.rt <- idx.rt[-which(nas)]
    }
  } else if (!missing(min_rt) | !missing(max_rt)){
    if (missing(min_rt))
      min_rt <- 0
    if (missing(max_rt))
      max_rt <- peak_table$pk_meta["rt",] + 1
    idx.rt <- which(peak_table$pk_meta["rt",] > min_rt &
                      peak_table$pk_meta["rt",] < max_rt)
  } else{idx.rt <- seq_along(peak_table$pk_meta)}
  if (!missing(min_value)){
    val <- apply(peak_table$tab, 2, eval(what))
    idx.val <- which(val >= min_value)
  } else (idx.val <- seq_along(peak_table$tab))
  if (!missing(comp)){
    idx.comp <- which(peak_table$pk_meta["component",] %in% comp)
  } else (idx.comp <- seq_along(peak_table$tab))
  idx <- Reduce(intersect, list(idx.rt, idx.val, idx.comp))
  peak_table$tab <- peak_table$tab[,idx, drop = FALSE]
  peak_table$pk_meta <- peak_table$pk_meta[,idx, drop = FALSE]
  if (inherits(peak_table$ref_spectra, c("data.frame", "matrix"))){
    peak_table$ref_spectra <- peak_table$ref_spectra[,idx, drop = FALSE]
  }
  peak_table
}
