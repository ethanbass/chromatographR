#' Filter peak lists
#' 
#' Utility function for filtering peaks from a `peak_list` object. Peaks can be
#' filtered according to height, area, standard deviation, and/or retention time.
#' 
#' @param peak_list A `peak_list` object, consisting of a nested list of peak
#' tables, where the first level is the sample, and the second level is the 
#' spectral component. Every component is described by a `data.frame` where 
#' every row represents a peak, and the columns contain information on 
#' retention time,  peak shape, height, and area.
#' @param min_height Minimum peak height.
#' @param min_area Minimum peak area.
#' @param min_sd Minimal standard deviation.
#' @param max_sd Maximum standard deviation.
#' @param min_rt Minimum retention time.
#' @param max_rt Maximum retention time.
#' @return A filtered `peak_list` object containing only peaks that satisfy the
#' specified criteria.
#' @author Ron Wehrens, Ethan Bass
#' @seealso [`get_peaks`], [`filter_peaktable`]
#' @examples
#' data(Sa_warp)
#' pks <- get_peaks(Sa_warp[1], lambda = "210")
#' filter_peaks(pks, min_height = 100)
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
      warning("The 'min_height' is less than minimum peak height.")}
  if (!missing(min_area)){
    if (min_area < min(x$area))
      warning("The 'min_area' is less than minimum peak area.")}
  if (!missing(min_sd)){
    if (min_sd < min(x$sd))
      warning("The 'min_sd' is less than minimum peak standard deviation.")}
  if (!missing(max_sd)){
    if (max_sd > max(x$sd))
      warning("The 'max_sd' is greater than maximum peak standard deviation,")}
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
                    comp[which(comp[,"sd"] <= max_sd &
                           comp[, "sd"] >= min_sd &
                           comp[, "height"] >= min_height &
                           comp[, "area"] >= min_area &
                           comp[, "rt"] >= as.numeric(min_rt) &
                           comp[, "rt"] <= as.numeric(max_rt)), , drop = FALSE]))
  att <- attributes(peak_list)
  transfer_metadata(result, peak_list, transfer_class = TRUE)
}

#' Filter peak table
#' 
#' Utility function to filter features (columns) in a `peak_table`. Filtering is
#' applied consistently across all components of the peak table, including the 
#' intensity matrix (`tab`), feature metadata (`pk_meta`), and reference spectra
#' (if present). Filtering can be based on retention time, wavelength, feature
#' intensity summarized across samples (mean, median, or max), and/or feature
#' sparsity (i.e., proportion of zero values across samples).
#' 
#' @param peak_table A peak_table object from [`get_peaktable`].
#' @param rts Vector of retention times used to select features. Values are
#' mapped to the closest retention times in the peak table. If no 
#' match is found within `tol`, the value is ignored and a warning is issued.
#' @param min_rt Minimum retention time for features to be retained.
#' @param max_rt Maximum retention time for features to be retained.
#' @param min_value Minimum threshold for feature intensity summarized across
#' samples (using the method specified by `what`).
#' @param max_zeros Maximum allowed feature sparsity (proportion of zero values 
#' across samples).
#' @param what Method used to summarize feature intensity across samples for
#' filtering. One of `"mean"`, `"median"` (default),  or `"max"`.
#' @param lambda Wavelength(s) used to select features in the peak table.
#' Only features matching the specified wavelength(s) are retained.
#' @param tol Tolerance for matching of retention times to `rts`. A
#' feature is considered a match if the absolute difference is ≤ `tol`.
#' @return A filtered `peak_table` containing only features (columns) that
#' satisfy the specified criteria. The same feature indices are applied to all
#' `peak_table` components.
#' @author Ethan Bass
#' @seealso [`get_peaktable`], [`filter_peaks`]
#' @examples
#' data(pk_tab)
#' pk_tab <- filter_peaktable(pk_tab, min_rt = 10, max_rt = 16)
#' @export filter_peaktable

filter_peaktable <- function(peak_table, rts, min_rt, max_rt, min_value, 
                             max_zeros, lambda,
                             what = c("median", "mean", "max"), tol = 0){
  check_peaktable(peak_table)
  if (missing(rts) & missing(min_rt) &
      missing(max_rt) & missing(min_value) & 
      missing(lambda) & missing(max_zeros)) {
    warning("Nothing to filter...")
    return(peak_table)
  }
  peak_table$pk_meta["rt",] <- as.numeric(peak_table$pk_meta["rt",])
  what <- match.arg(what, c("median", "mean", "max"))
  if (!missing(rts)){
    rts <- as.numeric(rts)
    if (!inherits(rts, c("numeric"))){
      stop("`rts` should be a vector of retention times.")
    }
    idx.rt <- as.numeric(sapply(rts, function(rt){
      diffs <- abs(rt - peak_table$pk_meta["rt", ])
      ifelse(min(diffs) <= tol, which.min(diffs), NA)
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
    idx.rt <- which(peak_table$pk_meta["rt",] >= min_rt &
                      peak_table$pk_meta["rt",] <= max_rt)
  } else{idx.rt <- seq_along(peak_table$pk_meta)}
  if (!missing(min_value)){
    val <- apply(peak_table$tab, 2, eval(what), na.rm = TRUE)
    idx.val <- which(val >= min_value)
  } else (idx.val <- seq_along(peak_table$tab))
  if (!missing(max_zeros)){
    if (!is.numeric(max_zeros) || max_zeros < 0 || max_zeros > 1)
      stop("`max_zeros` must be a numeric value between 0 and 1.")
    prop_zero <- colMeans(peak_table$tab == 0, na.rm = TRUE)
    idx.zeros <- which(prop_zero <= max_zeros)
  } else (idx.zeros <- seq_along(peak_table$tab))
  if (!missing(lambda)){
    idx.lambda <- which(peak_table$pk_meta["lambda",] %in% lambda)
  } else (idx.lambda <- seq_along(peak_table$tab))
  idx <- Reduce(intersect, list(idx.rt, idx.val, idx.lambda, idx.zeros))
  peak_table$tab <- peak_table$tab[,idx, drop = FALSE]
  peak_table$pk_meta <- peak_table$pk_meta[, idx, drop = FALSE]
  if (inherits(peak_table$ref_spectra, c("data.frame", "matrix"))){
    peak_table$ref_spectra <- peak_table$ref_spectra[, idx, drop = FALSE]
  }
  peak_table
}
