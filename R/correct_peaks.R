
#' Correct peak positions according to a ptw warping model
#' 
#' Corrects retention time differences using parametric time warping as 
#' implemented in \code{\link[ptw]{ptw}}.
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in a separate column
#' in the list of peak tables.
#' 
#' @importFrom stats predict
#' @param peak_list A `peak_list` object created by \code{\link{get_peaks}},
#' containing a nested list of peak tables where the first level is the sample,
#' and the second level is the spectral wavelength. Every component is described
#' by a matrix where every row is one peak, and the columns contain information on
#' retention time, peak width (FWHM), peak width, height, and
#' area.
#' @param mod_list A list of ptw models.
#' @param chrom_list List of chromatograms supplied to create ptw models.
#' @param match_names Logical. Whether to actively match the names of the 
#' \code{peak_list} to the list of models (\code{mod_list}). Defaults to 
#' \code{TRUE}.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention time.
#' @author Ron Wehrens, Ethan Bass
#' @note This function is adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/correctPeaks.R}{getPeakTable}
#' function in the alsace package by Ron Wehrens.
#' @seealso \code{\link{correct_rt}}
#' @export correct_peaks

correct_peaks <- function(peak_list, mod_list, chrom_list, match_names = TRUE){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(mod_list)
  }
  ref_times <- get_times(chrom_list, idx = attr(mod_list, "reference"))
  if (match_names){
    mod_list <- mod_list[match(names(peak_list), names(mod_list))]
  } else {
    if (length(peak_list) != length(mod_list)){
      stop("The length of the provided peaklist does not match the length of the model list.")
    }
  }
  if (length(ref_times) != length(mod_list[[1]]$warp.fun)){
    stop("Dimensions of the warping models and chromatograms do not match.")
  }
  corrected_pks <- mapply(function(samp, mod){
    lapply(samp, function(profile){
      if (nrow(profile) > 0){
        profile <- cbind(profile, rt.cor = c(predict.ptw(mod, profile[, "rt"],
                                                         what = "time",
                                                         RTref = ref_times)))
        
        if (all(c("start", "end") %in% colnames(profile))){
          profile <- cbind(profile,
                           start.cor = c(predict.ptw(mod, profile[, "start"],
                                                     what = "time",
                                                     RTref = ref_times)),
                           end.cor = c(predict.ptw(mod, profile[, "end"],
                                                   what = "time",
                                                   RTref = ref_times)))
        }
      } else {
        profile <- cbind(profile, rt.cor = rep(0, 0))
      }
      profile
    }
    )}, peak_list, mod_list, SIMPLIFY = FALSE)
  corrected_pks <- transfer_metadata(corrected_pks, peak_list, transfer_class = TRUE)
  corrected_pks
}

#' Plot PTW alignments
#' @importFrom graphics matplot
#' @param x A \code{ptw_list} object created by \code{\link{correct_rt}}.
#' @param lambdas Which lambdas to plot.
#' @param legend Logical. Whether to label the plots.
#' @param ... Additional arguments to \code{\link[graphics]{matplot}}.
#' @author Ethan Bass
#' @export

plot.ptw_list <- function(x, lambdas, legend = TRUE, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,1), mar=c(2,3,2,3))
  
  all.lambdas <- as.numeric(rownames(x[[1]]$warped.sample))
  ts <- as.numeric(colnames(x[[1]]$sample))
  
  if (missing(lambdas)){
    lambdas <- all.lambdas
  }
  if (any(!(lambdas %in% all.lambdas))){
    stop("Lambdas not found. Please check argument and try again")
  }
  
  lambda.idx <- which(lambdas %in% all.lambdas)
  # plot warped samples
  plot.new()
  plot.window(xlim=c(head(ts,1), tail(ts,1)),
              ylim=c(0, max(sapply(x, function(xx) xx$warped.sample), na.rm = TRUE)*1.2))
  for (i in seq_along(x)){
    matplot(ts, t(x[[i]]$warped.sample[lambda.idx, , drop = FALSE]),
            type = 'l', add = TRUE)
  }
  if (legend){
    legend("topright", legend = "ptw", bty = "n")
  }
  # plot reference
  plot.new()
  plot.window(xlim = c(head(ts,1), tail(ts,1)),
              ylim = c(0, max(x[[i]]$reference, na.rm = TRUE)*1.2))
  for (i in seq_along(x)){
    matplot(ts, t(x[[i]]$sample[lambda.idx, , drop = FALSE]),
            type = 'l', add = TRUE)
  }
  if (legend){
    legend("topright", legend = "queries", bty = "n")
  }
}

#' @note This is the function from the ptw package, reproduced here because it
#' isn't exported from ptw.
#' @noRd
predict.ptw <- function (object, newdata, what = c("response", "time"), RTref = NULL, 
                         ...) 
{
  what <- match.arg(what)
  switch(what, response = {
    if (missing(newdata)) return(object$warped.sample)
    if (!is.matrix(newdata)) newdata <- matrix(newdata, nrow = 1)
    if (object$warp.type == "individual" & nrow(newdata) > 
        1 & nrow(newdata) != nrow(object$warp.fun)) stop("Incorrect number of rows in newdata")
    if (object$warp.type == "individual") {
      WF <- object$warp.fun
    } else {
      WF <- matrix(object$warp.fun, nrow(object$sample), 
                   ncol(object$warp.fun), byrow = TRUE)
    }
    if (object$mode == "backward") {
      t(sapply(seq_len(nrow(newdata)), function(i) approx(x = seq_len(ncol(newdata)), 
                                                          y = newdata[i, ], xout = WF[i, ])$y))
    } else {
      t(sapply(seq_len(nrow(newdata)), function(i) approx(x = WF[i, 
      ], y = newdata[i, ], xout = seq_len(ncol(newdata)))$y))
    }
  }, time = {
    correctedTime <- switch(object$mode, backward = -sweep(object$warp.fun, 
                                                           2, 2 * (seq_len(ncol(object$ref))), FUN = "-"), object$warp.fun)
    if (is.null(RTref)) {
      if (is.null(colnames(object$ref))) {
        RTref <- seq_len(ncol(object$ref))
      } else {
        RTref <- as.numeric(colnames(object$ref))
      }
    }
    if (missing(newdata)) {
      newdata <- RTref
      newdataIndices <- seq_len(length(RTref))
    } else {
      newdataIndices <- round((newdata - min(RTref)) * 
                                (length(RTref) - 1)/diff(range(RTref)) + 1)
    }
    t(sapply(seq_len(nrow(correctedTime)), function(i) approx(RTref, 
                                                              NULL, correctedTime[i, newdataIndices])$y))
  })
}
