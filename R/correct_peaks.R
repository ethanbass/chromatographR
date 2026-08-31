#' Correct peak positions according to a PTW warping model
#' 
#' Corrects retention time differences in `peak_list` using parametric time
#' warping as implemented in the [`ptw`][ptw::ptw] package.
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in a separate column
#' in the list of peak tables.
#' 
#' @importFrom stats predict
#' @param peak_list A `peak_list` object created by [`get_peaks`], containing a 
#' nested list of peak tables where the first level is the sample,
#' and the second level is the spectral wavelength. Every wavelength is described
#' by a matrix where each row corresponds to a feature, and the columns contain 
#' information on that feature (e.g., retention time, peak width (FWHM), height,
#' area, etc.)
#' @param mod_list A list of `ptw` models.
#' @param chrom_list List of chromatograms from which the `ptw` models are
#' derived.
#' @param match_names Logical. Whether to actively match the names of the 
#' `peak_list` to the list of models (`mod_list`). Defaults to `TRUE`.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention times.
#' @author Ron Wehrens, Ethan Bass
#' @note This function is adapted from the `correctPeaks` function in the alsace
#' package by Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/correctPeaks.R>.
#' @seealso [`correct_rt`]
#' @examplesIf interactive()
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, lambdas = 210, what = "models")
#' pks <- get_peaks(Sa_pr, lambda = "210")
#' correct_peaks(pks, mod_list = warp, chrom_list = Sa_pr)
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
  corrected_pks <- transfer_metadata(corrected_pks, peak_list, 
                                     transfer_class = TRUE)
  corrected_pks
}

#' Predict PTW
#' 
#' Predict method for objects of class `ptw`. Reproduced from the [`ptw`]
#' package as the original is not exported.
#'
#' @param object an object of class `ptw`.
#' @param newdata an optional matrix of new data to predict. If missing,
#'   returns the warped sample.
#' @param what character string specifying whether to return the warped
#'   `"response"` or corrected `"time"`. Defaults to `"response"`.
#' @param RTref optional numeric vector of reference retention times.
#' @param ... additional arguments (currently ignored).
#'
#' @return a matrix of warped responses or corrected retention times.
#'
#' @author Ron Wehrens
#' @keywords internal
#' @exportS3Method
predict.ptw <- function (object, newdata, what = c("response", "time"), 
                         RTref = NULL, 
                         ...){
  what <- match.arg(what)
  switch(what, response = {
    if (missing(newdata)) return(object$warped.sample)
    if (!is.matrix(newdata)) newdata <- matrix(newdata, nrow = 1)
    if (object$warp.type == "individual" & nrow(newdata) > 
        1 & nrow(newdata) != nrow(object$warp.fun))
      stop("Incorrect number of rows in newdata")
    if (object$warp.type == "individual") {
      WF <- object$warp.fun
    } else{
      WF <- matrix(object$warp.fun, nrow(object$sample), 
                   ncol(object$warp.fun), byrow = TRUE)
    }
    if (object$mode == "backward") {
      t(sapply(seq_len(nrow(newdata)), function(i){
        approx(x = seq_len(ncol(newdata)), y = newdata[i, ], xout = WF[i, ])$y
        }))
    } else{
      t(sapply(seq_len(nrow(newdata)), function(i){
        approx(x = WF[i,], y = newdata[i, ], xout = seq_len(ncol(newdata)))$y
      }))
    }
  }, time = {
    correctedTime <- switch(object$mode, 
                            backward = -sweep(object$warp.fun, 2, 2 * 
                                                (seq_len(ncol(object$ref))), 
                                              FUN = "-"),
                            object$warp.fun)
    if (is.null(RTref)) {
      if (is.null(colnames(object$ref))) {
        RTref <- seq_len(ncol(object$ref))
      } else {
        RTref <- as.numeric(colnames(object$ref))
      }
    }
    if (missing(newdata)){
      newdata <- RTref
      newdataIndices <- seq_len(length(RTref))
    } else{
      newdataIndices <- round((newdata - min(RTref)) * 
                                (length(RTref) - 1)/diff(range(RTref)) + 1)
    }
    t(sapply(seq_len(nrow(correctedTime)), function(i){
      approx(RTref, NULL, correctedTime[i, newdataIndices])$y
      }))
  })
}
