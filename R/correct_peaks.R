#' Correct peak positions according to a warping model
#'
#' Corrects retention time differences in `peak_list` using warping models
#' fitted by [`correct_rt`] with any of its algorithms.
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in an `rt.cor`
#' column, alongside `start.cor` and `end.cor` when the peak tables have
#' `start` and `end` columns.
#'
#' The corrected times refer to the aligned time axis, so functions that later
#' read spectra or traces at those times need the aligned chromatograms. The
#' returned peak list stores the warping models for this purpose, and
#' [`get_peaktable`] carries them into the peak table, where they are
#' re-applied to the unaligned chromatograms when needed. To use a particular
#' aligned object instead, pass it as `chrom_list_aligned`.
#' 
#' @importFrom stats predict
#' @param peak_list A `peak_list` object created by [`get_peaks`], containing a 
#' nested list of peak tables where the first level is the sample,
#' and the second level is the spectral wavelength. Every wavelength is described
#' by a matrix where each row corresponds to a feature, and the columns contain 
#' information on that feature (e.g., retention time, peak width (FWHM), height,
#' area, etc.)
#' @param mod_list Warping models, as returned by [`correct_rt`] with
#' `what = "models"`: a `ptw_list`, an `lptw_list`, or a `VPdtw` object.
#' Corrected times are on the time axis of the chromatograms `correct_rt`
#' returns, which for `"vpdtw"` extends beyond the reference at either end.
#' @param chrom_list List of chromatograms from which the models are derived.
#' If omitted, it is retrieved using the name stored in `mod_list`, which only
#' `ptw_list` objects record.
#' @param match_names Logical. If `TRUE` (default), models are matched to the
#' samples of `peak_list` by name when both are named; otherwise by position.
#' @param chrom_list_aligned Optional aligned chromatograms, as returned by
#' [`correct_rt`] with `what = "corrected.values"`. Only its name is recorded,
#' so the object must still exist under that name when it is used.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention times.
#' @author Ron Wehrens, Ethan Bass
#' @note This function is adapted from the `correctPeaks` function in the alsace
#' package by Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/correctPeaks.R>.
#' @seealso [`correct_rt`]
#' @examples
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, lambdas = 210, what = "models")
#' pks <- get_peaks(Sa_pr, lambda = "210")
#' correct_peaks(pks, mod_list = warp, chrom_list = Sa_pr)
#' @export correct_peaks

correct_peaks <- function(peak_list, mod_list, chrom_list, match_names = TRUE,
                          chrom_list_aligned = NULL){
  if (missing(chrom_list)){
    if (!inherits(mod_list, "ptw_list"))
      stop("Please supply `chrom_list` to correct peaks with `lptw` or `vpdtw` models.")
    chrom_list <- get_chrom_list(mod_list)
  }
  if (any(vapply(peak_list, is.data.frame, logical(1))))
    stop("`correct_peaks` does not accept peak lists built with `collapse = TRUE`.")
  to_ref <- warp_maps(mod_list, chrom_list)
  if (match_names){
    to_ref <- to_ref[model_positions(to_ref, peak_list)]
  } else if (length(peak_list) != length(to_ref)){
    stop("The length of the provided peaklist does not match the length of the model list.")
  }
  corrected_pks <- mapply(function(samp, map_rt){
    transfer_metadata(lapply(samp, function(profile){
      old <- colnames(profile) %in% c("rt.cor", "start.cor", "end.cor")
      if (any(old)) profile <- profile[, !old, drop = FALSE]
      if (nrow(profile) > 0){
        profile <- cbind(profile, rt.cor = map_rt(profile[, "rt"]))

        if (all(c("start", "end") %in% colnames(profile))){
          profile <- cbind(profile,
                           start.cor = map_rt(profile[, "start"]),
                           end.cor = map_rt(profile[, "end"]))
        }
      } else {
        profile <- cbind(profile, rt.cor = rep(0, 0))
      }
      profile
    }), samp)}, peak_list, to_ref, SIMPLIFY = FALSE)
  corrected_pks <- transfer_metadata(corrected_pks, peak_list, 
                                     transfer_class = TRUE)
  attr(corrected_pks, "warp_models") <- strip_models(mod_list)
  attr(corrected_pks, "chrom_list_aligned") <- if (!is.null(chrom_list_aligned))
    deparse(substitute(chrom_list_aligned))
  corrected_pks
}

#' Per-sample functions mapping retention times onto the aligned axis
#' @noRd
warp_maps <- function(mod_list, chrom_list){
  if (inherits(mod_list, "VPdtw")){
    to_ref <- vpdtw_to_ref(mod_list, chrom_list, stored_reference(mod_list))
    names(to_ref) <- names(chrom_list)
    return(to_ref)
  }
  ref_times <- get_times(chrom_list, idx = reference_axis_idx(
    stored_reference(mod_list), chrom_list))
  lapply(mod_list, function(mod){
    if (is.null(mod$warp.fun))
      mod$warp.fun <- matrix(ptw::warp.time(seq_along(ref_times),
                                            mod$warp.coef[1, ]), nrow = 1)
    if (length(ref_times) != length(mod$warp.fun))
      stop("Dimensions of the warping models and chromatograms do not match.")
    function(rt) c(predict.ptw(mod, rt, what = "time", RTref = ref_times))
  })
}

#' Keep only what is needed to re-apply warping models
#'
#' `ptw` objects carry the reference, sample and warped sample, which can be
#' larger than the warp itself. `VPdtw` objects carry the query, reference and
#' warped query.
#' @noRd
strip_models <- function(mod_list){
  if (inherits(mod_list, "VPdtw"))
    return(structure(
      list(xVals = mod_list$xVals, shift = mod_list$shift,
           samples = colnames(mod_list$query),
           n_timepoints = nrow(mod_list$query)),
      class = "VPdtw", reference = attr(mod_list, "reference"),
      parameters = attr(mod_list, "parameters")))
  stripped <- lapply(mod_list, function(m) m[intersect(names(m), c("warp.coef", "mode"))])
  transfer_metadata(stripped, mod_list)
}

#' Predict PTW
#' 
#' Predict method for objects of class `ptw`. Reproduced from the ptw package,
#' which does not export the original.
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
    newdataIndices <- if (missing(newdata)) seq_along(RTref) else
      approx(RTref, seq_along(RTref), xout = newdata)$y
    t(sapply(seq_len(nrow(correctedTime)), function(i){
      warped <- approx(seq_len(ncol(correctedTime)), correctedTime[i, ],
                       xout = newdataIndices)$y
      approx(RTref, NULL, warped)$y
      }))
  })
}
