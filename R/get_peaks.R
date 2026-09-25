#' Get peak list
#' 
#' Detects chromatographic peaks, fits peak models, and extracts peak parameters
#' at one or more wavelengths.
#' 
#' Peaks are detected by [`find_peaks`] as zero-crossings in the smoothed first
#' derivative of each trace, then fitted by [`fit_peaks`] according to `fit`.
#' Peak areas are estimated by trapezoidal approximation.
#'
#' Pass additional arguments to `find_peaks` via `...` to adjust detection.
#' Increase `smooth_window` to reduce splitting of broad peaks; too wide a
#' window will merge closely-eluting peaks.
#'
#' `sd`, `FWHM`, `tau`, `tau_right`, `tau_left`, and `area` are in the units
#' set by `time_unit`, which defaults to minutes. To compare directly with
#' ChemStation integration results, set `time_unit = "s"`.
#'
#' `rt`, `start`, and `end` are *not* affected by `time_unit`: they are read
#' off the chromatogram's own time axis and are always in its units. With
#' `time_unit = "s"` against a chromatogram indexed in minutes, `rt` is in
#' minutes while `FWHM` is in seconds. `max_fwhm` is compared against `FWHM`
#' and is in the same units.
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of chromatograms in matrix format (timepoints ×
#' wavelengths).
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' find peaks at. Omit for one-dimensional chromatograms.
#' @param fit Peak model. Options are bidirectional exponentially modified
#' Gaussian (`"bemg"`), exponential-Gaussian hybrid (`"egh"`), `"gaussian"`,
#' or `"raw"`. `"raw"` integrates the signal directly by trapezoidal
#' approximation without fitting a model. Defaults to `"egh"`.
#' @param windows Optional two-column matrix or `data.frame` of start and end
#' retention times (or a single pair of times), giving manual integration
#' windows that are applied to every chromatogram. If provided, peak detection
#' is skipped entirely and one peak is returned for every window in every
#' sample, whether or not a peak is detectable there. This lets you
#' quantify a compound in samples where it is not detected, and guarantees that
#' every sample contributes to the same set of peak table columns. Two adjacent
#' windows sharing a boundary time act as a vertical drop line: their areas sum
#' exactly to the area of their union, which is how a peak is split at a chosen
#' time. `max_fwhm` is ignored, since the width of each window is chosen by
#' the user. Defaults to `NULL`.
#' @param max_fwhm Maximum peak width, as a full width at half maximum, in the
#' units set by `time_unit`. Peaks wider than `max_fwhm` are discarded.
#' Expressed in the same units as `max_fwhm` in [`filter_peaks`], so the same
#' value transfers between the two, and, unlike a standard deviation, it is
#' comparable across values of `fit`. Peaks whose `FWHM` could not be measured
#' (`NA`, see the `FWHM` entry under *Value*) are retained. Defaults to
#' `NULL`, which disables width filtering: a "typical" peak width varies too
#' much across instruments and methods for a fixed cutoff to be safe as a
#' default. See the section on calibrating these parameters from your own
#' data in `vignette("chromatographR")` for how to choose a value.
#' Ignored if `windows` is provided. A warning fires if filtering discards an
#' unusually large share of detected peaks.
#' @param sd_max The `sd_max` argument is deprecated, because it is expressed as
#' a standard deviation in timepoints, which is a peak width only for
#' `fit = "gaussian"`. It is also not expressed in the same units as
#' `max_fwhm` in [`filter_peaks`]. Please use `max_fwhm` instead. If supplied,
#' it is converted to a `FWHM` cutoff of `2.355 * sd_max * dt` in the units of
#' `time_unit`, where `dt` is the chromatogram's time step. Pass `sd_max = NULL`
#' to disable width filtering. Ignored if `windows` is provided.
#' @param sd.max The `sd.max` argument is deprecated. Please use `sd_max` instead.
#' @param max_iter Maximum number of iterations for non-linear least squares
#' in [`fit_peaks`]. Defaults to `1000`.
#' @param max.iter The `max.iter` argument is deprecated. Please use `max_iter` 
#' instead.
#' @param time_unit Specifies the units for `sd`, `FWHM`, `area`, and `tau` (if
#' applicable). Options are minutes (`"min"`), seconds (`"s"`), or 
#' milliseconds (`"ms"`). Defaults to `"min"`.
#' @param time.units The `time.units` argument is deprecated. Please use 
#' `time_unit` instead.
#' @param estimate_purity Logical. Whether to estimate purity. Defaults to
#' `FALSE`. When `TRUE`, calls [`get_purity`] for every peak in every sample,
#' at a cost that grows with the number of wavelengths.
#' @param noise_threshold Noise threshold. Input to `get_purity`.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed.
#' @param cl Number of cores (integer) or a cluster object from
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows,
#' integers are ignored.
#' @param baseline Whether to subtract a local baseline beneath each peak.
#' Options are `"none"` (no baseline, default), `"flat"` (constant offset per
#' peak), or `"sloped"` (linearly varying offset per peak, useful for gradient
#' elution). For the peak models, the offset is estimated by non-linear least
#' squares along with the shape parameters. When `fit = "raw"`, it is instead
#' read directly off the boundaries of the integration window: `"flat"` draws a
#' horizontal line at the lower of the two boundary intensities, while
#' `"sloped"` draws a line between them (valley-to-valley integration).
#' @param collapse Logical. Whether to collapse multiple peak lists per sample
#' into a single list when multiple wavelengths (`lambdas`) are provided.
#' @param ... Additional arguments passed to [`find_peaks`]. Increase
#' `smooth_window` if peaks are split across multiple detections. Other
#' useful arguments: `smooth_type`, `min_curvature`, `min_amplitude`.
#' @return An S3 object of class `peak_list`: a nested list
#' (sample → wavelength) of `data.frame`s, one row per peak, with columns:
#' * `rt`: The retention time of the peak maximum.
#' * `start`: The retention time where the peak is estimated to begin.
#' * `end`: The retention time where the peak is estimated to end.
#' * `sd`: The standard deviation (\eqn{\sigma}) of the Gaussian component of
#' the fitted peak shape. This is a parameter of the model rather than the width
#' of the peak: for `"egh"` and `"bemg"` the asymmetry of the peak is carried by
#' `tau` or by `tau_left`/`tau_right`, so \eqn{\sigma} describes only the core
#' of the peak and is not comparable across values of `fit`. Use `FWHM` to
#' compare widths. `NA` when `fit = "raw"`, which fits no model; the width of
#' the integration window is available as `end - start`.
#' * `tau`: Empirical skewness parameter (in units of time) determining peak
#' asymmetry when `fit == "egh"`.
#' * `tau_right`: Exponential rate constant controlling right-sided tailing
#' when `fit = "bemg"`. Note that unlike `tau` in the EGH model,
#' this parameter has units of 1/time.
#' * `tau_left`: Exponential rate constant controlling left-sided tailing
#' when `fit = "bemg"`. Note that unlike `tau` in the EGH model,
#' this parameter has units of 1/time.
#' * `h`: The amplitude parameter of the fitted BEMG shape function (only when
#' `fit = "bemg"`), i.e. the value of the pure shape at `center` before any
#' baseline is added. Not the same quantity as the `height` column.
#' * `center`: The center parameter of the fitted BEMG shape function (only
#' when `fit = "bemg"`). This is the location around which the exponential
#' tails are applied. Because BEMG peaks can be asymmetric
#' (`tau_right ≠ tau_left`) and a sloped baseline shifts the fitted curve,
#' the actual peak maximum (`rt`) may differ from `center`.
#' * `FWHM`: The full width at half maximum of the peak. For `"gaussian"` this
#' is \eqn{2\sqrt{2\ln 2}\,\sigma \approx 2.355 \sigma}; for `"egh"` it is the
#' closed form \eqn{\sqrt{(\tau\ln 2)^2 + 8\sigma^2\ln 2}} implied by the model
#' (Lan & Jorgenson, 2001); and for `"bemg"`, which has no closed form, it is
#' measured numerically from the fitted shape. For `"raw"` it is measured
#' directly from the signal above the local `baseline`, and is `NA` when the
#' integration window does not fall below half the height of the peak on both
#' sides, which happens for peaks that are fused with a neighbor. Unlike `sd`,
#' `FWHM` is comparable across values of `fit`.
#' * `height`: The height of the peak above the fitted baseline.
#' * `area`: The area of the peak above the fitted baseline, estimated by
#' trapezoidal approximation.
#' * `r.squared`: Coefficient of determination (\eqn{R^2}), calculated by
#' linear regression of fitted values to the raw data. Statistically informal
#' (the model was fit by nonlinear least squares, not OLS), but useful as a
#' fit diagnostic.
#' * `purity`: The peak purity as estimated by [`get_purity`].
#' * `window_limited`: `1` if a fitted shape parameter came to rest on a bound
#' derived from the integration window, meaning the window rather than the data
#' determined the reported shape, and `0` otherwise. Not returned when
#' `fit = "raw"`.
#' * `converged`: `1` if the nonlinear least squares optimizer reported
#' convergence for the final stage of the fit, and `0` otherwise. A `0` means
#' the reported parameters are wherever the optimizer stopped rather than a
#' converged fit, so treat the corresponding `height` and `area` with caution.
#' It also covers the case where no optimizer ran at all: a fitting stage whose
#' integration window holds fewer points than the stage has parameters is
#' skipped, since such a window cannot determine the model, and the reported
#' parameters are then the starting estimates rather than a measurement. Both
#' are common in peak lists that contain many narrow, low windows.
#' Not returned when `fit = "raw"`.
#' * `floor`: Constant baseline level (only when `baseline = "flat"`).
#' * `floor_start`: Baseline level at the left edge of the peak window (only
#' when `baseline = "sloped"`).
#' * `floor_end`: Baseline level at the right edge of the peak window (only
#' when `baseline = "sloped"`). Both `floor_start` and `floor_end` are
#' constrained to be non-negative.
#' @author Ethan Bass
#' @note Adapted from `getAllPeaks` by Ron Wehrens
#' (<https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R>); however,
#' the underlying algorithms for finding and integrating peaks differ.
#' @references 
#' * Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
#' functions as a simple model of asymmetric chromatographic peaks. *Journal of
#' Chromatography A* 915:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
#'
#' * Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions 
#' — A good model for chromatographic peaks in isocratic HPLC? *Chromatographia*,
#' 26: 285-296. \doi{10.1007/BF02268168}.
#' 
#' * Arase, Shuntaro, Kanta Horie, Takashi Kato, et al. 2016. Intelligent Peak 
#' Deconvolution through In-Depth Study of the Data Matrix from Liquid 
#' Chromatography Coupled with a Photo-Diode Array Detector Applied to 
#' Pharmaceutical Analysis. *Journal of Chromatography. A*, 1469: 35–47. 
#' \doi{doi.org/10.1016/j.chroma.2016.09.037}.
#'
#' * O'Haver, Tom. Pragmatic Introduction to Signal Processing:
#' Applications in scientific measurement.
#' <https://terpconnect.umd.edu/~toh/spectrum/> (Accessed January, 2022).
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. 
#' *Metabolomics* 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examples
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'), max_fwhm = 1, fit = "egh")
#' @seealso [`find_peaks`], [`fit_peaks`], [`print.peak_list`], 
#' \code{\link{[.peak_list}}
#' @export get_peaks
#' @md

get_peaks <- function(chrom_list, lambdas,
                      fit = c("egh", "bemg", "gaussian", "raw"),
                      windows = NULL,
                      max_fwhm = NULL, sd_max = NULL, max_iter = 1000,
                      time_unit = c("min", "s", "ms"),
                      baseline = c("none", "flat", "sloped"),
                      estimate_purity = FALSE,  noise_threshold = .001,
                      show_progress = NULL, cl = 2, collapse = FALSE,
                      time.units = NULL, sd.max = NULL, max.iter = NULL, ...){
  check_duplicated_names(names(chrom_list))
  max_iter <- resolve_deprecated(max.iter, max_iter)
  # `sd_max = NULL` was the documented way to switch width filtering off, so an
  # explicitly supplied `NULL` has to stay distinguishable from the default.
  sd_max_supplied <- !missing(sd_max) || !is.null(sd.max)
  sd_max <- resolve_deprecated(sd.max, sd_max)
  if (sd_max_supplied){
    if (!missing(max_fwhm))
      stop("Supply either `max_fwhm` or `sd_max`, not both. `max_fwhm` is a ",
           "full width at half maximum in the time units of the chromatogram, ",
           "whereas `sd_max` is a standard deviation in timepoints.")
    if (!is.null(sd_max))
      warning("The `sd_max` argument is deprecated. It is expressed as a ",
              "standard deviation in timepoints, which is not a peak width ",
              "for any model except `\"gaussian\"`. Please use `max_fwhm` ",
              "instead, which is expressed as a full width at half maximum in ",
              "the time units of the chromatogram, like `max_fwhm` in ",
              "`filter_peaks`.", immediate. = TRUE)
  }
  time_unit <- resolve_deprecated(time.units, time_unit)
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- time_factor(time_unit)
  fit_is_default <- missing(fit)
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  baseline <- match.arg(baseline, c("none", "flat", "sloped"))
  chrom_list_string <- deparse(substitute(chrom_list))
  if (class(chrom_list)[1] == "matrix")
    chrom_list <- list(chrom_list)
  if (missing(lambdas)){
    if (ncol(chrom_list[[1]]) == 1){
      lambdas <- colnames(chrom_list[[1]])
    } else stop("Wavelengths (`lambdas`) must be provided.")
  }
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  if (is.null(names(chrom_list))){
    warning("Sample names not found. It is recommended to include names for your samples.",
            immediate. = TRUE)
    names(chrom_list) <- seq_along(chrom_list)
  }
  if (!is.null(windows)){
    windows <- check_windows(windows, times = get_times(chrom_list))
    if (fit_is_default){
      message("`windows` was supplied with the default `fit = \"egh\"`, so a peak model will be fit within each window. Set `fit = \"raw\"` to integrate the raw signal against a zero baseline instead.")
    }
    if (estimate_purity){
      warning("Spectral purity is not meaningful across a manual integration window spanning more than one component.",
              immediate. = TRUE)
    }
  }
  peaks <- list()
  laplee <- choose_apply_fnc(show_progress, cl = cl)
  
  fit_sample <- function(sample){
    # `timepoints` and `tdiff` are properties of the chromatogram, not of the
    # wavelength, so they are derived once per sample
    timepoints <- get_times(chrom_list, idx = sample)
    tdiff <- get_time_resolution(chrom_list, idx = sample)
    cutoff <- if (sd_max_supplied){
      if (is.null(sd_max)) NULL else FWHM_CONST * sd_max * tdiff * tfac
    } else max_fwhm
    suppressWarnings(ptable <- lapply(lambdas, function(lambda){
      pks <- fit_peaks(chrom_list[[sample]], lambda = lambda, fit = fit,
                       windows = windows,
                       max_iter = max_iter, max_fwhm = NULL,
                       estimate_purity = estimate_purity,
                       noise_threshold = noise_threshold,
                       baseline = baseline, ...)
      if (nrow(pks) == 0){
        # No peaks in this sample/wavelength.
        # Return an empty table with the full schema so downstream `rbind`s and
        # column lookups behave as they do for non-empty samples.
        return(cbind(sample = character(0), lambda = character(0), pks))
      }
      pks <- cbind(sample = names(chrom_list)[sample], lambda, pks)
      pks <- remove_bad_peaks(pks)
      pks <- convert_indices_to_times(pks, timepoints = timepoints,
                                      tdiff = tdiff, tfac = tfac)
      if (is.null(windows)){
        n_width_before <- nrow(pks)
        pks <- pks[keep_by_width(pks$FWHM, cutoff), , drop = FALSE]
        attr(pks, "n_width_before") <- n_width_before
        attr(pks, "n_width_after") <- nrow(pks)
      }
      pks
    }))
    names(ptable) <- lambdas
    n_width_before <- sum(vapply(ptable, function(x) attr(x, "n_width_before") %||% 0L, numeric(1)))
    n_width_after <- sum(vapply(ptable, function(x) attr(x, "n_width_after") %||% 0L, numeric(1)))
    if (collapse){
      ptable <- do.call(rbind, ptable)
    }
    attr(ptable, "n_width_before") <- n_width_before
    attr(ptable, "n_width_after") <- n_width_after
    ptable <- transfer_metadata(ptable, chrom_list[[sample]],
                                exclude = c('names','row.names','dim',
                                              'dimnames', "format_out", 
                                              "data_format", "parser"))
    ptable
  }
  result <- laplee(seq_along(chrom_list), fit_sample)
  result <- refit_lost_samples(result, fit_sample, laplee, names(chrom_list))
  names(result) <- names(chrom_list)
  # `result` is either a data.frame per sample (when `collapse`) or a list of
  # data.frames, one per wavelength
  collect <- function(var){
    unlist(lapply(result, function(s){
      if (is.data.frame(s)) s[[var]] else lapply(s, function(d) d[[var]])
    }), use.names = FALSE)
  }
  if (!is.null(windows)){
    n_expected <- nrow(windows) * length(lambdas)
    dropped <- vapply(result, function(s){
      n_expected - ifelse(is.data.frame(s), nrow(s), sum(vapply(s, nrow, numeric(1))))
    }, numeric(1))
    if (any(dropped > 0)){
      warning(sprintf(
        "Some windows were dropped because their apex fell on the edge of the chromatogram, in the following sample%s: %s.",
        plural(sum(dropped > 0)),
        paste(names(dropped)[dropped > 0], collapse = ", ")), immediate. = TRUE)
    }
    n_negative <- sum(collect("area") < 0, na.rm = TRUE)
    if (n_negative > 0){
      warning(sprintf(
        "The local baseline lies above the signal in %d window%s, giving a negative area.\n     Negative areas are reported as such, but are coerced to zero by `get_peaktable`.",
        n_negative, plural(n_negative)), immediate. = TRUE)
    }
  }

  if (is.null(windows)){
    warn_if_width_filtered(result)
    warn_if_underdetermined(result, chrom_list, fit, baseline)
  }
  structure(result,
            chrom_list = chrom_list_string,
            lambdas = lambdas, fit = fit, sd_max = sd_max, max_fwhm = max_fwhm,
            max_iter = max_iter, windows = windows,
            time_unit = time_unit,
            baseline = baseline,
            intensity_unit = get_metadata_attribute(chrom_list[[1]], "detector_y_unit"),
            class = c("peak_list", "list"))
}

#' Refit samples whose worker process died
#'
#' A forked worker that dies takes its whole chunk of the peak list with it, and
#' `mclapply` reports the loss only as `NULL` elements, which would otherwise be
#' returned as samples containing no peaks. Losses are retried through the same
#' apply function, so a repeat crash stays inside a worker rather than reaching
#' the caller, and anything still missing is an error rather than silently
#' absent data.
#'
#' @param result The list returned by the parallel apply.
#' @param fit_sample The per-sample function that produced it.
#' @param laplee The apply function in use.
#' @param nms Sample names, for reporting.
#' @param tries Maximum number of refitting passes.
#' @return `result`, with lost elements refitted.
#' @author Ethan Bass
#' @noRd
refit_lost_samples <- function(result, fit_sample, laplee, nms, tries = 2L){
  lost <- function(x) vapply(x, function(e) is.null(e) || inherits(e, "try-error"),
                             logical(1))
  for (i in seq_len(tries)){
    bad <- which(lost(result))
    if (length(bad) == 0) return(result)
    warning(sprintf(
      "%d sample%s returned no result, which usually means a worker process died. Refitting.",
      length(bad), plural(length(bad))), immediate. = TRUE, call. = FALSE)
    result[bad] <- laplee(bad, fit_sample)
  }
  bad <- which(lost(result))
  if (length(bad) > 0){
    stop(sprintf(
      "%d of %d samples could not be fitted because their worker process died repeatedly: %s.\n  This is a crash inside the optimizer rather than a problem with the data. Try `fit = \"egh\"`, or `cl = 1` to see the underlying error.",
      length(bad), length(result),
      paste(head(nms[bad], 10), collapse = ", ")), call. = FALSE)
  }
  result
}

#' Report a peak list dominated by width-filtered peaks
#'
#' `max_fwhm`/`sd_max` discard peaks silently, which is invisible in the
#' returned peak list unless the count of detected peaks is compared against
#' the count before filtering. A large discarded share usually means the
#' cutoff was set too tight for this data, so it's surfaced as a warning
#' rather than left for the user to notice from a suspiciously short peak
#' list.
#'
#' @param result The assembled per-sample peak tables, carrying
#' `n_width_before`/`n_width_after` attributes set by `get_peaks`.
#' @return `NULL`, invisibly; called for its warning.
#' @author Ethan Bass
#' @noRd
warn_if_width_filtered <- function(result){
  n_before <- sum(vapply(result, function(s) attr(s, "n_width_before") %||% 0L, numeric(1)))
  n_after <- sum(vapply(result, function(s) attr(s, "n_width_after") %||% 0L, numeric(1)))
  if (n_before == 0) return(invisible(NULL))
  frac <- 1 - n_after / n_before
  if (frac > WIDTH_FILTER_WARN_FRAC){
    warning(sprintf(
      "%d of %d detected peaks (%.0f%%) were discarded for exceeding `max_fwhm`.\n     If this is unexpected, `max_fwhm` may be set too tight for this data; see `?get_peaks`.",
      n_before - n_after, n_before, 100 * frac), immediate. = TRUE, call. = FALSE)
  }
  invisible(NULL)
}

#' Report a peak list dominated by unfittable windows
#'
#' A window holding fewer timepoints than the model has parameters cannot
#' determine it, so the stage is skipped and the reported shape is a starting
#' estimate. A few such peaks are unremarkable, but a large share of them means
#' `find_peaks` is detecting noise, which is a problem upstream of the fit and
#' invisible in the peak list itself. Skipped when `windows` is supplied, since
#' the user chose those windows.
#'
#' @param result The assembled per-sample peak tables.
#' @param chrom_list The chromatograms, for their timepoints.
#' @param fit The peak model.
#' @param baseline The baseline model.
#' @return `NULL`, invisibly; called for its warning.
#' @author Ethan Bass
#' @noRd
warn_if_underdetermined <- function(result, chrom_list, fit, baseline){
  p <- n_fit_params(fit, baseline)
  if (is.na(p)) return(invisible(NULL))
  n_bad <- 0L
  n_all <- 0L
  for (i in seq_along(result)){
    s <- result[[i]]
    dfs <- if (is.data.frame(s)) list(s) else s
    tdiff <- get_time_resolution(chrom_list, idx = i)
    if (!is.finite(tdiff) || tdiff <= 0) next
    for (df in dfs){
      if (is.null(df) || nrow(df) == 0) next
      npts <- (df$end - df$start) / tdiff + 1
      n_bad <- n_bad + sum(npts < p, na.rm = TRUE)
      n_all <- n_all + nrow(df)
    }
  }
  if (n_all == 0) return(invisible(NULL))
  frac <- n_bad / n_all
  if (frac > UNDERDETERMINED_WARN_FRAC){
    warning(sprintf(
      "%d of %d detected peaks (%.0f%%) sit in an integration window with fewer timepoints than the `%s` model has parameters (%d), so their reported shape is a starting estimate rather than a fit.\n     This usually means noise is being detected as peaks. Consider raising `min_prominence_snr` or widening `smooth_window`, which are passed through to `find_peaks`.",
      n_bad, n_all, 100 * frac, fit, p), immediate. = TRUE, call. = FALSE)
  }
  invisible(NULL)
}

#' Remove bad peaks
#' This function is called internally by `get_peaks`.
#' @author Ethan Bass
#' @noRd
remove_bad_peaks <- function(pks){
  pks[which(apply(pks, 1, function(x)!all(is.na(x))) & 
              apply(pks[,c("rt","start","end")], 1, function(x)all(!is.na(x))) &
              pks[,"rt"] >= 1 
              # pks[,"rt"] > pks[,"start"] &
              # pks[,"rt"] < pks[,"end"]
            ), , drop = FALSE]
}

#' Convert indices to times
#' This function is called internally by `get_peaks`.
#' @author Ethan Bass
#' @noRd
convert_indices_to_times <- function(x, timepoints, tdiff, tfac){
  idx_cols <- intersect(c('rt', 'start', 'end', 'center'), colnames(x))
  x[, idx_cols] <- sapply(idx_cols, function(j){
    i <- x[, j]
    i[i < 1 | i > length(timepoints)] <- NA_integer_
    timepoints[i]
  })
  x[, c('sd', 'FWHM', 'area')] <- x[, c('sd', 'FWHM', 'area')] * tdiff * tfac
  if (!is.null(x$tau)){
    x[, c('tau')] <- x[, c('tau')] * tdiff * tfac
  } 
  if (!is.null(x$tau_right)){
    x[, c('tau_right', 'tau_left')] <- x[, c('tau_right', 'tau_left')] / (tdiff * tfac)
  } 
  x
}

#' Subset peak list
#'
#' Subsets a `peak_list` while preserving metadata attributes.
#'
#' @param x A `peak_list` object.
#' @param i Indices specifying elements to extract.
#' @param j Not used.
#' @param ... Additional arguments (ignored).
#' @param drop Not used.
#' @return A `peak_list` object.
#' @author Ethan Bass
#' @method [ peak_list
#' @export
`[.peak_list` <- function(x, i, j, ..., drop = FALSE) {
  # `exclude` is given explicitly because a `peak_list`'s `baseline` attribute
  # is the local-baseline *mode*, not the estimated baseline that
  # `transfer_metadata` drops by default
  transfer_metadata(.subset(x, i), x,
                    exclude = c('names', 'row.names', 'dim', 'dimnames'))
}

#' Print peak list
#'
#' Prints a summary of a `peak_list` object, including the number of
#' samples, wavelengths, fit method, and total number of peaks.
#'
#' @param x A `peak_list` object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @author Ethan Bass
#' @method print peak_list
#' @export
print.peak_list <- function(x, ...) {
  meta <- attributes(x)
  n_samples <- length(x)
  lambdas <- paste(meta$lambdas, collapse = ", ")
  n_peaks <- sum(sapply(x, function(s) sum(sapply(s, nrow))))
  mean_peaks <- round(n_peaks / n_samples)
  
  cat(sprintf("A peak_list with %d samples and %d wavelength(s) (%s nm)\n",
              n_samples, length(meta$lambdas), lambdas))
  cat(sprintf("Fit method: %s | Time unit: %s | max_fwhm: %s\n",
              meta$fit, meta$time_unit,
              if (!is.null(meta$max_fwhm)) meta$max_fwhm
              else if (!is.null(meta$sd_max))
                sprintf("%s * %s timepoints (from sd_max)",
                        signif(FWHM_CONST, 4), meta$sd_max)
              else "none"))
  cat(sprintf("Total peaks: %d (mean %s per sample)\n",
              n_peaks, mean_peaks))
  if (!is.na(meta$chrom_list))
    cat(sprintf("Source: %s\n", meta$chrom_list))
  invisible(x)
}

