#' Get peak list
#' 
#' Detects chromatographic peaks, fits peak models, and extracts peak parameters
#' at one or more wavelengths.
#' 
#' Peaks are detected by [`find_peaks`] as zero-crossings in the smoothed first
#' derivative of each trace, then fitted by [`fit_peaks`] according to `fit`.
#' Peak areas are the trapezoidal integral of the fitted curve (or, for
#' `fit = "raw"`, of the signal) above the local `baseline`.
#' 
#' @aliases get_peaks
#' @importFrom stats median
#' @param chrom_list A list of chromatograms in matrix format (timepoints ×
#' wavelengths).
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' find peaks at. Omit for one-dimensional chromatograms, or when `targets` has a
#' `lambda` column.
#' @param fit Peak model. Options are bidirectional exponentially modified
#' Gaussian (`"bemg"`), exponential-Gaussian hybrid (`"egh"`), `"gaussian"`,
#' or `"raw"`. `"raw"` integrates the signal directly by trapezoidal
#' approximation without fitting a model. Defaults to `"egh"`, or to `"raw"`
#' when `windows` is provided, which allows no other value.
#' @param windows Optional two-column matrix or `data.frame` of start and end
#' retention times (or a single pair of times), giving manual integration
#' windows applied to every chromatogram. If provided, peak detection is
#' skipped and each window is integrated in every sample. Requires
#' `fit = "raw"`, and `max_fwhm` is ignored. Cannot be combined with
#' `targets`. Defaults to `NULL`. See *Manual integration windows*.
#' @param targets Target retention times, as a numeric vector or a table of
#' targets (see [`targeted_analysis`]). If provided, only peaks with an apex
#' falling within `tol` of one of these retention times are fit, and [`get_peaktable`]
#' groups peaks by target instead of clustering them. A `lambda` column takes
#' the place of `lambdas`, and a `tol` column the place of the `tol` argument,
#' which must then be omitted. A target's peak whose model fails is kept, with `height` and `area` set to
#' `NA` and `converged` to `0`, whereas without `targets` such peaks are dropped.
#' Defaults to `NULL`, in which case all detected peaks are fit.
#' @param tol Numeric. Half-width of the window around each value of `targets`, in
#' minutes, whatever the time unit of `chrom_list`. Either a single value applied to
#' every target, or a vector giving a separate tolerance for each element of
#' `targets`. Defaults to `0.1`. `tol` is compared against the apex of the
#' detected peak rather than the returned `rt`, which is the apex of the fitted
#' model, so a returned `rt` may fall outside the window. Where the windows of
#' two targets overlap, a peak in the overlap is considered only for the nearer
#' target. Targets closer together than the drift between samples can
#' therefore be assigned each other's peaks; supply `mod_list` to read them on
#' the aligned time axis. Ignored if `targets` is `NULL`.
#' @param select_by Character. How to proceed when more than one peak falls
#' within the window around a target retention time. `"height"` (the default)
#' retains the most intense peak, `"proximity"` retains the peak closest to the
#' target retention time, and `"none"` retains every peak in the window.
#' Because selection precedes model fitting, `"height"` compares raw apex
#' intensities, so the retained peak does not always have the largest fitted
#' `height`. With `"none"`, [`get_peaktable`] keeps only the peak with the
#' largest `response` for each target in each sample, but every peak in the
#' window contributes to the target's `pk_meta`, including its `rt`. Ignored if
#' `targets` is `NULL`.
#' @param mod_list Optional warping models from [`correct_rt`] (with
#' `what = "models"`). If provided, the times in `targets` are read on the
#' aligned time axis, so they can be copied from the `rt` row of an aligned peak table. Peaks are
#' still detected and fitted on the unaligned data; use [`correct_peaks`]
#' afterwards to add corrected retention times. Requires `targets`.
#' @param max_fwhm Maximum peak width, as a full width at half maximum, in the
#' units set by `time_unit`. Peaks wider than `max_fwhm` are discarded, and a
#' warning is issued if that removes more than 20% of detected peaks. Peaks
#' whose `FWHM` could not be measured (`NA`, see the `FWHM` entry under
#' *Value*) are retained. The same value means the same thing as `max_fwhm` in
#' [`filter_peaks`] and, unlike a standard deviation, is comparable across
#' values of `fit`. Defaults to `NULL`, which disables width filtering. To
#' choose a value, see the section on calibrating these parameters from your
#' own data in `vignette("chromatographR")`. Ignored if `windows` is provided.
#' @param sd_max The `sd_max` argument is deprecated. Please use `max_fwhm`
#' instead. `sd_max` is a standard deviation in timepoints, which is a peak
#' width only for `fit = "gaussian"`. If supplied, it is converted to a `FWHM`
#' cutoff of `2.355 * sd_max * dt` in the units of `time_unit`, where `dt` is
#' the chromatogram's time step. Pass `sd_max = NULL` to disable width
#' filtering. Cannot be combined with `max_fwhm`. Ignored if `windows` is
#' provided.
#' @param sd.max The `sd.max` argument is deprecated. Please use `sd_max` instead.
#' @param max_iter Maximum number of iterations for non-linear least squares
#' in [`fit_peaks`]. Defaults to `1000`.
#' @param max.iter The `max.iter` argument is deprecated. Please use `max_iter` 
#' instead.
#' @param time_unit Units for `sd`, `FWHM`, `tau`, `tau_right`, `tau_left`,
#' `area` and `max_fwhm`. Options are minutes
#' (`"min"`), seconds (`"s"`), or milliseconds (`"ms"`). Defaults to the time
#' unit of `chrom_list` (minutes if none is recorded), so that every column
#' shares one unit. Set `"s"` to compare areas with ChemStation integration
#' results. `rt`, `start` and `end` are not converted: they are always in the
#' units of the chromatogram's own time axis.
#' @param time.units The `time.units` argument is deprecated. Please use 
#' `time_unit` instead.
#' @param estimate_purity Logical. Whether to estimate purity. Defaults to
#' `FALSE`. When `TRUE`, calls [`get_purity`] for every peak in every sample,
#' at a cost that grows with the number of wavelengths.
#' @param purity_args A named list of additional arguments passed to
#' [`get_purity`], e.g. `list(background_correction = "none", metric = "max")`.
#' A `lambdas` element here selects the wavelengths used to compare spectra,
#' not the wavelengths at which peaks are integrated.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed.
#' @param cl Number of cores (integer) or a cluster object from
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows, only a
#' cluster object runs in parallel; an integer runs serially.
#' @param baseline Whether to subtract a local baseline beneath each peak.
#' Options are `"none"` (no baseline, default), `"flat"` (constant offset per
#' peak), or `"sloped"` (linearly varying offset per peak, useful for gradient
#' elution). For the peak models, the offset is estimated by non-linear least
#' squares along with the shape parameters. When `fit = "raw"`, it is instead
#' read directly off the boundaries of the integration window: `"flat"` draws a
#' horizontal line at the lower of the two boundary intensities, while
#' `"sloped"` draws a line between them (valley-to-valley integration).
#' @param collapse Logical. Whether to combine the peaks from every wavelength
#' into a single `data.frame` per sample. Defaults to `FALSE`.
#' @param ... Additional arguments passed to [`find_peaks`]. Increase
#' `smooth_window` if peaks are split across multiple detections, but too wide
#' a window merges closely-eluting peaks. Other useful arguments:
#' `smooth_type`, `min_prominence_snr`, `min_curvature`, `min_amplitude`.
#' @return An S3 object of class `peak_list`: a nested list
#' (sample, then wavelength) of `data.frame`s, or a list of one `data.frame` per
#' sample when `collapse = TRUE`, with one row per peak and columns:
#' * `rt`: The retention time of the peak maximum. For the fitted models it
#' generally falls between timepoints, so match it to the time axis with a
#' tolerance rather than `==` or [match()]; [`get_spectra`] with `what = "rt"`
#' returns the spectrum at the nearest timepoint.
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
#' * `tau_right`: Time constant (in units of time) of the exponential tail on
#' the right side of the peak when `fit = "bemg"`. A larger value gives a
#' longer tail, and a value near 0 means no tail on that side.
#' * `tau_left`: Time constant (in units of time) of the exponential tail on
#' the left side of the peak when `fit = "bemg"`. A larger value gives a
#' longer tail, and a value near 0 means no tail on that side.
#' * `h`: The amplitude parameter of the fitted BEMG shape function (only when
#' `fit = "bemg"`). It scales the shape but is not its value at `center`, so
#' use the `height` column for the height of the peak.
#' * `center`: The center parameter of the fitted BEMG shape function (only
#' when `fit = "bemg"`). This is the location around which the exponential
#' tails are applied. Because BEMG peaks can be asymmetric
#' (`tau_right != tau_left`) and a sloped baseline shifts the fitted curve,
#' the actual peak maximum (`rt`) may differ from `center`.
#' * `FWHM`: The full width at half maximum of the peak. For `"gaussian"` this
#' is \eqn{2\sqrt{2\ln 2}\,\sigma \approx 2.355 \sigma}; for `"egh"` it is the
#' closed form \eqn{\sqrt{(\tau\ln 2)^2 + 8\sigma^2\ln 2}} implied by the model
#' (Lan & Jorgenson, 2001); and for `"bemg"`, which has no closed form, it is
#' measured numerically from the fitted shape. For `"raw"` it is measured
#' directly from the signal above the local `baseline`, and is `NA` when the
#' integration window does not fall below half the height of the peak on both
#' sides, as for a peak fused with a neighbor. Unlike `sd`,
#' `FWHM` is comparable across values of `fit`.
#' * `height`: The height of the peak above the fitted baseline.
#' * `area`: The area of the peak above the fitted baseline, estimated by
#' trapezoidal approximation.
#' * `r.squared`: Coefficient of determination (\eqn{R^2}) of the fitted
#' curve over the peak window, `1 - SS_res / SS_tot`. Negative when the fit is
#' worse than a flat line at the mean. Not returned when `fit = "raw"`.
#' * `purity`: The peak purity as estimated by [`get_purity`]. `NA` unless
#' `estimate_purity = TRUE`.
#' * `window_limited`: `1` if a fitted shape parameter hit a limit set by the
#' integration window, meaning the window rather than the data
#' determined the reported shape, and `0` otherwise. Not returned when
#' `fit = "raw"`.
#' * `converged`: `1` if the nonlinear least squares optimizer reported
#' convergence for the final stage of the fit, and `0` otherwise. A `0` means
#' the reported parameters are wherever the optimizer stopped rather than a
#' converged fit, so treat the corresponding `height` and `area` with caution.
#' It also covers the case where no optimizer ran at all: a fitting stage whose
#' integration window holds fewer points than the stage has parameters is
#' skipped, since such a window cannot determine the model, and the reported
#' parameters are then the starting estimates rather than a measurement.
#' `get_peaks` warns when more than 20% of peaks sit in such a window.
#' Not returned when `fit = "raw"`.
#' * `floor`: Constant baseline level (only when `baseline = "flat"`).
#' * `floor_start`: Baseline level at the left edge of the peak window (only
#' when `baseline = "sloped"`).
#' * `floor_end`: Baseline level at the right edge of the peak window (only
#' when `baseline = "sloped"`). For the peak models, `floor`, `floor_start`
#' and `floor_end` are constrained to be non-negative.
#' * `target`: The row of `targets` the peak was selected for (only when `targets` is
#' supplied).
#' @section Manual integration windows:
#' Manual windows integrate every sample over the same retention-time ranges,
#' so each compound is measured the same way in every sample, and every sample
#' gets the same set of peak table columns. The raw signal in each window is
#' integrated whether or not a peak is detected there, so a window with no peak
#' gets an area near zero. With `baseline = "none"`, two adjacent windows
#' sharing a boundary time act as a vertical drop line: their areas sum exactly
#' to the area of their union, so a peak can be split at a chosen time. The
#' `"flat"` and `"sloped"` baselines are drawn within each window separately, so
#' their areas do not sum to that of the union. To fit peak models at known
#' retention times, use `targets` instead.
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
#' \doi{10.1016/j.chroma.2016.09.037}.
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
#' \code{\link{[.peak_list}}, [`targeted_analysis`]
#' @export get_peaks
#' @md

get_peaks <- function(chrom_list, lambdas,
                      fit = c("egh", "bemg", "gaussian", "raw"),
                      windows = NULL, targets = NULL, tol = 0.1,
                      select_by = c("height", "proximity", "none"),
                      mod_list = NULL,
                      max_fwhm = NULL, sd_max = NULL, max_iter = 1000,
                      time_unit = NULL,
                      baseline = c("none", "flat", "sloped"),
                      estimate_purity = FALSE, purity_args = list(),
                      show_progress = NULL, cl = 2, collapse = FALSE,
                      time.units = NULL, sd.max = NULL, max.iter = NULL, ...){
  check_duplicated_names(names(chrom_list), what = "chrom_list")
  if (estimate_purity) check_purity_args(purity_args)
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
  axis_unit <- get_time_unit(chrom_list, na_value = "min")
  if (!axis_unit %in% c("min", "s", "ms")) axis_unit <- "min"
  time_unit <- match.arg(time_unit %||% axis_unit, c("min", "s", "ms"))
  tfac <- time_factor(time_unit) / time_factor(axis_unit)
  fit_is_default <- missing(fit)
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  baseline <- match.arg(baseline, c("none", "flat", "sloped"))
  chrom_list_string <- deparse(substitute(chrom_list))
  if (class(chrom_list)[1] == "matrix")
    chrom_list <- list(chrom_list)
  if (!missing(tol) && is.data.frame(targets) && "tol" %in% colnames(targets))
    stop("Supply tolerances either in `tol` or in the `tol` column of `targets`, not both.")
  if (!is.null(targets)) targets <- check_targets(targets, tol)
  if (!is.null(targets)) check_bounds_arg(...)
  if (!is.null(targets) && !anyNA(targets$lambda)){
    if (!missing(lambdas))
      stop("Supply wavelengths either in `lambdas` or in the `lambda` column of `targets`, not both.")
    lambdas <- unique(targets$lambda)
  }
  if (missing(lambdas)){
    if (ncol(chrom_list[[1]]) == 1){
      lambdas <- colnames(chrom_list[[1]])
    } else stop("Wavelengths (`lambdas`) must be provided.")
  }
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  if (!is.null(targets) && anyNA(targets$lambda) && any(!is.na(targets$name)) &&
      length(lambdas) > 1){
    stop("Named targets would be duplicated across wavelengths. Add a `lambda` column to `targets` giving the wavelength for each target.")
  }
  if (is.null(names(chrom_list))){
    warning("Sample names not found. It is recommended to include names for your samples.",
            immediate. = TRUE)
    names(chrom_list) <- seq_along(chrom_list)
  }
  select_by <- match.arg(select_by)
  if (!is.null(windows)){
    if (!is.null(targets)){
      stop("`windows` and `targets` cannot be used together. `targets` selects among detected peaks, whereas `windows` integrates fixed windows.")
    }
    windows <- check_windows(windows, times = get_times(chrom_list))
    if (!fit_is_default && fit != "raw"){
      stop("`windows` integrates the raw signal, so `fit` must be \"raw\". To fit peak models at known retention times, use `targets`.",
           call. = FALSE)
    }
    fit <- "raw"
    if (estimate_purity){
      warning("Spectral purity is not meaningful across a manual integration window spanning more than one component.",
              immediate. = TRUE)
    }
  }
  to_ref <- NULL
  if (!is.null(mod_list)){
    if (is.null(targets))
      stop("`mod_list` is only used to place the target retention times in `targets`.")
    to_ref <- warp_maps(mod_list, chrom_list)
    to_ref <- to_ref[model_positions(to_ref, chrom_list)]
  }
  if (!is.null(targets)){
    all_times <- lapply(seq_along(chrom_list), function(i){
      times <- get_times(chrom_list, idx = i)
      if (is.null(to_ref)) times else to_ref[[i]](times)
    })
    rts <- targets$rt
    ts_range <- range(unlist(all_times), na.rm = TRUE)
    oob <- which(rts < ts_range[1] | rts > ts_range[2])
    if (length(oob) > 0){
      stop(sprintf(
        "The following retention times fall outside the bounds of the chromatograms (%s-%s):\n%s",
        signif(ts_range[1], 4), signif(ts_range[2], 4),
        paste(sprintf("         %d (%s)", oob, sQuote(rts[oob])), collapse = ", \n")))
    }
    span <- diff(ts_range) / rt_factor(axis_unit)
    wide <- unique(targets$tol[targets$tol > span / 10])
    if (length(wide) > 0){
      warning(sprintf(
        "`tol` is always in minutes. A tolerance of %s min is more than a tenth of the run (%s min).",
        paste(wide, collapse = ", "), signif(span, 3)), immediate. = TRUE)
    }
  }
  peaks <- list()
  laplee <- choose_apply_fnc(show_progress, cl = cl)
  
  fit_sample <- function(sample){
    # `timepoints`, `tdiff` and the noise variance are properties of the
    # chromatogram, not of the wavelength, so they are derived once per sample
    timepoints <- get_times(chrom_list, idx = sample)
    tdiff <- get_time_resolution(chrom_list, idx = sample)
    noise_variance <- if (estimate_purity)
      hoist_noise_variance(chrom_list[[sample]], purity_args) else NULL
    cutoff <- if (sd_max_supplied){
      if (is.null(sd_max)) NULL else FWHM_CONST * sd_max * tdiff * tfac
    } else max_fwhm
    aligned_times <- if (!is.null(to_ref)) all_times[[sample]]
    n_width <- c(0L, 0L)
    ptable <- lapply(lambdas, function(lambda){
      y <- chrom_list[[sample]][, get_lambda_idx(lambda, get_lambdas(chrom_list[[sample]]),
                                                 allow_max = FALSE)]
      pos <- if (!is.null(windows)){
        windows_to_pos(windows, y = y, times = timepoints)
      } else if (!is.null(targets)){
        idx <- which(is.na(targets$lambda) | targets$lambda == lambda)
        p <- select_peaks_by_rt(find_peaks(y, ...), y = y,
                                times = aligned_times %||% timepoints,
                                rts = targets$rt[idx],
                                tol = targets$tol[idx] * rt_factor(axis_unit),
                                select_by = select_by)
        p$target <- idx[attr(p, "target")]
        p[!duplicated(p[, c("lower", "upper")]), , drop = FALSE]
      }
      pks <- fit_peaks(chrom_list[[sample]], lambda = lambda, fit = fit,
                       pos = pos, max_iter = max_iter, max_fwhm = NULL,
                       estimate_purity = estimate_purity,
                       purity_args = purity_args,
                       noise_variance = noise_variance,
                       baseline = baseline, ...)
      if (nrow(pks) == 0){
        # No peaks in this sample/wavelength (e.g. no peak near any of `targets`).
        # Return an empty table with the full schema so downstream `rbind`s and
        # column lookups behave as they do for non-empty samples.
        return(cbind(sample = character(0), lambda = character(0), pks))
      }
      pks <- cbind(sample = names(chrom_list)[sample], lambda, pks)
      pks <- remove_bad_peaks(pks)
      pks <- convert_indices_to_times(pks, timepoints = timepoints,
                                      tdiff = tdiff, tfac = tfac)
      if (is.null(windows)){
        n_before <- nrow(pks)
        pks <- pks[keep_by_width(pks$FWHM, cutoff), , drop = FALSE]
        n_width <<- n_width + c(n_before, nrow(pks))
      }
      pks
    })
    names(ptable) <- lambdas
    if (collapse){
      ptable <- do.call(rbind, ptable)
    }
    attr(ptable, "n_width") <- n_width
    ptable <- transfer_metadata(ptable, chrom_list[[sample]],
                                exclude = c('names','row.names','dim',
                                              'dimnames', 'estimated_baseline', "format_out",
                                              "data_format", "parser"))
    ptable
  }
  result <- laplee(setNames(seq_along(chrom_list), names(chrom_list)), fit_sample)
  n_width <- rowSums(vapply(result, attr, numeric(2), "n_width"))
  result[] <- lapply(result, `attr<-`, "n_width", NULL)
  tables <- lapply(result, function(s) if (is.data.frame(s)) list(s) else s)
  if (!is.null(targets)){
    found <- seq_len(nrow(targets)) %in%
      unlist(lapply(tables, lapply, `[[`, "target"))
    if (any(!found)){
      label <- ifelse(is.na(targets$name), targets$rt, targets$name)
      warning(sprintf(
        "No peak was assigned to the following targets in any sample: %s.\n     Either none was found within `tol` (try a larger `tol`, or check whether `max_fwhm` removed it), or its peak went to another target.",
        paste(label[!found], collapse = ", ")), immediate. = TRUE)
    }
  }

  if (is.null(windows)){
    warn_if_width_filtered(n_width[1], n_width[2])
    warn_if_underdetermined(tables, chrom_list, fit, baseline)
  }
  structure(result,
            chrom_list = chrom_list_string,
            lambdas = lambdas, fit = fit, sd_max = sd_max, max_fwhm = max_fwhm,
            max_iter = max_iter, windows = windows, targets = targets,
            time_unit = time_unit,
            baseline = baseline,
            intensity_unit = get_metadata_attribute(chrom_list[[1]], "detector_y_unit"),
            class = c("peak_list", "list"))
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
#' @param n_before,n_after Peak counts before and after the width filter,
#' summed over samples and wavelengths.
#' @return `NULL`, invisibly; called for its warning.
#' @author Ethan Bass
#' @noRd
warn_if_width_filtered <- function(n_before, n_after){
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
#' @param tables The assembled peak tables: per sample, a list of
#' data.frames.
#' @param chrom_list The chromatograms, for their timepoints.
#' @param fit The peak model.
#' @param baseline The baseline model.
#' @return `NULL`, invisibly; called for its warning.
#' @author Ethan Bass
#' @noRd
warn_if_underdetermined <- function(tables, chrom_list, fit, baseline){
  p <- n_fit_params(fit, baseline)
  if (is.na(p)) return(invisible(NULL))
  n_bad <- 0L
  n_all <- 0L
  for (i in seq_along(tables)){
    tdiff <- get_time_resolution(chrom_list, idx = i)
    if (!is.finite(tdiff) || tdiff <= 0) next
    for (df in tables[[i]]){
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
    if (length(timepoints) < 2) return(timepoints[i])
    approx(seq_along(timepoints), timepoints, xout = i)$y
  })
  x[, c('sd', 'FWHM', 'area')] <- x[, c('sd', 'FWHM', 'area')] * tdiff * tfac
  tau_cols <- intersect(c('tau', 'tau_right', 'tau_left'), colnames(x))
  x[, tau_cols] <- x[, tau_cols] * tdiff * tfac
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
  transfer_metadata(.subset(x, i), x)
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
  n_peaks <- sum(vapply(x, function(s) sum(vapply(s, NROW, numeric(1))), numeric(1)))
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

