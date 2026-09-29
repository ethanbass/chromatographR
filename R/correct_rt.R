#' Retention time alignment via time warping
#'
#' Chromatograms are aligned using one of three algorithms, according to the
#' value of `alg`: parametric time warping (`"ptw"`), as implemented in
#' [`ptw`][ptw::ptw]; variable penalty dynamic time warping (`"vpdtw"`), as
#' implemented in [`VPdtw`][VPdtw::VPdtw]; or landmark parametric time warping
#' (`"lptw"`), which fits the warping polynomial directly to matched peak
#' positions rather than full chromatographic signals.
#'
#' Some arguments are specific to particular warping functions. For example
#' `n_traces` applies only to `"ptw"` warping;
#' `penalty`, `max_shift`, and `dilation_span` apply only to `"vpdtw"` warping;
#' and `peaks`, `rt_window`, `spectral_cutoff`, and `regularization` apply only
#' to `"lptw"` warping.
#'
#' `trwdth`, `max_shift` and `dilation_span` are widths along the time axis,
#' given in minutes whatever the time unit of `chrom_list`, and rounded
#' internally to a whole number of timepoints. `trwdth` was a count of
#' timepoints up to chromatographR 0.8.0; multiply an old value by the sampling
#' interval in minutes to convert it. A width wider than the chromatogram itself
#' is rejected with an error, since it is almost certainly a stale timepoint
#' count. `max_shift` replaces `maxshift`, which is still accepted as a count of
#' timepoints.
#'
#' `regularization` (`"lptw"`) and `penalty` (`"vpdtw"`) are not comparable:
#' both are dimensionless and higher means less warping, but they are on
#' unrelated scales.
#'
#' @aliases correct_rt
#' @import ptw
#' @importFrom scales rescale
#' @importFrom stats approx cor
#' @param chrom_list List of chromatograms in matrix format, sharing a common
#' time axis (see [`preprocess`]).
#' @param lambdas A character or numeric vector specifying the wavelengths to
#' use for alignment. `"vpdtw"` and `"lptw"` accept only a single wavelength
#' and raise an error otherwise. For one-dimensional chromatograms, this
#' argument can be omitted.
#' @param models List of models to warp by. The models provided here (if any)
#' must match the algorithm selected in `alg`, and are matched to the samples of
#' `chrom_list` by name.
#' @param reference Index or name of the sample to be used as the reference.
#' Fitted models record it by name. If not specified, chosen by
#' [`bestref`][ptw::bestref] from a similarity matrix of the chromatograms.
#' @param alg Alignment algorithm to use: parametric time warping (`"ptw"`),
#' landmark parametric time warping (`"lptw"`), or variable penalty dynamic time
#' warping (`"vpdtw"`).
#' @param peaks A peak list object from [`get_peaks`] or a `matched_peak_list`
#' object from [`match_peaks`]. If a `matched_peak_list` is provided, peak
#' matching will be skipped. Only used when `alg = "lptw"`.
#' @param what Output type: either the `"corrected.values"` (useful for
#' visual inspection and downstream analysis) or the warping `"models"`
#' (for further programmatic use).
#' @param init.coef Starting values for the optimization. Named to match the
#' argument of the same name in [`ptw`][ptw::ptw], which receives it. For
#' `alg = "lptw"`, only its length is used, setting the degree of the warping
#' polynomial.
#' @param n_traces Number of wavelengths to warp by, chosen from `lambdas` (or
#' from all wavelengths if `lambdas` is missing) by the CODA criterion in
#' [`select.traces`][ptw::select.traces]. Defaults to `NULL`, which uses every
#' wavelength in `lambdas`.
#' @param n.traces The `n.traces` argument is deprecated. Please use `n_traces`
#' instead.
#' @param fill_zeros Logical. If `TRUE`, out-of-bounds regions produced by
#' warping are filled with zeros; otherwise they are returned as `NA` (default).
#' @param n_zeros Number of zeros to pad onto each end of the chromatograms
#' before fitting. Only used when `alg = "ptw"`.
#' @param n.zeros The `n.zeros` argument is deprecated. Please use `n_zeros`
#' instead.
#' @param scale Logical. If `TRUE`, rescale each chromatogram to the range
#' \[0, 1\] before warping. The returned corrected values stay rescaled.
#' @param trwdth Width of the triangle in the weighted cross-correlation (WCC)
#' criterion, passed to
#' [`ptw`][ptw::ptw]. Defaults to `0.2`, which matches the default of `20`
#' timepoints in [`ptw`][ptw::ptw] itself on a chromatogram sampled every
#' `0.01` minutes. Roughly a peak width is a reasonable choice.
#' @param plot_it Logical. Whether to plot alignment. Defaults to `FALSE`.
#' @param penalty Numeric >= 0. Strength of the penalty against warping
#' for `alg = "vpdtw"`: the constant of proportionality in
#' \eqn{\lambda \propto} [`dilation`][VPdtw::dilation]`(reference)` over a
#' `dilation_span`-wide window, as recommended by Clifford and Stone (2012).
#' **A higher
#' number allows less warping**; `0` applies no penalty at all. Defaults to
#' `0.2`. In chromatographR <= 0.8.0 `penalty` was a *divisor*, so use the
#' reciprocal of an old value. Values above `1` raise a warning, since they are
#' more likely to be old divisors than deliberate choices.
#' @param dilation_span Width of the moving-maximum window used to build the
#' `"vpdtw"` penalty. The penalty at each point is the largest reference
#' intensity within `dilation_span` either side of it, computed by
#' [`dilation`][VPdtw::dilation] after conversion to timepoints. Smaller
#' values track individual peaks more closely, allowing more warping in the
#' baseline between them. Larger values approach a constant penalty, which
#' becomes exactly constant once `dilation_span` exceeds the distance between
#' the tallest peak of the reference and the furthest end of the chromatogram.
#' Defaults to `0.3`, around a couple of peak widths on the example data.
#' @param max_shift Maximum allowable shift for `VPdtw` warping, in minutes.
#' Set it from the retention time drift you actually expect; the default of `1`
#' allows a minute of drift. Although `max_shift` is a time, `VPdtw` itself
#' works on the index axis, so the `shift` element of the returned models is
#' reported in timepoints.
#' @param maxshift The `maxshift` argument is deprecated. Please use
#' `max_shift`, which is in minutes. `maxshift` is read as a count of
#' timepoints, as in chromatographR 0.8.0, and converted.
#' @param rt_window Numeric. Argument to [`match_peaks`]. Maximum allowed
#' difference in retention time, in minutes, between a reference and sample
#' peak to be considered a candidate match. Defaults to `0.3`.
#' @param spectral_cutoff Numeric. Argument to [`match_peaks`]. Minimum
#' spectral correlation for a candidate peak pair to be retained. Defaults to
#' `0.95`.
#' @param regularization Numeric >= 0. Strength of the penalty against warping
#' for `alg = "lptw"`: a dimensionless weight pulling the warp toward the
#' identity warp, passed as `penalty` to `fit_penalized_warp`. **A higher
#' number allows less warping**; `0` (the default) applies no penalty at all.
#' `1` weights a unit of squared warp displacement equally against a unit of
#' squared peak-position mismatch, and around `2.4` halves the applied warp, so
#' useful values run roughly `0`-`10`. Unbounded above, approaching the identity
#' warp asymptotically.
#' @param verbose Logical. Whether to print a progress message as the warping
#' models are fit and applied. Defaults to `getOption("verbose")`.
#' @param show_progress Logical. Whether to show progress bar. Defaults to
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed. Has no effect when
#' `alg = "vpdtw"`.
#' @param cl Either an integer specifying the number of cores to use for
#' parallel processing or a cluster object created by
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows systems,
#' integers are ignored. Has no effect when `alg = "vpdtw"`.
#' @param ... For `alg = "ptw"`, optional additional arguments passed to
#' [`ptw`][ptw::ptw]. The only argument that cannot be changed is `warp.type`,
#' which is hard-coded to `"global"` to permit warping on multiple wavelengths.
#' For `alg = "lptw"`, optional arguments passed to `fit_penalized_warp`:
#' `penalty_scale` (per-coefficient penalty multipliers) and `min_peaks`
#' (minimum matched peaks required to fit the polynomial). Its `penalty`
#' argument is set by `regularization` rather than through `...`.
#' @return A list of warping models or a list of warped absorbance profiles,
#' according to the value of the `what` argument. With `alg = "vpdtw"`, the
#' warped chromatograms share a new time axis, which can be longer than the
#' input. `lptw_list` models returned when `alg = "lptw"` are not compatible
#' with [`correct_peaks`].
#' @author Ethan Bass
#' @note Adapted from the `correctRT` function in the alsace package by Ron 
#' Wehrens (<https://github.com/rwehrens/alsace/blob/master/R/correctRT.R>).
#' @seealso [`correct_rt_group`], [`correct_peaks`], [`ptw`][ptw::ptw], 
#' [`VPdtw`][VPdtw::VPdtw]
#' @references 
#' * Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy, P.,
#' Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty dynamic time
#' warping. *Analytical chemistry*, 81(3):1000-1007. \doi{10.1021/ac802041e}.
#'
#' * Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping Code
#' for Aligning Mass Spectrometry Chromatograms in R. *Journal of
#' Statistical Software*, 47(8):1-17. \doi{10.18637/jss.v047.i08}.
#' 
#' * Eilers, P.H.C. 2004. Parametric Time Warping.
#' *Analytical Chemistry*, 76:404-411. \doi{10.1021/ac034800e}.
#' 
#' * Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. *Bioinformatics*,
#' 31:3063-3065. \doi{10.1093/bioinformatics/btv299}.
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R.
#' *Metabolomics*, 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examples
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, lambdas=210)
#' @md
#' @export correct_rt
correct_rt <- function(chrom_list, lambdas, models = NULL, reference = 'best',
                       alg = c("ptw", "lptw", "vpdtw"), peaks = NULL,
                       what = c("corrected.values", "models"),
                       init.coef = c(0, 1, 0), n_traces = NULL,
                       fill_zeros = FALSE, n_zeros = 0,
                       scale = FALSE, trwdth = 0.2, plot_it = FALSE,
                       penalty = 0.2, max_shift = 1, dilation_span = 0.3,
                       rt_window = 0.3, spectral_cutoff = 0.95,
                       regularization = 0,
                       verbose = getOption("verbose"), show_progress = NULL,
                       cl = 2, ..., n.traces = NULL, n.zeros = NULL,
                       maxshift = NULL){
  chrom_list_name <- deparse(substitute(chrom_list))
  n_traces <- resolve_deprecated(n.traces, n_traces)
  n_zeros <- resolve_deprecated(n.zeros, n_zeros)
  default_shift <- missing(max_shift) && is.null(maxshift)
  if (!is.null(maxshift))
    max_shift <- resolve_maxshift(maxshift, chrom_list)
  what <- match.arg(what, c("corrected.values", "models"))
  alg <- match.arg(tolower(alg), c("ptw", "lptw", "vpdtw"))
  if (length(regularization) != 1 || !is.numeric(regularization) ||
      is.na(regularization) || regularization < 0)
    stop("`regularization` must be a single number >= 0.")
  check_duplicated_names(names(chrom_list), what = "chrom_list")
  check_common_time_axis(chrom_list)
  
  if (!is.null(models)){
    model_class <- switch(alg, ptw = "ptw_list", lptw = "lptw_list",
                          vpdtw = "VPdtw")
    if (!inherits(models, model_class)){
      stop("The supplied models do not match the selected algorithm. Please check arguments
           and try again.")
    }
    stored_ref <- stored_reference(models)
    if (!is.null(stored_ref))
      reference <- reconcile_reference(
        reference, stored_ref, names(chrom_list), "the supplied `models`",
        "Re-fit the models with the desired reference, or omit `reference` to use the stored one.")
  }
  has_precomputed_matches <- alg == "lptw" && inherits(peaks, "matched_peak_list")
  if (missing(lambdas)) {
    if (ncol(chrom_list[[1]]) != 1 & is.null(models) & is.null(n_traces) &
        !has_precomputed_matches){
      stop("Must specify wavelengths ('lambdas') or number of traces ('n_traces')
           to use for alignment.")
    } else lambdas <- colnames(chrom_list[[1]])
  }
  lambdas <- as.character(lambdas)
  if (!all(lambdas %in% colnames(chrom_list[[1]]))){
    stop("Lambdas not found!")
  }
  if (scale){
    chrom_list <- lapply(chrom_list, scales::rescale)
  }
  chrom_list_og <- chrom_list
  if (n_zeros > 0){
    chrom_list <- lapply(chrom_list, function(x){
      x_p <- pad_zeros(x, n_zeros = n_zeros, side = "both")
      transfer_metadata(x_p, x)
    })
  }
  needs_reference <- (is.null(models) || alg == "vpdtw") &&
    !has_precomputed_matches
  allmats.t <- NULL
  if (identical(reference, 'best') && needs_reference){
    allmats.t <- stack_chroms(chrom_list, lambdas)
    reference <- resolve_best_reference(allmats.t, verbose)
  }
  if (is.null(models) && !has_precomputed_matches)
    check_reference(reference, chrom_list)
  args <- substitute(list(lambdas = lambdas, models = models, 
                          reference = reference, alg = alg, 
                          init.coef = init.coef, n_traces = n_traces,
                          n_zeros = n_zeros, scale = scale, trwdth = trwdth,
                          penalty = penalty, max_shift = max_shift,
                          dilation_span = dilation_span))
  if (alg %in% c("ptw", "lptw")){
    laplee <- choose_apply_fnc(show_progress, cl = cl)
    if (is.null(models)){
      if (alg == "ptw"){
        if (verbose) message("Fitting PTW warping models.")
        times <- get_times(chrom_list)
        trwdth_pts <- time_to_points(trwdth, "trwdth", chrom_list_og,
                                     clamp = missing(trwdth))
        allmats.t <- allmats.t %||% stack_chroms(chrom_list, lambdas)
        if (is.null(n_traces)){
          traces <- ifelse(length(lambdas) == 1, 1, list(lambdas))[[1]]
        } else {
          traces <- ptw::select.traces(X = allmats.t, criterion = 'coda')
          traces <- traces$trace.nrs[seq_len(n_traces)]
        }
        models <- laplee(seq_len(dim(allmats.t)[3]), function(ii){
          pw <- ptw::ptw(allmats.t[,, reference],
              allmats.t[,, ii], selected.traces = traces, init.coef = init.coef,
              warp.type = "global", trwdth = trwdth_pts, ...)
          if (nrow(pw$reference) == 1){
            colnames(pw$reference) <- times
            colnames(pw$sample) <- times
            colnames(pw$warped.sample) <- times
            rownames(pw$reference) <- lambdas
            rownames(pw$sample) <- lambdas
            rownames(pw$warped.sample) <- lambdas
          }
          if (n_zeros > 0){
            pw$warp.coef[1, ] <- shift_warp_coef(pw$warp.coef[1, ], n_zeros)
            keep <- (n_zeros + 1):(ncol(pw$reference) - n_zeros)
            pw$warp.fun <- pw$warp.fun[, keep, drop = FALSE] - n_zeros
            for (el in c("reference", "sample", "warped.sample"))
              pw[[el]] <- pw[[el]][, keep, drop = FALSE]
          }
          pw
        })
        names(models) <- names(chrom_list)
        models <- structure(models, chrom_list = chrom_list_name,
                            reference = reference_label(reference, names(chrom_list)),
                            init.coef = init.coef,
                            n_traces = n_traces, n_zeros = n_zeros, scale = scale,
                            trwdth = trwdth, penalty = penalty,
                            max_shift = max_shift, class = "ptw_list")
        if (plot_it){
          plot(models)
        }
      } else if (alg == "lptw"){
        if (inherits(peaks, "matched_peak_list")) {
          matched <- align_matches(peaks, chrom_list)
          stored_ref <- reference_position(matched)
          if (is.null(stored_ref)){
            if (identical(reference, "best"))
              reference <- resolve_best_reference(stack_chroms(chrom_list, lambdas),
                                                  verbose)
          } else {
            reference <- reconcile_reference(
              reference, stored_ref, names(chrom_list), "the `matched_peak_list`",
              "Re-run `match_peaks` with the desired reference.",
              label = attr(matched, "reference"))
          }
        } else {
          if (length(lambdas) > 1)
            stop("LPTW only supports warping by a single wavelength")
          if (is.null(peaks)){
            peaks <- get_peaks(chrom_list_og, lambdas = lambdas)
          }
          matched <- match_peaks(peaks, chrom_list = chrom_list_og,
                                 reference = reference, lambda = lambdas,
                                 rt_window = rt_window,
                                 spectral_cutoff = spectral_cutoff,
                                 verbose = verbose)
        }
        identity_coef <- identity_warp(length(init.coef))
        models <- laplee(matched, function(m) {
          if (is.null(m))
            return(list(warp.coef = matrix(identity_coef, nrow = 1),
                        mode = "forward"))
          w <- (m$cor_scores %||% 1) * pmin(m$ref_i[, "I"], m$samp_i[, "I"])
          warnings <- character(0)
          coef <- withCallingHandlers(
            fit_penalized_warp(
              ref_i = m$ref_i, samp_i = m$samp_i,
              weights = w, init.coef = init.coef, penalty = regularization, ...
            ),
            warning = function(cond){
              warnings <<- c(warnings, conditionMessage(cond))
              invokeRestart("muffleWarning")
            })
          list(warp.coef = matrix(coef, nrow = 1), mode = "forward",
               warnings = warnings)
        })
        for (i in seq_along(models)){
          tp <- seq_len(nrow(chrom_list_og[[i]]))
          if (is.unsorted(ptw::warp.time(tp, models[[i]]$warp.coef[1, ]))){
            models[[i]]$warp.coef[1, ] <- identity_coef
            models[[i]]$warnings <- c(models[[i]]$warnings, paste(
              "The fitted warp reverses direction within the chromatogram.",
              "Returning identity warp. Try a lower-degree `init.coef` or a",
              "higher `regularization`."))
          }
          for (msg in models[[i]]$warnings)
            warning(sprintf("Sample %s: %s", names(chrom_list)[i] %||% i, msg),
                    call. = FALSE)
          models[[i]]$warnings <- NULL
        }
        names(models) <- names(chrom_list)
        models <- structure(models, class = "lptw_list",
                            reference = reference_label(reference, names(chrom_list)))
      }
    }
    if (what == "corrected.values"){
      if (verbose) message(sprintf("Applying %s warping models to chromatograms.", toupper(alg)))
      pos <- model_positions(models, chrom_list_og)
      result <- laplee(seq_along(chrom_list_og), function(ii){
        mod <- models[[pos[ii]]]
        x <- warp_sample_ptw(chrom_list_og[[ii]], mod$warp.coef[1, ],
                             fill_zeros = fill_zeros, mode = mod$mode)
        rownames(x) <- rownames(chrom_list_og[[ii]])
        colnames(x) <- colnames(chrom_list_og[[ii]])
        transfer_metadata(x, chrom_list_og[[ii]])
      })
      result <- structure(result, warped = TRUE, warping_args = args)
      names(result) <- names(chrom_list)
      result
    } else {
      return(models)
    }
  } else if (alg == "vpdtw"){
    if (is.null(models)){
      if (length(lambdas) > 1)
        stop("VPdtw only supports warping by a single wavelength")
      check_vpdtw_penalty(penalty)
      if (verbose) message("Fitting VPdtw warping models.")
      span_pts <- time_to_points(dilation_span, "dilation_span", chrom_list_og,
                                 idx = reference, clamp = missing(dilation_span))
      check_dilation_span(span_pts, chrom_list_og)
      maxshift_pts <- time_to_points(max_shift, "max_shift", chrom_list_og,
                                     idx = reference, clamp = default_shift)
      allmats <- sapply(chrom_list_og, function(x) x[, lambdas, drop = FALSE])
      penalty <- VPdtw::dilation(allmats[,reference], span_pts) * penalty
      models <- VPdtw::VPdtw(query = allmats, reference = allmats[, reference],
                             penalty = penalty, maxshift = maxshift_pts)
      attr(models, "reference") <- reference_label(reference, names(chrom_list_og))
      attr(models, "parameters") <- list(reference = attr(models, "reference"),
                                         max_shift = max_shift)
    }
    if (plot_it){
      VPdtw::plot.VPdtw(models)
    }
    if (what == "corrected.values"){
      if (verbose) message("Applying VPdtw warping models to chromatograms.")
      apply_vpdtw(chrom_list_og, models, reference, args = args,
                  fill_zeros = fill_zeros)
    } else {
      return(models)
    }
  }
}

#' Check that chromatograms share a time axis
#'
#' Samples may differ from the first by less than one sampling interval.
#' @noRd
check_common_time_axis <- function(chrom_list){
  dt <- get_time_resolution(chrom_list, idx = 1)
  if (length(chrom_list) < 2 || !is.finite(dt)) return(invisible(chrom_list))
  t1 <- get_times(chrom_list, idx = 1)
  off <- vapply(seq_along(chrom_list), function(i){
    t <- get_times(chrom_list, idx = i)
    if (length(t) != length(t1)) Inf else max(abs(t - t1))
  }, numeric(1))
  bad <- which(!(off < dt))
  if (length(bad)){
    ids <- names(chrom_list) %||% as.character(seq_along(chrom_list))
    stop(sprintf(paste0(
      "All chromatograms must share a time axis for retention time ",
      "correction, but %s %s %s from %s by a sampling interval or more. ",
      "Interpolate them onto a common grid with `preprocess(dim1 = ...)`."),
      if (length(bad) == 1) "sample" else "samples",
      paste(sQuote(ids[bad]), collapse = ", "),
      if (length(bad) == 1) "differs" else "differ",
      sQuote(ids[1])), call. = FALSE)
  }
  invisible(chrom_list)
}

#' Check the VPdtw penalty
#' @noRd
check_vpdtw_penalty <- function(penalty){
  if (length(penalty) != 1 || !is.numeric(penalty) || is.na(penalty))
    stop("`penalty` must be a single non-missing number.")
  if (penalty < 0)
    stop("`penalty` must be >= 0.")
  if (penalty > 1)
    warning(sprintf(paste0(
      "`penalty` is above 1: you supplied %s.\n",
      "  `penalty` changed meaning in chromatographR 0.9.0: it is now a ",
      "multiplier applied to the dilation of the reference,\n",
      "  so a higher value ",
      "allows less warping. Up to 0.8.0 it was a divisor.\n",
      "  If your %s was meant as the old divisor, use its reciprocal: ",
      "`penalty = %s`."),
      format(penalty), format(penalty), format(signif(1 / penalty, 4))))
  invisible(penalty)
}

#' Check the VPdtw dilation span
#' @param span Width of the dilation window, in timepoints.
#' @noRd
check_dilation_span <- function(span, chrom_list){
  n <- nrow(chrom_list[[1]])
  window <- 2 * span - 1
  if (window >= n)
    warning(sprintf(paste0(
      "`dilation_span` spans %s points, which is wider than the %s-point ",
      "chromatogram, so the VPdtw penalty will be nearly constant rather ",
      "than variable. Consider a value below %s minutes."),
      format(window), format(n),
      format(signif(ceiling((n + 1) / 2) *
                      get_time_resolution(chrom_list, unit = "min"), 3))))
  invisible(span)
}

#' Convert a width in minutes to a number of timepoints
#'
#' @param value Width in minutes.
#' @param name Name of the argument, for messages.
#' @param chrom_list List of chromatograms supplying the time axis.
#' @param idx Index of the chromatogram to read the time axis from.
#' @param clamp If `TRUE`, a width wider than the chromatogram is narrowed to
#'   its span instead of raising an error.
#' @return The width in timepoints, a whole number >= 1.
#' @noRd
time_to_points <- function(value, name, chrom_list, idx = 1, clamp = FALSE){
  if (length(value) != 1 || !is.numeric(value) || !is.finite(value) ||
      value <= 0)
    stop(sprintf("`%s` must be a single positive number.", name))
  tfac <- time_factor(get_time_unit(chrom_list, idx = idx, na_value = "min"))
  times <- get_times(chrom_list, idx = idx) / tfac
  dt <- get_time_resolution(chrom_list, idx = idx) / tfac
  if (length(times) < 2 || !is.finite(dt) || dt <= 0)
    stop(sprintf(paste0(
      "`%s` is a width in minutes, but the time axis could not be read. ",
      "Chromatograms must carry numeric retention times as their rownames."),
      name))
  span <- diff(range(times))
  if (value > span && clamp)
    value <- span
  if (value > span)
    stop(sprintf(paste0(
      "`%s = %s` is wider than the chromatogram itself, which spans %s ",
      "minutes.\n",
      "  `%s` changed units in chromatographR 0.9.0: it is now a width in ",
      "minutes, not a number of timepoints.\n",
      "  To convert a timepoint count, multiply it by the sampling interval ",
      "of %s minutes: `%s = %s`."),
      name, format(value), format(signif(span, 4)), name, format(dt), name,
      format(signif(value * dt, 4))))
  gaps <- diff(times)
  if (max(abs(gaps - dt)) > 0.01 * dt)
    warning(sprintf(paste0(
      "The time axis is not evenly spaced, so `%s` has been converted to ",
      "timepoints using the median sampling interval of %s minutes."),
      name, format(dt)))
  pts <- round(value / dt)
  if (pts < 1){
    warning(sprintf(paste0(
      "`%s = %s` is narrower than the sampling interval of %s minutes, so it ",
      "has been widened to a single timepoint."), name, format(value),
      format(dt)))
    pts <- 1
  }
  pts
}

#' Resolve the deprecated `maxshift` argument, a count of timepoints
#' @noRd
resolve_maxshift <- function(maxshift, chrom_list){
  if (length(maxshift) != 1 || !is.numeric(maxshift) || !is.finite(maxshift) ||
      maxshift <= 0)
    stop("`maxshift` must be a single positive number.")
  max_shift <- maxshift * get_time_resolution(chrom_list, unit = "min")
  warning(sprintf(paste0(
    "The `maxshift` argument is deprecated. Please use `max_shift`, which is ",
    "in minutes. `maxshift = %s` was read as a count of timepoints, ",
    "i.e. `max_shift = %s`."), format(maxshift), format(signif(max_shift, 4))),
    immediate. = TRUE)
  max_shift
}

#' Coefficients of an identity warp
#' @noRd
identity_warp <- function(n_coef){
  c(0, 1, rep(0, n_coef - 2))
}

#' Reconcile a requested reference against one stored on a fitted object
#' @noRd
reconcile_reference <- function(reference, stored_ref, sample_names, source,
                                remedy, label = stored_ref){
  as_idx <- function(r) if (is.character(r) && r %in% sample_names)
    match(r, sample_names) else suppressWarnings(as.numeric(r))
  if (!identical(reference, "best") &&
      !isTRUE(all.equal(as_idx(reference), as_idx(stored_ref))))
    stop("The `reference` argument (", reference, ") differs from the ",
         "reference stored in ", source, " (", label, "). ", remedy)
  stored_ref
}

#' Check that a reference names or indexes one chromatogram
#' @noRd
check_reference <- function(reference, chrom_list){
  ok <- length(reference) == 1 && !is.na(reference) &&
    if (is.character(reference)) reference %in% names(chrom_list)
    else is.numeric(reference) && reference == round(reference) &&
      reference >= 1 && reference <= length(chrom_list)
  if (!ok)
    stop(sprintf(paste0(
      "`reference` must be \"best\", the name of a sample in `chrom_list`, ",
      "or an index between 1 and %d."), length(chrom_list)), call. = FALSE)
  invisible(reference)
}

#' Reference stored on fitted warping models
#'
#' `VPdtw` models from chromatographR <= 0.8.0 record it only in their
#' `parameters` attribute.
#' @noRd
stored_reference <- function(models){
  attr(models, "reference") %||% attr(models, "parameters")$reference
}

#' Name of a reference sample, for storing on fitted models
#' @return The sample name, or `reference` itself if the samples are unnamed.
#' @noRd
reference_label <- function(reference, sample_names){
  if (is.numeric(reference) && !is.null(sample_names))
    sample_names[reference]
  else reference
}

#' Position of the chromatogram whose time axis a stored reference points to
#'
#' Accepts a sample name or an index. Falls back to the first chromatogram when
#' the reference sample is not in `chrom_list`.
#' @noRd
reference_axis_idx <- function(reference, chrom_list){
  i <- if (is.character(reference)) match(reference, names(chrom_list))
       else reference
  if (length(i) != 1 || is.na(i) || i < 1 || i > length(chrom_list)) 1L
  else as.integer(i)
}

#' Position of each sample's warping model
#'
#' Matches by sample name when both the models and `chrom_list` are named, and
#' by position otherwise.
#' @return An integer vector with one element per chromatogram in `chrom_list`.
#' @noRd
model_positions <- function(models, chrom_list){
  vpdtw <- inherits(models, "VPdtw")
  model_names <- if (vpdtw) models$samples %||% colnames(models$query)
                 else names(models)
  n_models <- if (vpdtw) ncol(models$shift) else length(models)
  if (!is.null(model_names) && !is.null(names(chrom_list))){
    pos <- match(names(chrom_list), model_names)
    if (anyNA(pos))
      stop("The supplied `models` have no warp for sample(s) ",
           paste(sQuote(names(chrom_list)[is.na(pos)]), collapse = ", "), ".",
           call. = FALSE)
    return(pos)
  }
  if (n_models != length(chrom_list))
    stop("The supplied `models` warp ", n_models, " samples, but `chrom_list` ",
         "has ", length(chrom_list), ".", call. = FALSE)
  seq_along(chrom_list)
}

#' Check that chromatograms have as many timepoints as a VPdtw fit
#' @noRd
check_vpdtw_length <- function(models, chrom_list){
  n <- models$n_timepoints %||% nrow(models$query)
  bad <- vapply(chrom_list, nrow, integer(1)) != n
  if (any(bad)){
    ids <- names(chrom_list) %||% as.character(seq_along(chrom_list))
    stop(sprintf(paste0(
      "The supplied `models` were fitted to chromatograms of %d timepoints, ",
      "but %s %s a different number."), n,
      if (sum(bad) == 1) "sample" else "samples",
      paste(sQuote(ids[bad]), collapse = ", ")), call. = FALSE)
  }
  invisible(chrom_list)
}

#' Apply VPdtw
#' @noRd
apply_vpdtw <- function(chrom_list, models, reference, args = NULL, 
                        fill_zeros = FALSE) {
  check_vpdtw_length(models, chrom_list)
  src <- models$xVals - models$shift
  pos <- model_positions(models, chrom_list)
  new_rownames <- format_axis(vpdtw_times(
    chrom_list, models$xVals, reference_axis_idx(reference, chrom_list)))
  result <- lapply(seq_along(chrom_list), function(ii) {
    x <- warp_sample_vpdtw(chrom_list[[ii]], src[, pos[ii]])
    rownames(x) <- new_rownames
    transfer_metadata(x, chrom_list[[ii]])
  })
  names(result) <- names(chrom_list)
  
  lapply(result, function(xx) {
    if (fill_zeros && any(is.na(xx))) xx[which(is.na(xx))] <- 0
    structure(xx, warped = TRUE, warping_args = args)
  })
}

#' Time axis of the VPdtw target grid
#' @param x Positions on the target grid, in the index units of `xVals`, where
#'   `1` to `n` are the rows of the reference.
#' @noRd
vpdtw_times <- function(chrom_list, x, reference){
  times <- get_times(chrom_list, idx = reference)
  dt <- get_time_resolution(chrom_list, idx = reference)
  n <- length(times)
  inside <- pmin(pmax(x, 1), n)
  approx(seq_len(n), times, xout = inside)$y + (x - inside) * dt
}

#' Map sample retention times onto the VPdtw target axis
#' @return One function per sample, taking retention times on that sample's
#'   own axis and returning them on the target axis of `apply_vpdtw`.
#' @noRd
vpdtw_to_ref <- function(models, chrom_list, reference){
  check_vpdtw_length(models, chrom_list)
  src <- models$xVals - models$shift
  pos <- model_positions(models, chrom_list)
  ref_idx <- reference_axis_idx(reference, chrom_list)
  lapply(seq_along(chrom_list), function(ii){
    times <- get_times(chrom_list, idx = ii)
    ok <- !is.na(src[, pos[ii]])
    function(rt){
      j <- approx(times, seq_along(times), xout = rt)$y
      x <- approx(src[ok, pos[ii]], models$xVals[ok], xout = j, ties = mean)$y
      vpdtw_times(chrom_list, x, ref_idx)
    }
  })
}

#' Warp one sample onto the VPdtw target axis
#' @param src Sample row feeding each row of the target grid, as in
#'   `VPdtw`'s `warpedQuery`: `xVals - shift`.
#' @noRd
warp_sample_vpdtw <- function(mat, src) {
  ok <- !is.na(src) & src >= 1 & src <= nrow(mat)
  out <- matrix(NA_real_, length(src), ncol(mat),
                dimnames = list(NULL, colnames(mat)))
  out[ok, ] <- mat[src[ok], , drop = FALSE]
  out
}

#' Align a matched_peak_list to a chrom_list
#' @param matched A `matched_peak_list`.
#' @param chrom_list The chromatograms to be warped.
#' @return `matched`, ordered to match `chrom_list`.
#' @noRd
align_matches <- function(matched, chrom_list){
  if (!is.null(names(matched)) && !is.null(names(chrom_list))){
    missing_names <- setdiff(names(chrom_list), names(matched))
    if (length(missing_names))
      stop("The `matched_peak_list` has no entry for sample(s) ",
           paste(missing_names, collapse = ", "), ".")
    ref <- attr(matched, "reference")
    if (!is.null(ref) && !ref %in% names(chrom_list))
      stop("The reference sample of the `matched_peak_list` (", ref,
           ") is not in `chrom_list`.")
    return(matched[names(chrom_list)])
  }
  if (length(matched) != length(chrom_list))
    stop("The `matched_peak_list` has ", length(matched), " samples but ",
         "`chrom_list` has ", length(chrom_list), ".")
  matched
}

#' Re-express a PTW warp fitted on zero-padded indices in unpadded indices
#' @param coef Warp polynomial coefficients.
#' @param shift Number of points padded onto the start of the chromatogram.
#' @return Coefficients `c` such that `warp.time(t, c) = warp.time(t + shift, coef) - shift`.
#' @noRd
shift_warp_coef <- function(coef, shift){
  k <- seq_along(coef) - 1
  m <- outer(k, k, function(j, kk) ifelse(kk >= j, choose(kk, j) * shift^pmax(kk - j, 0), 0))
  out <- c(m %*% coef)
  out[1] <- out[1] - shift
  out
}

#' Warp sample PTW
#' @noRd
warp_sample_ptw <- function(mat, coef, fill_zeros = FALSE,
                            mode = c("forward", "backward")){
  mode <- match.arg(mode)
  tp <- seq_len(nrow(mat))
  w <- ptw::warp.time(tp, coef)
  if (is.unsorted(w))
    stop("The warping function reverses direction within the chromatogram. ",
         "Try a lower-degree `init.coef`, or for `alg = \"lptw\"` a higher ",
         "`regularization`.", call. = FALSE)
  if (mode == "forward"){
    lo <- findInterval(tp, w)
    hi <- lo + 1
    oob <- lo < 1 | hi > length(w)
    lo_c <- pmax(1, lo)
    hi_c <- pmin(length(w), hi)
    frac <- (tp - w[lo_c]) / (w[hi_c] - w[lo_c])
  } else {
    stop("Backward mode is not yet supported.")
  }
  frac[oob] <- 0
  frac_mat <- matrix(frac, nrow = length(tp), ncol = ncol(mat))
  x <- mat[lo_c, ] + frac_mat * (mat[hi_c, ] - mat[lo_c, ])
  x[which(tp < min(w) | tp > max(w)), ] <- ifelse(fill_zeros, 0, NA)
  x
}


#' Fit a penalized warp polynomial to matched peak pairs
#' @param ref_i,samp_i Matrices with columns `"rt"` (time-point index) and `"I"`
#'   for the matched reference and sample peaks, in the same row order.
#' @param weights Non-negative weights, one per pair. Defaults to uniform.
#' @param init.coef Its length sets the polynomial degree.
#' @param penalty Dimensionless strength of the pull toward the identity warp.
#' @param penalty_scale Per-coefficient multipliers on `penalty`.
#' @param min_peaks Minimum number of pairs. Defaults to `length(init.coef)`.
#' @return Warp coefficients, or the identity warp with a warning if the fit is
#'   impossible.
#' @noRd
fit_penalized_warp <- function(ref_i, samp_i,
                               weights       = NULL,
                               init.coef     = c(0, 1, 0),
                               penalty       = 0,
                               penalty_scale = NULL,
                               min_peaks     = NULL) {
  n_coef <- length(init.coef)
  identity_coef <- identity_warp(n_coef)
  if (is.null(min_peaks)) min_peaks <- n_coef

  if (is.null(weights)) weights <- rep(1.0, nrow(ref_i))
  usable <- is.finite(weights) & is.finite(ref_i[, "rt"]) &
    is.finite(samp_i[, "rt"])
  ref_i <- ref_i[usable, , drop = FALSE]
  samp_i <- samp_i[usable, , drop = FALSE]
  weights <- weights[usable]

  n_peaks <- nrow(ref_i)
  if (n_peaks < min_peaks) {
    warning(sprintf(
      "fit_penalized_warp: %d peak(s) insufficient for degree-%d polynomial (need >= %d). Returning identity warp.",
      n_peaks, n_coef - 1L, min_peaks))
    return(identity_coef)
  }

  weights <- pmax(weights, 0)
  if (sum(weights) == 0) {
    warning("fit_penalized_warp: all weights are zero. Returning identity warp.")
    return(identity_coef)
  }
  weights <- weights / sum(weights)

  if (is.null(penalty_scale)) {
    penalty_scale <- rep(1.0, n_coef)
  } else {
    if (length(penalty_scale) != n_coef)
      stop("'penalty_scale' must have the same length as 'init.coef'.")
  }

  ref_rt  <- ref_i[, "rt"]
  samp_rt <- samp_i[, "rt"]

  X <- outer(samp_rt, seq_len(n_coef) - 1, "^")
  coef_scale <- sqrt(colSums(weights * X^2))
  Xs <- sweep(X, 2, coef_scale, "/")
  D <- diag(penalty * penalty_scale, n_coef)
  A <- crossprod(Xs, weights * Xs) + D
  if (rcond(A) < 1e-12) {
    warning("fit_penalized_warp: matched peaks do not determine a degree-",
            n_coef - 1L, " polynomial. Returning identity warp.")
    return(identity_coef)
  }
  b <- solve(A, crossprod(Xs, weights * ref_rt) + D %*% (identity_coef * coef_scale))
  c(b) / coef_scale
}

#' Correct retention time with group-based warping
#'
#' Aligns chromatograms using parametric time warping (`"ptw"`), as
#' implemented in [`ptw`][ptw::ptw], with warping functions estimated
#' from within-group averages. This is useful when samples fall into batches
#' with shared retention time shifts, or when individual samples contain peaks
#' absent in other groups that would otherwise confound alignment. One warp is
#' fit per group rather than one per sample, so the cost of warping falls
#' roughly in proportion to the number of groups.
#'
#' @aliases correct_rt_group
#' @import ptw
#' @param chrom_list A list of chromatograms in matrix format, sharing a common
#' time axis (see [`preprocess`]).
#' @param lambdas A character or numeric vector specifying the wavelengths to
#' use for alignment.
#' @param groups A vector of group assignments for each chromatogram, or a
#' single string naming a metadata attribute from which to extract the group
#' assignments (e.g. `"batch"`). If a vector is provided, it must be the same 
#' length as `chrom_list` and may be named (matching `names(chrom_list)`) or
#' positional. All samples must have a group assignment; `NA` values will
#' trigger an error.
#' @param reference Index or name of a sample in `chrom_list` whose group will
#' be used as the alignment reference. Defaults to `"best"`, which selects the
#' reference group automatically using [`bestref`][ptw::bestref] applied to the
#' group averages.
#' @param reference_group Name of the group to use as the reference group. If
#' supplied, overrides `reference`. Defaults to `NULL`.
#' @param penalty Unused. `correct_rt_group` warps with `"ptw"` only, so the
#' `"vpdtw"` penalty has no effect here.
#' @param max_shift Unused. `correct_rt_group` warps with `"ptw"` only, so the
#' `"vpdtw"` shift limit has no effect here.
#' @param maxshift Deprecated and unused. See `max_shift`.
#' @inheritParams correct_rt
#' @return A list of warped chromatogram matrices in the same order as
#' `chrom_list`, with each sample warped using the warp coefficients
#' estimated from its group average.
#' @author Ethan Bass
#' @seealso [`correct_rt`], [`ptw`][ptw::ptw]
#' @examples
#' data(Sa_pr)
#' groups <- c("a", "a", "b", "b")
#' warp <- correct_rt_group(chrom_list = Sa_pr, lambdas = 210, groups = groups)
#' @md
#' @export

correct_rt_group <- function(chrom_list, lambdas, groups, reference = 'best',
                             reference_group = NULL,
                             init.coef = c(0, 1, 0), n_traces = NULL,
                             fill_zeros = FALSE, n_zeros = 0, scale = FALSE, 
                             trwdth = 0.2, plot_it = FALSE, penalty = NULL, 
                             max_shift = 1, verbose = getOption("verbose"),
                             show_progress = NULL, cl = 2,
                             ..., n.traces = NULL, n.zeros = NULL,
                             maxshift = NULL) {
  n_traces <- resolve_deprecated(n.traces, n_traces)
  n_zeros <- resolve_deprecated(n.zeros, n_zeros)
  max_shift <- resolve_deprecated(maxshift, max_shift)
  check_duplicated_names(names(chrom_list), what = "chrom_list")
  check_common_time_axis(chrom_list)
  if (length(groups) == 1 && is.character(groups)) {
    groups_vec <- lapply(chrom_list, function(x) attr(x, groups))
    absent <- vapply(groups_vec,
                     function(g) is.null(g) || length(g) != 1 || is.na(g),
                     logical(1))
    if (any(absent)) {
      stop(sprintf(
        "Metadata attribute '%s' not found on all chromatograms. Missing for:\n%s",
        groups,
        paste(sprintf("         %d (%s)", which(absent),
                      sQuote(names(chrom_list)[absent])), collapse = ", \n")))
    }
    groups_vec <- unlist(groups_vec)
  } else {
    groups_vec <- groups
    if (is.null(names(groups_vec))) {
      if (length(groups_vec) != length(chrom_list)){
        stop("Length of `groups` must match the number of samples in `chrom_list`.")
      }
      names(groups_vec) <- names(chrom_list)
    } else {
      groups_vec <- groups_vec[names(chrom_list)] 
    }
  }
  groups_vec <- as.factor(groups_vec)
  if (any(is.na(groups_vec))) {
    na_idx <- which(is.na(groups_vec))
    na_samples <- names(groups_vec)[na_idx]
    stop(sprintf("All samples must have a group assignment. NA%s found for:\n%s",
                 plural(length(na_idx)),
                 paste(
                   sprintf("         %d (%s)", na_idx, sQuote(na_samples)),
                   collapse = ", \n")))
  }
  
  group_levels <- levels(groups_vec)
  
  if (verbose) message("Computing within-group averages.")
  group_avg_list <- lapply(setNames(group_levels, group_levels), function(g) {
    members <- names(groups_vec)[which(groups_vec == g)]
    mats <- lapply(chrom_list[members], as.matrix)
    avg <- Reduce("+", mats) / length(mats)
    dimnames(avg) <- dimnames(mats[[1]])
    avg
  })
  
  if (is.null(reference_group)) {
    if (reference == 'best') {
      reference_group <- NULL  # passed through, correct_rt will resolve
    } else {
      ref_name <- if (is.numeric(reference)) names(chrom_list)[reference] else reference
      if (is.na(ref_name) || !ref_name %in% names(groups_vec)){
        stop(sprintf(
          "`reference` (%s) does not match any sample in `chrom_list`.",
          sQuote(if (is.numeric(reference)) reference else ref_name)))
      }
      reference_group <- groups_vec[ref_name]
      reference <- which(names(group_avg_list) == reference_group)
    }
  } else {
    reference <- which(names(group_avg_list) == reference_group)
    if (length(reference) == 0){
      stop(sprintf(
        "`reference_group` (%s) does not match any group. Valid groups are: %s.",
        sQuote(reference_group),
        paste(sQuote(group_levels), collapse = ", ")))
    }
  }
  
  if (verbose) message("Fitting PTW warp models on group averages.")
  group_models <- correct_rt(
    chrom_list  = group_avg_list,
    lambdas     = lambdas,
    reference   = reference,
    alg         = "ptw",
    what        = "models",
    init.coef   = init.coef,
    n_traces    = n_traces,
    n_zeros     = n_zeros,
    scale       = scale,
    trwdth      = trwdth,
    plot_it     = plot_it,
    penalty     = penalty,
    max_shift   = max_shift,
    verbose     = verbose,
    show_progress = show_progress,
    cl          = cl,
    ...
  )
  laplee <- choose_apply_fnc(show_progress, cl = cl)
  if (verbose) message("Applying group warp models to individual chromatograms.")
  result <- laplee(group_levels, function(g) {
    members <- names(groups_vec)[which(groups_vec == g)]
    member_list <- chrom_list[members]
    correct_rt(
      chrom_list  = member_list,
      lambdas     = lambdas,
      models      = structure(setNames(rep(group_models[g],
                                                  length(member_list)),
                                              names(member_list)),
                              class = "ptw_list"),
      alg         = "ptw",
      what        = "corrected.values",
      init.coef   = init.coef,
      n_zeros     = n_zeros,
      fill_zeros = fill_zeros,
      scale       = scale,
      trwdth      = trwdth,
      verbose     = verbose,
      show_progress = FALSE,
      cl          = 1,
      ...
    )
  })
  
  result_flat <- do.call(c, result)
  result_flat[names(chrom_list)]
}

#' Pad zeros
#' Helper for [`correct_rt`].
#' @noRd
pad_zeros <- function(x, n_zeros, side = "both") {
  if (is.vector(x)) x <- as.matrix(x)
  zeros <- matrix(0, nrow = n_zeros, ncol = ncol(x))
  result <- switch(side,
         both = rbind(zeros, x, zeros),
         left = rbind(zeros, x),
         right = rbind(x, zeros)
  )
  if (!is.null(rownames(x))) {
    times <- as.numeric(rownames(x))
    dt <- median(diff(times))
    pre <- if (side %in% c("both", "left"))
      seq(times[1] - n_zeros * dt, by = dt, length.out = n_zeros) else numeric(0)
    post <- if (side %in% c("both", "right"))
      seq(times[length(times)] + dt, by = dt, length.out = n_zeros) else numeric(0)
    rownames(result) <- format_axis(c(pre, times, post))
  }
  result
}

#' Strip zeros
#' @noRd
strip_zeros <- function(x, n_zeros, side = "both") {
  n <- nrow(x)
  switch(side,
         both = x[(n_zeros + 1):(n - n_zeros), ],
         left = x[(n_zeros + 1):n, ],
         right = x[1:(n - n_zeros), ]
  )
}

#' Match peaks across chromatograms
#'
#' Matches peaks across chromatograms by spectral correlation and retention
#' time proximity. Peak shape similarity can be added, either as a hard filter
#' on candidate pairs (`shape_cutoff`) or as a penalty on their ranking
#' (`shape_weight`); both are off by default. Each reference peak is matched to
#' at most one sample peak, taking candidate pairs in order of score, highest
#' first. The resulting `matched_peak_list` object can be passed to
#' [`correct_rt`] to perform landmark parametric time warping
#' (`alg = "lptw"`).
#'
#' @param peak_list A peak list object from [`get_peaks`].
#' @param chrom_list List of chromatograms in matrix format. If `NULL`, will be
#'   retrieved from the `peak_list` object.
#' @param reference Index or name of the reference sample, or `"best"` to select the
#'   reference algorithmically using [`bestref`][ptw::bestref]. Defaults to `1`.
#' @param lambda A character or numeric value specifying the wavelength to use
#'   for peak matching. Required if peaks were detected at multiple wavelengths.
#' @param rt_window Numeric. Maximum allowed difference in retention time, in
#'   minutes, between a reference and sample peak to be considered a candidate
#'   match. Defaults to `0.3`.
#' @param spectral_cutoff Numeric. Minimum spectral correlation for a candidate
#'   peak pair to be retained. Ignored, with a warning, for one-dimensional
#'   chromatograms. Defaults to `0.95`.
#' @param shape_cutoff Numeric. Maximum shape distance for a candidate peak pair
#'   to be retained. Raises an error when `fit = "raw"`, since peak shape
#'   parameters are not fit. Defaults to `NULL` (disabled).
#' @param rt_weight Numeric in \[0, 1\]. Penalizes retention time offset when
#'   ranking candidate pairs: each pair's spectral correlation is multiplied by
#'   `1 - rt_weight * d / rt_window`, where `d` is the pair's retention time
#'   difference. At `rt_weight = 1`, a pair at the edge of the window scores
#'   zero. Defaults to `0.15`.
#' @param shape_weight Numeric in \[0, 1\]. Applies a multiplicative
#'   shape-distance penalty analogous to `rt_weight`. Shape distances are
#'   normalized by `shape_cutoff` (if provided) or by the maximum observed shape
#'   distance. Silently ignored when `fit = "raw"`. Defaults to `0`.
#' @param plot_it Logical. Whether to plot the matches with
#'   [`plot.matched_peak_list`]. Defaults to `FALSE`.
#' @param fit Character. The type of fit used to construct the peak list, one
#'   of `"egh"`, `"bemg"`, `"gaussian"` or `"raw"`. Normally this does not need
#'   to be specified, since it is detected from the `fit` attribute stored in
#'   the `peak_list` object. Supply it explicitly only if that attribute is
#'   missing, or to override it (e.g. to use the simpler gaussian shape
#'   distance, which depends only on `sd`, for an `"egh"` peak list).
#' @param verbose Whether to print messages to console.
#' @return A `matched_peak_list` object: a list with one element per sample
#'   containing the matched peak pairs, per-pair spectral correlation scores,
#'   and peak data in index space for use with [`correct_rt`]. A sample with no
#'   peaks or no matches is `NULL`, with a warning.
#' @author Ethan Bass
#' @seealso [`correct_rt`], [`get_peaks`]
#' @examples
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = 210, fit = "egh")
#' match_peaks(pks, lambda = 210)
#' @export match_peaks
match_peaks <- function(peak_list, chrom_list = NULL, reference = 1,
                        lambda = NULL, rt_window = 0.3, spectral_cutoff = 0.95,
                        shape_cutoff = NULL, rt_weight = 0.15, shape_weight = 0,
                        plot_it = FALSE, fit = NULL,
                        verbose = getOption("verbose")) {
  if (is.null(chrom_list)){
    chrom_list <- get_chrom_list(peak_list)
  }
  rt_window_rt <- rt_window * rt_factor(get_time_unit(peak_list))
  if (is.null(lambda)){
    if (length(peak_list[[1]]) == 1){
      lambda_idx <- names(peak_list[[1]])[1]
    } else{
      stop("Please specify the wavelength to match peaks at")
    }
  } else{
    lambda_idx <- as.character(lambda)
  }
  if (identical(reference, "best")) {
    reference <- resolve_best_reference(stack_chroms(chrom_list, lambda_idx),
                                        isTRUE(verbose))
  } else if (is.character(reference)) {
    ref_idx <- match(reference, names(peak_list))
    if (length(reference) != 1 || is.na(ref_idx))
      stop(sprintf("`reference` %s is not a sample in `peak_list`.",
                   sQuote(paste(reference, collapse = ", "))))
    reference <- ref_idx
  }
  check_chrom_list(peak_list, chrom_list)
  if (is.null(names(peak_list)) || anyDuplicated(names(peak_list)))
    stop("The samples in `peak_list` must have unique names.")
  if (is.null(fit)){
    fit <- attr(peak_list, "fit")
    if (is.null(fit)){
      stop("Please specify `fit` to proceed.")
    }
  }
  fit <- match.arg(fit, c("egh", "bemg", "gaussian", "raw"))
  pks_flat <- lapply(peak_list, function(x) x[[lambda_idx]])
  pks_ref <- pks_flat[[reference]]
  pks_ref <- pks_ref[!is.na(pks_ref$rt),]
  rt_ref <- pks_ref$rt
  if (length(rt_ref) == 0){
    stop(sprintf(
      "No peaks were found in the reference sample (%s) at wavelength %s. Please choose a different `reference` or relax the criteria used to detect peaks.",
      sQuote(names(peak_list)[reference]), sQuote(lambda_idx)))
  }
  use_shape <- !is.null(shape_cutoff) || (shape_weight > 0 && fit != "raw")
  if (use_shape){
    if (fit == "raw"){
      stop("Peak shape parameters are not fit when `fit = \"raw\"`. Please set ",
           "`shape_cutoff = NULL` to proceed.")
    }
    shape_cols <- check_shape_cols(fit, pks_ref)
  }
  time_res <- get_time_resolution(chrom_list, reference)
  t_start <- get_times(chrom_list, reference)[1]

  spectra_at <- function(rts, i){
    m <- t(chrom_list[[i]][get_retention_idx(rts, get_times(chrom_list, i)), ,
                           drop = FALSE])
    colnames(m) <- rts
    m
  }
  spectra_ref <- spectra_at(rt_ref, reference)

  use_spectra <- nrow(spectra_ref) > 1
  if (!use_spectra)
    warning("Spectral matching is unavailable for one-dimensional ",
            "chromatograms, so peaks are matched on retention time",
            if (use_shape) " and peak shape" else "", " alone. ",
            "`spectral_cutoff` is ignored.", call. = FALSE)

  sample_indices <- seq_along(chrom_list)
  names_peaks <- names(peak_list)
  matches_list <- lapply(sample_indices, function(idx) {
    if (verbose){
      message(sprintf("Comparing %s against reference (%s)",
                      sQuote(names_peaks[idx]), sQuote(names_peaks[reference])))
    }
    pks_samp <- pks_flat[[idx]]
    pks_samp <- pks_samp[!is.na(pks_samp$rt),]
    rt_samp <- pks_samp$rt

    if (length(rt_samp) == 0) {
      warning(sprintf("No peaks found for sample %d ('%s'). Skipping.",
                      idx, names(chrom_list)[idx]), call. = FALSE)
      return(NULL)
    }

    if (use_spectra){
      cor_matrix <- cor(spectra_at(rt_samp, idx), spectra_ref)
      spectral_mask <- cor_matrix >= spectral_cutoff
    } else {
      cor_matrix <- matrix(1, nrow = length(rt_samp), ncol = length(rt_ref),
                           dimnames = list(rt_samp, rt_ref))
      spectral_mask <- TRUE
    }

    rt_dist_matrix <- abs(outer(rt_samp, rt_ref, "-"))
    rt_mask <- rt_dist_matrix <= rt_window_rt

    shape_dist_matrix <- NULL
    if (use_shape) {
      shape_dist_matrix <- shape_distance(pks_samp, pks_ref, shape_cols)
      shape_mask <- if (!is.null(shape_cutoff)) shape_dist_matrix <= shape_cutoff else TRUE
    } else {
      shape_mask <- TRUE
    }

    combined_mask <- spectral_mask & shape_mask & rt_mask
    cor_matrix[!combined_mask] <- NA

    if (sum(!is.na(cor_matrix)) == 0) {
      warning(sprintf("No peaks matched for sample %d ('%s'). Skipping.",
                      idx, names(chrom_list)[idx]), call. = FALSE)
      return(NULL)
    }

    score_matrix <- cor_matrix
    if (rt_weight > 0)
      score_matrix <- score_matrix * (1 - rt_weight * rt_dist_matrix / rt_window_rt)
    if (shape_weight > 0 && !is.null(shape_dist_matrix)) {
      norm <- if (!is.null(shape_cutoff)) shape_cutoff else max(shape_dist_matrix, na.rm = TRUE)
      if (norm > 0)
        score_matrix <- score_matrix * (1 - shape_weight * pmin(shape_dist_matrix / norm, 1))
    }

    pairs <- which(!is.na(cor_matrix), arr.ind = TRUE)
    pairs <- pairs[order(score_matrix[pairs], decreasing = TRUE), , drop = FALSE]
    used_ref <- logical(length(rt_ref))
    used_samp <- logical(length(rt_samp))
    mr <- ms <- integer(nrow(pairs))
    n <- 0L

    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k, 2]  # ref (col)
      j <- pairs[k, 1]  # samp (row)
      if (!used_ref[i] && !used_samp[j]) {
        n <- n + 1L
        mr[n] <- i
        ms[n] <- j
        used_ref[i] <- TRUE
        used_samp[j] <- TRUE
      }
    }
    matched_rows <- mr[seq_len(n)]
    matched_cols <- ms[seq_len(n)]

    matches <- cbind(rt_ref = unname(rt_ref[matched_rows]),
                     rt_samp = unname(rt_samp[matched_cols]))
    pair_idx      <- cbind(matched_cols, matched_rows)
    cor_scores    <- cor_matrix[pair_idx]
    combined_scores <- if (rt_weight > 0 || shape_weight > 0)
      score_matrix[pair_idx] else NULL

    to_idx <- function(rt) (rt - t_start) / time_res + 1
    ref_i  <- cbind(rt = to_idx(matches[, "rt_ref"]),
                    I  = pks_ref$height[matched_rows])
    samp_i <- cbind(rt = to_idx(matches[, "rt_samp"]),
                    I  = pks_samp$height[matched_cols])

    list(matches = matches, ref_i = ref_i, samp_i = samp_i,
         cor_scores = cor_scores, combined_scores = combined_scores,
         rt_samp = rt_samp, cor_matrix = cor_matrix,
         matched_rows = matched_rows, matched_cols = matched_cols)
  })

  names(matches_list) <- names_peaks
  result <- structure(matches_list, reference = names_peaks[reference],
                      rt_ref = rt_ref, lambda = lambda_idx,
                      spectral_cutoff = spectral_cutoff,
                      rt_window = rt_window_rt,
                      class = "matched_peak_list")
  if (plot_it) plot(result)
  result
}

#' Stack chromatograms into a wavelengths x time x samples array for ptw
#' @noRd
stack_chroms <- function(chrom_list, lambdas){
  sapply(chrom_list, function(x) t(x[, lambdas, drop = FALSE]), simplify = "array")
}

#' Position of the reference sample in a `matched_peak_list`
#' @noRd
reference_position <- function(matched){
  ref <- attr(matched, "reference")
  if (is.null(ref)) return(NULL)
  pos <- match(ref, names(matched))
  if (is.na(pos))
    stop("The reference sample (", ref, ") is not among the samples of the ",
         "`matched_peak_list`. Every match is relative to it, so keep it when ",
         "subsetting.")
  pos
}

#' Choose the best reference chromatogram
#' @noRd
resolve_best_reference <- function(allmats.t, verbose = FALSE){
  best <- ptw::bestref(allmats.t)
  reference <- as.numeric(names(sort(table(best$best.ref),
                                     decreasing = TRUE))[1])
  if (verbose)
    message(paste("Selected chromatogram", reference, "as best reference."))
  reference
}

#' Euclidean distance between the shape parameters of sample and reference peaks
#' @param p1,p2 Peak tables for the sample and reference peaks.
#' @param cols Shape-parameter columns to compare.
#' @return A matrix of relative distances, sample peaks by reference peaks.
#' @noRd
shape_distance <- function(p1, p2, cols){
  rel2 <- function(a, b)
    (outer(a, b, "-") / outer(a, b, function(u, v) (u + v) / 2))^2
  sqrt(Reduce(`+`, lapply(cols, function(k) rel2(p1[[k]], p2[[k]]))))
}

#' Subset matched peak list
#'
#' Subsets a `matched_peak_list` while preserving metadata attributes.
#'
#' @param x A `matched_peak_list` object.
#' @param i Indices specifying elements to extract.
#' @param j Not used.
#' @param ... Additional arguments (ignored).
#' @param drop Not used.
#' @return A `matched_peak_list` object.
#' @author Ethan Bass
#' @method [ matched_peak_list
#' @export
`[.matched_peak_list` <- function(x, i, j, ..., drop = FALSE) {
  transfer_metadata(.subset(x, i), x,
                    exclude = c('names', 'row.names', 'dim', 'dimnames'))
}
