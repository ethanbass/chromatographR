#' Calculate peak purity
#' 
#' Estimates peak purity by assessing spectral dissimilarity across the peak
#' region, following Stahl (2003).
#' 
#' At each timepoint in the peak, the spectrum is compared to a reference
#' spectrum, by default the mean of the assessed spectra. The resulting
#' dissimilarity (\eqn{1 - r^2}) is divided by the dissimilarity that noise
#' alone could explain, one minus the threshold returned by
#' [`get_agilent_threshold`]. A ratio below one means the spectrum is as similar
#' to the reference as the noise permits. With the default
#' `metric = "fraction"`, a timepoint only counts as impure if it belongs to a
#' run of at least `min_run` consecutive failing timepoints, so isolated
#' failures are ignored. Stahl requires a run of three.
#' 
#' Spectra are background-corrected before comparison by subtracting the
#' spectra at the peak boundaries, interpolated across the peak
#' (`background = "linear"`) or averaged (`background = "constant"`). Anything
#' that changes linearly across the peak is removed with the background,
#' including a co-eluting impurity much broader than the peak. If
#' [`preprocess`] has removed the baseline completely, `background = "none"`
#' avoids this, but any remaining baseline then enters the comparison as part of
#' each spectrum.
#'
#' With the default `metric = "fraction"`, pure peaks do not always score
#' exactly one, so treat purity below about 0.9, rather than any value below
#' one, as a sign of impurity. Wide, tailing and weak peaks fall below 0.9 more
#' often than others, even when pure. So do peaks above about 2 AU, whose
#' spectra near the apex are distorted by the detector's non-linear response.
#' The purity calculation is still being refined, and its values may change in
#' future releases.
#'
#' @param x A chromatographic matrix (timepoints × wavelengths).
#' @param pos A length-3 integer vector `c(apex, start, end)` giving the row
#' indices of the peak apex, left boundary, and right boundary in `x`.
#' @param noise_inflation Multiplier applied to the noise variance, whether
#' supplied via `noise_variance` or estimated. Values above `1` treat the noise
#' as larger than measured, which loosens the threshold. Defaults to `1`.
#' @param cutoff Timepoints at or below `cutoff` times the background-corrected
#' apex height are excluded, measured at the wavelength of maximum absorbance at
#' the apex. Defaults to `0.05`.
#' @param noise_variance Variance of noise. If provided, it replaces the
#' estimate, and the `noise_*` arguments that control estimation are ignored.
#' Estimating the noise variance requires a pass over the whole chromatogram,
#' so supplying it is much faster when many peaks are assessed.
#' @param noise_method How to select the timepoints used to estimate the noise
#' variance. Either `"quantile"` (the default), which takes the timepoints in
#' the lowest `noise_quantile` of maximum absorbance, or `"threshold"`, the
#' legacy behavior, which takes timepoints falling below `noise_threshold`
#' times the largest absorbance in the chromatogram.
#' @param noise_quantile Quantile of maximum absorbance below which timepoints
#' are treated as noise. Defaults to `0.1`.
#' @param noise_window Optional range of row indices to use as the noise
#' region. Takes precedence over `noise_method`.
#' @param noise_threshold Threshold to define noise, as a proportion of the
#' largest absorbance in the chromatogram. Used only when
#' `noise_method == "threshold"`, where `NULL` (the default) means `0.005`.
#' @param max_threshold Ceiling on the purity threshold. Without it, the
#' threshold approaches one for high signal-to-noise peaks, so the tolerated
#' dissimilarity approaches zero, and spectral variation that noise cannot
#' explain, such as detector non-linearity or a changing solvent composition,
#' makes such peaks impure. Corresponds to the user-defined similarity
#' threshold in the Agilent ChemStation, where `0.99` is a similarity of 990.
#' Defaults to `0.99`. Set to `1` to disable.
#' @param background How to correct the peak spectra for background
#' absorbance. Either `"linear"` (the default), `"constant"`, or `"none"`.
#' @param reference Which spectrum to compare the peak spectra against.
#' Either `"average"` (the default), the mean of the assessed spectra, or
#' `"apex"`, the spectrum at the peak apex. The apex spectrum correlates
#' perfectly with itself, so it forces a purity ratio of exactly zero at the
#' apex of every peak. Stahl recommends the average spectrum for unknown
#' impurities.
#' @param metric Statistic to summarize the purity ratios. Either
#' `"fraction"` (the default), the proportion of in-peak timepoints that are
#' pure, or `"max"`, the largest purity ratio. The run-length rule applies to
#' `"fraction"` only.
#' @param min_run Minimum number of consecutive failing timepoints required to
#' count as impure. Defaults to `3`.
#' @param lambdas Wavelengths to include in calculations. Character values are
#' matched to the column names of `x`; numeric values are taken as column
#' indices. Defaults to all wavelengths.
#' @param try Logical. If `FALSE`, returns `NA` without estimating purity.
#' Defaults to `TRUE`.
#' @param weight The `weight` argument is deprecated. Please use `noise_inflation`
#' instead.
#' @return Returns the purity of the peak specified by `pos`: with
#' `metric = "fraction"`, the proportion of timepoints above `cutoff` not
#' counted as impure; with `metric = "max"`, the largest purity ratio. Returns
#' `NA_real_` if the peak cannot be assessed, including when fewer than three
#' timepoints lie above `cutoff`. An error during the calculation is converted
#' to a warning and also returns `NA_real_`.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' [5988-8647EN.pdf](https://web.archive.org/web/20220615145615/https://www.agilent.com/cs/library/applications/5988-8647EN.pdf)
#' @author Ethan Bass
#' @keywords internal
#' @examples
#' data(Sa_warp)
#' x <- Sa_warp[[1]]
#' pos <- c(37, 22, 51)
#' get_purity(x, pos)
#' @export

get_purity <- function(x, pos, noise_inflation = 1, cutoff = 0.05,
                       noise_variance = NULL,
                       noise_method = c("quantile", "threshold"),
                       noise_quantile = 0.1, noise_window = NULL,
                       noise_threshold = NULL, max_threshold = 0.99,
                       background = c("linear", "constant", "none"),
                       reference = c("average", "apex"),
                       metric = c("fraction", "max"), min_run = 3,
                       lambdas, try = TRUE, weight = NULL){
  if (!try) return(NA_real_)
  metric <- match.arg(metric, c("fraction", "max"))
  noise_inflation <- resolve_deprecated(weight, noise_inflation)
  p <- tryCatch(
    get_purity_values(x, pos, noise_inflation = noise_inflation, cutoff = cutoff,
                      noise_variance = noise_variance,
                      noise_method = noise_method,
                      noise_quantile = noise_quantile,
                      noise_window = noise_window,
                      noise_threshold = noise_threshold,
                      max_threshold = max_threshold,
                      background = background, reference = reference,
                      lambdas = if (missing(lambdas)) NULL else lambdas),
    error = function(e){
      warning(sprintf("Could not estimate purity for the peak at index %s: %s",
                      pos[[1]], conditionMessage(e)), call. = FALSE)
      NULL
    })
  summarize_purity(p, metric = metric, min_run = min_run)
}

#' Summarize purity ratios
#' @param p Purity ratios, as returned by `get_purity_values`.
#' @param metric Either `"fraction"` or `"max"`.
#' @param min_run Minimum number of consecutive failing timepoints required to
#' count as impure.
#' @return Returns a single purity value, or `NA_real_` if fewer than three
#' timepoints can be assessed.
#' @author Ethan Bass
#' @noRd

summarize_purity <- function(p, metric = "fraction", min_run = 3){
  if (is.null(p)) return(NA_real_)
  v <- p[attr(p, "keep")]
  v <- v[!is.na(v)]
  if (length(v) < 3) return(NA_real_)
  if (metric == "max") return(max(v))
  fail <- v >= 1
  if (min_run > 1) fail <- fail & in_run(fail, min_run)
  mean(!fail)
}

#' Flag elements belonging to a run
#' @param v A logical vector.
#' @param k Minimum run length.
#' @return Returns a logical vector marking elements of `v` that belong to a
#' `TRUE` run of at least `k` elements.
#' @author Ethan Bass
#' @noRd

in_run <- function(v, k){
  r <- rle(v)
  rep(r$values & r$lengths >= k, r$lengths)
}

#' Prepare the spectra within a peak
#' 
#' Slices the peak window out of a chromatogram, drops wavelengths that are
#' incomplete within the window, and applies background correction.
#' 
#' @param x A chromatographic matrix (timepoints × wavelengths).
#' @param pos A length-3 integer vector `c(apex, start, end)`.
#' @param lambdas Wavelengths to include.
#' @param background Either `"linear"`, `"constant"`, or `"none"`.
#' @param cutoff Proportion of the corrected apex height below which timepoints
#' are excluded.
#' @param reference Either `"average"` or `"apex"`.
#' @return Returns a list containing the background-corrected window (`m`), the
#' reference spectrum (`reference`), the trace at the wavelength of maximum apex
#' absorbance (`trace`), the index of the apex within the window (`apex_i`), the
#' indices of the timepoints above `cutoff` (`keep`), and the correction
#' applied (`background`). Returns `NULL` if the peak cannot be
#' assessed.
#' @author Ethan Bass
#' @noRd

prep_purity_window <- function(x, pos, lambdas = NULL, background = "none",
                               cutoff = 0.05, reference = "average"){
  pos <- as.numeric(pos)
  if (length(pos) < 3 || anyNA(pos)) return(NULL)
  if (pos[2] < 1 || pos[3] > nrow(x) || pos[2] >= pos[3]) return(NULL)
  if (pos[1] < pos[2] || pos[1] > pos[3]) return(NULL)
  lambdas <- resolve_lambdas(x, lambdas)
  if (length(lambdas) < 3) return(NULL)
  m <- x[seq.int(pos[2], pos[3]), lambdas, drop = FALSE]
  m <- m[, colSums(!is.finite(m)) == 0, drop = FALSE]
  if (ncol(m) < 3 || nrow(m) < 3) return(NULL)
  apex_i <- as.integer(pos[1] - pos[2] + 1)
  if (background != "none"){
    edges <- m[c(1L, nrow(m)), , drop = FALSE]
    w <- (seq_len(nrow(m)) - 1) / (nrow(m) - 1)
    m <- m - if (background == "linear"){
      outer(1 - w, edges[1, ]) + outer(w, edges[2, ])
    } else {
      matrix(colMeans(edges), nrow(m), ncol(m), byrow = TRUE)
    }
  }
  apex <- m[apex_i, ]
  trace <- m[, which.max(apex)]
  keep <- which(trace > cutoff * trace[apex_i])
  ref <- if (reference == "apex" || length(keep) < 3){
    apex
  } else {
    colMeans(m[keep, , drop = FALSE])
  }
  list(m = m, reference = ref, trace = trace, apex_i = apex_i, keep = keep,
       background = background)
}

#' Resolve wavelengths to column indices
#' @param x A chromatographic matrix.
#' @param lambdas Wavelengths, as indices or as names. `NULL` selects all.
#' @return Returns a vector of column indices.
#' @author Ethan Bass
#' @noRd

resolve_lambdas <- function(x, lambdas = NULL){
  if (is.null(lambdas) || length(lambdas) == 0) return(seq_len(ncol(x)))
  if (is.character(lambdas)){
    return(which(get_lambdas(x) %in% as.numeric(lambdas)))
  }
  lambdas
}

#' Inflation of the noise variance due to background correction
#' @param background Either `"linear"`, `"constant"`, or `"none"`.
#' @return Returns the factor by which the noise variance of a spectrum is
#' increased by subtracting the boundary spectra.
#' @author Ethan Bass
#' @noRd

background_inflation <- function(background){
  switch(background, "linear" = 5/3, "constant" = 3/2, 1)
}

#' Calculate variance of noise regions
#' 
#' Estimates the across-wavelength variance of spectral noise regions used in
#' peak purity calculations.
#' 
#' @param x A chromatogram in matrix format
#' @param noise_method Either `"quantile"` or `"threshold"`.
#' @param noise_quantile Quantile of maximum absorbance below which timepoints
#' are treated as noise.
#' @param noise_window Optional range of row indices to use as the noise
#' region. Takes precedence over `noise_method`.
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance. Used only when `noise_method == "threshold"`.
#' @param lambdas Wavelengths to include.
#' @importFrom stats var quantile
#' @return Returns the average across-wavelength variance of the signal over
#' the retention times treated as noise, or `NA_real_` if no timepoints
#' qualify.
#' @note The variance is measured across wavelengths, not across time, since
#' that is the quantity the purity threshold compares against the variance of
#' each spectrum. It is therefore not directly comparable to the temporal noise
#' estimate returned by `estimate_noise`, and the two diverge on data that has
#' been smoothed along the wavelength axis by [`preprocess`].
#' @author Ethan Bass
#' @keywords internal

get_noise_variance <- function(x, noise_method = c("quantile", "threshold"),
                               noise_quantile = 0.1, noise_window = NULL,
                               noise_threshold = 0.005, lambdas){
  noise_method <- match.arg(noise_method, c("quantile", "threshold"))
  lambdas <- resolve_lambdas(x, if (missing(lambdas)) NULL else lambdas)
  xx <- x[, lambdas, drop = FALSE]
  if (min(xx, na.rm = TRUE) == 0 && mean(xx == 0, na.rm = TRUE) > 0.01){
    warning("The chromatogram appears to be clipped at zero, which biases the ",
            "noise variance low and inflates impurity. Preprocess with ",
            "`zero_floor = FALSE` for purity estimation.", call. = FALSE)
  }
  noise_idx <- if (!is.null(noise_window)){
    which(seq_len(nrow(xx)) >= min(noise_window) &
            seq_len(nrow(xx)) <= max(noise_window))
  } else if (noise_method == "threshold"){
    find_noise(xx, noise_threshold = noise_threshold,
               lambdas = seq_len(ncol(xx)))
  } else {
    find_noise_quantile(xx, noise_quantile = noise_quantile)
  }
  if (length(noise_idx) == 0) return(NA_real_)
  mean(apply(xx[noise_idx, , drop = FALSE], 1, var, na.rm = TRUE), na.rm = TRUE)
}

#' Define noise spectra by quantile of maximum absorbance
#' 
#' Unlike `find_noise`, the selection does not depend on the height of the
#' largest peak in the chromatogram and always returns at least one timepoint.
#' 
#' @param x A chromatogram in matrix format
#' @param noise_quantile Quantile of maximum absorbance below which timepoints
#' are treated as noise.
#' @return Returns indices of the lowest-intensity retention times.
#' @author Ethan Bass
#' @noRd

find_noise_quantile <- function(x, noise_quantile = 0.1){
  max_abs <- suppressWarnings(apply(x, 1, max, na.rm = TRUE))
  max_abs[!is.finite(max_abs)] <- NA
  if (all(is.na(max_abs))) return(integer(0))
  which(max_abs <= quantile(max_abs, noise_quantile, na.rm = TRUE))
}

#' Define noise spectra based on specified threshold
#' @param x A chromatogram in matrix format 
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance.
#' @param lambdas Wavelengths to include.
#' @return Returns indices of retention times where the signal falls below the
#' specified noise threshold.
#' @author Ethan Bass
#' @noRd

find_noise <- function(x, noise_threshold = 0.01, lambdas){
  max_abs <- apply(x[,lambdas, drop = FALSE], 1, max)
  which(max_abs < max(max_abs, na.rm = TRUE) * noise_threshold)
}

#' Calculate row variances
#' @param m A numeric matrix with no missing values.
#' @return Returns the variance of each row.
#' @author Ethan Bass
#' @noRd

row_var <- function(m){
  rowSums((m - rowMeans(m))^2) / (ncol(m) - 1)
}

#' Calculate purity thresholds
#' 
#' Computes the threshold used by [`get_purity`]: the similarity expected
#' between a spectrum and the reference if they differed only by noise, the
#' product of their reliabilities
#' \eqn{(1 - \sigma^2_{noise} / \sigma^2_i)(1 - \sigma^2_{noise} / \sigma^2_{ref})},
#' capped at `max_threshold`.
#' 
#' @param x A chromatographic matrix (timepoints × wavelengths).
#' @param pos A length-3 integer vector `c(apex, start, end)` giving the row
#' indices of the peak apex, left boundary, and right boundary in `x`.
#' @param noise_inflation Multiplier applied to the noise variance. Defaults to
#' `1`.
#' @param noise_variance Variance of noise.
#' @param noise_method How to select the timepoints used to estimate the noise
#' variance. Either `"quantile"` or `"threshold"`.
#' @param noise_quantile Quantile of maximum absorbance below which timepoints
#' are treated as noise.
#' @param noise_window Optional range of row indices to use as the noise
#' region.
#' @param noise_threshold Threshold to define noise. Used only when
#' `noise_method == "threshold"`, where `NULL` (the default) means `0.005`.
#' @param max_threshold Ceiling on the threshold. Defaults to `0.99`.
#' @param background How to correct the peak spectra for background
#' absorbance. Either `"linear"`, `"constant"`, or `"none"`.
#' @param reference Which spectrum to compare against. Either `"average"` (the
#' default) or `"apex"`. See [`get_purity`].
#' @param cutoff Proportion of the corrected apex height below which
#' timepoints are excluded from the average reference spectrum. Defaults to
#' `0.05`.
#' @param lambdas Wavelengths to include. Character values are matched to the
#' column names of `x`; numeric values are taken as column indices.
#' @param weight The `weight` argument is deprecated. Please use `noise_inflation`
#' instead.
#' @return Returns a vector of purity thresholds, one for each timepoint from
#' the start to the end of the peak specified by `pos`, or `NULL` if the peak
#' cannot be assessed.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' [5988-8647EN.pdf](https://web.archive.org/web/20220615145615/https://www.agilent.com/cs/library/applications/5988-8647EN.pdf)
#' @author Ethan Bass
#' @keywords internal
#' @examples
#' data(Sa_warp)
#' x <- Sa_warp[[1]]
#' pos <- c(37, 22, 51)
#' get_agilent_threshold(x, pos)
#' @export

get_agilent_threshold <- function(x, pos, noise_inflation = 1,
                                  noise_variance = NULL,
                                  noise_method = c("quantile", "threshold"),
                                  noise_quantile = 0.1, noise_window = NULL,
                                  noise_threshold = NULL, max_threshold = 0.99,
                                  background = c("linear", "constant", "none"),
                                  reference = c("average", "apex"),
                                  cutoff = 0.05, lambdas, weight = NULL){
  p <- prep_purity_inputs(x, pos, noise_inflation = noise_inflation, cutoff = cutoff,
                          noise_variance = noise_variance,
                          noise_method = noise_method,
                          noise_quantile = noise_quantile,
                          noise_window = noise_window,
                          noise_threshold = noise_threshold,
                          background = background,
                          reference = reference,
                          lambdas = if (missing(lambdas)) NULL else lambdas,
                          weight = weight)
  if (is.null(p)) return(NULL)
  agilent_threshold(p$w$m, p$w$reference, noise_inflation = p$noise_inflation,
                    var_noise = p$var_noise, max_threshold = max_threshold)
}

#' Resolve the noise variance for a purity calculation
#' @param x A chromatogram in matrix format.
#' @param noise_variance Variance of noise, or `NULL` to estimate it.
#' @param noise_method How to select the noise timepoints.
#' @param noise_quantile Quantile of maximum absorbance defining noise.
#' @param noise_window Optional range of row indices to use as noise.
#' @param noise_threshold Threshold to define noise.
#' @param lambdas Wavelengths to include.
#' @param background The background correction applied to the peak spectra.
#' @return Returns the noise variance, scaled to account for background
#' correction.
#' @author Ethan Bass
#' @noRd

resolve_noise_variance <- function(x, noise_variance = NULL,
                                   noise_method = c("quantile", "threshold"),
                                   noise_quantile = 0.1, noise_window = NULL,
                                   noise_threshold = NULL, lambdas = NULL,
                                   background = "none"){
  noise_method <- match.arg(noise_method, c("quantile", "threshold"))
  if (is.null(noise_variance)){
    if (!is.null(noise_threshold) && noise_method != "threshold" &&
        is.null(noise_window)){
      warning("`noise_threshold` is ignored unless `noise_method == \"threshold\"`.",
              call. = FALSE)
    }
    noise_variance <- get_noise_variance(
      x, noise_method = noise_method, noise_quantile = noise_quantile,
      noise_window = noise_window,
      noise_threshold = if (is.null(noise_threshold)) 0.005 else noise_threshold,
      lambdas = lambdas)
  }
  noise_variance * background_inflation(background)
}

#' Calculate the Agilent purity threshold for a prepared window
#' 
#' The threshold is the similarity expected between two noisy measurements of
#' the same underlying spectrum, which is the product of the two reliabilities
#' \eqn{(1 - \sigma^2_{noise} / \sigma^2_i)(1 - \sigma^2_{noise} / \sigma^2_{apex})}.
#' 
#' @param m A background-corrected peak window (timepoints × wavelengths).
#' @param reference The reference spectrum.
#' @param noise_inflation Multiplier applied to `var_noise`.
#' @param var_noise Variance of noise.
#' @param max_threshold Ceiling on the threshold.
#' @return Returns a vector of purity thresholds, one per timepoint in `m`.
#' @author Ethan Bass
#' @noRd

agilent_threshold <- function(m, reference, noise_inflation = 1, var_noise,
                              max_threshold = 0.99){
  reliability <- function(v){
    if (var_noise == 0) rep(1, length(v)) else pmax(0, 1 - noise_inflation * var_noise / v)
  }
  pmin(reliability(row_var(m)) * reliability(var(reference)), max_threshold)
}

#' Prepare the window and noise variance shared by the purity calculations
#'
#' Returns `NULL` when the window cannot be prepared or the noise variance is
#' not finite, so callers can return `NULL` in turn.
#' @noRd
prep_purity_inputs <- function(x, pos, noise_inflation, cutoff, noise_variance,
                               noise_method, noise_quantile, noise_window,
                               noise_threshold, background, reference,
                               lambdas, weight){
  noise_inflation <- resolve_deprecated(weight, noise_inflation)
  w <- prep_purity_window(x, pos, lambdas = lambdas,
                          background = match.arg(background,
                                                 c("linear", "constant", "none")),
                          cutoff = cutoff,
                          reference = match.arg(reference,
                                                c("average", "apex")))
  if (is.null(w)) return(NULL)
  var_noise <- resolve_noise_variance(
    x, noise_variance = noise_variance, noise_method = noise_method,
    noise_quantile = noise_quantile, noise_window = noise_window,
    noise_threshold = noise_threshold, lambdas = lambdas,
    background = w$background)
  if (!is.finite(var_noise)) return(NULL)
  list(w = w, var_noise = var_noise, noise_inflation = noise_inflation)
}

#' Calculate peak purity values
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param noise_inflation Multiplier applied to the noise variance. Defaults to
#' `1`.
#' @param cutoff Proportion of the background-corrected apex height below which
#' timepoints are excluded.
#' @param noise_variance Variance of noise.
#' @param noise_method How to select the timepoints used to estimate the noise
#' variance. Either `"quantile"` or `"threshold"`.
#' @param noise_quantile Quantile of maximum absorbance below which timepoints
#' are treated as noise.
#' @param noise_window Optional range of row indices to use as the noise
#' region.
#' @param noise_threshold Threshold to define noise. Used only when
#' `noise_method == "threshold"`.
#' @param max_threshold Ceiling on the threshold. Defaults to `0.99`.
#' @param background How to correct the peak spectra for background
#' absorbance. Either `"linear"`, `"constant"`, or `"none"`.
#' @param reference Which spectrum to compare against. Either `"average"` (the
#' default) or `"apex"`. See [`get_purity`].
#' @param lambdas Wavelengths to include in calculations.
#' @param weight The `weight` argument is deprecated. Please use `noise_inflation`
#' instead.
#' @importFrom stats cor var
#' @return Returns a vector of peak purity values at each timepoint within the
#' peak specified by `pos`, carrying a `keep` attribute giving the indices of
#' the timepoints above `cutoff`. Returns `NULL` if the peak cannot be
#' assessed.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' [5988-8647EN.pdf](https://web.archive.org/web/20220615145615/https://www.agilent.com/cs/library/applications/5988-8647EN.pdf)
#' @author Ethan Bass
#' @noRd

get_purity_values <- function(x, pos, noise_inflation = 1, cutoff = 0.05,
                              noise_variance = NULL,
                              noise_method = c("quantile", "threshold"),
                              noise_quantile = 0.1, noise_window = NULL,
                              noise_threshold = NULL, max_threshold = 0.99,
                              background = c("linear", "constant", "none"),
                              reference = c("average", "apex"),
                              lambdas = NULL, weight = NULL){
  p <- prep_purity_inputs(x, pos, noise_inflation = noise_inflation, cutoff = cutoff,
                          noise_variance = noise_variance,
                          noise_method = noise_method,
                          noise_quantile = noise_quantile,
                          noise_window = noise_window,
                          noise_threshold = noise_threshold,
                          background = background,
                          reference = reference, lambdas = lambdas,
                          weight = weight)
  if (is.null(p)) return(NULL)
  w <- p$w
  r <- suppressWarnings(as.vector(cor(w$reference, t(w$m))))
  num <- 1 - r^2
  den <- 1 - agilent_threshold(w$m, w$reference, noise_inflation = p$noise_inflation,
                               var_noise = p$var_noise,
                               max_threshold = max_threshold)
  out <- num / den
  degenerate <- which(den == 0)
  if (length(degenerate) > 0){
    out[degenerate] <- ifelse(num[degenerate] < sqrt(.Machine$double.eps), 0, Inf)
  }
  attr(out, "similarity") <- 1 - num
  attr(out, "threshold")  <- 1 - den
  attr(out, "keep") <- w$keep
  out
}

#' Trim peak
#' @param x A numeric vector of signal intensities (a single wavelength trace).
#' @param pos A length-3 integer vector `c(apex, start, end)` giving the row
#' indices of the peak apex, left boundary, and right boundary respectively.
#' @param cutoff Proportion of the apex intensity to use as cutoff. Defaults to
#' `0.05`.
#' @return Returns the positions within `start:end` (so `1` is `start`) where
#' the signal exceeds `cutoff` times the apex intensity.
#' @author Ethan Bass
#' @keywords internal
#' @examples
#' data(Sa_warp)
#' x <- Sa_warp[[1]]
#' signal <- x[, "210"]
#' pos <- c(37, 22, 51)
#' trim_peak(signal, pos)
#' @export

trim_peak <- function(x, pos, cutoff = 0.05){
  idx <- pos[2]:pos[3]
  which(x[idx] > cutoff*x[pos[1]])
}

#' Estimate the noise variance once for a whole chromatogram
#' 
#' Called by [`fit_peaks`] so that the noise variance is estimated once per
#' chromatogram rather than once per peak.
#' 
#' @param x A chromatogram in matrix format.
#' @param purity_args A list of arguments destined for [`get_purity`].
#' @return Returns the noise variance, before any adjustment for background
#' correction.
#' @author Ethan Bass
#' @noRd

hoist_noise_variance <- function(x, purity_args = list()){
  if (!is.null(purity_args$noise_variance)) return(purity_args$noise_variance)
  do.call(get_noise_variance,
          c(list(x), purity_args[intersect(names(purity_args),
                                           names(formals(get_noise_variance)))]))
}
