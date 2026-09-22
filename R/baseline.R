#' Baseline estimators
#'
#' Reference for the baseline estimators available through the `baseline`
#' argument of [`preprocess`] and [`correct_baseline`].
#'
#' @section Baseline estimators:
#' Four estimators are available. Which performs best on a given chromatogram is an empirical question that
#' depends on the noise, the drift, and how crowded the trace is, and is not
#' one this page can settle for your data.
#'
#' * `"als"` (the default) is asymmetric least squares, as implemented by
#' [`asysm`][ptw::asysm]. A Whittaker smoother is fit to the signal with
#' residuals above the fit weighted by `p` and residuals below it by `1 - p`.
#' `p` is fixed for the whole trace, so the asymmetry is set in advance rather
#' than taken from the data.

#' * `"arpls"` is asymmetrically reweighted penalized least squares (Baek et
#' al. 2015). It uses the same Whittaker solve, but replaces the fixed `p` with
#' weights recomputed at each iteration from the mean and standard deviation of
#' the negative residuals, so the asymmetry follows the trace. That also gives
#' it a failure mode on very clean traces: as the noise falls, the standard
#' deviation the weights are scaled by collapses, the weighting becomes a hard
#' one-sided clamp, and the fit can no longer be pulled up to follow a steep
#' drift, detaching from the trace where the drift is steepest. This is a
#' property of the published method rather than of this implementation —
#' `pybaselines` reproduces it — and `"aspls"` does not share it.

#' * `"aspls"` is adaptive smoothness parameter penalized least squares (Zhang
#' et al. 2020). It keeps the re-weighting and additionally scales the penalty
#' at each timepoint by the normalized absolute residual, so the fit is
#' stiffest where residuals are largest and most flexible where they are
#' smallest. Since that scaling is at most one, it lowers the effective
#' penalty: the useful `stiffness` range sits about 100x above the other two,
#' and a value carried over from `"als"` or `"arpls"` will under-smooth.

#' * `"snip"` is Statistics-sensitive Non-linear Iterative Peak-clipping
#' (Ryan et al. 1988), an iterated minimum filter applied in log-log space.
#' Two of its properties follow from its construction rather than from any
#' measurement: the estimate can never exceed the signal, so the corrected
#' trace is non-negative; and each point depends only on the `2m + 1` window
#' around it, so a distortion in one part of the trace cannot propagate along
#' it the way it can through a global solve. Unlike the others it has no
#' default width and requires `m`.
#'
#' One difference between them is directional rather than quantitative, and
#' follows from the loss functions. `"als"` and `"snip"` are both one-sided:
#' the first weights positive residuals down by a fixed factor, the second
#' cannot exceed the signal at all. Both therefore tend to fit below the center
#' of the noise, which leaves the corrected trace biased high. That offset is
#' absolute, so it matters in proportion to the size of the peak it sits under.
#' `"arpls"` and `"aspls"` take their asymmetry from the residuals instead, so
#' they have no fixed offset of this kind. How large any of this is in practice
#' is a property of the chromatogram, not of the estimator.
#'
#' All four depend on their smoothing parameter, and a value that suits one
#' chromatogram need not suit another. Rather than choosing an
#' estimator on general grounds, fit a baseline and look at it:
#' `attach_baseline = TRUE` in [`preprocess`] retains the estimate for exactly
#' that purpose.
#'
#' A baseline is estimated separately for each wavelength, since the drift
#' differs across the spectrum.
#'
#' @section Acknowledgment:
#' The `"aspls"` implementation follows `pybaselines` (Erb 2025) at two points
#' that Zhang et al. leave open: the default asymmetric coefficient `k` is
#' `0.5` rather than the `2` printed in the paper, and the adaptive penalty is
#' applied by left-multiplying the second-difference matrix by `diag(alpha)`,
#' an ordering the paper's notation does not settle.
#'
#' @name baseline_estimators
#' @author Ethan Bass
#' @seealso [`preprocess`], [`correct_baseline`], [`get_baseline`]
#' @references
#' * Baek, S.-J., Park, A., Ahn, Y.-J., and Choo, J. 2015. Baseline correction
#' using asymmetrically reweighted penalized least squares smoothing.
#' *Analyst* 140:250-257. \doi{10.1039/C4AN01061B}.
#' * Zhang, F., Tang, X., Tong, A., Wang, B., Wang, J., Lv, Y., Tang, C., and
#' Wang, J. 2020. Baseline correction for infrared spectra using adaptive
#' smoothness parameter penalized least squares method. *Spectroscopy Letters*
#' 53:222-233. \doi{10.1080/00387010.2020.1730908}.
#' * Ryan, C.G., Clayton, E., Griffin, W.L., Sie, S.H., and Cousens, D.R. 1988.
#' SNIP, a statistics-sensitive background treatment for the quantitative
#' analysis of PIXE spectra in geoscience applications. *Nuclear Instruments
#' and Methods in Physics Research B* 34:396-402.
#' \doi{10.1016/0168-583X(88)90063-8}.
#' * Erb, D. 2025. *pybaselines: A Python library of algorithms for the
#' baseline correction of experimental data*.
#' \doi{10.5281/zenodo.5608581}.
#' @md
NULL

#' Smoothing parameter an estimator will actually use
#'
#' @param baseline Estimator name.
#' @param dots Arguments passed to the estimator, which may name the parameter
#' directly.
#' @param dt Sampling interval the parameter was applied against.
#' @return A length-one named numeric carrying the estimator name as a
#' `"method"` attribute and the sampling interval as a `"dt"` attribute, or
#' `NULL` if the estimator has no such parameter.
#' @author Ethan Bass
#' @noRd
effective_baseline_param <- function(baseline, dots = list(), dt = NA_real_){
  nm <- if (baseline == "snip") "m" else "stiffness"
  param <- if (nm %in% names(dots)){
    setNames(as.numeric(dots[[nm]]), nm)
  } else if (baseline != "snip"){
    setNames(default_stiffness(baseline), nm)
  } else NULL
  if (is.null(param)) return(NULL)
  structure(param, method = baseline, dt = dt)
}

#' Default smoothing penalty
#'
#' @param baseline Estimator name.
#' @return A single numeric.
#' @author Ethan Bass
#' @noRd
default_stiffness <- function(baseline){
  if (baseline == "aspls") ASPLS_DEFAULT_STIFFNESS else DEFAULT_STIFFNESS
}

#' @noRd
baseline_fn <- function(baseline){
  switch(baseline, als = als_baseline, arpls = arpls_baseline,
         aspls = aspls_baseline, snip = snip_baseline)
}

#' @noRd
als_baseline <- function(y, stiffness = 1e7, p = 0.001, eps = 1e-8,
                         maxit = 25){
  ptw::asysm(y, lambda = stiffness, p = p, eps = eps, maxit = maxit)
}

#' @noRd
baseline_params <- function(baseline){
  setdiff(names(formals(baseline_fn(baseline))), "y")
}

#' @noRd
check_baseline_dots <- function(dots, baseline){
  valid <- baseline_params(baseline)
  unknown <- setdiff(names(dots), valid)
  if (length(unknown) == 0) return(invisible(NULL))
  estimators <- c("als", "arpls", "aspls", "snip")
  hints <- vapply(unknown, function(nm){
    owners <- estimators[vapply(estimators,
                                function(b) nm %in% baseline_params(b),
                                logical(1))]
    if (length(owners) == 0) return(sQuote(nm))
    sprintf("%s (an argument of %s)", sQuote(nm),
            paste(sprintf("`baseline = \"%s\"`", owners), collapse = ", "))
  }, character(1))
  stop(sprintf("`baseline = \"%s\"` does not accept %s.\n  It accepts: %s.",
               baseline, paste(hints, collapse = ", "),
               paste(sQuote(valid), collapse = ", ")))
}

#' Estimate baseline
#'
#' Estimates a chromatographic baseline along the time axis using one of four
#' estimators. Called by `preprocess` and `correct_baseline`. The user-facing
#' description of the estimators is the `baseline_estimators` topic.
#'
#' @param y A numeric vector (a chromatographic signal).
#' @param baseline Estimator to use: `"als"`, `"arpls"`, `"aspls"`, or
#' `"snip"`.
#' @param dt Sampling interval, in the time units of the chromatogram.
#' @param ... Additional arguments to the individual estimators
#' ([`asysm`][ptw::asysm] for `"als"`, `arpls_baseline` for `"arpls"`,
#' `aspls_baseline` for `"aspls"`, or `snip_baseline` for `"snip"`).
#' @return A numeric vector of the same length as `y` containing the estimated
#' baseline.
#' @author Ethan Bass
#' @noRd
estimate_baseline <- function(y, baseline = c("als", "arpls", "aspls", "snip"),
                              ..., dt = NULL){
  baseline <- match.arg(baseline, c("als", "arpls", "aspls", "snip"))
  dots <- list(...)
  check_baseline_dots(dots, baseline)
  if (length(y) < MIN_BASELINE_POINTS)
    stop(sprintf(paste(
      "Baseline estimation needs at least %d timepoints, but this trace has",
      "%d.\n  A second-difference penalty has nothing to constrain it over so",
      "short a span, and `ptw::asysm` does not survive a single point."),
      MIN_BASELINE_POINTS, length(y)))
  if (!all(is.finite(y)))
    stop("Baseline estimation requires a complete trace, but `y` contains ",
         "missing or non-finite values. Interpolate across the gap or trim it ",
         "before calling `baseline = \"", baseline, "\"`.")
  if (baseline == "snip"){
    if (is.null(dots$m))
      stop("`baseline = \"snip\"` requires `m`, the filter half-width in ",
           "timepoints. It has no default: the window has to be wider than ",
           "the peaks it clips, which depends on your chromatography.")
  } else if (is.null(dots$stiffness)){
    dots$stiffness <- default_stiffness(baseline)
  }
  do.call(baseline_fn(baseline), c(list(y), dots))
}

#' Asymmetrically reweighted penalized least squares (arPLS) baseline
#'
#' Fits a baseline by iteratively reweighted Whittaker smoothing, where the
#' weights are derived from the distribution of the negative residuals rather
#' than from a fixed asymmetry parameter.
#'
#' @details
#' At each iteration a Whittaker smoother ([`whit2`][ptw::whit2]) is fit to `y`
#' with the current weights, and the residuals `d <- y - z` are computed. The
#' mean `m` and standard deviation `s` of the *negative* residuals describe the
#' noise below the current fit. New weights are then assigned by a logistic
#' function centered at `-m + 2 * s`, so that points well above the fit (likely
#' peak) receive weights near zero while points scattered symmetrically about it
#' (likely noise) receive weights near one. Iteration stops when the relative
#' change in the weight vector falls below `tol`.
#'
#' @param y A numeric vector (a chromatographic signal).
#' @param stiffness Smoothing penalty (the \eqn{\lambda} of the Whittaker
#' smoother). Larger values give a stiffer baseline.
#' Defaults to `1e7`. The fit depends on it and it is worth setting
#' deliberately.
#' @param tol Relative convergence tolerance on the weight vector. Defaults to
#' `1e-6`.
#' @param maxit Maximum number of iterations. Defaults to `50`.
#' @return A numeric vector of the same length as `y` containing the estimated
#' baseline.
#' @author Ethan Bass
#' @references Baek, S.-J., Park, A., Ahn, Y.-J., and Choo, J. 2015. Baseline
#' correction using asymmetrically reweighted penalized least squares smoothing.
#' *Analyst* 140:250-257. \doi{10.1039/C4AN01061B}.
#' @noRd
arpls_baseline <- function(y, stiffness = 1e7, tol = 1e-6, maxit = 50){
  w <- rep(1, length(y))
  z <- y
  for (i in seq_len(maxit)){
    z <- ptw::whit2(y, lambda = stiffness, w = w)
    d <- y - z
    dn <- d[d < 0]
    if (length(dn) < 2) break
    m <- mean(dn)
    s <- sd(dn)
    if (!is.finite(s) || s == 0) break
    wt <- 1 / (1 + exp(2 * (d - (-m + 2 * s)) / s))
    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < tol) break
    w <- wt
  }
  z
}

#' Adaptive smoothness parameter penalized least squares (asPLS) baseline
#'
#' Fits a baseline by iteratively reweighted Whittaker smoothing in which the
#' smoothing penalty itself varies along the trace, so the baseline is stiff
#' across peaks and flexible where the trace is peak-free.
#'
#' @details
#' asPLS augments the arPLS system with a per-timepoint smoothness vector
#' `alpha`, solving `(W + stiffness * diag(alpha) %*% t(D) %*% D) z = W y`,
#' where `alpha <- abs(d) / max(abs(d))` is the normalized absolute residual from the
#' previous iteration. Because `alpha` is near one at peak apices and near zero
#' on peak-free stretches, the penalty is strongest exactly where the baseline
#' must not bend into a peak.
#'
#' Weights use the asPLS logistic `w <- 1 / (1 + exp(k * (d - s) / s))`, where
#' `s` is the standard deviation of the negative residuals. Note that this is
#' centered on `s` rather than arPLS's `-m + 2 * s`, so it is a distinct
#' weighting, not just a reparameterization.
#'
#' Dividing that system through by `stiffness * alpha` gives
#' `(diag(w / (stiffness * alpha)) + t(D) %*% D) z = diag(w / (stiffness *
#' alpha)) y`, which is exactly the system [`whit2`][ptw::whit2] solves at
#' `lambda = 1` with weights `w / (stiffness * alpha)`. Using that identity
#' keeps the whole solve in `ptw`'s compiled banded code, which avoids a sparse
#' solve and an additional dependency. `alpha` is floored at `alpha_min` so the
#' reciprocal
#' stays finite where a residual is exactly zero.
#'
#' @param y A numeric vector (a chromatographic signal).
#' @param stiffness Smoothing penalty (the \eqn{\lambda} of the Whittaker
#' smoother). Defaults to `1e9`, 100x the `"arpls"` default, because
#' `alpha <= 1` scales the penalty down, shifting the useful range up
#' accordingly — an `"arpls"` `stiffness` will under-smooth here.
#' @param k Asymmetric coefficient of the weighting logistic. Larger values
#' give a more step-like weight curve. Defaults to `0.5` rather than the `2`
#' printed in Zhang et al., following `pybaselines`, which found `0.5` both
#' reproduces the paper's own Table 2 / Figure 5 more closely and behaves
#' better on noisy data.
#' @param tol Relative convergence tolerance on the weight vector. Defaults to
#' `1e-3`.
#' @param maxit Maximum number of iterations. Defaults to `100`.
#' @param alpha_min Floor on the adaptive smoothness vector, guarding the
#' reciprocal used in the `whit2` reformulation. Defaults to `1e-9`.
#' @return A numeric vector of the same length as `y` containing the estimated
#' baseline.
#' @author Ethan Bass
#' @references Zhang, F., Tang, X., Tong, A., Wang, B., Wang, J., Lv, Y., Tang,
#' C., and Wang, J. 2020. Baseline correction for infrared spectra using
#' adaptive smoothness parameter penalized least squares method. *Spectroscopy
#' Letters* 53:222-233. \doi{10.1080/00387010.2020.1730908}.
#' @noRd
aspls_baseline <- function(y, stiffness = 1e9, k = 0.5, tol = 1e-3, maxit = 100,
                           alpha_min = 1e-9){
  n <- length(y)
  w <- rep(1, n)
  alpha <- rep(1, n)
  z <- y
  for (i in seq_len(maxit)){
    z <- ptw::whit2(y, lambda = 1,
                    w = w / (stiffness * pmax(alpha, alpha_min)))
    d <- y - z
    dn <- d[d < 0]
    if (length(dn) < 2) break
    s <- sd(dn)
    if (!is.finite(s) || s == 0) break
    wt <- 1 / (1 + exp(k * (d - s) / s))
    ad <- abs(d)
    mx <- max(ad)
    if (mx > 0) alpha <- ad / mx
    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < tol) break
    w <- wt
  }
  z
}

#' SNIP baseline
#'
#' Estimates a baseline by Statistics-sensitive Non-linear Iterative
#' Peak-clipping: an iterated minimum filter of increasing width applied to a
#' log-log transform of the signal.
#'
#' @details
#' The signal is transformed with the LLS operator
#' `log(log(sqrt(s + 1) + 1) + 1)`, filtered by replacing each point with the
#' minimum of itself and the mean of the points `i` positions to either side
#' for `i` in `1:m`, and transformed back with the exact inverse
#' `(exp(exp(t) - 1) - 1)^2 - 1`.
#'
#' The signal is padded with `m` copies of its first and last value before
#' filtering, and the padding discarded afterwards. Without padding the first
#' and last `m` timepoints are never filtered and the estimate is wrong there.
#' Edge-value replication was chosen over extrapolating the trace, which is
#' unstable at the ends.
#'
#' Because each step takes a minimum, the estimate can never exceed the signal,
#' so `max(baseline - y) <= 0` always holds and the corrected trace requires no
#' clipping.
#'
#' @param y A numeric vector (a chromatographic signal).
#' @param m Half-width of the filter, in timepoints. Required; there is no
#' default, since a usable width depends on the chromatography.
#' @return A numeric vector of the same length as `y` containing the estimated
#' baseline.
#' @author Ethan Bass
#' @references Ryan, C.G., Clayton, E., Griffin, W.L., Sie, S.H., and Cousens,
#' D.R. 1988. SNIP, a statistics-sensitive background treatment for the
#' quantitative analysis of PIXE spectra in geoscience applications. *Nuclear
#' Instruments and Methods in Physics Research B* 34:396-402.
#' \doi{10.1016/0168-583X(88)90063-8}.
#' @noRd
snip_baseline <- function(y, m){
  n <- length(y)
  m <- as.integer(max(1L, min(m, floor((n - 1) / 2))))

  ## pad with edge values so the first and last `m` points are also filtered
  s <- c(rep(y[1], m), y, rep(y[n], m))
  ## The LLS transform requires a non-negative argument. Shift rather than
  ## clamp: clamping would put the estimate above the signal wherever the
  ## signal is negative, breaking the non-negativity guarantee that is the
  ## whole point of this estimator. A constant shift is undone exactly below,
  ## and the transform is monotone, so `baseline <= signal` still holds.
  shift <- min(0, min(s))
  s <- s - shift
  ## LLS transform
  tf <- log(log(sqrt(s + 1) + 1) + 1)
  np <- length(tf)
  for (i in seq_len(m)){
    idx <- (i + 1):(np - i)
    tf[idx] <- pmin(tf[idx], (tf[idx - i] + tf[idx + i]) / 2)
  }
  ## exact inverse LLS transform
  bl <- (exp(exp(tf) - 1) - 1)^2 - 1 + shift
  bl[(m + 1):(m + n)]
}

#' Standard deviation by iterative sigma clipping
#' @param z A numeric vector.
#' @param k Multiple of the standard deviation at which to clip.
#' @param tol Relative convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param min_sd Value returned when the clipped standard deviation collapses
#' to (or below) it.
#' @return The clipped standard deviation.
#' @author Ethan Bass
#' @noRd
sigma_clipped_sd <- function(z, k = 4, tol = 1e-5, max_iter = 50, min_sd = 0){
  z <- z[is.finite(z)]
  if (length(z) < 2) return(min_sd)
  # `sd` is measured about the mean, so the clip has to be centered too:
  # clipping about zero collapses to `min_sd` on any trace with an offset
  z <- z - median(z)
  old <- 0
  for (i in seq_len(max_iter)){
    s <- sd(z)
    if (!is.finite(s) || s <= min_sd) return(min_sd)
    if (abs(old - s) / s < tol) return(s)
    z <- pmin(pmax(z, -k * s), k * s)
    old <- s
  }
  sd(z)
}

## These defaults are grid-dependent: the penalty applies on the index axis, so
## the same value spans a different amount of time on a differently sampled
## run. They are starting points, not calibrated optima.
##
## asPLS multiplies the penalty by `alpha = |d| / max|d| <= 1`, so its useful
## range sits ~2 decades higher than the uniform-penalty estimators'.
DEFAULT_STIFFNESS <- 1e7
ASPLS_DEFAULT_STIFFNESS <- 1e9

## Below four points the Whittaker estimators' second-difference penalty is
## unconstrained (`asysm` returns 0 for two points), and `ptw::asysm` corrupts
## memory outright on a single point rather than raising an error.
MIN_BASELINE_POINTS <- 4L

