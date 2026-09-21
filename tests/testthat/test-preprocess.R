# Tests for the baseline estimators and the `preprocess` baseline plumbing.

sim_chrom <- function(dt = 0.01, front = FALSE, early = FALSE, seed = 1,
                      sd_pk = 0.12, drift_rate = 1,
                      step = 0, step_width = 1){
  set.seed(seed)
  tt <- seq(0, 30, by = dt)
  b <- 50 + 30 * exp(tt * drift_rate / 12) +
    step * plogis((tt - 15) / step_width)
  pk <- 1000 * dnorm(tt, 10, sd_pk) + 20 * dnorm(tt, 15, sd_pk) +
    5 * dnorm(tt, 20, sd_pk)
  if (front) pk <- pk + 20000 * dnorm(tt, 1.5, 0.05)
  if (early) pk <- pk + 500 * dnorm(tt, 0.4, 0.12)
  list(t = tt, b = b, y = b + pk + rnorm(length(tt)), dt = dt)
}
# FWHM of the simulated peaks, and a SNIP half-width that spans it
PW <- 2.355 * 0.12
SNIP_M <- floor(((2.25 * PW / 0.01) - 1) / 2)
rmse <- function(a, b) sqrt(mean((a - b)^2))


# ── Estimators ────────────────────────────────────────────────────────────────

test_that("estimate_baseline dispatches on the `baseline` argument", {
  s <- sim_chrom()
  expect_equal(estimate_baseline(s$y, "als"), ptw::asysm(s$y))
  expect_length(estimate_baseline(s$y, "arpls"), length(s$y))
  expect_length(estimate_baseline(s$y, "aspls"), length(s$y))
  expect_length(estimate_baseline(s$y, "snip", m = SNIP_M),
                length(s$y))
  # the snip width is chosen once per chromatogram by the matrix-level caller,
  # so this vector-level entry point requires it
  expect_error(estimate_baseline(s$y, "snip"), regexp = "requires `m`")
})

test_that("every estimator refuses an incomplete trace, rather than splicing it", {
  s <- sim_chrom()
  yn <- s$y
  yn[c(5, 500, 1000)] <- NA
  for (m in c("als", "arpls", "aspls"))
    expect_error(estimate_baseline(yn, m), regexp = "complete trace", info = m)
  expect_error(estimate_baseline(yn, "snip", m = SNIP_M),
               regexp = "complete trace")
  expect_error(estimate_baseline(rep(NA_real_, 50), "arpls"),
               regexp = "complete trace")
  expect_error(estimate_baseline(c(s$y[-1], Inf), "arpls"),
               regexp = "complete trace")
})

test_that("arPLS is far more accurate than AsLS on a known baseline", {
  skip_on_cran()
  s <- sim_chrom()
  r_als <- rmse(estimate_baseline(s$y, "als"), s$b)
  r_arp <- rmse(estimate_baseline(s$y, "arpls"), s$b)
  expect_lt(r_arp, r_als / 5)
})

test_that("arPLS removes the area bias that AsLS introduces on a small peak", {
  # AsLS tracks the lower envelope of the noise rather than its centre, so the
  # baseline sits low and every integrated area comes out high. The effect is
  # ruinous for a small peak.
  skip_on_cran()
  s <- sim_chrom()
  idx <- s$t > 19 & s$t < 21
  area <- function(m) sum((s$y - estimate_baseline(s$y, m))[idx]) * s$dt
  expect_gt(area("als"), 1.5 * 5)          # >50% high
  expect_equal(area("arpls"), 5, tolerance = 0.1)   # within 10%
})

test_that("SNIP's estimate never exceeds the signal", {
  skip_on_cran()
  s <- sim_chrom()
  expect_lte(max(estimate_baseline(s$y, "snip", m = SNIP_M) - s$y), 1e-8)
  # including where the signal is negative, which a naive zero-clamp breaks
  yn <- s$y - 200
  expect_lte(max(estimate_baseline(yn, "snip", m = SNIP_M) - yn), 1e-8)
})

test_that("SNIP's non-negativity guarantee holds on real data", {
  skip_on_cran()
  data(Sa)
  pr <- preprocess(Sa, dim1 = seq(10, 18.66, by = .01),
                   dim2 = seq(200, 318, by = 2),
                   baseline = "none", show_progress = FALSE, cl = 1)
  worst <- max(vapply(pr, function(x){
    max(apply(x, 2, function(y)
      max(estimate_baseline(y, "snip", m = 11) - y)))
  }, numeric(1)))
  expect_lte(worst, 1e-8)
})

test_that("SNIP is robust to a dominant solvent front", {
  # AsLS solves globally, so a saturated front distorts the baseline
  # everywhere; SNIP's minimum filter has bounded support and cannot be pulled.
  skip_on_cran()
  s <- sim_chrom(front = TRUE)
  away <- s$t > 3
  r_snip <- rmse(estimate_baseline(s$y, "snip", m = SNIP_M)[away], s$b[away])
  r_als <- rmse(estimate_baseline(s$y, "als")[away], s$b[away])
  expect_lt(r_snip, r_als)
})

test_that("SNIP pads with edge values so the trace start is usable", {
  # Without padding the first `m` timepoints are never filtered and the
  # first-minute RMSE blows up (5.5 vs ~1 with padding).
  skip_on_cran()
  s <- sim_chrom(early = TRUE)
  first <- s$t <= 1
  r <- rmse(estimate_baseline(s$y, "snip", m = SNIP_M)[first], s$b[first])
  expect_lt(r, 5)
})


test_that("each estimator falls back to its own default stiffness", {
  expect_equal(default_stiffness("als"), DEFAULT_STIFFNESS)
  expect_equal(default_stiffness("arpls"), DEFAULT_STIFFNESS)
  expect_equal(default_stiffness("aspls"), ASPLS_DEFAULT_STIFFNESS)
  s <- sim_chrom()
  expect_equal(estimate_baseline(s$y, "als"),
               als_baseline(s$y, stiffness = DEFAULT_STIFFNESS))
})

# ── preprocess plumbing ───────────────────────────────────────────────────────

test_that("preprocess accepts every baseline estimator", {
  skip_on_cran()
  data(Sa)
  args <- list(Sa[1], dim1 = seq(10, 18.66, by = .01),
               dim2 = seq(200, 318, by = 2), show_progress = FALSE, cl = 1)
  for (m in c("als", "arpls", "aspls")){
    x <- do.call(preprocess, c(args, baseline = m))[[1]]
    expect_equal(dim(x), c(867L, 60L), info = m)
    expect_true(all(is.finite(x)), info = m)
  }
  x <- do.call(preprocess, c(args, baseline = "snip", m = 11))[[1]]
  expect_equal(dim(x), c(867L, 60L))
  # snip cannot overshoot, so the corrected data needs no clipping
  expect_gte(min(x), -1e-8)
})

test_that("zero_floor clamps to zero and reproduces the historical output", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  new <- do.call(preprocess, c(list(Sa), args, zero_floor = TRUE))
  # the historical implementation: baseline.corr then an unconditional clip
  old <- lapply(Sa, function(x){
    m <- do.call(preprocess, c(list(x), args, baseline = "none"))
    z <- apply(m, 2, ptw::baseline.corr)
    z[z < 0] <- 0
    dimnames(z) <- dimnames(m)
    z
  })
  expect_equal(lapply(new, as.vector), lapply(old, as.vector),
               ignore_attr = TRUE)
  expect_gte(min(new[[1]]), 0)
  # negatives are retained by default
  expect_lt(min(do.call(preprocess, c(list(Sa[1]), args))[[1]]), 0)
})

test_that("attach_baseline retains a baseline that reconstructs the input", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  withbl <- do.call(preprocess, c(list(Sa), args, attach_baseline = TRUE))
  nobl <- do.call(preprocess, c(list(Sa), args, baseline = "none"))
  for (i in seq_along(withbl)){
    bl <- attr(withbl[[i]], "baseline", exact = TRUE)
    expect_equal(dim(bl), dim(withbl[[i]]))
    # `as.vector` because arithmetic carries the "baseline" attribute forward
    expect_equal(as.vector(withbl[[i]] + bl), as.vector(nobl[[i]]),
                 tolerance = 1e-8)
  }
  # the corrected data is unaffected by keeping the estimate
  plain <- do.call(preprocess, c(list(Sa), args))
  expect_equal(lapply(withbl, as.vector), lapply(plain, as.vector))
  # ...but clipping is applied after the subtraction, so it breaks the identity
  clipped <- do.call(preprocess, c(list(Sa[1]), args, attach_baseline = TRUE,
                                   zero_floor = TRUE))[[1]]
  expect_gt(max(abs(as.vector(clipped + attr(clipped, "baseline", exact = TRUE)) -
                      as.vector(nobl[[1]]))), 1)
})

test_that("attach_baseline warns and is ignored without baseline correction", {
  data(Sa)
  expect_warning(
    x <- preprocess(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
                    dim2 = seq(200, 318, by = 2), attach_baseline = TRUE,
                    baseline = "none", show_progress = FALSE),
    regexp = "ignored")
  expect_null(attr(x, "baseline", exact = TRUE))
})

test_that("the deprecated `remove.time.baseline` argument maps to `baseline`", {
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  expect_warning(off <- do.call(preprocess, c(list(Sa[[1]]), args,
                                              remove.time.baseline = FALSE)),
                 regexp = "deprecated")
  expect_equal(as.vector(off),
               as.vector(do.call(preprocess, c(list(Sa[[1]]), args,
                                               baseline = "none"))))
  expect_warning(on <- do.call(preprocess, c(list(Sa[[1]]), args,
                                             remove.time.baseline = TRUE)),
                 regexp = "deprecated")
  expect_equal(as.vector(on),
               as.vector(do.call(preprocess, c(list(Sa[[1]]), args))))
})

test_that("`remove.time.baseline` errors when it contradicts `baseline`", {
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  expect_error(suppressWarnings(
    do.call(preprocess, c(list(Sa[[1]]), args, remove.time.baseline = FALSE,
                          baseline = "arpls"))), regexp = "conflicts")
  expect_error(suppressWarnings(
    do.call(preprocess, c(list(Sa[[1]]), args, remove.time.baseline = TRUE,
                          baseline = "none"))), regexp = "conflicts")
  expect_error(
    do.call(preprocess, c(list(Sa[[1]]), args, remove.time.baseline = "yes")),
    regexp = "TRUE")
})


test_that("the smoothing parameter the estimator used is recorded", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 12),
               show_progress = FALSE, cl = 1)
  snipped <- do.call(preprocess, c(list(Sa[1]), args, baseline = "snip",
                                   m = 10))[[1]]
  expect_equal(attr(snipped, "baseline_param"),
               effective_baseline_param("snip", list(m = 10), dt = 0.01))
  smoothed <- do.call(preprocess, c(list(Sa[1]), args, baseline = "arpls"))[[1]]
  expect_equal(attr(smoothed, "baseline_param"),
               effective_baseline_param("arpls", dt = 0.01))

  raw <- do.call(preprocess, c(list(Sa[1]), args, baseline = "arpls",
                               stiffness = 1e7))[[1]]
  expect_equal(attr(raw, "baseline_param"),
               structure(c(stiffness = 1e7), method = "arpls", dt = 0.01))

  none_default <- do.call(preprocess, c(list(Sa[1]), args,
                                        baseline = "als"))[[1]]
  expect_equal(attr(none_default, "baseline_param"),
               structure(c(stiffness = 1e7), method = "als", dt = 0.01))

  none <- do.call(preprocess, c(list(Sa[1]), args, baseline = "none"))[[1]]
  expect_null(attr(none, "baseline_param"))

  from_cb <- correct_baseline(none, baseline = "snip", m = 10,
                              show_progress = FALSE, cl = 1)
  expect_equal(attr(from_cb, "baseline_param"),
               effective_baseline_param("snip", list(m = 10), dt = 0.01))
  expect_null(attr(preprocess(from_cb, dim1 = seq(11, 18, by = 0.02),
                              dim2 = seq(200, 318, by = 12), baseline = "none",
                              show_progress = FALSE, cl = 1),
                   "baseline_param"))
})

test_that("get_baseline retrieves the attached baseline", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  l <- do.call(preprocess, c(list(Sa[1:2]), args, attach_baseline = TRUE))
  bl <- get_baseline(l[[1]])
  expect_true(is.matrix(bl))
  expect_equal(dim(bl), dim(l[[1]]))
  expect_equal(dimnames(bl), dimnames(l[[1]]))
  # a list in, a list of the same class out, so it can go straight to plotting
  bls <- get_baseline(l)
  expect_length(bls, 2L)
  expect_equal(class(bls), class(l))
  expect_equal(bls[[1]], bl)
  expect_error(get_baseline(do.call(preprocess, c(list(Sa[[1]]), args))),
               regexp = "attach_baseline")
})

test_that("`p` and `m` reach the estimator, not `maxI`", {
  skip_on_cran()
  data(Sa)
  args <- list(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
               dim2 = seq(200, 318, by = 2), show_progress = FALSE, cl = 1)
  mild <- do.call(preprocess, c(args, p = .001, stiffness = 1e5))
  strong <- do.call(preprocess, c(args, p = .1, stiffness = 1e5))
  expect_gt(max(abs(mild - strong)), 1)
  snipped <- do.call(preprocess, c(args, baseline = "snip", m = 5))
  expect_gt(max(snipped), 5)
  expect_equal(max(do.call(preprocess, c(args, maxI = 100))), 100)
})

test_that("preprocess warns that it is bridging missing values", {
  skip_on_cran()
  data(Sa)
  X <- Sa[[1]]
  X[500:519, 5] <- NA
  expect_warning(
    p <- preprocess(X, dim1 = seq(10, 18.66, by = .01),
                    dim2 = seq(200, 318, by = 2), show_progress = FALSE),
    regexp = "longest run 20 timepoints")
  expect_false(anyNA(p))
  expect_silent(preprocess(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
                           dim2 = seq(200, 318, by = 2),
                           show_progress = FALSE))
})

test_that("the bridging warning only claims what will actually happen", {
  skip_on_cran()
  data(Sa)
  X <- Sa[[1]]
  X[500:519, 5] <- NA
  expect_warning(
    preprocess(X, dim1 = seq(10, 18.66, by = .01),
               dim2 = seq(200, 318, by = 2), baseline = "none",
               show_progress = FALSE),
    regexp = "^Missing values will be bridged by linear interpolation\\.")
  expect_error(
    expect_no_warning(
      preprocess(X, interpolate_rows = FALSE, interpolate_cols = FALSE,
                 show_progress = FALSE)))
})

test_that("warnings raised per-column survive the default parallel path", {
  # `asysm` warns on non-convergence from inside the workers that `preprocess`
  # dispatches to, and `parLapply` does not propagate warnings.
  skip_on_cran()
  data(Sa)
  for (cl in c(1, 2))
    expect_warning(
      preprocess(Sa[[1]], dim1 = seq(10, 18.6, by = 0.05),
                 dim2 = seq(200, 318, by = 6), maxit = 1,
                 show_progress = FALSE, cl = cl),
      regexp = "did not reach convergence", info = cl)
  p <- preprocess(Sa[1], dim1 = seq(10, 18.66, by = .01),
                  dim2 = seq(200, 318, by = 2), show_progress = FALSE)
  expect_null(attr(p[[1]], "baseline_warnings"))
})

test_that("warping drops the baseline attribute and keeps fixed-width times", {
  # The baseline is indexed by the time axis, so any transform that resamples
  # the signal invalidates it. Under `ptw` the shapes even still agree, which is
  # worse -- it would look valid while being silently misaligned.
  skip_on_cran()
  data(Sa)
  pr <- preprocess(Sa, dim1 = seq(10, 18.66, by = .01),
                   dim2 = seq(200, 318, by = 2), attach_baseline = TRUE,
                   show_progress = FALSE, cl = 1)
  has_bl <- function(l) any(vapply(l, function(x) !is.null(attr(x, "baseline", exact = TRUE)),
                                   logical(1)))
  ndec <- function(x) length(unique(nchar(sub("^.*\\.", "", rownames(x)))))
  for (alg in c("ptw", "lptw", "vpdtw")){
    w <- suppressWarnings(suppressMessages(
      correct_rt(pr, alg = alg, lambdas = "210", show_progress = FALSE,
                 cl = 1)))
    expect_false(has_bl(w), info = alg)
    expect_true(all(vapply(w, ndec, integer(1)) == 1L), info = alg)
  }
  for (m in c("mean", "nearest", "preceding"))
    expect_false(has_bl(subtract_blanks(pr, blank_id = 1, method = m)),
                 info = m)
})

# ── Dimname formatting ────────────────────────────────────────────────────────

test_that("format_axis gives every value the same number of decimals", {
  # `as.character` formats each element independently and drops trailing zeros
  expect_equal(as.character(c(10.09, 10.10, 10.11)), c("10.09", "10.1", "10.11"))
  expect_equal(format_axis(c(10.09, 10.10, 10.11)),
               c("10.09", "10.10", "10.11"))
  # integers stay integers
  expect_equal(format_axis(seq(200, 206, by = 2)), c("200", "202", "204", "206"))
  # enough decimals are used to keep every value distinct
  expect_equal(format_axis(c(1, 1.5)), c("1.0", "1.5"))
  expect_equal(format_axis(c(0.001, 0.002)), c("0.001", "0.002"))
  # non-numeric input passes through
  expect_equal(format_axis(c("a", "b")), c("a", "b"))
  expect_equal(format_axis(c(10, Inf)), c("10", "Inf"))
  expect_equal(format_axis(c(10.1, NA)), c("10.1", NA))
  expect_equal(format_axis(numeric(0)), numeric(0))
  # and a real axis round-trips exactly, with no collisions
  for (ax in list(seq(10, 18.66, by = .01), seq(600, 1119.6, by = .6),
                  seq(200, 318, by = 2), c(0.05, 0.1, 0.15))){
    out <- format_axis(ax)
    expect_equal(as.numeric(out), ax)
    expect_false(anyDuplicated(out) > 0)
  }
})

test_that("format_axis stays fixed-width when no width resolves the axis", {
  expect_equal(format_axis(c(1, 2, 2 + 4 * .Machine$double.eps)),
               c("1.000000000000", "2.000000000000", "2.000000000000"))
})

test_that("preprocess emits fixed-width retention times", {
  skip_on_cran()
  data(Sa)
  ts <- seq(10, 18.66, by = .01)
  x <- preprocess(Sa[[1]], dim1 = ts, dim2 = seq(200, 318, by = 2),
                  show_progress = FALSE, cl = 1)
  decimals <- nchar(sub("^.*\\.", "", rownames(x)))
  expect_length(unique(decimals), 1L)
  expect_equal(as.numeric(rownames(x)), ts)
  expect_true("10.10" %in% rownames(x))
  # wavelengths are integral, so they stay unchanged
  expect_equal(colnames(x), as.character(seq(200, 318, by = 2)))
})

# ── asPLS ─────────────────────────────────────────────────────────────────────

test_that("asPLS's whit2 reformulation matches an explicit solve of its system", {
  # asPLS's system is `(W + stiffness * diag(alpha) %*% t(D) %*% D) z = W y`.
  # Dividing through by `stiffness * alpha` turns it into the system `whit2`
  # solves at `lambda = 1` with weights `w / (stiffness * alpha)`, which is what
  # `aspls_baseline` actually calls. Guard that identity against a dense solve,
  # since a regression here would be silent -- the result would still look like
  # a plausible baseline.
  skip_on_cran()
  set.seed(3)
  n <- 120
  y <- 10 + seq_len(n) / 12 + 30 * dnorm(seq_len(n), 60, 6) + rnorm(n, 0, 0.2)
  lambda <- 1e5
  k <- 0.5
  DtD <- crossprod(diff(diag(n), differences = 2))
  w <- rep(1, n)
  alpha <- rep(1, n)
  for (i in seq_len(100)){
    z <- solve(diag(w) + lambda * (diag(alpha) %*% DtD), w * y)
    d <- y - z
    dn <- d[d < 0]
    if (length(dn) < 2) break
    s <- sd(dn)
    if (!is.finite(s) || s == 0) break
    wt <- 1 / (1 + exp(k * (d - s) / s))
    ad <- abs(d)
    if (max(ad) > 0) alpha <- ad / max(ad)
    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < 1e-3) break
    w <- wt
  }
  expect_equal(aspls_baseline(y, stiffness = lambda), as.numeric(z),
               tolerance = 1e-4)
})

test_that("arPLS matches an explicit dense solve of its own system", {
  skip_on_cran()
  set.seed(3)
  n <- 150
  y <- 10 + seq_len(n) / 12 + 30 * dnorm(seq_len(n), 70, 6) + rnorm(n, 0, 0.2)
  lambda <- 1e5
  DtD <- crossprod(diff(diag(n), differences = 2))
  w <- rep(1, n)
  for (i in seq_len(50)){
    z <- solve(diag(w) + lambda * DtD, w * y)
    d <- y - z
    dn <- d[d < 0]
    if (length(dn) < 2) break
    m <- mean(dn)
    sdev <- sd(dn)
    if (!is.finite(sdev) || sdev == 0) break
    wt <- 1 / (1 + exp(2 * (d - (-m + 2 * sdev)) / sdev))
    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < 1e-6) break
    w <- wt
  }
  expect_equal(arpls_baseline(y, stiffness = lambda), as.numeric(z),
               tolerance = 1e-6)
})

test_that("asPLS's baseline degrades more slowly than arPLS's at extreme stiffness", {
  # asPLS's documented selling point, measured: the advantage is about 2x, and
  # only at the stiff end of the range. An earlier 10x figure came from a
  # sinusoidal baseline term that has since been removed as unrealistic.
  skip_on_cran()
  skip_on_ci()
  s <- sim_chrom()
  grid <- 10^(6:12)
  worst <- function(fn) max(vapply(grid,
    function(l) rmse(fn(s$y, stiffness = l), s$b), numeric(1)))
  expect_lt(worst(aspls_baseline), worst(arpls_baseline) / 1.5)
})

test_that("arPLS recovers small-peak areas better than asPLS or AsLS", {
  # The metric that propagates to results. Unlike baseline RMSE, this ordering
  # the estimator to recommend for trace-level peaks.
  skip_on_cran()
  skip_on_ci()
  dt <- 0.0067; pw <- 0.15; sd_pk <- pw / 2.355
  tt <- seq(10, 18.67, by = dt)
  amps <- c(300, 30, 5, 1.5); ctr <- c(11, 13, 15, 17)
  b <- 20 * exp(-(tt - 10) / 4) + 2
  set.seed(1)
  y <- b + rowSums(vapply(seq_along(amps),
    function(i) amps[i] * dnorm(tt, ctr[i], sd_pk), numeric(length(tt)))) +
    rnorm(length(tt), sd = 0.5)
  small_err <- function(bl) max(vapply(3:4, function(i){
    w <- which(abs(tt - ctr[i]) < 3 * sd_pk)
    abs(sum((y - bl)[w]) * dt - amps[i]) / amps[i] * 100
  }, numeric(1)))
  for (k in c(1, 2, 4)){
    lam <- (k * 1.3 * pw / dt)^4
    lam_aspls <- (k * 6.3 * pw / dt)^4
    e_arpls <- small_err(arpls_baseline(y, stiffness = lam))
    expect_lt(e_arpls, small_err(aspls_baseline(y, stiffness = lam_aspls)))
    expect_lt(e_arpls, small_err(ptw::asysm(y, lambda = lam)))
    expect_lt(e_arpls, 2)
  }
})

test_that("asPLS gets its own default stiffness, not arPLS's", {
  expect_gt(ASPLS_DEFAULT_STIFFNESS, DEFAULT_STIFFNESS)
  s <- sim_chrom()
  expect_false(isTRUE(all.equal(
    estimate_baseline(s$y, "aspls", dt = s$dt),
    estimate_baseline(s$y, "arpls", dt = s$dt))))
})

test_that("correct_baseline reproduces the baseline step of preprocess", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 12),
               show_progress = FALSE, cl = 1)
  for (b in c("als", "arpls", "aspls", "snip")){
    # snip has no default width, so it needs one either way
    extra <- if (b == "snip") list(m = 10) else list()
    full <- do.call(preprocess, c(list(Sa[1:2]), args, baseline = b, extra))
    none <- do.call(preprocess, c(list(Sa[1:2]), args, baseline = "none"))
    comp <- do.call(correct_baseline, c(list(none, baseline = b),
                                        list(show_progress = FALSE, cl = 1), extra))
    expect_equal(comp, full, tolerance = 0, info = b)
  }
})

test_that("correct_baseline preserves shape, names, dimnames and metadata", {
  data(Sa_pr)
  out <- correct_baseline(Sa_pr[1:2], show_progress = FALSE, cl = 1)
  expect_length(out, 2)
  expect_equal(names(out), names(Sa_pr[1:2]))
  expect_equal(class(out), class(Sa_pr[1:2]))
  expect_equal(dimnames(out[[1]]), dimnames(Sa_pr[[1]]))

  m <- Sa_pr[[1]]
  attr(m, "time_unit") <- "min"
  attr(m, "sample_name") <- "s1"
  attr(m, "parser") <- "chromConverter"
  attr(m, "data_format") <- "wide"
  cb <- correct_baseline(m, show_progress = FALSE, cl = 1)
  expect_true(is.matrix(cb))
  for (a in c("time_unit", "sample_name", "parser", "data_format"))
    expect_equal(attr(cb, a), attr(m, a), info = a)
})

test_that("correct_baseline attaches and floors on request", {
  data(Sa_pr)
  a <- correct_baseline(Sa_pr[[1]], attach_baseline = TRUE,
                        show_progress = FALSE, cl = 1)
  bl <- get_baseline(a)
  expect_equal(dim(bl), dim(a))
  expect_equal(dimnames(bl), dimnames(a))
  expect_equal(a + bl, Sa_pr[[1]], ignore_attr = TRUE)

  expect_error(get_baseline(correct_baseline(Sa_pr[[1]], show_progress = FALSE,
                                             cl = 1)),
               regexp = "No baseline found")
  expect_lt(min(correct_baseline(Sa_pr[[1]], show_progress = FALSE, cl = 1)), 0)
  expect_gte(min(correct_baseline(Sa_pr[[1]], zero_floor = TRUE,
                                  show_progress = FALSE, cl = 1)), 0)
})

test_that("correct_baseline surfaces per-column warnings like preprocess", {
  skip_on_cran()
  data(Sa_pr)
  x <- Sa_pr[1]
  for (cl in c(1, 2))
    expect_warning(correct_baseline(x, baseline = "als", maxit = 1,
                                    show_progress = FALSE, cl = cl),
                   regexp = "did not reach convergence", info = cl)
  expect_null(attr(suppressWarnings(
    correct_baseline(x, baseline = "als", maxit = 1, show_progress = FALSE,
                     cl = 1))[[1]], "baseline_warnings"))
  expect_null(attr(correct_baseline(x, show_progress = FALSE, cl = 1)[[1]],
                   "baseline_warnings"))
})

test_that("correct_baseline refuses an incomplete trace and offers no `none`", {
  data(Sa_pr)
  x <- Sa_pr[[1]]
  x[100:110, 3] <- NA
  expect_error(correct_baseline(x, show_progress = FALSE, cl = 1),
               regexp = "complete trace")
  expect_error(correct_baseline(Sa_pr[[1]], baseline = "none"),
               regexp = "should be one of")
  expect_error(correct_baseline("foo"),
               regexp = "x should be a matrix or a list of matrices")
})

test_that("crop_chroms returns an exact submatrix, dimnames included", {
  data(Sa)
  x <- Sa[[1]]
  ts <- as.numeric(rownames(x))
  lambdas <- as.numeric(colnames(x))

  expect_identical(crop_chroms(x, rt_range = c(12, 16)),
                   x[ts >= 12 & ts <= 16, , drop = FALSE])
  expect_identical(crop_chroms(x, lambda_range = c(220, 300)),
                   x[, lambdas >= 220 & lambdas <= 300, drop = FALSE])
  expect_identical(crop_chroms(x, rt_range = c(12, 16),
                               lambda_range = c(220, 300)),
                   x[ts >= 12 & ts <= 16, lambdas >= 220 & lambdas <= 300,
                     drop = FALSE])
  # bounds are inclusive, and an NA leaves that end open
  expect_identical(crop_chroms(x, rt_range = c(NA, 12)),
                   x[ts <= 12, , drop = FALSE])
  expect_identical(crop_chroms(x, rt_range = c(12, NA)),
                   x[ts >= 12, , drop = FALSE])
})

test_that("crop_chroms never resamples, unlike cropping via preprocess", {
  data(Sa)
  x <- Sa[[1]]
  ts <- as.numeric(rownames(x))
  keep <- ts[ts >= 12 & ts <= 16]
  cropped <- crop_chroms(x, rt_range = c(12, 16))
  interpolated <- preprocess(x, dim1 = keep, interpolate_cols = FALSE,
                             baseline = "none", smooth_spectra = FALSE,
                             show_progress = FALSE, cl = 1)
  # the values coincide, but `preprocess` reformats the retention times
  expect_equal(unname(interpolated), unname(cropped), tolerance = 0)
  expect_false(identical(rownames(interpolated), rownames(cropped)))
  expect_identical(rownames(cropped), rownames(x)[ts >= 12 & ts <= 16])
})

test_that("crop_chroms preserves list shape, names and metadata", {
  data(Sa)
  out <- crop_chroms(Sa[1:2], rt_range = c(12, 16))
  expect_length(out, 2)
  expect_equal(names(out), names(Sa[1:2]))
  expect_equal(class(out), class(Sa[1:2]))

  m <- Sa[[1]]
  attr(m, "time_unit") <- "min"
  attr(m, "sample_name") <- "s1"
  attr(m, "parser") <- "chromConverter"
  cm <- crop_chroms(m, rt_range = c(12, 16))
  for (a in c("time_unit", "sample_name", "parser"))
    expect_equal(attr(cm, a), attr(m, a), info = a)
})

test_that("crop_chroms crops an attached baseline alongside the data", {
  data(Sa_pr)
  a <- correct_baseline(Sa_pr[[1]], attach_baseline = TRUE,
                        show_progress = FALSE, cl = 1)
  cropped <- crop_chroms(a, rt_range = c(12, 16), lambda_range = c(220, 300))
  bl <- get_baseline(cropped)
  expect_equal(dim(bl), dim(cropped))
  expect_equal(dimnames(bl), dimnames(cropped))
  expect_identical(bl, get_baseline(a)[rownames(cropped), colnames(cropped)])
})

test_that("crop_chroms rejects unusable ranges", {
  data(Sa)
  x <- Sa[[1]]
  expect_error(crop_chroms(x), regexp = "at least one of")
  expect_error(crop_chroms(x, rt_range = c(16, 12)), regexp = "must be increasing")
  expect_error(crop_chroms(x, rt_range = 12), regexp = "length 2")
  expect_error(crop_chroms(x, rt_range = c("a", "b")), regexp = "length 2")
  expect_error(crop_chroms(x, rt_range = c(NA, NA)), regexp = "at least one bound")
  expect_error(crop_chroms(x, rt_range = c(50, 60)), regexp = "does not overlap")
  expect_error(crop_chroms(x, lambda_range = c(900, 1000)), regexp = "does not overlap")
  expect_error(crop_chroms("foo", rt_range = c(12, 16)),
               regexp = "x should be a matrix or a list of matrices")
})

test_that("preprocess smooths only when there are enough wavelengths", {
  data(Sa)
  ts <- seq(10, 18.66, by = .01)
  for (k in 1:4){
    m <- Sa[[1]][, seq_len(k), drop = FALSE]
    out <- preprocess(m, dim1 = ts, interpolate_cols = FALSE,
                      show_progress = FALSE, cl = 1)
    expect_equal(ncol(out), k, info = k)
  }
  # smoothing is skipped below four wavelengths, so it matches an explicit skip
  m3 <- Sa[[1]][, 1:3, drop = FALSE]
  expect_equal(
    preprocess(m3, dim1 = ts, interpolate_cols = FALSE, show_progress = FALSE,
               cl = 1),
    preprocess(m3, dim1 = ts, interpolate_cols = FALSE, smooth_spectra = FALSE,
               show_progress = FALSE, cl = 1),
    tolerance = 0)
})

test_that("estimator arguments are matched exactly, not by prefix", {
  data(Sa)
  y <- Sa[[1]][, 20]
  # `m` is a snip argument; it used to partial-match arpls's `maxit`
  expect_error(estimate_baseline(y, "arpls", m = 20), regexp = "does not accept")
  expect_error(estimate_baseline(y, "arpls", m = 20), regexp = "snip")
  expect_error(estimate_baseline(y, "arpls", to = 1e-6), regexp = "does not accept")
  expect_error(estimate_baseline(y, "als", k = 0.5), regexp = "aspls")
  # documented per-estimator arguments still work
  expect_length(estimate_baseline(y, "arpls", maxit = 20), length(y))
  expect_length(estimate_baseline(y, "als", p = .001), length(y))
  expect_length(estimate_baseline(y, "snip", m = 20), length(y))
  expect_error(preprocess(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
                          dim2 = seq(200, 318, by = 12), baseline = "arpls",
                          m = 20, show_progress = FALSE, cl = 1),
               regexp = "does not accept")
})

test_that("crop_chroms warns when a range reaches past the data", {
  data(Sa)
  x <- Sa[[1]]
  expect_silent(crop_chroms(x, rt_range = c(12, 16)))
  expect_silent(crop_chroms(x, rt_range = c(12, NA)))
  expect_silent(crop_chroms(x, lambda_range = c(NA, 300)))
  expect_warning(crop_chroms(x, rt_range = c(5, 12)),
                 regexp = "extends beyond the available data")
  expect_warning(crop_chroms(x, lambda_range = c(300, 400)),
                 regexp = "extends beyond the available data")
  # the warning reports what was actually kept
  expect_warning(crop_chroms(x, rt_range = c(5, 12)), regexp = "kept")
  # named chromatograms are identified
  expect_warning(crop_chroms(Sa[1:2], rt_range = c(5, 25)), regexp = "119")
})

test_that("crop_chroms warns when offset time axes give unequal lengths", {
  a <- matrix(1, 11, 4,
              dimnames = list(seq(0, 1, by = .1), c(200, 210, 220, 230)))
  b <- matrix(1, 11, 4,
              dimnames = list(seq(0.05, 1.05, by = .1), c(200, 210, 220, 230)))
  # rectangular in, ragged out: the requested range lies inside both spans,
  # but the samples do not sample it identically
  expect_warning(out <- crop_chroms(list(a = a, b = b), rt_range = c(0.1, 0.5)),
                 regexp = "unequal dimensions")
  expect_equal(unname(sapply(out, nrow)), c(5L, 4L))
  # a shared axis crops to equal lengths without complaint
  expect_silent(crop_chroms(list(a = a, a2 = a), rt_range = c(0.1, 0.5)))
})

# ── interpolation grids ───────────────────────────────────────────────────────

test_that("c() written where seq() was meant is named as such", {
  skip_on_cran()
  data(Sa)
  lam <- seq(200, 318, by = 2)
  # `c(from, to, by = x)` silently returns a three-element vector
  expect_error(preprocess(Sa[1], dim1 = c(0.2, 20, by = 0.005), dim2 = lam,
                          cl = 1, show_progress = FALSE),
               regexp = "Did you mean `seq\\(0.2, 20, by = 0.005\\)`")
  # caught even when the values are inside the data range
  expect_error(preprocess(Sa[1], dim1 = c(10, 18, by = 0.005), dim2 = lam,
                          cl = 1, show_progress = FALSE),
               regexp = "3-element vector")
  expect_error(preprocess(Sa[1], dim1 = c(10, 18, length.out = 100),
                          dim2 = lam, cl = 1, show_progress = FALSE),
               regexp = "length.out = 100")
  # and on the spectral axis
  expect_error(preprocess(Sa[1], dim1 = seq(10, 18, by = 0.01),
                          dim2 = c(200, 318, by = 2), cl = 1,
                          show_progress = FALSE),
               regexp = "`dim2` was given")
})

test_that("an interpolation grid must be strictly increasing", {
  skip_on_cran()
  data(Sa)
  lam <- seq(200, 318, by = 2)
  expect_error(preprocess(Sa[1], dim1 = c(10, 18, 12), dim2 = lam, cl = 1,
                          show_progress = FALSE),
               regexp = "strictly increasing")
  expect_error(preprocess(Sa[1], dim1 = c(10, 11, 11, 12), dim2 = lam, cl = 1,
                          show_progress = FALSE),
               regexp = "element 2 \\(11\\) does not precede element 3 \\(11\\)")
  expect_error(preprocess(Sa[1], dim1 = c(10, NA, 12), dim2 = lam, cl = 1,
                          show_progress = FALSE),
               regexp = "finite values")
  # a scalar `dim2` stays legal when the spectral axis is not interpolated
  expect_silent(preprocess(Sa[[1]][, 1, drop = FALSE],
                           dim1 = seq(10, 18, by = 0.01), dim2 = 210,
                           interpolate_cols = FALSE, baseline = "none",
                           cl = 1, show_progress = FALSE))
})
