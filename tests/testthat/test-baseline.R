# Regression tests for the baseline estimators and their `preprocess` /
# `correct_baseline` plumbing. See also test-preprocess.R, which covers the
# SNIP width criterion itself.

sim_offset_chrom <- function(dt = 0.01, seed = 1){
  set.seed(seed)
  tt <- seq(0, 30, by = dt)
  b <- 50 + 30 * exp(tt / 12)
  pk <- 1000 * dnorm(tt, 10, 0.12) + 200 * dnorm(tt, 15, 0.12) +
    60 * dnorm(tt, 20, 0.12)
  list(t = tt, b = b, y = b + pk + rnorm(length(tt)))
}

test_that("sigma_clipped_sd is invariant to a DC offset", {
  set.seed(1)
  z <- rnorm(1000)
  expect_equal(sigma_clipped_sd(z + 10), sigma_clipped_sd(z))
  expect_equal(sigma_clipped_sd(z + 1e4), sigma_clipped_sd(z))
  # and it still measures the scale it is given
  expect_equal(sigma_clipped_sd(z * 5), sigma_clipped_sd(z) * 5)
})

test_that("estimate_noise survives a detector offset", {
  y <- sim_offset_chrom()$y
  expect_equal(estimate_noise(y + 100), estimate_noise(y))
  expect_gt(estimate_noise(y), 0)
})

test_that("the prominence SNR filter is not disabled by an offset", {
  # clipping about zero collapsed the noise estimate, which silently switched
  # `min_prominence_snr` off on any trace with a baseline offset
  y <- sim_offset_chrom()$y
  expect_equal(nrow(find_peaks(y + 100, bounds = TRUE, min_prominence_snr = 5)),
               nrow(find_peaks(y, bounds = TRUE, min_prominence_snr = 5)))
})

test_that("a stale baseline attribute does not survive re-correction", {
  skip_on_cran()
  data(Sa)
  x1 <- preprocess(Sa[1], dim1 = seq(10, 18.6, by = 0.01),
                   dim2 = seq(200, 318, by = 2), attach_baseline = TRUE,
                   cl = 1, show_progress = FALSE)
  expect_equal(dim(get_baseline(x1)[[1]]), dim(x1[[1]]))

  # re-preprocessing onto a different grid must not carry the old baseline
  x2 <- preprocess(x1, dim1 = seq(11, 18, by = 0.02),
                   dim2 = seq(200, 318, by = 4), attach_baseline = FALSE,
                   cl = 1, show_progress = FALSE)
  expect_null(attr(x2[[1]], "baseline", exact = TRUE))
  expect_error(get_baseline(x2[[1]]), regexp = "No baseline found")
  # `attr()` partial-matches, so a bare read would return `baseline_param`
  expect_type(attr(x2[[1]], "baseline"), "double")

  # and neither must `correct_baseline`
  x3 <- correct_baseline(x1, attach_baseline = FALSE, cl = 1,
                         show_progress = FALSE)
  expect_null(attr(x3[[1]], "baseline", exact = TRUE))
})

test_that("a trace too short to baseline is rejected, not silently mangled", {
  # `ptw::asysm` corrupts memory on a single point rather than erroring, and
  # the R-level `dimnames` failure that used to mask it was only a symptom
  m <- matrix(c(10, 20, 30), nrow = 1,
              dimnames = list("10.00", c("210", "220", "230")))
  expect_error(correct_baseline(list(a = m), attach_baseline = TRUE, cl = 1,
                                show_progress = FALSE),
               regexp = "at least 4 timepoints")
  expect_error(estimate_baseline(c(1, 2, 3), "als"),
               regexp = "at least 4 timepoints")
  # four points is the shortest trace that is accepted
  expect_length(estimate_baseline(c(1, 2, 3, 4), "als"), 4L)
})

test_that("subtract_baseline always returns a baseline shaped like the data", {
  # `apply(X, 2, ...)` drops to a bare vector for a single-row matrix, so the
  # return shape is normalized rather than assumed
  X <- matrix(as.numeric(1:20), ncol = 2,
              dimnames = list(seq(10, 10.9, by = 0.1), c("210", "220")))
  res <- subtract_baseline(X, baseline = "als", attach_baseline = TRUE)
  expect_equal(dim(res$baseline), dim(X))
  expect_equal(dim(res$x), dim(X))
  expect_null(attr(res$x, "baseline", exact = TRUE))
})

test_that("correct_baseline accepts a bare numeric vector", {
  sim <- sim_offset_chrom()
  y <- setNames(sim$y, sim$t)
  v <- correct_baseline(y, cl = 1, show_progress = FALSE,
                        attach_baseline = TRUE)
  m <- correct_baseline(as.matrix(y), cl = 1, show_progress = FALSE,
                        attach_baseline = TRUE)
  expect_null(dim(v))
  expect_equal(names(v), names(y))
  expect_equal(as.numeric(v), as.numeric(m))
  expect_equal(as.numeric(get_baseline(v)), as.numeric(get_baseline(m)))
  expect_null(dim(get_baseline(v)))
  expect_equal(names(get_baseline(v)), names(y))
  # the time axis is still read off the names
  expect_equal(attr(attr(v, "baseline_param"), "dt"), 0.01)
  expect_null(names(correct_baseline(unname(y), cl = 1,
                                     show_progress = FALSE)))
})
