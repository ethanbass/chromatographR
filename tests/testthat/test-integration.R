# Tests for peak integration accuracy using synthetic peaks with known areas.

# Build a single-column chromatogram matrix containing one Gaussian peak.
# All positions are in index space (1-indexed integers).
make_gaussian_chrom <- function(n = 501, center = 251, width = 20, height = 1000,
                                baseline = 0) {
  x_idx <- seq_len(n)
  y <- baseline + height * exp(-(x_idx - center)^2 / (2 * width^2))
  mat <- matrix(y, ncol = 1,
                dimnames = list(as.character(x_idx), "220"))
  list(mat = mat, n = n, center = center, width = width, height = height,
       baseline = baseline,
       true_area = height * width * sqrt(2 * pi))  # analytical Gaussian area
}

# Build a single-column chromatogram matrix from an arbitrary signal vector.
# true_area is the trapezoidal integral (ground truth for discretely-sampled
# signals without a closed-form expression).
make_chrom_from_y <- function(y) {
  n <- length(y)
  mat <- matrix(y, ncol = 1,
                dimnames = list(as.character(seq_len(n)), "220"))
  list(mat = mat, n = n, center = which.max(y),
       height = max(y),
       true_area = sum((y[-length(y)] + y[-1]) / 2))
}

# ── Symmetric Gaussian peak ───────────────────────────────────────────────────

test_that("gaussian fit recovers area, height, center, and sd of a Gaussian peak", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "gaussian",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 1e-3)
  expect_equal(pks$height, p$height,    tolerance = 1e-3)
  expect_equal(pks$rt,     p$center,    tolerance = 1e-3)
  expect_equal(pks$sd,     p$width,     tolerance = 1e-3)
})

test_that("egh fit recovers area and height of a symmetric Gaussian peak", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "egh",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 1e-3)
})

test_that("bemg fit recovers area and height of a symmetric Gaussian peak", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "bemg",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 0.01)
})

# ── Asymmetric EGH peaks ─────────────────────────────────────────────────────

test_that("egh fit recovers area and height of a right-tailed EGH peak", {
  skip_on_cran()
  x_idx <- seq_len(601)
  y <- egh(x_idx, center = 251, width = 20, 
                            height = 1000, tau = 40)
  p <- make_chrom_from_y(y)
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "egh",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 0.01)
})

#is this test really useful?
test_that("gaussian fit underestimates area of a strongly tailed EGH peak", {
  skip_on_cran()
  x_idx <- seq_len(601)
  y <- egh(x_idx, center = 251, width = 20, height = 1000, tau = 60)
  p <- make_chrom_from_y(y)
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks_egh   <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "egh",
                         estimate_purity = FALSE)
  pks_gauss <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "gaussian",
                         estimate_purity = FALSE)
  expect_lt(pks_gauss$area, pks_egh$area)
})

# ── Asymmetric BEMG peaks ─────────────────────────────────────────────────────
# fit_bemg starts at a = b = 1.0. With width = 20, exp(a * width^2 / 2) = exp(200)
# overflows on the first NLS iteration. width = 5 keeps exp(a * 12.5) ≈ 3e5, safe.

test_that("bemg fit recovers area and height of a right-tailed BEMG peak", {
  skip_on_cran()
  x_idx <- seq_len(201)
  y <- bemg(x_idx, center = 101, width = 5, height = 1000,
                              a = 0.1, b = 0.5)
  p <- make_chrom_from_y(y)
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "bemg",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 0.01)
})

test_that("bemg fit recovers area and height of a left-tailed BEMG peak", {
  skip_on_cran()
  x_idx <- seq_len(201)
  y <- bemg(x_idx, center = 101, width = 5, height = 1000,
                              a = 0.5, b = 0.1)
  p <- make_chrom_from_y(y)
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "bemg",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 0.01)
})

# ── Raw integration ───────────────────────────────────────────────────────────

test_that("raw integration matches trapezoidal sum of the signal", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  y <- p$mat[, 1]
  expected_area <- sum((y[-length(y)] + y[-1]) / 2)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "raw",
                   max_fwhm = NULL, estimate_purity = FALSE)
  expect_equal(pks$area, expected_area, tolerance = 1e-6)
})

# ── Baseline subtraction ──────────────────────────────────────────────────────
# A symmetric Gaussian is a valid target for all three peak models (see the
# tests above), so every fit is held to the same area recovery here. The
# baseline branches of `fit_egh` and `fit_bemg` are exercised only through this
# loop; testing `fit = "gaussian"` alone leaves them entirely unrun.
BASELINE_FITS <- c("gaussian", "egh", "bemg")

test_that("flat baseline fit subtracts the offset and recovers the Gaussian area", {
  skip_on_cran()
  p_floor <- make_gaussian_chrom(baseline = 50)
  p_clean <- make_gaussian_chrom(baseline = 0)
  pos <- data.frame(pos = p_floor$center, lower = 1, upper = p_floor$n)

  for (f in BASELINE_FITS){
    pks_baseline <- fit_peaks(p_floor$mat, lambda = "220", pos = pos, fit = f,
                              baseline = "flat", estimate_purity = FALSE)
    pks_none     <- fit_peaks(p_floor$mat, lambda = "220", pos = pos, fit = f,
                              baseline = "none", estimate_purity = FALSE)

    expect_equal(pks_baseline$area, p_clean$true_area, tolerance = 0.01,
                 info = f)
    expect_lt(pks_baseline$area, pks_none$area)
    # the fitted floor recovers the offset that was added to the signal
    expect_equal(pks_baseline$floor, p_floor$baseline, tolerance = 0.01,
                 info = f)
  }
})

# ── Sloped baseline ───────────────────────────────────────────────────────────

test_that("sloped baseline fit subtracts a linear ramp and recovers the Gaussian area", {
  skip_on_cran()
  # Signal = Gaussian + linear ramp from 20 to 80 over the window
  n <- 501; center <- 251; width <- 20; height <- 1000
  x_idx <- seq_len(n)
  ramp <- seq(20, 80, length.out = n)
  y <- height * exp(-(x_idx - center)^2 / (2 * width^2)) + ramp
  mat <- matrix(y, ncol = 1, dimnames = list(as.character(x_idx), "220"))
  true_area <- height * width * sqrt(2 * pi)
  pos <- data.frame(pos = center, lower = 1, upper = n)

  for (f in BASELINE_FITS){
    pks_sloped <- fit_peaks(mat, lambda = "220", pos = pos, fit = f,
                            baseline = "sloped", estimate_purity = FALSE)
    pks_none   <- fit_peaks(mat, lambda = "220", pos = pos, fit = f,
                            baseline = "none",   estimate_purity = FALSE)

    expect_equal(pks_sloped$area, true_area, tolerance = 0.02, info = f)
    expect_lt(pks_sloped$area, pks_none$area)
    expect_gt(pks_sloped$floor_start, 0)
    expect_gt(pks_sloped$floor_end,   pks_sloped$floor_start)
    # and the recovered ramp matches the one that was added
    expect_equal(pks_sloped$floor_start, ramp[1], tolerance = 0.02, info = f)
    expect_equal(pks_sloped$floor_end,   ramp[n], tolerance = 0.02, info = f)
  }
})

# ── nlsLM failure fallbacks ───────────────────────────────────────────────────

test_that("egh and bemg return start values for an unfittable window", {
  skip_on_cran()
  # A three-point window has fewer observations than either model has
  # parameters, so `nlsLM` cannot move off its starting point and the reported
  # peak is the start-value estimate. Asserting the exact values guards against
  # a regression silently substituting plausible-looking wrong numbers.
  y <- c(1, 20, 1)
  p <- make_chrom_from_y(y)
  # the start values `fit_egh`/`fit_bemg` compute from the window
  start_center <- which.max(y)
  start_height <- max(y)
  start_width  <- sum(y > (start_height / 2)) / 2
  w <- rbind(c(1, length(y)))

  for (b in c("none", "flat", "sloped")){
    pks <- suppressWarnings(
      fit_peaks(p$mat, lambda = "220", windows = w, fit = "egh", baseline = b,
                max_fwhm = NULL, estimate_purity = FALSE))
    expect_equal(pks$rt,     start_center, info = b)
    expect_equal(pks$sd,     start_width,  info = b)
    expect_equal(pks$height, start_height, info = b)
    # tau starts at 0, so the fallback peak is a plain Gaussian
    expect_equal(pks$tau, 0, info = b)

    pks <- suppressWarnings(
      fit_peaks(p$mat, lambda = "220", windows = w, fit = "bemg", baseline = b,
                max_fwhm = NULL, estimate_purity = FALSE))
    expect_equal(pks$rt,     start_center, info = b)
    expect_equal(pks$sd,     start_width,  info = b)
    # `h` is the amplitude parameter, which is what starts at the window
    # maximum; `height` is read off the evaluated curve and so is lower
    expect_equal(pks$h, start_height, info = b)
    expect_lt(pks$height, pks$h)
    # `a` and `b` both start at 1
    expect_equal(pks$tau_left,  1, info = b)
    expect_equal(pks$tau_right, 1, info = b)
    # no optimizer ran on a three-point window, so nothing converged
    expect_equal(pks$converged, 0, info = b)
  }
  # and the floor terms are the start value of 0, not a fitted value
  pks_flat <- suppressWarnings(
    fit_peaks(p$mat, lambda = "220", windows = w, fit = "egh",
              baseline = "flat", max_fwhm = NULL, estimate_purity = FALSE))
  expect_equal(pks_flat$floor, 0)
  pks_sloped <- suppressWarnings(
    fit_peaks(p$mat, lambda = "220", windows = w, fit = "egh",
              baseline = "sloped", max_fwhm = NULL, estimate_purity = FALSE))
  expect_equal(pks_sloped$floor_start, 0)
  expect_equal(pks_sloped$floor_end, 0)
})

test_that("egh and bemg fall back to start values when nlsLM errors", {
  skip_on_cran()
  # `nls.control(warnOnly = TRUE)` turns most convergence failures into
  # warnings, so the `try-error` fallback is reached only when `nlsLM` raises a
  # hard error. An all-zero window does it: the gradient is singular at the
  # initial estimates. `fit_peaks` filters such a window out before fitting, so
  # this is tested against the fitters directly.
  y <- rep(0, 9)
  x <- seq_along(y)
  for (b in c("none", "flat", "sloped")){
    fit <- suppressWarnings(fit_egh(x, y, baseline = b))
    # the start values: the window max is 0, at the first index
    expect_equal(unname(fit$center), 1, info = b)
    expect_equal(unname(fit$height), 0, info = b)
    expect_equal(unname(fit$width),  0, info = b)
    expect_equal(unname(fit$tau),    0, info = b)
    # the floor terms are start values too, never a fitted value, and are
    # reported only for the baseline modes that fit them
    expect_equal(unname(fit$floor), if (b == "none") NULL else 0, info = b)
    expect_equal(unname(fit$slope), if (b == "sloped") 0 else NULL, info = b)
    # and the fitted curve is all zeros, so the residual is the input
    expect_equal(unname(fit$residual), y, info = b)

    fit <- suppressWarnings(fit_bemg(x, y, baseline = b))
    expect_equal(unname(fit$center), 1, info = b)
    expect_equal(unname(fit$height), 0, info = b)
    expect_equal(unname(fit$width),  0, info = b)
    # `a` and `b` start at 1
    expect_equal(unname(fit$a), 1, info = b)
    expect_equal(unname(fit$b), 1, info = b)
    expect_equal(unname(fit$floor), if (b == "none") NULL else 0, info = b)
    expect_equal(unname(fit$slope), if (b == "sloped") 0 else NULL, info = b)
  }
})

# ── Automatic peak detection ───────────────────────────────────────────────────────

test_that("get_peaks detects and integrates a synthetic Gaussian peak automatically", {
  skip_on_cran()
  # No pos supplied: find_peaks must locate the peak from the signal itself.
  p <- make_gaussian_chrom(n = 501, center = 251, width = 20, height = 1000)
  pks <- get_peaks(list(s = p$mat), lambdas = "220", fit = "gaussian",
                   smooth_type = "none", estimate_purity = FALSE,
                   show_progress = FALSE)[[1]][[1]]
  expect_equal(nrow(pks), 1L)
  expect_equal(as.numeric(pks$rt), p$center, tolerance = 0.5)
  expect_equal(pks$area, p$true_area, tolerance = 0.01)
})

# ── get_peaks end-to-end: time-unit scaling ───────────────────────────────────

test_that("get_peaks area scales with dt and matches the analytical value", {
  skip_on_cran()
  n <- 501; center_idx <- 251; width_idx <- 20; height <- 1000
  dt <- 0.01  # minutes per index step
  y  <- height * exp(-(seq_len(n) - center_idx)^2 / (2 * width_idx^2))

  ts_coarse <- seq(0, (n - 1) * dt,      by = dt)        # dt = 0.01 min
  ts_fine   <- seq(0, (n - 1) * dt / 10, by = dt / 10)   # dt = 0.001 min

  make_mat <- function(ts) matrix(y, ncol = 1, dimnames = list(as.character(ts), "220"))

  fit <- function(mat) get_peaks(list(s = mat), lambdas = "220", fit = "gaussian",
                                 smooth_type = "none", estimate_purity = FALSE,
                                 show_progress = FALSE)[[1]][[1]]$area

  area_coarse <- fit(make_mat(ts_coarse))
  area_fine   <- fit(make_mat(ts_fine))

  # area should scale proportionally to dt
  expect_equal(area_coarse / area_fine, 10, tolerance = 0.01)

  # absolute area should match height * sigma_time * sqrt(2*pi)
  expect_equal(area_coarse, height * (width_idx * dt) * sqrt(2 * pi), tolerance = 0.01)
})

# should there be tests that time units are translated correctly inside get_peaks?
# ── Manual integration windows ────────────────────────────────────────────────
# `make_gaussian_chrom`/`make_chrom_from_y` label rows 1..n, so retention times
# and index positions coincide, which keeps the expected areas easy to state.

test_that("`windows` integrates the supplied window and recovers the area", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pks <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, p$n)),
                   fit = "raw", estimate_purity = FALSE)
  expect_equal(nrow(pks), 1L)
  expect_equal(pks$start, 1)
  expect_equal(pks$end, p$n)
  expect_equal(pks$area, p$true_area, tolerance = 0.01)
})

test_that("`windows` ignores `max_fwhm`", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pks <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, p$n)),
                   fit = "raw", max_fwhm = 1e-9, estimate_purity = FALSE)
  expect_equal(nrow(pks), 1L)
  expect_true(is.na(pks$sd))
  expect_equal(pks$end - pks$start, p$n - 1)
})

test_that("windows sharing a boundary act as a drop line and are exactly additive", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  whole <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, p$n)),
                     fit = "raw", estimate_purity = FALSE)
  split <- fit_peaks(p$mat, lambda = "220",
                     windows = rbind(c(1, 200), c(200, 350), c(350, p$n)),
                     fit = "raw", estimate_purity = FALSE)
  expect_equal(nrow(split), 3L)
  expect_equal(sum(split$area), whole$area)
})

test_that("`windows` returns a peak for a window containing no peak", {
  skip_on_cran()
  # forced integration: this is what `rts` cannot do, since it filters detected apexes
  p <- make_gaussian_chrom()
  pks <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, 20)),
                   fit = "raw", estimate_purity = FALSE)
  expect_equal(nrow(pks), 1L)
  expect_equal(pks$area, 0, tolerance = 1e-6)
})

test_that("`rt` and `height` come from the apex within the window", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  # a window whose midpoint (251) is the apex, and one whose midpoint is not
  off <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, 400)),
                   fit = "raw", estimate_purity = FALSE)
  expect_equal(off$rt, p$center)
  expect_equal(off$height, p$height, tolerance = 1e-6)
  # midpoint of c(1, 400) is 200, where the signal is well below the maximum
  expect_gt(off$height, p$mat[200, 1])
})



test_that("`windows` is validated", {
  p <- make_gaussian_chrom()
  fp <- function(w) fit_peaks(p$mat, lambda = "220", windows = w, fit = "raw")
  expect_error(fp(rbind(c(300, 100))), "must precede its end")
  expect_error(fp(rbind(c(1, 2, 3))), "exactly two columns")
  expect_error(fp(c(1, 2, 3)), "two-column matrix")
  expect_error(fp(rbind(c("a", "b"))), "numeric retention times")
  expect_error(fp(rbind(c(1, p$n + 100))), "outside the bounds")
})

# ── Raw integration with a local baseline ─────────────────────────────────────

test_that("raw flat baseline sits at the lower of the two window boundaries", {
  skip_on_cran()
  # the window minimum (10) lies below both boundaries (50 and 60), so the two
  # candidate definitions of a flat baseline give different answers here
  y <- c(50, 30, 10, 200, 10, 30, 60)
  p <- make_chrom_from_y(y)
  w <- rbind(c(1, length(y)))
  none <- fit_peaks(p$mat, lambda = "220", windows = w, fit = "raw",
                    baseline = "none", estimate_purity = FALSE)
  flat <- fit_peaks(p$mat, lambda = "220", windows = w, fit = "raw",
                    baseline = "flat", estimate_purity = FALSE)
  expect_equal(flat$floor, min(y[1], y[length(y)]))
  expect_equal(flat$area, none$area - min(y[1], y[length(y)]) * (length(y) - 1))
  # guards against defining the flat baseline as min(y), which would be a no-op
  # on baseline-corrected data (`preprocess` clips negatives to zero)
  expect_lt(flat$area, none$area)
  expect_equal(flat$height, max(y) - min(y[1], y[length(y)]))
})

test_that("raw sloped baseline recovers a Gaussian sitting on a linear ramp", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  idx <- seq_len(p$n)
  ramp <- 10 + 0.5 * (idx - 1)
  p_ramp <- make_chrom_from_y(p$mat[, 1] + ramp)
  pks <- fit_peaks(p_ramp$mat, lambda = "220", windows = rbind(c(1, p$n)),
                   fit = "raw", baseline = "sloped", estimate_purity = FALSE)
  expect_equal(pks$floor_start, ramp[1], tolerance = 1e-6)
  expect_equal(pks$floor_end, ramp[p$n], tolerance = 1e-6)
  expect_equal(pks$area, p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height, tolerance = 0.01)
})

test_that("a baseline above the signal gives a negative area, with a warning", {
  skip_on_cran()
  # both boundaries are high and the interior is low, so a flat baseline drawn
  # at the boundaries sits above the trace
  y <- c(1, 1, 100, 5, 4, 3, 2, 1, 100, 1, 1)
  p <- make_chrom_from_y(y)
  expect_warning(
    pks <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(3, 9)),
                     fit = "raw", baseline = "flat", estimate_purity = FALSE),
    "negative area")
  expect_lt(pks$area, 0)
})

# ── Manual windows through get_peaks and get_peaktable ────────────────────────

test_that("get_peaks applies `windows` to every sample and feeds get_peaktable", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  q <- make_gaussian_chrom(height = 500)
  cl <- list(a = p$mat, b = q$mat)
  w <- rbind(c(1, 200), c(200, 350), c(350, p$n))
  pks <- get_peaks(cl, lambdas = "220", windows = w, fit = "raw", cl = 1,
                   show_progress = FALSE)
  expect_s3_class(pks, "peak_list")
  # every sample contributes a row for every window, regardless of detection
  expect_equal(vapply(pks, function(s) nrow(s[["220"]]), numeric(1)),
               c(a = 3, b = 3))
  expect_equal(pks$a[["220"]]$start, w[, 1])
  expect_equal(pks$a[["220"]]$end, w[, 2])

  pt <- get_peaktable(pks)
  expect_equal(ncol(pt$tab), nrow(w))
  expect_true(all(pt$tab > 0))
  expect_equal(unname(unlist(pt$pk_meta["start", ])), w[, 1])
  expect_equal(unname(unlist(pt$pk_meta["end", ])), w[, 2])
})

test_that("get_peaks points at `fit = \"raw\"` when `windows` is used with the default fit", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  cl <- list(a = p$mat)
  w <- rbind(c(1, p$n))
  expect_message(get_peaks(cl, lambdas = "220", windows = w, cl = 1,
                           show_progress = FALSE),
                 'fit = "raw"')
  expect_message(get_peaks(cl, lambdas = "220", windows = w, fit = "raw", cl = 1,
                           show_progress = FALSE),
                 regexp = NA)
})


test_that("`windows` works with a peak model as well as raw integration", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  pks <- fit_peaks(p$mat, lambda = "220", windows = rbind(c(1, p$n)),
                   fit = "gaussian", estimate_purity = FALSE)
  expect_equal(nrow(pks), 1L)
  expect_equal(pks$area, p$true_area, tolerance = 0.01)
  expect_equal(pks$rt, p$center, tolerance = 0.01)
})

test_that("get_peaktable groups by window rather than by apex retention time", {
  skip_on_cran()
  # The apex within a window drifts between samples, so clustering on retention
  # time fragments a single window across several columns even though the
  # windows are identical. Here the apexes are 200 apart in one sample and 260
  # apart in the other, well beyond the default `hmax`.
  a <- make_chrom_from_y(c(rep(0, 99), 1000, rep(0, 100), 500, rep(0, 100)))
  b <- make_chrom_from_y(c(rep(0, 39), 1000, rep(0, 120), 500, rep(0, 140)))
  cl <- list(a = a$mat, b = b$mat)
  w <- rbind(c(1, 150), c(150, 301))
  pks <- get_peaks(cl, lambdas = "220", windows = w, fit = "raw", cl = 1,
                   show_progress = FALSE)
  pt <- get_peaktable(pks)
  expect_equal(ncol(pt$tab), nrow(w))
  # no sample is zero-filled, and the boundaries are the windows themselves
  expect_true(all(pt$tab > 0))
  expect_equal(unname(unlist(pt$pk_meta["start", ])), w[, 1])
  expect_equal(unname(unlist(pt$pk_meta["end", ])), w[, 2])
})

test_that("get_peaktable falls back to retention time clustering without `windows`", {
  skip_on_cran()
  p <- make_gaussian_chrom()
  cl <- list(a = p$mat, b = p$mat)
  pks <- get_peaks(cl, lambdas = "220", fit = "raw", max_fwhm = NULL, cl = 1,
                   show_progress = FALSE)
  expect_null(attr(pks, "windows"))
  expect_s3_class(get_peaktable(pks), "peak_table")
})

test_that("get_peaktable still groups by window after retention time correction", {
  skip_on_cran()
  a <- make_chrom_from_y(c(rep(0, 99), 1000, rep(0, 100), 500, rep(0, 100)))
  b <- make_chrom_from_y(c(rep(0, 39), 1000, rep(0, 120), 500, rep(0, 140)))
  w <- rbind(c(1, 150), c(150, 301))
  pks <- get_peaks(list(a = a$mat, b = b$mat), lambdas = "220", windows = w,
                   fit = "raw", cl = 1, show_progress = FALSE)
  for (s in c("a", "b")){
    for (col in c("rt", "start", "end")){
      pks[[s]][["220"]][[paste0(col, ".cor")]] <- pks[[s]][["220"]][[col]] +
        if (s == "a") 0.5 else -0.5
    }
  }
  pt <- get_peaktable(pks)
  expect_true(pt$args$use.cor)
  expect_equal(ncol(pt$tab), nrow(w))
  expect_true(all(pt$tab > 0))
})
