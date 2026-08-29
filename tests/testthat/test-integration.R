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
  y <- chromatographR:::egh(x_idx, center = 251, width = 20, height = 1000, tau = 40)
  p <- make_chrom_from_y(y)
  pos <- data.frame(pos = p$center, lower = 1, upper = p$n)
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "egh",
                   estimate_purity = FALSE)
  expect_equal(pks$area,   p$true_area, tolerance = 0.01)
  expect_equal(pks$height, p$height,    tolerance = 0.01)
})

test_that("gaussian fit underestimates area of a strongly tailed EGH peak", {
  skip_on_cran()
  x_idx <- seq_len(601)
  y <- chromatographR:::egh(x_idx, center = 251, width = 20, height = 1000, tau = 60)
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
  y <- chromatographR:::bemg(x_idx, center = 101, width = 5, height = 1000,
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
  y <- chromatographR:::bemg(x_idx, center = 101, width = 5, height = 1000,
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
  # sd_max = NULL: raw sd = window width (500) would be filtered by default sd_max = 50
  pks <- fit_peaks(p$mat, lambda = "220", pos = pos, fit = "raw",
                   sd_max = NULL, estimate_purity = FALSE)
  expect_equal(pks$area, expected_area, tolerance = 1e-6)
})

# ── Baseline subtraction ──────────────────────────────────────────────────────

test_that("flat baseline fit subtracts the offset and recovers the Gaussian area", {
  skip_on_cran()
  p_floor <- make_gaussian_chrom(baseline = 50)
  p_clean <- make_gaussian_chrom(baseline = 0)
  pos <- data.frame(pos = p_floor$center, lower = 1, upper = p_floor$n)

  pks_baseline <- fit_peaks(p_floor$mat, lambda = "220", pos = pos, fit = "gaussian",
                            baseline = "flat", estimate_purity = FALSE)
  pks_none     <- fit_peaks(p_floor$mat, lambda = "220", pos = pos, fit = "gaussian",
                            baseline = "none", estimate_purity = FALSE)

  expect_equal(pks_baseline$area, p_clean$true_area, tolerance = 0.01)
  expect_lt(pks_baseline$area, pks_none$area)
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
