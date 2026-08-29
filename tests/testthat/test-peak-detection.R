# Tests for find_peaks correctness using synthetic signals with known peak positions.

make_gaussian_signal <- function(n, centers, widths, heights) {
  x <- seq_len(n)
  y <- numeric(n)
  for (i in seq_along(centers))
    y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * widths[i]^2))
  y
}

# ── Single peak ───────────────────────────────────────────────────────────────

test_that("find_peaks detects a single Gaussian peak at the correct position", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = 151, widths = 10, heights = 1000)
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 1L)
  expect_equal(p$pos, 151, tolerance = 1)
  expect_true(p$lower < p$pos)
  expect_true(p$upper > p$pos)
})

# ── Two peaks ─────────────────────────────────────────────────────────────────

test_that("find_peaks detects two well-separated Gaussian peaks in the right order", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(101, 201), widths = c(10, 10),
                    heights = c(1000, 800))
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 2L)
  expect_true(p$pos[1] < p$pos[2])
  expect_equal(p$pos[1], 101, tolerance = 1)
  expect_equal(p$pos[2], 201, tolerance = 1)
})

# ── amp_thresh filtering ──────────────────────────────────────────────────────

test_that("find_peaks amp_thresh removes peaks below the amplitude threshold", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(101, 201), widths = c(10, 10),
                    heights = c(1000, 10))
  p_all  <- find_peaks(y, smooth_type = "none", amp_thresh = 0)
  p_tall <- find_peaks(y, smooth_type = "none", amp_thresh = 50)
  expect_equal(nrow(p_all),  2L)
  expect_equal(nrow(p_tall), 1L)
  expect_equal(p_tall$pos, 101, tolerance = 1)
})

# ── bounds = FALSE ────────────────────────────────────────────────────────────

test_that("find_peaks with bounds = FALSE returns a numeric vector", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = 151, widths = 10, heights = 1000)
  p <- find_peaks(y, smooth_type = "none", bounds = FALSE)
  expect_true(is.numeric(p))
  expect_equal(length(p), 1L)
  expect_equal(p, 151, tolerance = 1)
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("find_peaks errors when y is not a vector", {
  expect_error(find_peaks(matrix(1:9, 3)), regexp = "Please provide a vector")
})
