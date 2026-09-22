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

# ── min_amplitude filtering ──────────────────────────────────────────────────────

test_that("find_peaks min_amplitude removes peaks below the amplitude threshold", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(101, 201), widths = c(10, 10),
                    heights = c(1000, 10))
  p_all  <- find_peaks(y, smooth_type = "none", min_amplitude = 0)
  p_tall <- find_peaks(y, smooth_type = "none", min_amplitude = 50)
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

# ── Apex alignment (C0) ───────────────────────────────────────────────────────

test_that("find_peaks reports the true apex, not the timepoint before it", {
  # `d[i] = y[i+1] - y[i]`, so a downward zero-crossing between `d[i]` and
  # `d[i+1]` puts the apex at `y[i+1]`. Reporting `i` made every `pos` one
  # sampling interval early.
  expect_equal(find_peaks(c(0, 1, 2, 5, 20, 100, 20, 5, 2, 1, 0),
                          smooth_type = "none")$pos, 6)
})

test_that("find_peaks pos equals which.max for every smoothing method", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  for (st in c("none", "gaussian", "box", "mva", "tmva")){
    p <- find_peaks(y, smooth_type = st)
    expect_equal(p$pos[which.max(p$prominence)], which.max(y),
                 info = st)
  }
})

test_that("min_amplitude is tested at the apex, not on the rising flank", {
  # `preprocess` used to clip negatives to exactly zero, and `y > min_amplitude`
  # is false at zero, so a peak whose preceding timepoint had been clipped was
  # dropped even at the inert `min_amplitude = 0` default.
  mk <- function(preceding){
    y <- rep(0, 50); y[20] <- preceding; y[21] <- 8; y[22] <- 3; y
  }
  expect_equal(nrow(find_peaks(mk(3), smooth_type = "none")), 1L)
  expect_equal(nrow(find_peaks(mk(0), smooth_type = "none")), 1L)
})

test_that("the deprecated `slope_thresh` and `amp_thresh` still work", {
  y <- make_gaussian_signal(n = 601, centers = c(150, 301, 450),
                            widths = c(4, 20, 40), heights = c(20, 1000, 500))
  expect_warning(old <- find_peaks(y, smooth_type = "none", amp_thresh = 50),
                 regexp = "`amp_thresh` argument is deprecated")
  expect_equal(old, find_peaks(y, smooth_type = "none", min_amplitude = 50))
  expect_warning(old <- find_peaks(y, smooth_type = "none", slope_thresh = 1),
                 regexp = "`slope_thresh` argument is deprecated")
  expect_equal(old, find_peaks(y, smooth_type = "none", min_curvature = 1))
})

test_that("min_curvature filters on curvature, excluding broad peaks first", {
  y <- make_gaussian_signal(n = 601, centers = c(150, 450),
                            widths = c(4, 40), heights = c(500, 500))
  expect_equal(nrow(find_peaks(y, smooth_type = "none", min_curvature = 0)), 2L)
  broad_gone <- find_peaks(y, smooth_type = "none", min_curvature = 1)
  expect_equal(broad_gone$pos, 150L)
})

test_that("savgol rejects a window too narrow to smooth", {
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  expect_error(find_peaks(y, smooth_type = "savgol"), regexp = "too narrow")
  expect_error(find_peaks(y, smooth_type = "savgol", smooth_window = 5),
               regexp = "returns the signal unchanged")
  p <- find_peaks(y, smooth_type = "savgol", smooth_window = 21)
  expect_equal(p$pos[which.max(p$prominence)], 301L)
  expect_equal(nrow(find_peaks(y, smooth_type = "savgol", smooth_window = 20)),
               nrow(find_peaks(y, smooth_type = "savgol", smooth_window = 21)))
})

test_that("even-width moving averages do not shift the reported apex", {
  # An even-width running mean is not centered on a sample: it creates a
  # plateau at the extremum, giving two consecutive sign changes.
  skip_on_cran()
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  for (st in c("mva", "tmva")){
    p <- find_peaks(y, smooth_type = st, smooth_window = 4)
    expect_equal(p$pos[which.max(p$prominence)], 301L, info = st)
  }
})

# ── Noise-referenced prominence filter (C) ────────────────────────────────────

test_that("min_prominence_snr = 0 leaves find_peaks output unchanged", {
  skip_on_cran()
  data(Sa_pr)
  for (nm in names(Sa_pr)) for (lam in c("210", "260", "310")){
    y <- Sa_pr[[nm]][, lam]
    expect_identical(find_peaks(y), find_peaks(y, min_prominence_snr = 0))
  }
})

test_that("find_peaks returns a prominence column", {
  data(Sa_pr)
  p <- find_peaks(Sa_pr[[1]][, "210"])
  expect_true("prominence" %in% names(p))
  expect_equal(nrow(p), length(p$prominence))
  expect_true(all(p$prominence >= 0))
})

test_that("raising min_prominence_snr drops the lowest-prominence apices first", {
  skip_on_cran()
  data(Sa_pr)
  y <- Sa_pr[[1]][, "210"]
  levels <- c(0, 1, 2, 4, 8, 16)
  sets <- lapply(levels, function(k) find_peaks(y, min_prominence_snr = k))
  # monotonically fewer peaks, and each set nested in the looser one
  expect_true(all(diff(vapply(sets, nrow, integer(1))) <= 0))
  for (i in seq_len(length(sets) - 1))
    expect_true(all(sets[[i + 1]]$pos %in% sets[[i]]$pos))
  # and the retained minimum prominence rises
  mins <- vapply(sets, function(p) min(p$prominence), numeric(1))
  expect_true(all(diff(mins) > 0))
})

test_that("a noise-referenced threshold survives a dominant solvent front", {
  # A single large peak turns any threshold expressed as a fraction of the
  # largest feature into a hundreds-of-sigma filter. Here 1% of full scale is
  # 200 units against a noise sd of ~1.
  skip_on_cran()
  set.seed(42)
  tt <- seq(0, 30, by = 0.01)
  g <- function(h, m, s) h * exp(-(tt - m)^2 / (2 * s^2))
  y <- g(20000, 2, 0.05) + g(30, 12, 0.15) + g(8, 20, 0.15) +
    rnorm(length(tt))
  sw <- 0.01  # ~30 timepoints, comparable to the analyte peak width
  near <- function(p, target)
    any(abs(p$pos - which.min(abs(tt - target))) <= 5)

  all_pks <- find_peaks(y, smooth_window = sw)
  rel <- all_pks[all_pks$prominence >
                   0.01 * max(all_pks$prominence, na.rm = TRUE), ]
  snr <- find_peaks(y, smooth_window = sw, min_prominence_snr = 4)

  # a max-referenced threshold loses both analytes ...
  expect_false(near(rel, 12))
  expect_false(near(rel, 20))
  # ... while a noise-referenced one keeps them
  expect_true(near(snr, 12))
  expect_true(near(snr, 20))
})

test_that("min_prominence_snr requires bounds = TRUE", {
  expect_error(find_peaks(c(0, 1, 5, 1, 0), min_prominence_snr = 4,
                          bounds = FALSE),
               regexp = "requires `bounds = TRUE`")
  expect_type(find_peaks(c(0, 1, 5, 1, 0), bounds = FALSE), "integer")
})
