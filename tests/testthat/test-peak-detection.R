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

# ── Restricting peak fitting to target retention times ────────────────────────

data(Sa_pr)
lam_rt   <- "210"
pks_full <- get_peaks(Sa_pr, lambdas = lam_rt, show_progress = FALSE)
# Two apex retention times known to be present in sample 1.
rts_2 <- pks_full[[1]][[lam_rt]]$rt[c(3, 8)]











# ── `bounds = FALSE` is not usable for integration ────────────────────────────

test_that("`bounds = FALSE` is rejected by fit_peaks and get_peaks", {
  data(Sa_pr)
  # `find_peaks(bounds = FALSE)` returns a bare vector of apex positions, which
  # cannot be integrated because every peak model needs a window. Previously
  # this reached `fit_peaks` through `...` and failed with the opaque message
  # "argument is of length zero".
  expect_error(fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw", bounds = FALSE),
               "cannot be used with `fit_peaks`")
  expect_error(get_peaks(Sa_pr[1], lambdas = 210, fit = "raw", bounds = FALSE,
                         cl = 1, show_progress = FALSE),
               "cannot be used with `fit_peaks`")
})

test_that("`bounds = TRUE` remains a harmless no-op and find_peaks is unaffected", {
  skip_on_cran()
  data(Sa_pr)
  expect_equal(fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw", bounds = TRUE),
               fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw"))
  # `find_peaks` still honors it when called directly, as `plot_spectrum` does
  expect_null(dim(find_peaks(Sa_pr[[1]][, "210"], bounds = FALSE)))
  expect_s3_class(find_peaks(Sa_pr[[1]][, "210"], bounds = TRUE), "data.frame")
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


# ── exponentially modified Gaussians must not overflow ───────────────────────

test_that("log_erfc agrees with log(erfc()) and survives the far upper tail", {
  x <- seq(-3, 5, by = 0.25)
  expect_equal(log_erfc(x), log(erfc(x)), tolerance = 1e-12)
  expect_true(all(is.finite(log_erfc(seq(10, 200, by = 10)))))
  expect_true(all(is.infinite(log(erfc(seq(30, 200, by = 10))))))
})

test_that("bemg stays finite where the naive exp(A) * erfc(z) form overflows", {
  x <- seq(2.79, 3.185, by = 0.005)
  naive <- function(center, width, height, a, b){
    tm <- x - center; as2 <- a * width^2; bs2 <- b * width^2
    s <- sqrt(2 * width^2)
    height * 0.5 * (exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / s) +
                    exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / s))
  }
  expect_true(any(!is.finite(naive(2.845, 40, 960, 1, 1))))
  expect_true(all(is.finite(bemg(x, 2.845, 40, 960, 1, 1))))
  expect_true(all(is.finite(bemg(x, 2.9, 25, 500, 2, 2))))
  expect_true(all(bemg(x, 2.845, 40, 960, 1, 1) >= 0))
})

test_that("bemg is unchanged where the naive form was already finite", {
  x <- seq(-1, 1, length.out = 101)
  naive <- function(center, width, height, a, b){
    tm <- x - center; as2 <- a * width^2; bs2 <- b * width^2
    s <- sqrt(2 * width^2)
    height * 0.5 * (exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / s) +
                    exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / s))
  }
  for (p in list(c(0.5, 0.05, 500, 2, 3), c(-0.3, 0.1, 50, 0.5, 0.8),
                 c(0, 0.2, 1000, 0.3, 4))){
    ref <- naive(p[1], p[2], p[3], p[4], p[5])
    skip_if(any(!is.finite(ref)))
    expect_equal(bemg(x, p[1], p[2], p[3], p[4], p[5]), ref, tolerance = 1e-10)
  }
})

test_that("a peak whose model yields a non-finite intensity is dropped, with a warning", {
  d <- data.frame(rt = c(10, 20, 30), sd = c(1, 2, 3),
                  height = c(100, NaN, 300), area = c(5, 6, NA))
  expect_warning(out <- drop_unfittable_peaks(d), "Discarded 2 peaks")
  expect_equal(nrow(out), 1)
  expect_equal(out$rt, 10)
  expect_silent(drop_unfittable_peaks(d[1, ]))
})

test_that("an out-of-range index becomes NA instead of aborting the run", {
  cl <- list(a = matrix(1, nrow = 5, ncol = 1,
                        dimnames = list(as.character(seq(0.1, 0.5, by = 0.1)), "220")))
  x <- data.frame(rt = 2, start = 1, end = 5, center = -1726,
                  sd = 1, FWHM = 2, area = 3)
  out <- convert_indices_to_times(x, timepoints = get_times(cl),
                                  tdiff = get_time_resolution(cl), tfac = 1)
  expect_true(is.na(out$center))
  expect_equal(out$rt, 0.2)
})

# ── fitted width is bounded by the integration window ───────────────────────

test_that("fit_window_bounds scales with the window and allows a margin", {
  b <- fit_window_bounds(seq(0, 10, by = 1), mult = 2)
  expect_equal(b$span, 10)
  expect_equal(b$width_max, 2 * 10 / FWHM_CONST)
  expect_equal(b$center_min, -5)
  expect_equal(b$center_max, 15)
  expect_equal(fit_window_bounds(c(4, 4))$span, 1)
  expect_equal(fit_window_bounds(seq(0, 10, by = 1), mult = 1)$center_min, 0)
})

test_that("clamp_start moves a start value inside its bounds", {
  expect_equal(clamp_start(40, 0, 4.2), 4.2)
  expect_equal(clamp_start(-3, 0, 4.2), 0)
  expect_equal(clamp_start(2, 0, 4.2), 2)
  expect_equal(clamp_start(NaN, 0, 4), 2)
})

test_that("no peak model fits wider than twice its own window", {
  set.seed(42)
  x <- seq(0, 2, length.out = 401)
  # a broad feature with a narrow window carved out of its flank: the case that
  # previously produced a fit many times wider than the data supporting it
  y <- 500 * exp(-(x - 1)^2 / (2 * 0.4^2)) + rnorm(length(x), sd = 1)
  for (w in list(c(1.30, 1.36), c(0.50, 0.58), c(1.00, 1.30))){
    i <- which(x >= w[1] & x <= w[2])
    span <- x[max(i)] - x[min(i)]
    for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
      m <- suppressWarnings(get(f)(x[i], y[i], baseline = "sloped"))
      expect_lte(FWHM_CONST * abs(unname(m$width)), FIT_WIDTH_MAX_MULT * span * 1.001)
      expect_gte(unname(m$width), 0)
      expect_true(all(is.finite(m$y)))
    }
  }
})


test_that("window_limited marks fits that came to rest on a window limit", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_false(m$window_limited, info = f)
  }

  set.seed(42)
  xb <- seq(0, 2, length.out = 401)
  broad <- 500 * exp(-(xb - 1)^2 / (2 * 0.4^2)) + rnorm(length(xb), sd = 1)
  n_pinned <- 0L
  for (w in list(c(1.30, 1.36), c(0.50, 0.58), c(1.00, 1.30), c(1.35, 1.45))){
    i <- which(xb >= w[1] & xb <= w[2])
    width_max <- FIT_WIDTH_MAX_MULT * (xb[max(i)] - xb[min(i)]) / FWHM_CONST
    for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
      m <- suppressWarnings(get(f)(xb[i], broad[i], baseline = "sloped"))
      info <- paste(f, "window", w[1], w[2])
      if (abs(unname(m$width)[1] - width_max) < 1e-6 * width_max){
        n_pinned <- n_pinned + 1L
        expect_true(m$window_limited, info = info)
      }
    }
  }
  expect_gt(n_pinned, 0)
})

test_that("converged reports the optimizer's own convergence status", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_true(m$converged, info = f)
  }

  flat <- rep(0, 20)
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(seq_along(flat), flat, baseline = "none"))
    expect_false(m$converged, info = f)
  }
})

test_that("the bound does not disturb a peak that fits comfortably inside it", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_equal(unname(m$width), 0.02, tolerance = 0.02)
    expect_equal(unname(m$center), 0.5, tolerance = 0.01)
    expect_lt(FWHM_CONST * abs(unname(m$width)), FIT_WIDTH_MAX_MULT * 1)
  }
})

# ── max_fwhm is in time units and agrees with filter_peaks ──────────────────

.mk_chrom <- function(dt){
  ts <- seq(0, 3, by = dt)
  y <- 1000 * exp(-(ts - 1)^2/(2*0.02^2)) + 600 * exp(-(ts - 2)^2/(2*0.06^2))
  list(s1 = matrix(y, ncol = 1, dimnames = list(as.character(ts), "220")))
}

test_that("get_peaks(max_fwhm) and filter_peaks(max_fwhm) select the same peaks", {
  skip_on_cran()
  for (dt in c(0.005, 0.01)){
    cl <- .mk_chrom(dt)
    all <- get_peaks(cl, lambdas = "220", max_fwhm = Inf, estimate_purity = FALSE,
                     show_progress = FALSE)
    for (v in c(0.03, 0.1, 1)){
      a <- suppressWarnings(get_peaks(cl, lambdas = "220", max_fwhm = v,
                                      estimate_purity = FALSE,
                                      show_progress = FALSE))[[1]][[1]]
      b <- suppressWarnings(filter_peaks(all, max_fwhm = v))[[1]][[1]]
      expect_equal(nrow(a), nrow(b))
      expect_equal(a$rt, b$rt)
    }
  }
})

test_that("the max_fwhm default applies no width filtering", {
  skip_on_cran()
  dt <- 0.005
  cl <- .mk_chrom(dt)
  d <- get_peaks(cl, lambdas = "220", estimate_purity = FALSE, show_progress = FALSE)[[1]][[1]]
  e <- get_peaks(cl, lambdas = "220", max_fwhm = Inf,
                 estimate_purity = FALSE, show_progress = FALSE)[[1]][[1]]
  expect_equal(d$rt, e$rt)
  expect_equal(d$sd, e$sd)
})

test_that("sd_max is deprecated but keeps its timepoint semantics", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  expect_warning(old <- get_peaks(cl, lambdas = "220", sd_max = 4,
                                  estimate_purity = FALSE, show_progress = FALSE),
                 "deprecated")
  new <- suppressWarnings(get_peaks(cl, lambdas = "220", max_fwhm = FWHM_CONST * 4 * 0.005,
                   estimate_purity = FALSE, show_progress = FALSE))
  expect_equal(old[[1]][[1]]$rt, new[[1]][[1]]$rt)
  expect_error(get_peaks(cl, lambdas = "220", max_fwhm = 1, sd_max = 50),
               "not both")
})

# ── FWHM is a measured width, not a rescaled sigma ──────────────────────────

.dense_fwhm <- function(f, half_span){
  g <- seq(-half_span, half_span, length.out = 1000001)
  v <- f(g)
  apex <- which.max(v)
  half <- v[[apex]] / 2
  g[apex - 1 + min(which(v[apex:length(v)] < half))] -
    g[max(which(v[seq_len(apex)] < half))]
}

test_that("egh_fwhm is the half-height width of the egh shape", {
  skip_on_cran()
  for (tau in c(0, 0.1, 1, 5, 20, 50, -1, -5, -20)){
    shape <- function(t) egh(t, 0, 1, 1, tau)
    expect_equal(egh_fwhm(1, tau), .dense_fwhm(shape, 300 + 50 * abs(tau)),
                 tolerance = 1e-3, info = paste("tau =", tau))
  }
  # and it reduces to the gaussian identity when the peak is symmetric
  expect_equal(egh_fwhm(2.5, 0), FWHM_CONST * 2.5)
})

test_that("bemg_fwhm is the half-height width of the bemg shape", {
  skip_on_cran()
  for (ab in list(c(1, 1), c(0.3, 2), c(5, 0.2), c(0.05, 0.05))){
    m <- list(center = 0, width = 1, height = 1, a = ab[1], b = ab[2])
    shape <- function(t) bemg(t, 0, 1, 1, ab[1], ab[2])
    expect_equal(bemg_fwhm(m), .dense_fwhm(shape, 500 + 100 / min(ab)),
                 tolerance = 1e-3, info = paste(ab, collapse = ", "))
  }
})

test_that("bemg_fwhm measures a peak whose gaussian core is degenerate", {
  skip_on_cran()
  m <- list(center = 0, width = 0, height = 10, a = 38.9, b = 3.125)
  expect_equal(bemg_fwhm(m),
               .dense_fwhm(function(t) bemg(t, 0, 0, 10, 38.9, 3.125), 20),
               tolerance = 1e-3)
})

test_that("FWHM is comparable across fits, unlike sd", {
  skip_on_cran()
  # each model divides a peak between core and tail differently, so sigma is a
  # property of the parameterisation while the half-height width is a property
  # of the peak
  data(Sa_pr, envir = environment())
  fits <- c("gaussian", "egh", "bemg")
  pks <- lapply(fits, function(f)
    suppressWarnings(get_peaks(Sa_pr[1], lambdas = "210", fit = f,
                               estimate_purity = FALSE, show_progress = FALSE,
                               cl = 1))[[1]][[1]])
  names(pks) <- fits
  # the same peaks, so the comparison is peak by peak
  for (f in fits[-1]) expect_equal(pks[[f]]$start, pks[[1]]$start)
  spread <- function(col){
    m <- sapply(pks, `[[`, col)
    apply(m, 1, function(v) max(v) / min(v) - 1)
  }
  expect_gt(median(spread("sd")), 0.4)
  expect_lt(median(spread("FWHM")), 0.3)
})

test_that("a raw peak has no fitted sigma and a measured FWHM", {
  skip_on_cran()
  # two peaks 0.08 apart with sigma 0.03: the valley between them sits above
  # half height, so there is no crossing to measure unless a local baseline
  # brings the window edges down
  ts <- seq(0, 3, by = 0.005)
  y <- 1000 * exp(-(ts - 1.4)^2 / (2 * 0.03^2)) +
    900 * exp(-(ts - 1.48)^2 / (2 * 0.03^2))
  cl <- list(s1 = matrix(y, ncol = 1, dimnames = list(as.character(ts), "220")))
  raw <- function(b) suppressWarnings(
    get_peaks(cl, lambdas = "220", fit = "raw", baseline = b,
              estimate_purity = FALSE, show_progress = FALSE))[[1]][[1]]
  none <- raw("none")
  expect_true(all(is.na(none$sd)))
  expect_true(all(is.na(none$FWHM)))
  sloped <- raw("sloped")
  expect_true(all(is.na(sloped$sd)))
  expect_true(all(is.finite(sloped$FWHM)))
  # an unmeasurable width is not grounds for discarding the peak
  expect_equal(nrow(none), nrow(sloped))
  none_list <- suppressWarnings(
    get_peaks(cl, lambdas = "220", fit = "raw", baseline = "none",
              estimate_purity = FALSE, show_progress = FALSE))
  expect_equal(nrow(filter_peaks(none_list, max_fwhm = 1e-9)[[1]][[1]]),
               nrow(none))
})

test_that("one max_fwhm transfers across every fit", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  n <- vapply(c("gaussian", "egh", "bemg", "raw"), function(f)
    nrow(suppressWarnings(
      get_peaks(cl, lambdas = "220", fit = f, max_fwhm = 0.2,
                estimate_purity = FALSE, show_progress = FALSE))[[1]][[1]]),
    numeric(1))
  expect_true(all(n == n[["gaussian"]]))
})


test_that("max_fwhm is ignored when windows are supplied", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  k <- suppressMessages(get_peaks(cl, lambdas = "220",
                                  windows = matrix(c(0.9, 1.1), nrow = 1),
                                  fit = "raw", max_fwhm = 1e-9,
                                  estimate_purity = FALSE, show_progress = FALSE))
  expect_equal(nrow(k[[1]][[1]]), 1)
})

test_that("filter_peaks can filter on converged and window_limited", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa_pr, lambdas = 210, fit = "egh",
                                    baseline = "sloped", sd_max = NULL,
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  n <- function(p) sum(vapply(p, function(s) nrow(s[["210"]]), numeric(1)))
  for (col in c("converged", "window_limited")){
    kept <- do.call(filter_peaks, c(list(pks), setNames(list(TRUE), col)))
    drop <- do.call(filter_peaks, c(list(pks), setNames(list(FALSE), col)))
    expect_equal(n(kept) + n(drop), n(pks), info = col)
    expect_true(all(unlist(lapply(kept, function(s) s[["210"]][[col]])) == 1),
                info = col)
    expect_true(all(unlist(lapply(drop, function(s) s[["210"]][[col]])) == 0),
                info = col)
  }
  both <- filter_peaks(pks, converged = TRUE, window_limited = FALSE)
  expect_lte(n(both), n(filter_peaks(pks, converged = TRUE)))
  expect_s3_class(both, "peak_list")

  expect_error(filter_peaks(pks, converged = "yes"), "TRUE or FALSE")
  expect_error(filter_peaks(pks, window_limited = NA), "TRUE or FALSE")
  expect_error(filter_peaks(pks, converged = c(TRUE, FALSE)), "TRUE or FALSE")
})

test_that("filter_peaks rejects fit diagnostics for a raw peak list", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  raw <- suppressWarnings(get_peaks(Sa_pr[1], lambdas = 210, fit = "raw",
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  expect_error(filter_peaks(raw, converged = TRUE), "no `converged` column")
  expect_error(filter_peaks(raw, window_limited = FALSE),
               "no `window_limited` column")
})

test_that("bemg rejects a fit whose tail grows away from the apex", {
  skip_on_cran()
  # a rising flank with no apex of its own: the bi-exponential can fit it with
  # a negative rate, which is a wall rather than a peak
  set.seed(1)
  x <- seq(0, 2, length.out = 60)
  y <- 200 * exp(-(x - 2.4)^2 / (2 * 0.5^2)) + rnorm(60, sd = 0.5)
  m <- suppressWarnings(fit_bemg(x, y, baseline = "none"))
  expect_gte(unname(m$a), 0)
  expect_gte(unname(m$b), 0)
  expect_false(m$converged)

  # the rejection is what makes the width measurable
  expect_true(is.finite(bemg_fwhm(m)))
})

test_that("reject_growing_tail passes through a well behaved fit", {
  skip_on_cran()
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  m <- suppressWarnings(fit_bemg(x, y, baseline = "none"))
  expect_gte(unname(m$a), 0)
  expect_gte(unname(m$b), 0)
  expect_true(m$converged)
})

test_that("find_peaks keeps its column schema when nothing is detected", {
  # `sapply` over an empty vector yields `list()`, which `data.frame` drops
  p <- find_peaks(rep(0, 100), bounds = TRUE)
  expect_equal(nrow(p), 0L)
  expect_named(p, c("pos", "lower", "upper", "prominence"))
})

test_that("fit_peaks rejects a `pos` that is not a find_peaks table", {
  x <- Sa_pr[[1]]
  expect_error(fit_peaks(x, lambda = "210", pos = c(50, 100)),
               regexp = "must be a data.frame with columns")
})

test_that("fit_peaks refuses `lambda = \"max\"` with a clear message", {
  expect_error(fit_peaks(Sa_pr[[1]], lambda = "max"),
               regexp = "must be a specific wavelength")
})

test_that("get_lambda_idx checks the wavelength of a single-column chromatogram", {
  expect_error(get_lambda_idx("999", lambdas = 210, allow_max = FALSE),
               regexp = "could not be found")
  expect_equal(get_lambda_idx("210", lambdas = 210, allow_max = FALSE), 1)
  expect_equal(get_lambda_idx("999", lambdas = NA_real_, allow_max = FALSE), 1)
})

test_that("print.peak_list reports no width filter by default", {
  pks <- get_peaks(Sa_pr[1], lambdas = "210", show_progress = FALSE)
  expect_output(print(pks), "max_fwhm: none")
})

test_that("find_peaks keeps peak windows clear of missing values at the ends", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(40, 270), widths = c(8, 8),
                            heights = c(1000, 1000))
  y[1:20] <- NA
  y[281:301] <- NA
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 2L)
  expect_equal(p$lower[1], 21)
  expect_equal(p$upper[2], 280)
  x <- matrix(y, ncol = 1, dimnames = list(seq_along(y) / 10, "210"))
  expect_no_warning(pks <- fit_peaks(x, lambda = "210"))
  expect_false(anyNA(pks$area))
})

test_that("get_peaks reports widths and areas in `time_unit` whatever the axis unit", {
  skip_on_cran()
  data(Sa_pr)
  s <- convert_time_units(Sa_pr[1], time_unit = "s")
  first <- function(x, u) get_peaks(x, lambdas = "210", time_unit = u,
                                    show_progress = FALSE, cl = 1)[[1]][["210"]][1, ]
  ref <- first(Sa_pr[1], "min")
  for (u in c("min", "s")){
    p <- first(s, u)
    f <- time_factor(u)
    expect_equal(p$area, ref$area * f, tolerance = 1e-3, info = u)
    expect_equal(p$sd, ref$sd * f, tolerance = 1e-3, info = u)
  }
  expect_equal(attr(get_peaks(s, lambdas = "210", show_progress = FALSE, cl = 1),
                    "time_unit"), "s")
})

test_that("fitted retention times are not snapped to the time grid", {
  skip_on_cran()
  ts <- seq(10, 12, by = 0.02)
  y <- 100 * exp(-(ts - 11.013)^2 / (2 * 0.05^2))
  chrom <- list(s = matrix(y, ncol = 1,
                           dimnames = list(format(ts, nsmall = 2), "210")))
  for (f in c("gaussian", "egh", "bemg")){
    rt <- suppressWarnings(get_peaks(chrom, lambdas = 210, fit = f, cl = 1,
                                     show_progress = FALSE))[[1]][["210"]]$rt
    expect_equal(rt, 11.013, tolerance = 1e-5, info = f)
  }
})
