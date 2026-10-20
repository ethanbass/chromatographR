### Peak purity ################################################################

gat <- get_agilent_threshold
gpv <- get_purity_values
gnv <- get_noise_variance
fnoise <- find_noise

wl_grid <- seq(200, 300, by = 10)
tt_grid <- 1:40
elution <- function(centre, width) exp(-(tt_grid - centre)^2 / (2 * width^2))
spectrum <- function(centre) dnorm(wl_grid, centre, 25) * 100

pure_peak <- function(){
  d <- outer(elution(20, 4), spectrum(250))
  dimnames(d) <- list(tt_grid, wl_grid)
  d
}

# a second component on the tail, with a clearly different spectrum
mixed_peak <- function(){
  d <- outer(elution(20, 4), spectrum(250)) +
       outer(elution(24, 3) * 0.35, spectrum(215))
  dimnames(d) <- list(tt_grid, wl_grid)
  d
}
peak_pos <- c(20, 12, 28)   # apex, start, end

### trim_peak ##################################################################

test_that("trim_peak returns window-relative indices above the cutoff", {
  # Indices are relative to the window `pos[2]:pos[3]`, not to `x`, because
  # `get_purity` uses them to subset purity values that are themselves computed
  # over the window. Changing this to absolute indices would silently misalign
  # the two.
  x <- c(0, 1, 2, 10, 2, 1, 0)          # apex at index 4, value 10
  expect_equal(trim_peak(x, c(4, 1, 7), cutoff = 0.05), c(2, 3, 4, 5, 6))
  expect_equal(trim_peak(x, c(4, 1, 7), cutoff = 0.25), 4)
  # window starting at 3 renumbers from 1
  expect_equal(trim_peak(x, c(4, 3, 7), cutoff = 0.05), c(1, 2, 3, 4))
  # nothing exceeds the apex itself
  expect_length(trim_peak(x, c(4, 1, 7), cutoff = 1), 0)
})

### find_noise / get_noise_variance ############################################

test_that("find_noise selects rows below the threshold times the global maximum", {
  d <- rbind(c(1, 1), c(50, 50), c(100, 100), c(50, 50), c(1, 1))
  dimnames(d) <- list(1:5, c("200", "210"))
  # row maxima are 1, 50, 100, 50, 1; global maximum is 100. `which()` keeps the
  # rownames, so compare unnamed.
  expect_equal(unname(fnoise(d, noise_threshold = 0.02)), c(1L, 5L))
  expect_equal(unname(fnoise(d, noise_threshold = 0.60)),
               c(1L, 2L, 4L, 5L))
  expect_length(fnoise(d, noise_threshold = 0.005), 0)
})

test_that("get_noise_variance averages the variance of the noise rows", {
  # `method = "threshold"` is the legacy selection rule, retained but no longer
  # the default.
  d <- rbind(c(1, 1), c(50, 50), c(100, 100), c(50, 50), c(1, 1))
  dimnames(d) <- list(1:5, c("200", "210"))
  # the selected rows are constant, so their variance is exactly zero
  expect_equal(gnv(d, noise_method = "threshold", noise_threshold = 0.02, lambdas = 1:2), 0)
  d[1, ] <- c(1, 4)                      # row maximum 4, variance 4.5
  expect_equal(gnv(d, noise_method = "threshold", noise_threshold = 0.05, lambdas = 1:2), (4.5 + 0) / 2)
})

test_that("get_noise_variance handles a single noise row", {
  # `x[noise_idx, ]` collapses to a vector when one row qualifies, which made
  # `apply` fail. `get_purity` swallows errors into NA, so purity silently
  # became NA for such chromatograms.
  d <- rbind(c(1, 1), c(50, 50), c(100, 100), c(50, 50), c(9, 9))
  dimnames(d) <- list(1:5, c("200", "210"))
  # only row 1 falls below 2% of the maximum
  expect_equal(unname(fnoise(d, noise_threshold = 0.02)), 1L)
  expect_equal(gnv(d, noise_method = "threshold", noise_threshold = 0.02, lambdas = 1:2), 0)
  expect_silent(gnv(d, noise_method = "threshold", noise_threshold = 0.02, lambdas = 1:2))
})

### get_agilent_threshold ######################################################

test_that("the Agilent threshold spans [0, 1] and falls as noise rises", {
  d <- pure_peak()
  # with no noise the threshold demands perfect similarity, up to the ceiling
  expect_true(all(gat(d, peak_pos, noise_variance = 0, background_correction = "none",
                      max_threshold = 1) == 1))
  expect_true(all(gat(d, peak_pos, noise_variance = 0, background_correction = "none") == 0.99))
  # the formula is squared after a max(0, .), so it cannot go negative
  expect_true(all(gat(d, peak_pos, noise_variance = 1e12, background_correction = "none") == 0))
  expect_true(all(gat(d, peak_pos, noise_variance = 1, background_correction = "none") >=
                  gat(d, peak_pos, noise_variance = 10, background_correction = "none")))
})

### get_purity_values ##########################################################

test_that("purity values are finite when the threshold reaches 1", {
  # A noise variance of zero drives the threshold to 1, so the denominator
  # `1 - threshold` vanishes. Timepoints matching the apex spectrum are pure
  # regardless, and previously produced NaN, which `get_purity` dropped and so
  # counted as impure.
  pv <- gpv(pure_peak(), peak_pos, noise_variance = 0,
            background_correction = "none", max_threshold = 1)
  expect_false(any(is.nan(pv)))
  expect_true(all(pv == 0))
})

test_that("infinite purity ratios count as impure", {
  expect_equal(get_purity(mixed_peak(), peak_pos, noise_variance = 0,
                          background_correction = "none", max_threshold = 1), 0)
  expect_equal(get_purity(mixed_peak(), peak_pos, noise_variance = 0,
                          background_correction = "none", max_threshold = 1,
                          metric = "max"), Inf)
})

### get_purity #################################################################

test_that("get_purity reports a rank-1 peak as pure and a mixed peak as less pure", {
  # purity is the proportion of in-peak timepoints with a ratio below 1, so 1
  # means every timepoint is pure
  expect_equal(get_purity(pure_peak(), peak_pos), 1)
  expect_lt(get_purity(mixed_peak(), peak_pos), 1)
})

test_that("get_purity returns NA when purity is not requested", {
  # `try` is the switch `fit_peaks` uses to skip purity entirely
  expect_true(is.na(get_purity(pure_peak(), peak_pos, try = FALSE)))
})

test_that("get_purity does not depend on the order of the wavelength columns", {
  # `trim_peak` selects the in-peak timepoints and expects a single trace, but
  # was given the whole matrix, which linear-indexes column-major. The cutoff
  # therefore always came from the first wavelength, so purity changed when the
  # columns were reordered. It now uses the wavelength of maximum absorbance at
  # the apex, which is a property of the peak rather than of column order.
  data(Sa_pr)
  x <- Sa_pr[[1]]
  fp <- find_peaks(x[, "210"])
  for (k in seq_len(min(8, nrow(fp)))){
    pos <- as.numeric(fp[k, ])
    expect_equal(get_purity(x, pos),
                 get_purity(x[, rev(colnames(x))], pos),
                 info = paste("peak", k))
  }
})

test_that("get_purity works as intended", {
  data(Sa_pr)
  pos <- as.numeric(find_peaks(Sa_pr[[1]][, "210"])[1, ])
  expect_equal(class(get_purity(Sa_pr[[1]], pos)), "numeric")
})

### noise region selection #####################################################

noisy_chrom <- function(sigma = 0.5, amp = 100, seed = 1){
  set.seed(seed)
  lams <- seq(200, 400, by = 2)
  tt <- 1:200
  d <- amp * outer(exp(-(tt - 100)^2 / (2 * 8^2)),
                   exp(-(lams - 260)^2 / (2 * 25^2))) +
       matrix(rnorm(length(tt) * length(lams), sd = sigma),
              length(tt), length(lams))
  dimnames(d) <- list(tt, lams)
  d
}

test_that("find_noise_quantile always returns timepoints", {
  # `find_noise` selects rows below a fraction of the largest absorbance in the
  # chromatogram, which can select nothing at all. The variance of an empty set
  # is NaN, which is not an error, so purity became NaN silently.
  d <- noisy_chrom()
  expect_length(find_noise(d, noise_threshold = 0.005), 0)
  expect_gt(length(find_noise_quantile(d)), 0)
  expect_true(is.na(get_noise_variance(d, noise_method = "threshold",
                                       noise_threshold = 0.005)))
  expect_true(is.finite(get_noise_variance(d)))
})

test_that("the noise variance does not track the height of the tallest peak", {
  data(Sa_warp)
  d <- Sa_warp[[3]]
  boosted <- d
  apex <- which.max(d[, "210"])
  boosted[apex + -20:20, ] <- boosted[apex + -20:20, ] * 10
  expect_lt(get_noise_variance(boosted) / get_noise_variance(d), 3)
  legacy <- function(x) get_noise_variance(x, noise_method = "threshold",
                                           noise_threshold = 0.001)
  expect_gt(legacy(boosted) / legacy(d), 3)
})

test_that("a single-peak chromatogram yields a finite purity", {
  d <- noisy_chrom()
  expect_false(is.na(get_purity(d, c(100, 75, 125))))
  expect_equal(get_purity(d, c(100, 75, 125)), 1)
})

test_that("noise_window overrides the noise region", {
  d <- noisy_chrom()
  expect_equal(get_noise_variance(d, noise_window = 1:20),
               mean(apply(d[1:20, ], 1, var)))
  expect_true(is.finite(get_noise_variance(d, noise_window = c(-5, 20))))
  expect_identical(get_noise_variance(d, noise_window = c(500, 600)), NA_real_)
})

test_that("a chromatogram clipped at zero warns", {
  d <- noisy_chrom()
  expect_no_warning(get_noise_variance(d))
  d[d < 0] <- 0
  expect_warning(get_noise_variance(d), "clipped at zero")
})

test_that("noise_threshold warns rather than being silently ignored", {
  # the plumbed value used to be dropped on the way to `get_purity_values`
  d <- noisy_chrom()
  expect_warning(get_purity(d, c(100, 75, 125), noise_threshold = 0.001),
                 "noise_threshold")
  expect_silent(get_purity(d, c(100, 75, 125), noise_method = "threshold",
                           noise_threshold = 0.05))
})

test_that("noise_threshold is honored end to end on the legacy path", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  pos <- c(37, 22, 51)
  vals <- vapply(c(0.001, 0.005, 0.01), function(nt){
    get_purity(x, pos, noise_method = "threshold", noise_threshold = nt)
  }, numeric(1))
  expect_false(anyDuplicated(vals) > 0)
  expect_true(all(diff(vals) > 0))
})

### row_var ####################################################################

test_that("row_var matches apply(., 1, var), including for a large offset", {
  data(Sa_warp)
  m <- Sa_warp[[3]][20:60, ]
  expect_equal(unname(row_var(m)), unname(apply(m, 1, var)), tolerance = 1e-12)
  # the one-pass formula loses precision for a large mean with a small variance
  expect_equal(unname(row_var(m + 1e6)), unname(apply(m + 1e6, 1, var)),
               tolerance = 1e-12)
})

### lambdas ####################################################################

test_that("lambdas restricts every part of the calculation", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  pos <- c(37, 22, 51)
  sub <- c("210", "240", "260", "280", "300", "318")
  expect_equal(get_purity(x, pos, lambdas = sub, noise_variance = 0.1),
               get_purity(x[, sub], pos, noise_variance = 0.1))
})

test_that("column order does not matter for a subset of wavelengths", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  pos <- c(37, 22, 51)
  sub <- c("210", "240", "260", "280", "300", "318")
  expect_equal(get_purity(x, pos, lambdas = sub),
               get_purity(x[, rev(colnames(x))], pos, lambdas = sub))
})

### background correction ######################################################

sloped_peak <- function(base = 0, quad = 0, frac = 0, seed = 1){
  set.seed(seed)
  lams <- seq(200, 400, by = 2)
  tt <- 1:200
  sp <- function(mu, sd) exp(-(lams - mu)^2 / (2 * sd^2))
  el <- function(centre, w) exp(-(tt - centre)^2 / (2 * w^2))
  d <- 100 * outer(el(100, 8), sp(260, 25)) +
       100 * frac * outer(el(106, 8), sp(320, 25)) +
       base * outer(seq(1, 2, length.out = length(tt)), sp(230, 60)) +
       quad * outer(seq(0, 1, length.out = length(tt))^2, sp(230, 60)) +
       matrix(rnorm(length(tt) * length(lams), sd = 0.5),
              length(tt), length(lams))
  dimnames(d) <- list(tt, lams)
  d
}

test_that("background correction keeps a pure peak on a curved baseline pure", {
  pos <- c(100, 75, 125)
  expect_lt(get_purity(sloped_peak(quad = 40), pos, background_correction = "none"), 0.75)
  expect_equal(get_purity(sloped_peak(quad = 40), pos), 1)
  expect_gt(get_purity(sloped_peak(quad = 40), pos, background_correction = "constant"),
            get_purity(sloped_peak(quad = 40), pos, background_correction = "none"))
  expect_lt(get_purity(sloped_peak(quad = 40), pos, background_correction = "constant"),
            get_purity(sloped_peak(quad = 40), pos))
})

test_that("background correction does not condemn a narrow peak", {
  # the boundary spectra of `pure_peak()` sit at 13.5% of the apex, so
  # subtracting them leaves a rank-1 residual with a negative coefficient and a
  # correlation of exactly -1 across the wings unless those timepoints are
  # excluded by the cutoff.
  expect_equal(get_purity(pure_peak(), peak_pos), 1)
  expect_equal(get_purity(mixed_peak(), peak_pos, background_correction = "constant") < 1, TRUE)
})

test_that("background correction is applied when a boundary sits on a neighbor", {
  d <- pure_peak()
  d[12, ] <- d[20, ] * 0.8
  for (bg in c("linear", "constant")){
    w <- prep_purity_window(d, peak_pos, background_correction = bg)
    expect_equal(w$background_correction, bg)
  }
  expect_equal(unname(prep_purity_window(d, peak_pos, background_correction = "linear")$m[1, ]),
               rep(0, ncol(d)))
})

### run-length rule ############################################################

test_that("an impure timepoint only counts as part of a run", {
  p <- structure(c(0.2, 0.3, 5, 0.3, 0.2, 0.1, 0.2, 0.3, 0.1),
                 keep = 1:9)
  expect_equal(summarize_purity(p), 1)
  p[6:8] <- 5
  expect_equal(summarize_purity(p), 6 / 9)
  expect_equal(summarize_purity(p, min_run = 1), 5 / 9)
  p <- structure(c(0.2, 5, 5, 0.3, 0.2, 0.1, 0.2, 0.3, 0.1), keep = 1:9)
  expect_equal(summarize_purity(p), 1)
  expect_equal(summarize_purity(p, min_run = 2), 7 / 9)
})

test_that("min_run is validated", {
  for (bad in list(NA, c(2, 3), 0, 2.5, "3"))
    expect_error(get_purity(pure_peak(), peak_pos, min_run = bad),
                 regexp = "`min_run` must be a single whole number")
})

test_that("a misspelled option stops get_purity rather than returning NA", {
  expect_error(get_purity(pure_peak(), peak_pos, background_correction = "bogus"),
               regexp = "should be one of")
  expect_error(get_purity(pure_peak(), peak_pos, reference = "bogus"),
               regexp = "should be one of")
  expect_error(get_purity(pure_peak(), peak_pos, noise_method = "bogus"),
               regexp = "should be one of")
})

test_that("metric = 'max' reports the largest ratio", {
  p <- structure(c(0.2, 0.3, 5, 0.3), keep = 1:4)
  expect_equal(summarize_purity(p, metric = "max"), 5)
})

### sensitivity ################################################################

test_that("purity decreases monotonically with impurity", {
  pos <- c(100, 75, 125)
  fracs <- c(0, 0.05, 0.1, 0.25, 0.5)
  vals <- vapply(fracs, function(f) get_purity(sloped_peak(frac = f), pos),
                 numeric(1))
  expect_equal(vals[1], 1)
  expect_true(all(diff(vals) < 0))
  maxes <- vapply(fracs, function(f){
    get_purity(sloped_peak(frac = f), pos, metric = "max")
  }, numeric(1))
  expect_true(all(diff(maxes[-1]) > 0))
})

### guards #####################################################################

test_that("unassessable peaks return NA_real_ rather than NaN or an error", {
  d <- pure_peak()
  expect_identical(get_purity(d, peak_pos, try = FALSE), NA_real_)
  expect_identical(get_purity(d[, 1, drop = FALSE], peak_pos), NA_real_)
  expect_identical(get_purity(d[, 1:2], peak_pos), NA_real_)
  expect_identical(get_purity(d, c(5, 12, 28)), NA_real_)
  expect_identical(get_purity(d, c(20, 20, 20)), NA_real_)
  expect_identical(get_purity(d, c(20, 0, 28)), NA_real_)
  expect_identical(get_purity(d, c(20, 12, 500)), NA_real_)
  expect_identical(get_purity(d, peak_pos, cutoff = 1), NA_real_)
  dd <- d
  dd[12:28, ] <- NA
  expect_identical(get_purity(dd, peak_pos), NA_real_)
})

test_that("purity tolerates NA-containing chromatograms", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  expect_true(anyNA(x))
  expect_false(is.na(get_purity(x, c(37, 22, 51))))
})

### hoisting ###################################################################

test_that("hoisting the noise variance does not change any purity value", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  fp <- find_peaks(x[, "210"], bounds = TRUE)
  nv <- hoist_noise_variance(x, list())
  pos <- lapply(seq_nrow(fp), function(j) as.numeric(fp[j, c("pos", "lower", "upper")]))
  expect_equal(vapply(pos, function(p) get_purity(x, p, noise_variance = nv), numeric(1)),
               vapply(pos, function(p) get_purity(x, p), numeric(1)))
})

test_that("a noise_variance in purity_args is used rather than estimated", {
  data(Sa_pr)
  expect_identical(hoist_noise_variance(Sa_pr[[1]], list(noise_variance = 1e-6)), 1e-6)
  purity <- function(nv){
    get_peaks(Sa_pr[1], lambdas = "210", estimate_purity = TRUE,
              purity_args = list(noise_variance = nv))[[1]][[1]]$purity
  }
  expect_true(all(purity(1e6) == 1))
  expect_false(all(purity(1e-6) == 1))
})

test_that("an unexpected failure gives NA rather than aborting a peak-fitting run", {
  # purity is a secondary annotation, so a failure should not take down the
  # integration
  d <- pure_peak()
  expect_no_warning(p <- get_purity(d, peak_pos, lambdas = 1:500))
  expect_identical(p, NA_real_)
})

test_that("the purity band fits inside the plot region", {
  # the band is drawn below the trace, so its depth has to be scaled to the
  # data and included in `ylim`, or it falls off the bottom of the plot
  data(Sa_pr)
  pks <- suppressWarnings(get_peaks(Sa_pr[1], lambdas = 210,
                                    estimate_purity = TRUE,
                                    show_progress = FALSE, cl = 1))
  pdf(NULL)
  on.exit(dev.off())
  plot(pks, plot_purity = TRUE)
  usr <- par("usr")
  y <- Sa_pr[[1]][, "210"]
  expect_lt(usr[3], -0.1 * diff(range(y, na.rm = TRUE)) * 2)
})

test_that("the purity band keeps the zoom when the window sits above zero", {
  skip_on_cran()
  data(Sa_pr)
  raised <- lapply(Sa_pr[1], function(x) transfer_metadata(x + 500, x))
  pks <- get_peaks(raised, lambdas = 210, baseline = "flat",
                   show_progress = FALSE, cl = 1)
  pdf(NULL)
  on.exit(dev.off())
  plot(pks, chrom_list = raised, targets = 11.36, xlim = c(11.26, 11.5))
  zoomed <- par("usr")
  plot(pks, chrom_list = raised, targets = 11.36, xlim = c(11.26, 11.5),
       plot_purity = "ratio")
  with_band <- par("usr")
  expect_gt(with_band[3], 400)
  expect_lt(with_band[3], zoomed[3])
})

### threshold formula ##########################################################

test_that("the threshold is the product of the two reliabilities", {
  # the threshold is the similarity expected between two noisy measurements of
  # the same spectrum, which is the attenuation product. The previous
  # `(1 - (a + b))^2` form matched this nowhere and hit zero once a + b >= 1,
  # at which point every peak passed unconditionally.
  d <- pure_peak()
  nv <- 0.5
  w <- prep_purity_window(d, peak_pos, background_correction = "none")
  a <- nv / row_var(w$m)
  b <- nv / var(w$reference)
  expect_equal(unname(gat(d, peak_pos, noise_variance = nv, background_correction = "none",
                          max_threshold = 1)),
               unname(pmax(0, 1 - a) * pmax(0, 1 - b)))
})

test_that("the threshold recovers the expected similarity of two noisy spectra", {
  skip_on_cran()
  set.seed(42)
  nl <- 101
  for (sig in c(10, 3, 1.5)){
    emp <- mean(replicate(2000, {
      true <- rnorm(nl, sd = sig)
      cor(true + rnorm(nl), true + rnorm(nl))^2
    }))
    a <- 1 / (sig^2 + 1)
    expect_equal((1 - a) * (1 - a), emp, tolerance = 0.02)
  }
})

test_that("max_threshold bounds the tolerated dissimilarity", {
  # without a ceiling the threshold approaches 1 for high signal-to-noise
  # peaks, so the denominator approaches 0 and purity tracks peak height
  # rather than spectral dissimilarity
  d <- pure_peak() * 1000
  expect_lt(max(gat(d, peak_pos, background_correction = "none")), 1)
  expect_equal(max(gat(d, peak_pos, background_correction = "none")), 0.99)
  expect_gt(max(gat(d, peak_pos, background_correction = "none", max_threshold = 1)), 0.999)
  expect_true(all(gat(d, peak_pos, background_correction = "none", max_threshold = 0.5) <= 0.5))
})

test_that("purity is less dependent on peak height than it was", {
  data(Sa_pr)
  x <- Sa_pr[[1]]
  fp <- find_peaks(x[, "210"], bounds = TRUE, min_prominence_snr = 0)
  pos <- lapply(seq_nrow(fp), function(j){
    as.numeric(fp[j, c("pos", "lower", "upper")])
  })
  h <- vapply(pos, function(p) x[p[1], "210"], numeric(1))
  capped <- vapply(pos, function(p) get_purity(x, p), numeric(1))
  uncapped <- vapply(pos, function(p) get_purity(x, p, max_threshold = 1),
                     numeric(1))
  expect_lt(abs(cor(log10(h), capped, use = "complete.obs")),
            abs(cor(log10(h), uncapped, use = "complete.obs")))
})

### noise_inflation ###############################################################

test_that("noise_inflation scales the noise variance", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  pos <- c(37, 22, 51)
  expect_equal(get_purity(x, pos, noise_inflation = 3, noise_variance = 0.1),
               get_purity(x, pos, noise_variance = 0.3))
  # a larger noise allowance cannot make a peak look less pure
  vals <- vapply(c(0.5, 1, 2, 4, 8), function(f){
    get_purity(x, pos, noise_inflation = f)
  }, numeric(1))
  expect_true(all(diff(vals) >= 0))
})

test_that("weight is deprecated in favor of noise_inflation", {
  data(Sa_warp)
  x <- Sa_warp[[1]]
  pos <- c(37, 22, 51)
  expect_warning(w <- get_purity(x, pos, weight = 4), "deprecated")
  expect_equal(w, get_purity(x, pos, noise_inflation = 4))
  expect_warning(gat(x, pos, weight = 4), "deprecated")
})

test_that("purity ratios are mapped onto a bounded log scale", {
  # the ratio is unbounded above, so a linear scale had to be capped, which
  # drew 29% of the timepoints on Sa_pr at exactly the same depth and read as
  # if the band were clipped. Stahl displays the purity ratio as a logarithm.
  lim <- -10
  expect_equal(scale_purity_ratio(1, lim), lim)
  expect_equal(scale_purity_ratio(100, lim), 2 * lim)
  expect_equal(scale_purity_ratio(0.01, lim), 0)
  # nothing can escape the band, however extreme the ratio
  extreme <- c(0, 1e-300, 1e300, Inf)
  y <- scale_purity_ratio(extreme, lim)
  expect_true(all(y >= 2 * lim & y <= 0))
  expect_true(is.na(scale_purity_ratio(NA_real_, lim)))
})

### reference spectrum #########################################################

test_that("the apex reference forces a purity ratio of zero at the apex", {
  # cor(apex, apex) is 1 by construction, so the apex spectrum declares every
  # peak perfectly pure at its own apex and the drawn purity trace forks
  data(Sa_pr)
  x <- Sa_pr[[1]]
  fp <- find_peaks(x[, "210"], bounds = TRUE, min_prominence_snr = 0)
  for (j in seq_len(min(6, nrow(fp)))){
    p <- as.numeric(fp[j, c("pos", "lower", "upper")])
    v <- get_purity_values(x, p, reference = "apex")
    if (is.null(v)) next
    expect_equal(v[[p[1] - p[2] + 1]], 0, info = paste("peak", j))
  }
})

test_that("the average reference does not privilege any timepoint", {
  data(Sa_pr)
  x <- Sa_pr[[1]]
  fp <- find_peaks(x[, "210"], bounds = TRUE, min_prominence_snr = 0)
  zeros <- vapply(seq_len(min(6, nrow(fp))), function(j){
    p <- as.numeric(fp[j, c("pos", "lower", "upper")])
    v <- get_purity_values(x, p)
    if (is.null(v)) return(NA_real_)
    v[[p[1] - p[2] + 1]]
  }, numeric(1))
  expect_true(all(zeros[!is.na(zeros)] > 0))
})

test_that("the average reference falls back to the apex for tiny windows", {
  d <- pure_peak()
  w <- prep_purity_window(d, peak_pos, cutoff = 0.99)
  expect_lt(length(w$keep), 3)
  expect_equal(w$reference, w$m[w$apex_i, ])
})

test_that("the average reference flags fewer pure peaks", {
  skip_on_cran()
  set.seed(3)
  nl <- 60; npt <- 21
  lams <- seq(200, 400, length.out = nl)
  sp <- exp(-(lams - 260)^2 / (2 * 25^2))
  el <- exp(-(seq_len(npt) - 11)^2 / (2 * 4^2))
  flagged <- function(ref){
    mean(replicate(150, {
      m <- outer(100 * el, sp) + matrix(rnorm(npt * nl), npt, nl)
      dimnames(m) <- list(seq_len(npt), lams)
      get_purity(m, c(11, 1, npt), background_correction = "none",
                 noise_variance = 1, reference = ref) < 1
    }))
  }
  expect_lt(flagged("average"), flagged("apex"))
})

### simulation-based calibration ###############################################
# `purity_sims()` is defined in helper-purity-sims.R

test_that("pure simulated peak has purity = 1 given the true noise variance", {
  s <- purity_sims()
  expect_equal(get_purity(s$pure$mat, s$pure$pos,
                          noise_variance = s$pure$noise_sd^2), 1)
})

test_that("a medium-SNR pure peak is judged below the threshold cap", {
  s <- purity_sims()
  th <- get_agilent_threshold(s$lo_snr$mat, s$lo_snr$pos,
                              noise_variance = s$lo_snr$noise_sd^2)
  expect_true(all(th < 0.99))
  expect_equal(get_purity(s$lo_snr$mat, s$lo_snr$pos,
                          noise_variance = s$lo_snr$noise_sd^2), 1)
})

test_that("noise variance estimate is within 35% of the true variance", {
  s  <- purity_sims()
  nv <- get_noise_variance(s$pure$mat)
  expect_equal(nv, s$pure$noise_sd^2, tolerance = 0.35)
})

test_that("co-eluting impurity lowers purity vs pure peak", {
  s <- purity_sims()
  expect_lt(
    get_purity(s$mixed$mat, s$mixed$pos),
    get_purity(s$pure$mat,  s$pure$pos)
  )
})

test_that("linear background correction restores purity on sloped baseline", {
  s <- purity_sims()
  purity_none   <- get_purity(s$bg$mat, s$bg$pos, background_correction = "none")
  purity_linear <- get_purity(s$bg$mat, s$bg$pos, background_correction = "linear")
  expect_gt(purity_linear, purity_none)
})

test_that("fit_peaks estimates purity only when asked", {
  data(Sa_pr)
  expect_true(all(is.na(fit_peaks(Sa_pr[[1]], lambda = "210", fit = "gaussian")$purity)))
  pur <- suppressWarnings(fit_peaks(Sa_pr[[1]], lambda = "210", fit = "gaussian",
                                    estimate_purity = TRUE)$purity)
  expect_false(all(is.na(pur)))
})
