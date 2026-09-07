### Peak purity ################################################################

gat <- get_agilent_threshold
gss <- get_spectral_similarity
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
  expect_equal(unname(fnoise(d, noise_threshold = 0.02, lambdas = 1:2)), c(1L, 5L))
  expect_equal(unname(fnoise(d, noise_threshold = 0.60, lambdas = 1:2)),
               c(1L, 2L, 4L, 5L))
  expect_length(fnoise(d, noise_threshold = 0.005, lambdas = 1:2), 0)
})

test_that("get_noise_variance averages the variance of the noise rows", {
  d <- rbind(c(1, 1), c(50, 50), c(100, 100), c(50, 50), c(1, 1))
  dimnames(d) <- list(1:5, c("200", "210"))
  # the selected rows are constant, so their variance is exactly zero
  expect_equal(gnv(d, noise_threshold = 0.02, lambdas = 1:2), 0)
  d[1, ] <- c(0, 3)                      # row maximum 3, variance 4.5
  expect_equal(gnv(d, noise_threshold = 0.05, lambdas = 1:2), (4.5 + 0) / 2)
})

test_that("get_noise_variance handles a single noise row", {
  # `x[noise_idx, ]` collapses to a vector when one row qualifies, which made
  # `apply` fail. `get_purity` swallows errors into NA, so purity silently
  # became NA for such chromatograms.
  d <- rbind(c(1, 1), c(50, 50), c(100, 100), c(50, 50), c(9, 9))
  dimnames(d) <- list(1:5, c("200", "210"))
  # only row 1 falls below 2% of the maximum
  expect_equal(unname(fnoise(d, noise_threshold = 0.02, lambdas = 1:2)), 1L)
  expect_equal(gnv(d, noise_threshold = 0.02, lambdas = 1:2), 0)
  expect_silent(gnv(d, noise_threshold = 0.02, lambdas = 1:2))
})

### get_spectral_similarity ####################################################

test_that("a rank-1 peak is spectrally identical at every timepoint", {
  expect_equal(as.numeric(gss(pure_peak(), peak_pos)), rep(1, 17))
})

test_that("a co-eluting second component reduces spectral similarity", {
  expect_lt(min(gss(mixed_peak(), peak_pos)), 1)
})

### get_agilent_threshold ######################################################

test_that("the Agilent threshold spans [0, 1] and falls as noise rises", {
  d <- pure_peak()
  # with no noise the threshold demands perfect similarity
  expect_true(all(gat(d, peak_pos, noise_variance = 0) == 1))
  # the formula is squared after a max(0, .), so it cannot go negative
  expect_true(all(gat(d, peak_pos, noise_variance = 1e12) == 0))
  expect_true(all(gat(d, peak_pos, noise_variance = 1) >=
                  gat(d, peak_pos, noise_variance = 10)))
})

### get_purity_values ##########################################################

test_that("purity values are finite when the threshold reaches 1", {
  # A noise variance of zero drives the threshold to 1, so the denominator
  # `1 - threshold` vanishes. Timepoints matching the apex spectrum are pure
  # regardless, and previously produced NaN, which `get_purity` dropped and so
  # counted as impure.
  pv <- gpv(pure_peak(), peak_pos, noise_variance = 0)
  expect_false(any(is.nan(pv)))
  expect_true(all(pv == 0))
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
