# Tests for get_peaktable clustering correctness using synthetic chromatograms.
#
# Two samples, each with two Gaussian peaks.
# Peak A: same RT in both samples (gap = 0).
# Peak B: 50-index gap between samples; on this time grid (dt ≈ 0.0167 min/index)
#         that is ~0.83 min.

.n  <- 301
.ts <- seq(0, 5, length.out = .n)   # dt ≈ 0.0167 min/index

.gauss <- function(centers, heights) {
  x <- seq_len(.n); y <- numeric(.n)
  for (i in seq_along(centers))
    y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * 10^2))
  matrix(y, ncol = 1, dimnames = list(as.character(.ts), "220"))
}

.cl <- list(
  s1 = .gauss(centers = c(101, 201), heights = c(1000, 800)),
  s2 = .gauss(centers = c(101, 251), heights = c(1000, 800))
)

.pks <- get_peaks(.cl, lambdas = "220", fit = "gaussian", smooth_type = "none",
                  estimate_purity = FALSE, show_progress = FALSE)

# ── hmax = 0.2: nearby peaks merge, distant peaks stay separate ───────────────

test_that("get_peaktable clusters nearby peaks into one column and separates distant peaks", {
  skip_on_cran()
  pt <- get_peaktable(.pks, hmax = 0.2)

  # 3 distinct peaks total → 3 columns
  expect_equal(ncol(pt$tab), 3L)

  # The shared peak A column: both samples are non-zero
  rt_A  <- as.numeric(.ts[101])
  col_A <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_A))
  expect_gt(pt$tab["s1", col_A], 0)
  expect_gt(pt$tab["s2", col_A], 0)

  # Peak B lands in two separate columns, one per sample
  rt_B_s1  <- as.numeric(.ts[201])
  rt_B_s2  <- as.numeric(.ts[251])
  col_B_s1 <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_B_s1))
  col_B_s2 <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_B_s2))
  expect_true(col_B_s1 != col_B_s2)
  expect_equal(pt$tab["s2", col_B_s1], 0)
  expect_equal(pt$tab["s1", col_B_s2], 0)
})

# ── hmax = 1.0: large enough to merge even the distant peaks ─────────────────

test_that("get_peaktable merges distant peaks when hmax is large enough", {
  skip_on_cran()
  pt <- get_peaktable(.pks, hmax = 1.0)

  # hmax = 1.0 min >> 0.83 min gap → all peaks collapse to 2 columns
  expect_equal(ncol(pt$tab), 2L)
  expect_true(all(pt$tab > 0))
})

# ── metadata reported by only some peaks must survive summarization ──────────

test_that("weighted.mean summarization keeps metadata that only some peaks report", {
  skip_on_cran()
  data(Sa, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa, lambdas = "210", fit = "raw", cl = 1,
                                    show_progress = FALSE))
  fwhm <- unlist(lapply(pks, function(s) s[["210"]]$FWHM))
  expect_true(any(is.finite(fwhm)))
  expect_true(any(!is.finite(fwhm)))

  m <- get_peaktable(pks, summarize_by = "weighted.mean")$pk_meta
  summarized <- as.numeric(m["FWHM", ])
  kept <- summarized[is.finite(summarized)]
  expect_gt(length(kept), sum(!is.finite(summarized)))
  expect_true(all(kept >= min(fwhm, na.rm = TRUE)))
  expect_true(all(kept <= max(fwhm, na.rm = TRUE)))

  expect_true(all(!is.finite(as.numeric(m["sd", ]))))
})

test_that("weighted_mean_finite ignores unusable weights and values", {
  expect_equal(weighted_mean_finite(c(1, 3), c(1, 1)), 2)
  expect_equal(weighted_mean_finite(c(1, 3), c(1, NaN)), 1)
  expect_equal(weighted_mean_finite(c(1, NaN), c(1, 5)), 1)
  expect_equal(weighted_mean_finite(c(2, 4), c(0, 0)), 3)
  expect_equal(weighted_mean_finite(c(2, 4), c(-1, -1)), 3)
  expect_identical(weighted_mean_finite(c(NA_real_, NaN), c(1, 1)), NA_real_)
})

test_that("filter_peaktable rejects min_converged without a converged row", {
  skip_on_cran()
  data(pk_tab, envir = environment())
  if (!"converged" %in% rownames(pk_tab$pk_meta)){
    expect_error(filter_peaktable(pk_tab, min_converged = 0.5),
                 "no `converged` row")
  } else {
    skip("stored pk_tab already carries a converged row")
  }
})

test_that("filter_peaktable matches one peak per requested retention time", {
  data(pk_tab, envir = environment())
  rt <- as.numeric(pk_tab$pk_meta["rt", ])
  one <- filter_peaktable(pk_tab, rts = rt[2], tol = 5)
  expect_equal(ncol(one$tab), 1)
  expect_equal(as.numeric(one$pk_meta["rt", ]), rt[2])

  two <- filter_peaktable(pk_tab, rts = rt[c(2, 10)], tol = 5)
  expect_equal(ncol(two$tab), 2)
  expect_setequal(as.numeric(two$pk_meta["rt", ]), rt[c(2, 10)])

  expect_warning(none <- filter_peaktable(pk_tab, rts = 99, tol = 0.01),
                 regexp = "not identified in the peak table")
  expect_equal(ncol(none$tab), 0)
})
