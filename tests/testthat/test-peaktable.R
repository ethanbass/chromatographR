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
