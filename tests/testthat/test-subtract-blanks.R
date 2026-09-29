# Tests for subtract_blanks.

data("Sa_pr")

# ── Correctness ───────────────────────────────────────────────────────────────

test_that("subtract_blanks (method = 'mean') subtracts correctly", {
  res <- subtract_blanks(Sa_pr, blank_id = 1, method = "mean", zero_floor = FALSE)
  expect_equal(length(res), length(Sa_pr) - 1)
  expect_equal(names(res), names(Sa_pr)[-1])
  expect_equal(unname(as.matrix(res[[1]])),
               unname(as.matrix(Sa_pr[[2]]) - as.matrix(Sa_pr[[1]])))
})

test_that("subtract_blanks zero_floor clamps negatives to zero", {
  expect_warning(res_floor <- subtract_blanks(Sa_pr, blank_id = 1,
                                              zero_floor = TRUE),
                 "more than 5% of values")
  res_nofloor <- subtract_blanks(Sa_pr, blank_id = 1, zero_floor = FALSE)
  expect_true(all(res_floor[[1]] >= 0))
  expect_true(any(res_nofloor[[1]] < 0))
  expect_equal(subtract_blanks(Sa_pr, blank_id = 1), res_nofloor)
})

test_that("subtract_blanks zero_floor does not warn when little is clipped", {
  cl <- Sa_pr
  cl[[1]][] <- min(sapply(Sa_pr, min)) - 1
  expect_no_warning(subtract_blanks(cl, blank_id = 1, zero_floor = TRUE))
})

test_that("subtract_blanks preserves matrix attributes", {
  attr(Sa_pr[[2]], "test_attr") <- "preserved"
  res <- subtract_blanks(Sa_pr, blank_id = 1)
  expect_equal(attr(res[[1]], "test_attr"), "preserved")
  attr(Sa_pr[[2]], "test_attr") <- NULL
})

test_that("subtract_blanks (method = 'nearest') picks closest blank", {
  res <- subtract_blanks(Sa_pr, blank_id = 1, method = "nearest",
                         zero_floor = FALSE)
  expect_equal(unname(as.matrix(res[[1]])),
               unname(as.matrix(Sa_pr[[2]]) - as.matrix(Sa_pr[[1]])))
})

test_that("subtract_blanks all methods select blanks correctly with multiple blanks", {
  # Layout: blank1(1), s1(2), s2(3), blank2(4), s3(5)
  # s2 is closer to blank2 by nearest, but only blank1 precedes it.
  # mean of blank1 (10) and blank2 (20) = 15.
  ts  <- as.character(seq(0, 1, length.out = 10))
  wls <- c("210", "254")
  mk  <- function(v) matrix(v, nrow = 10, ncol = 2, dimnames = list(ts, wls))
  cl  <- list(blank1 = mk(10), s1 = mk(100), s2 = mk(200), blank2 = mk(20), s3 = mk(300))

  res_mean <- subtract_blanks(cl, pattern = "blank", method = "mean",      zero_floor = FALSE)
  res_near <- subtract_blanks(cl, pattern = "blank", method = "nearest",   zero_floor = FALSE)
  res_prec <- subtract_blanks(cl, pattern = "blank", method = "preceding", zero_floor = FALSE)

  # mean subtracts the average of both blanks (15) from every sample
  expect_true(all(res_mean[["s1"]] == 85))    # 100 - 15
  expect_true(all(res_mean[["s2"]] == 185))   # 200 - 15
  expect_true(all(res_mean[["s3"]] == 285))   # 300 - 15

  # s1: nearest and preceding both pick blank1
  expect_true(all(res_near[["s1"]] == 90))    # 100 - 10
  expect_true(all(res_prec[["s1"]] == 90))

  # s2: nearest picks blank2 (dist 1 < dist 2), preceding picks blank1 (blank2 not yet seen)
  expect_true(all(res_near[["s2"]] == 180))   # 200 - 20
  expect_true(all(res_prec[["s2"]] == 190))   # 200 - 10

  # s3: nearest and preceding both pick blank2
  expect_true(all(res_near[["s3"]] == 280))   # 300 - 20
  expect_true(all(res_prec[["s3"]] == 280))
})

test_that("subtract_blanks interpolates blanks recorded on a different time grid", {
  wls <- c("210", "254")
  blank <- matrix(rep(seq(0, 10, by = 0.5), 2), ncol = 2,
                  dimnames = list(seq(0, 10, by = 0.5), wls))
  smp_t <- seq(0.25, 10.25, by = 0.5)
  smp <- matrix(100, nrow = length(smp_t), ncol = 2,
                dimnames = list(smp_t, wls))
  res <- subtract_blanks(list(blank = blank, smp = smp), blank_id = 1,
                         zero_floor = FALSE)
  expect_equal(dim(res$smp), dim(smp))
  expect_equal(rownames(res$smp), rownames(smp))
  expect_equal(unname(res$smp[, 1]), 100 - pmin(smp_t, 10))
})

test_that("subtract_blanks re-grids the mean blank for each sample's time axis", {
  wls <- c("210", "254")
  mk <- function(t, v) matrix(v, nrow = length(t), ncol = 2,
                              dimnames = list(t, wls))
  cl <- list(b1 = mk(0:10, 0:10), b2 = mk(0:10, 2 * 0:10),
             s1 = mk(0:10, 100), s2 = mk(0:10 + 0.5, 100))
  res <- subtract_blanks(cl, blank_id = 1:2, zero_floor = FALSE)
  expect_equal(unname(res$s1[, 1]), 100 - 1.5 * 0:10)
  expect_equal(unname(res$s2[, 1]), 100 - 1.5 * pmin(0:10 + 0.5, 10))
})

test_that("subtract_blanks warns when a blank is extrapolated past 5 points", {
  wls <- c("210", "254")
  mk <- function(t) matrix(1, nrow = length(t), ncol = 2, dimnames = list(t, wls))
  expect_no_warning(subtract_blanks(list(b = mk(0:10), s = mk(0:15)), blank_id = 1))
  expect_warning(subtract_blanks(list(b = mk(0:10), s = mk(0:16)), blank_id = 1),
                 "for: s")
})

test_that("subtract_blanks errors when wavelengths differ", {
  a <- matrix(1, 3, 2, dimnames = list(1:3, c("210", "254")))
  b <- matrix(1, 3, 2, dimnames = list(1:3, c("210", "280")))
  expect_error(subtract_blanks(list(a, b), blank_id = 1), "same wavelengths")
})

# ── Pattern matching ──────────────────────────────────────────────────────────

test_that("subtract_blanks pattern argument identifies blanks by name", {
  named <- Sa_pr[1:2]
  names(named) <- c("Blank_01", "Sample_01")
  res <- subtract_blanks(named, pattern = "^Blank")
  expect_equal(length(res), 1)
  expect_equal(names(res), "Sample_01")
})

# ── Error handling ────────────────────────────────────────────────────────────

test_that("subtract_blanks (method = 'preceding') errors when no preceding blank", {
  expect_error(
    subtract_blanks(Sa_pr, blank_id = 2, method = "preceding"),
    "No preceding blank found"
  )
})

test_that("subtract_blanks errors on out-of-range blank_id", {
  expect_error(subtract_blanks(Sa_pr, blank_id = 99), "out of range")
})

test_that("subtract_blanks errors when pattern matches nothing", {
  named <- Sa_pr[1:2]
  names(named) <- c("Sample_01", "Sample_02")
  expect_error(subtract_blanks(named, pattern = "^Blank"), "matched pattern")
})

test_that("subtract_blanks errors when blank_id name not found", {
  named <- Sa_pr[1:2]
  names(named) <- c("Sample_01", "Sample_02")
  expect_error(subtract_blanks(named, blank_id = "Blank_01"),
               "not found in")
})

test_that("subtract_blanks breaks nearest-blank ties toward the earlier blank", {
  mk <- function(v) matrix(v, nrow = 3, ncol = 1,
                           dimnames = list(c("1", "2", "3"), "210"))
  cl <- list(b1 = mk(1), s2 = mk(10), s3 = mk(10), s4 = mk(10), b5 = mk(5))
  res <- subtract_blanks(cl, blank_id = c(5, 1), method = "nearest")
  expect_equal(res$s3[, 1], rep(9, 3), ignore_attr = TRUE)
})

test_that("subtract_blanks subtracts row by row when times are not numeric", {
  mk <- function(v, t) matrix(v, nrow = length(t), ncol = 1,
                              dimnames = list(t, "210"))
  same <- list(b = mk(1, paste0("t", 1:3)), s = mk(10, c("1", "2", "3")))
  expect_equal(subtract_blanks(same, blank_id = 1)$s[, 1], rep(9, 3),
               ignore_attr = TRUE)
  longer <- list(b = mk(1, paste0("t", 1:4)), s = mk(10, c("1", "2", "3")))
  expect_error(subtract_blanks(longer, blank_id = 1),
               "need numeric retention times")
  one <- list(b = mk(1, "2"), s = mk(10, c("1", "2", "3")))
  expect_error(subtract_blanks(one, blank_id = 1), "at least two time points")
})
