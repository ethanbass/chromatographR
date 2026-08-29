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
  res_floor   <- subtract_blanks(Sa_pr, blank_id = 1, zero_floor = TRUE)
  res_nofloor <- subtract_blanks(Sa_pr, blank_id = 1, zero_floor = FALSE)
  expect_true(all(res_floor[[1]] >= 0))
  expect_true(any(res_nofloor[[1]] < 0))
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
