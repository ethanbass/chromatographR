test_that("warp_sample_ptw matches approx for forward mode", {
  data(Sa_pr)
  mat <- Sa_pr[[1]]
  coef <- c(0.5, 0.98, 0.0001)
  tp <- seq_len(nrow(mat))
  w <- ptw::warp.time(tp, coef)
  
  expected <- ptw:::warp.sample(t(mat), w, mode = "forward")
  result <- warp_sample_ptw(mat, coef, mode = "forward")
  
  expect_equal(result, expected, tolerance = 1e-6, ignore_attr=TRUE)
})

test_that("warp_sample_ptw returns NA for out-of-bounds regions", {
  data(Sa_pr)
  mat <- Sa_pr[[1]]
  coef <- c(50, 1, 0)
  result <- warp_sample_ptw(mat, coef, mode = "forward")
  expect_true(all(is.na(head(result, n = 50))))
})

test_that("warp_sample_ptw returns NA for out-of-bounds regions", {
  data(Sa_pr)
  mat <- Sa_pr[[1]]
  coef <- c(-50, 1, 0)
  result <- warp_sample_ptw(mat, coef, mode = "forward")
  expect_true(all(is.na(tail(result, n = 50))))
})

test_that("warp_sample_ptw handles identity warp", {
  data(Sa_pr)
  mat <- Sa_pr[[1]]
  coef <- c(0, 1, 0)
  result <- warp_sample_ptw(mat, coef, mode = "forward")
  expect_equal(result, mat, tolerance = 1e-6)
})
