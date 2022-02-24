data("Sa")
new.ts <- seq(1,38,by=.01) # choose time-points
new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
dat.pr <- preprocess(X = Sa, dim1 = new.ts, dim2 = new.lambdas)
lam <- "260"
out_egh <- get_peaks(dat.pr,lambdas = lam, fit = "egh")
out_gaussian <- get_peaks(dat.pr,lambdas = lam, fit = "gaussian")

test_that("get_peaks works", {
  expect_equal(names(out_egh), names(dat.pr))
  expect_equal(names(out_gaussian), names(dat.pr))
  expect_equal(names(out_egh[[1]]), lam)
})
