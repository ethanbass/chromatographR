data("Sa")
new.ts <- seq(1,38,by=.01) # choose time-points
new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
dat.pr <- preprocess(X = Sa, dim1 = new.ts, dim2 = new.lambdas)

warping.models <- correct_rt(dat.pr, what = "models", lambdas=c('210','260','360'), n.zeros = 250)
warp <- correct_rt(chrom_list=dat.pr, models=warping.models)

test_that("correct_rt works", {
  equals(length(warping.models), length(warp), length(Sa))
  expect_equal(names(warp), names(Sa))
})
