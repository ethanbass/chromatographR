data("Sa")
new.ts <- seq(1,38,by=.01) # choose time-points
new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
dat.pr <- preprocess(X = Sa, dim1 = new.ts, dim2 = new.lambdas)
lam <- "260"
pks_egh <- get_peaks(dat.pr,lambdas = lam, fit = "egh")
out <- get_peaktable(pks_egh)

test_that("get_peaktable works", {
  expect_equal(rownames(out[-c(1:3),]), names(dat.pr))
})
