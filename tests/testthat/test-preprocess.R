data("Sa")
new.ts <- seq(1,38,by=.01) # choose time-points
new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
out <- preprocess(X = Sa, dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works", {
  expect_equal(length(out),length(Sa))
  expect_equal(rownames(out[[1]]), as.character(new.ts))
  expect_equal(colnames(out[[1]]), as.character(new.lambdas))
  expect_equal(names(out), names(Sa))
  equals(sapply(out, dim)[1,])
  equals(sapply(out, dim)[2,])
})
