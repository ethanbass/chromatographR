# library(chromatographR)

### test preprocess ###

data("Sa")
new.ts <- seq(1,38,by=.01) # choose time-points
new.lambdas <- seq(200, 400, by = 2) # choose wavelengths

out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works on matrix", {
  expect_equal(class(out)[1],"matrix")
  expect_equal(rownames(out), as.character(new.ts))
  expect_equal(colnames(out), as.character(new.lambdas))
})

test_that("preprocess returns correct errors", {
  expect_error(preprocess(X=as.data.frame(Sa[[1]])))
})

dat.pr <- preprocess(X = Sa[1:2], dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works on a list", {
  expect_equal(length(dat.pr),length(Sa[1:2]))
  expect_equal(rownames(dat.pr[[1]]), as.character(new.ts))
  expect_equal(colnames(dat.pr[[1]]), as.character(new.lambdas))
  expect_equal(names(dat.pr), names(Sa[1:2]))
})

### test correct_rt ###

warping.models <- correct_rt(dat.pr, what = "models", lambdas=c('210','260','360'), n.zeros = 250)
warp <- correct_rt(chrom_list=dat.pr, models=warping.models)

test_that("correct_rt works", {
  equals(length(warping.models), length(warp), length(Sa[1:2]))
  expect_equal(names(warp), names(Sa[1:2]))
})

### test get_peaks ###
lam <- "260"
pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh")
pks_gaussian <- get_peaks(dat.pr, lambdas = lam, fit = "gaussian")

test_that("get_peaks works", {
  expect_equal(names(pks_egh), names(dat.pr))
  expect_equal(names(pks_gaussian), names(dat.pr))
  expect_equal(names(pks_egh[[1]]), lam)
  expect_equal(names(pks_gaussian[[1]]), lam)
  expect_equal(class(pks_egh), "peak_list")
})

### test get_peaktable ###

pk_tab <- suppressWarnings(get_peaktable(pks_egh, dat.pr))
test_that("get_peaktable works", {
  expect_equal(rownames(pk_tab$tab), names(dat.pr))
  expect_equal(class(pk_tab), "peak_table")
})

