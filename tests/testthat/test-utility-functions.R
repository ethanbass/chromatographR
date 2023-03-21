#### test utilities
data(Sa_pr)

test_that("get_times works as expected", {
  ts <- get_times(Sa_pr)
  expect_equal(ts, as.numeric(rownames(Sa_pr[[1]])))
})

test_that("get_lambdas works as expected", {
  lambdas <- get_lambdas(Sa_pr)
  expect_equal(lambdas, as.numeric(colnames(Sa_pr[[1]])))
})

test_that("choose_apply_fnc works as expected", {
  skip_if_not_installed("pbapply")
  fn <- choose_apply_fnc(progress_bar = TRUE)
  expect_equal(fn, pbapply::pblapply)
  fn <- choose_apply_fnc(progress_bar = FALSE)
  expect_equal(fn, lapply)
})

test_that("check_peaktable works as expected", {
  data(pk_tab)
  expect_null(check_peaktable(pk_tab))
  expect_error(check_peaktable(pks_egh))
})


# test_that("get_chromlist works as expected", {
#   expect_equal(get_chrom_list(x = pk_tab)[[1]], dat.pr[[1]])
# })

test_that("get_retention_idx works as expected", {
  idx <- get_retention_idx(RT = 18.49, times = get_times(Sa_pr))
  expect_equal(idx, 425)
  expect_error(get_retention_idx(RT = 9, times = get_times(Sa_pr)))
  expect_error(get_retention_idx(RT = 20, times = get_times(Sa_pr)))
})

test_that("check_idx works as expected", {
  expect_error(check_idx(idx = 1000, chrom_list = Sa_pr))
})

test_that("get_lambda_idx works as expected", {
  lam <- get_lambda_idx(lambda = 200, lambdas = get_lambdas(chrom_list = Sa_pr))
  expect_equal(lam, 1)
  expect_error(get_lambda_idx(lambda=190, lambdas = get_lambdas(chrom_list = Sa_pr)))
  expect_error(get_lambda_idx(lambda=400, lambdas = get_lambdas(chrom_list = Sa_pr)))
  expect_error(get_lambda_idx(lambda="max", lambdas = get_lambdas(chrom_list = Sa_pr),
                              allow_max = FALSE))
  # lam <- get_lambda_idx(lambda="max", lambdas = get_lambdas(chrom_list = dat.pr),
  #                       y = unlist(dat.pr[[1]][500, , drop = TRUE]),
  #                       allow_max = TRUE)
})

test_that("check_for_pkg functions as expected",{
  expect_error(check_for_pkg(pkg = "nonsense-package-that-doesn't-exist"))
})

test_that("reshape_peaktable works as expected",{
  data(pk_tab)
  pktab_long<-reshape_peaktable(pk_tab)
  expect_equal(ncol(pktab_long),4)
  expect_equal(nrow(pktab_long), nrow(pk_tab)*ncol(pk_tab))
})

test_that("reshape_chroms works as expected", {
  chrom_list_long <- reshape_chroms(Sa_pr,lambdas = "210")
  expect_equal(ncol(chrom_list_long), 4)
  expect_equal(nrow(chrom_list_long), length(Sa_pr)*nrow(Sa_pr[[1]]))
  expect_equal(unique(chrom_list_long$lambda), 210)
  expect_equal(unique(chrom_list_long$sample), names(Sa_pr))
  chrom_list_long_subset <- reshape_chroms(Sa_pr,lambdas = "210",idx = c(1:2))
  expect_equal(ncol(chrom_list_long_subset), 4)
  expect_equal(nrow(chrom_list_long_subset), 2*nrow(Sa_pr[[1]]))
  expect_equal(unique(chrom_list_long_subset$sample), names(Sa_pr)[1:2])
})