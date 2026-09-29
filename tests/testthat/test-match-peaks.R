# Regression tests for `match_peaks` and the `alg = "lptw"` path in
# `correct_rt`, which had no coverage.

test_that("a matched_peak_list does not conflict with the default reference", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)
  expect_equal(attr(mp, "reference"), "119")

  # `reference` defaults to "best", which used to be resolved via
  # `ptw::bestref` before the stored reference was consulted -- so the stored
  # reference was reported as a conflict with an argument nobody supplied
  expect_no_error(correct_rt(Sa_pr, lambdas = 210, alg = "lptw", peaks = mp,
                             cl = 1, show_progress = FALSE))
})

test_that("the stored reference check tolerates integer/double and still catches real conflicts", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)
  # 1L against a stored 1 is not a conflict
  expect_no_error(correct_rt(Sa_pr, lambdas = 210, alg = "lptw", peaks = mp,
                             reference = 1L, cl = 1, show_progress = FALSE))
  # but a genuinely different reference still is
  expect_error(correct_rt(Sa_pr, lambdas = 210, alg = "lptw", peaks = mp,
                          reference = 2, cl = 1, show_progress = FALSE),
               regexp = "differs from the")
})

test_that("match_peaks accepts the reference by name", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = "121")
  expect_equal(attr(mp, "reference"), "121")
  expect_s3_class(get_peaktable(mp, pks = pks), "peak_table")
  expect_error(get_peaktable(mp, pks = pks, hmax = 0.1), regexp = "unused argument.*hmax")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(mp, xlab = "RT", ylim = c(0, 5)))
  expect_error(match_peaks(pks, chrom_list = Sa_pr, reference = "nope"),
               regexp = "not a sample")
  expect_error(match_peaks(pks, chrom_list = Sa_pr, reference = 9),
               regexp = "not a sample")
})

test_that("subsetting a matched_peak_list keeps its reference sample", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = "121")
  sub <- mp[c(4, 2, 1)]
  expect_equal(names(sub)[reference_position(sub)], "121")
  expect_equal(get_peaktable(sub, pks = pks[c(4, 2, 1)])$tab,
               get_peaktable(mp, pks = pks)$tab[c(4, 2, 1), ])
  expect_error(plot(mp[c(1, 3)]), regexp = "reference sample \\(121\\)")
  names(pks) <- c("a", "a", "b", "c")
  expect_error(match_peaks(pks, chrom_list = Sa_pr), regexp = "unique names")
})

test_that("lptw aligns one-dimensional chromatograms instead of silently giving up", {
  skip_on_cran()
  data(Sa_pr)
  x1 <- structure(lapply(Sa_pr, function(m) m[, "210", drop = FALSE]),
                  names = names(Sa_pr))
  # a 1-column "spectrum" is a single number, so `cor` returned all-NA and
  # every sample fell through to an identity warp with only a warning
  msgs <- character(0)
  models <- withCallingHandlers(
    correct_rt(x1, lambdas = 210, alg = "lptw", what = "models", reference = 3,
               cl = 1, show_progress = FALSE),
    warning = function(w){
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # warned once, not once per sample
  expect_length(grep("Spectral matching is unavailable", msgs), 1)

  coefs <- t(vapply(models, function(m) as.numeric(m$warp.coef), numeric(3)))
  identity_coef <- c(0, 1, 0)
  is_identity <- apply(coefs, 1, function(r) isTRUE(all.equal(r, identity_coef)))
  # the reference warps to itself; every other sample must actually move
  expect_equal(sum(is_identity), 1L)
})

test_that("the stored reference follows the samples when chrom_list is reordered", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)
  rev_pr <- rev(Sa_pr)
  mods <- correct_rt(rev_pr, lambdas = 210, alg = "lptw", peaks = mp,
                     reference = names(Sa_pr)[1], cl = 1,
                     show_progress = FALSE, what = "models")
  expect_equal(attr(mods, "reference"), names(Sa_pr)[1])
  mods <- correct_rt(rev_pr, lambdas = 210, alg = "lptw", peaks = mp, cl = 1,
                     show_progress = FALSE, what = "models")
  expect_equal(attr(mods, "reference"), names(Sa_pr)[1])
  expect_error(correct_rt(Sa_pr[2:4], lambdas = 210, alg = "lptw", peaks = mp,
                          cl = 1, show_progress = FALSE),
               regexp = "is not in `chrom_list`")
})

test_that("peaks are matched by sample name, not position", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  shuffled <- pks[c(2, 1, 3, 4)]
  expect_equal(match_peaks(shuffled, chrom_list = Sa_pr, reference = 1)[names(Sa_pr)],
               match_peaks(pks, chrom_list = Sa_pr, reference = 2)[names(Sa_pr)],
               ignore_attr = TRUE)
  fit <- function(peaks)
    correct_rt(Sa_pr, lambdas = 210, alg = "lptw", peaks = peaks, reference = 1,
               cl = 1, show_progress = FALSE, what = "models")
  expect_equal(fit(shuffled), fit(pks))
})

test_that("an unusable reference value is filled with `na.fill`", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)
  pks[[1]][["210"]]$area[1] <- NA
  pt <- get_peaktable(mp, pks = pks, na.fill = NA)
  expect_true(is.na(pt$tab[1, 1]))
})

test_that("lptw drops matched pairs with a missing height", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)
  mp[[2]]$samp_i[1, "I"] <- NA
  expect_no_error(correct_rt(Sa_pr, lambdas = "210", alg = "lptw", peaks = mp,
                             cl = 1, show_progress = FALSE, what = "models"))
  mp[[2]]$samp_i[, "I"] <- NA
  expect_warning(mods <- correct_rt(Sa_pr, lambdas = "210", alg = "lptw",
                                    peaks = mp, cl = 1, show_progress = FALSE,
                                    what = "models"),
                 regexp = "insufficient")
  expect_equal(c(mods[[2]]$warp.coef), c(0, 1, 0))
})

test_that("rt_window is minutes whatever the time unit", {
  skip_on_cran()
  data(Sa_pr)
  pm <- get_peaks(Sa_pr, lambdas = "210", show_progress = FALSE, cl = 1)
  ps <- pm
  for (i in seq_along(ps)){
    ps[[i]][["210"]][c("rt", "start", "end")] <-
      ps[[i]][["210"]][c("rt", "start", "end")] * 60
    attr(ps[[i]], "time_unit") <- "s"
  }
  n_matched <- function(p, x) vapply(
    match_peaks(p, chrom_list = x, reference = 1, spectral_cutoff = 0,
                verbose = FALSE)[-1],
    function(m) nrow(m$matches), numeric(1))
  expect_equal(n_matched(ps, convert_time_units(Sa_pr, time_unit = "s")),
               n_matched(pm, Sa_pr))
})

test_that("matched retention times map to rows on an unevenly spaced axis", {
  skip_on_cran()
  data(Sa_pr)
  gap <- lapply(Sa_pr, function(x) transfer_metadata(x[-(150:160), , drop = FALSE], x))
  pks <- suppressWarnings(get_peaks(gap, lambdas = 210, cl = 1,
                                    show_progress = FALSE))
  m <- suppressWarnings(match_peaks(pks, chrom_list = gap, reference = 1,
                                    lambda = 210))[[2]]
  times <- get_times(gap, 2)
  after <- m$matches[, "rt_samp"] > times[150]
  expect_gt(sum(after), 0)
  expect_lt(max(abs(m$samp_i[after, "rt"] -
                      get_retention_idx(m$matches[after, "rt_samp"], times))),
            0.5)
})
