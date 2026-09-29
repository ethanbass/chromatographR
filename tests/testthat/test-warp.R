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
  expect_equal(result, mat, tolerance = 1e-6, ignore_attr = "baseline_param")
})

test_that("VPdtw models can be re-applied to reproduce the direct warp", {
  skip_on_cran()
  data(Sa_pr)
  direct <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                       show_progress = FALSE)
  models <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                       what = "models", show_progress = FALSE)
  expect_s3_class(models, "VPdtw")
  # the reference chosen when fitting is recorded on the models
  expect_equal(attr(models, "reference"), names(Sa_pr)[3])
  expect_equal(attr(models, "parameters")$reference, names(Sa_pr)[3])
  expect_equal(attr(models, "parameters")$max_shift, 1)
  old <- models
  attr(old, "reference") <- NULL
  expect_error(correct_rt(Sa_pr, alg = "vpdtw", models = old, reference = 1,
                          show_progress = FALSE),
               regexp = "differs from the reference stored")

  reused <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", models = models,
                       show_progress = FALSE)
  expect_equal(names(reused), names(Sa_pr))
  # `warping_args` records how each call was made, so it differs by design;
  # everything else must match
  for (i in seq_along(reused)){
    expect_equal(reused[[i]], direct[[i]],
                 ignore_attr = "warping_args")
  }
})

test_that("VPdtw `fill_zeros` controls out-of-bounds regions", {
  skip_on_cran()
  data(Sa_pr)
  models <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                       what = "models", show_progress = FALSE)
  na_warp <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", models = models,
                        fill_zeros = FALSE, show_progress = FALSE)
  zero_warp <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          models = models, fill_zeros = TRUE,
                          show_progress = FALSE)
  # warping shifts samples off the ends of the reference grid, so at least one
  # sample must have out-of-bounds cells for this to be a real test
  expect_true(any(vapply(na_warp, function(x) any(is.na(x)), logical(1))))
  expect_false(any(vapply(zero_warp, function(x) any(is.na(x)), logical(1))))
  for (i in seq_along(na_warp)){
    oob <- is.na(na_warp[[i]])
    expect_true(all(zero_warp[[i]][oob] == 0))
    expect_equal(zero_warp[[i]][!oob], na_warp[[i]][!oob])
  }
})

test_that("the pre-0.9.0 VPdtw default is reproduced by explicit arguments", {
  skip_on_cran()
  data(Sa_pr)
  allmats <- sapply(Sa_pr, function(x) x[, "210", drop = FALSE])
  old <- VPdtw::VPdtw(query = allmats, reference = allmats[, 3],
                      penalty = VPdtw::dilation(allmats[, 3], 350) / 5,
                      maxshift = 50)
  dt <- get_time_resolution(Sa_pr)
  expect_warning(new <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                                   reference = 3, penalty = 0.2,
                                   dilation_span = 350 * dt,
                                   max_shift = 50 * dt, what = "models",
                                   show_progress = FALSE),
                 regexp = "nearly constant")
  expect_equal(new$xVals, old$xVals)
  expect_equal(new$shift, old$shift)
})

test_that("a higher VPdtw penalty allows less warping", {
  skip_on_cran()
  data(Sa_pr)
  allmats <- sapply(Sa_pr, function(x) x[, "210", drop = FALSE])
  dil <- VPdtw::dilation(allmats[, 3], 16)
  moves <- function(p){
    m <- VPdtw::VPdtw(query = allmats, reference = allmats[, 3],
                      penalty = dil * p, maxshift = 50)
    sum(diff(m$shift) != 0, na.rm = TRUE)
  }
  expect_gt(moves(0.01), moves(0.5))
  expect_gt(moves(0.5),  moves(1))
  expect_gt(moves(1),    moves(5))
})

test_that("a divisor-style VPdtw penalty warns but is still used", {
  skip_on_cran()
  data(Sa_pr)
  # the warning echoes the supplied value and names its reciprocal, so it is
  # clear that both numbers refer to the caller's own value
  expect_warning(m <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                                 penalty = 20, what = "models",
                                 show_progress = FALSE),
                 regexp = "(?s)you supplied 20.*reciprocal: `penalty = 0\\.05`",
                 perl = TRUE)
  expect_s3_class(m, "VPdtw")
  # a value needing rounding is reported to 4 significant figures
  expect_warning(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", penalty = 3,
                            what = "models", show_progress = FALSE),
                 regexp = "`penalty = 0\\.3333`")
  expect_no_warning(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                               penalty = 1, what = "models",
                               show_progress = FALSE))
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          reference = "nope", show_progress = FALSE),
               regexp = "`reference` must be")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          reference = "2", show_progress = FALSE),
               regexp = "`reference` must be")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", penalty = -1,
                          show_progress = FALSE),
               regexp = "must be >= 0")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", penalty = NA,
                          show_progress = FALSE),
               regexp = "single non-missing number")
  # a zero penalty is legitimate: it is the paper's unpenalized DTW alignment
  expect_s3_class(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                             penalty = 0, what = "models",
                             show_progress = FALSE), "VPdtw")
})

test_that("`regularization` reaches fit_penalized_warp for alg = 'lptw'", {
  skip_on_cran()
  data(Sa_pr)
  as_coefs <- function(m) do.call(rbind, lapply(m, function(x) x$warp.coef[1, ]))
  coefs <- function(regularization, penalty = 0.2)
    as_coefs(correct_rt(Sa_pr, lambdas = "210", alg = "lptw", what = "models",
                        regularization = regularization, penalty = penalty,
                        show_progress = FALSE))
  expect_false(identical(coefs(regularization = 0), coefs(regularization = 1e6)))
  default <- as_coefs(correct_rt(Sa_pr, lambdas = "210", alg = "lptw",
                                 what = "models", show_progress = FALSE))
  expect_equal(default, coefs(regularization = 0))
  # the vpdtw `penalty` must not affect an lptw alignment
  expect_equal(coefs(0, penalty = 0.2), coefs(0, penalty = 1))

  # a larger value pulls the warp toward the identity warp
  dist_from_identity <- function(p){
    co <- coefs(regularization = p)
    mean(sqrt(rowSums((co - matrix(c(0, 1, 0), nrow(co), 3, byrow = TRUE))^2)))
  }
  d <- vapply(c(0, 1, 100, 1e6), dist_from_identity, numeric(1))
  expect_true(all(diff(d) < 0))
  expect_lt(d[4], 0.01)

  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "lptw",
                          regularization = -1, show_progress = FALSE),
               regexp = "`regularization` must be a single number >= 0")
  # unbounded above
  expect_s3_class(correct_rt(Sa_pr, lambdas = "210", alg = "lptw",
                             what = "models", regularization = 100,
                             show_progress = FALSE), "lptw_list")
})

test_that("`dilation_span` controls the VPdtw penalty window", {
  skip_on_cran()
  data(Sa_pr)
  models <- function(dilation_span)
    correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", what = "models",
               dilation_span = dilation_span, show_progress = FALSE)
  default <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", what = "models",
                        show_progress = FALSE)
  # Sa_pr is sampled every 0.02 minutes, so the default of 0.3 is 15 timepoints
  expect_equal(default$shift, models(dilation_span = 0.3)$shift)
  expect_false(identical(default$shift, models(dilation_span = 1)$shift))
  # narrower windows permit more warping, because the penalty falls toward
  # zero in the baseline between peaks
  max_shift <- function(sp) max(abs(models(dilation_span = sp)$shift), na.rm = TRUE)
  expect_gt(max_shift(0.3), max_shift(2))
})

test_that("the same `dilation_span` gives the same penalty at any sampling rate", {
  skip_on_cran()
  data(Sa)
  rng <- range(as.numeric(rownames(Sa[[1]])))
  grid <- function(by) suppressWarnings(
    preprocess(Sa, dim1 = seq(rng[1] + 0.1, rng[2] - 0.1, by = by),
               dim2 = c(200, 210, 220), cl = 1, show_progress = FALSE))
  times <- seq(11, 27, by = 0.02)
  lambda <- function(pr, span){
    m <- suppressWarnings(
      correct_rt(pr, lambdas = "210", alg = "vpdtw", what = "models",
                 reference = 3, dilation_span = span, show_progress = FALSE))
    approx(as.numeric(rownames(pr[[3]])), m$penalty, times)$y
  }
  rel_diff <- function(a, b) max(abs(a - b), na.rm = TRUE) / max(a, na.rm = TRUE)
  coarse <- grid(0.02)
  fine <- grid(0.005)
  # a width in time units is the same physical window on either grid, so the
  # penalty agrees; the same count of timepoints would not
  expect_lt(rel_diff(lambda(coarse, 0.3), lambda(fine, 0.3)), 0.2)
  expect_gt(rel_diff(lambda(coarse, 16 * 0.02), lambda(fine, 16 * 0.005)), 0.5)
})

test_that("a degenerate `dilation_span` warns", {
  skip_on_cran()
  data(Sa_pr)
  # Sa_pr is 434 points at 0.02 minutes, so a span of 5 minutes covers the
  # whole trace and flattens the penalty
  expect_silent(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                           what = "models", show_progress = FALSE))
  expect_warning(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                            what = "models", dilation_span = 5,
                            show_progress = FALSE),
                 regexp = "wider than the 434-point chromatogram")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          dilation_span = 0, show_progress = FALSE),
               regexp = "must be a single positive number")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          dilation_span = c(10, 20), show_progress = FALSE),
               regexp = "must be a single positive number")
})

test_that("retention time correction requires a common time axis", {
  skip_on_cran()
  data(Sa_pr)
  dt <- get_time_resolution(Sa_pr)
  shift_axis <- function(x, by){
    rownames(x) <- as.numeric(rownames(x)) + by
    x
  }
  shifted <- Sa_pr
  shifted[[2]] <- shift_axis(shifted[[2]], 5 * dt)
  for (alg in c("ptw", "lptw", "vpdtw"))
    expect_error(correct_rt(shifted, lambdas = "210", alg = alg, reference = 1,
                            cl = 1, show_progress = FALSE),
                 regexp = "share a time axis", info = alg)
  expect_error(correct_rt_group(shifted, lambdas = "210",
                                groups = c("a", "a", "b", "b"), cl = 1,
                                show_progress = FALSE),
               regexp = "share a time axis")
  shorter <- Sa_pr
  shorter[[3]] <- shorter[[3]][-1, , drop = FALSE]
  expect_error(correct_rt(shorter, lambdas = "210", reference = 1, cl = 1,
                          show_progress = FALSE),
               regexp = "share a time axis")
  jittered <- Sa_pr
  jittered[[2]] <- shift_axis(jittered[[2]], 0.5 * dt)
  expect_no_error(correct_rt(jittered, lambdas = "210", reference = 1, cl = 1,
                             show_progress = FALSE))
})

test_that("stale timepoint counts are rejected rather than read as times", {
  skip_on_cran()
  data(Sa_pr)
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          dilation_span = 350, show_progress = FALSE),
               regexp = "changed units in chromatographR 0.9.0")
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw",
                          max_shift = 50, show_progress = FALSE),
               regexp = "`max_shift = 50` is wider than the chromatogram")
  expect_error(correct_rt(Sa_pr, lambdas = "210", trwdth = 200,
                          show_progress = FALSE, cl = 1),
               regexp = "`trwdth = 200` is wider than the chromatogram")
})

test_that("the default max_shift is narrowed to fit a short chromatogram", {
  skip_on_cran()
  data(Sa_pr)
  short <- lapply(Sa_pr, function(x){
    t <- as.numeric(rownames(x))
    transfer_metadata(x[t >= 10 & t < 10.8, , drop = FALSE], x)
  })
  expect_no_error(withCallingHandlers(
    warped <- correct_rt(short, lambdas = "210", alg = "vpdtw", reference = 1),
    warning = function(w){
      if (grepl("three quarters of maxshift", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }))
  expect_equal(lapply(warped, colnames), lapply(short, colnames))
  expect_error(correct_rt(short, lambdas = "210", alg = "vpdtw", reference = 1,
                          max_shift = 1),
               regexp = "`max_shift = 1` is wider than the chromatogram")
})

test_that("the deprecated `maxshift` is read as a count of timepoints", {
  skip_on_cran()
  data(Sa_pr)
  fit <- function(...) suppressWarnings(
    correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", reference = 1,
               what = "models", show_progress = FALSE, ...))
  expect_warning(
    correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", reference = 1,
               what = "models", show_progress = FALSE, maxshift = 25),
    "`maxshift = 25` was read as a count of timepoints, i.e. `max_shift = 0.5`")
  expect_equal(fit(maxshift = 25)$shift, fit(max_shift = 0.5)$shift)
})

test_that("warping widths are minutes whatever the time unit of the axis", {
  skip_on_cran()
  data(Sa_pr)
  fit <- function(x) suppressWarnings(
    correct_rt(x, lambdas = "210", alg = "vpdtw", reference = 1,
               what = "models", show_progress = FALSE))
  expect_equal(fit(convert_time_units(Sa_pr, time_unit = "s"))$shift,
               fit(Sa_pr)$shift)
  expect_equal(time_to_points(0.2, "trwdth",
                              convert_time_units(Sa_pr, time_unit = "s")),
               time_to_points(0.2, "trwdth", Sa_pr))
})

test_that("`trwdth` reaches `ptw` rather than being silently dropped", {
  skip_on_cran()
  data(Sa_pr)
  coefs <- function(tw) lapply(
    correct_rt(Sa_pr, lambdas = "210", what = "models", trwdth = tw,
               show_progress = FALSE, cl = 1),
    function(x) x$warp.coef)
  expect_false(identical(coefs(0.05), coefs(2)))
})

test_that("pre-fit models are applied without re-deriving the reference", {
  skip_on_cran()
  data(Sa_pr)
  models <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw", what = "models",
                       show_progress = FALSE)
  one <- correct_rt(Sa_pr[1], lambdas = "210", alg = "ptw",
                    models = models[1], show_progress = FALSE)
  expect_length(one, 1)
  expect_true(is.matrix(one[[1]]))
  expect_equal(dim(one[[1]]), dim(Sa_pr[[1]]))
  # the reference stored on the models is used rather than recomputed, and a
  # conflicting explicit reference is refused rather than silently ignored
  expect_equal(attr(models, "reference"), names(Sa_pr)[3])
  expect_error(correct_rt(Sa_pr, lambdas = "210", alg = "ptw", models = models,
                          reference = 1, show_progress = FALSE),
               regexp = "differs from the reference stored in the supplied `models`")
  # supplying the stored reference explicitly is fine
  expect_equal(correct_rt(Sa_pr, lambdas = "210", alg = "ptw", models = models,
                          reference = 3, show_progress = FALSE),
               correct_rt(Sa_pr, lambdas = "210", alg = "ptw", models = models,
                          show_progress = FALSE))
})

test_that("correct_rt_group handles groups containing a single sample", {
  skip_on_cran()
  data(Sa_pr)
  # each group's warp is applied via `models =`, so a one-member group means a
  # one-sample `correct_rt` call
  expect_type(correct_rt_group(Sa_pr, lambdas = 210,
                               groups = c("a", "a", "a", "b"),
                               show_progress = FALSE, cl = 1), "list")
  # with every sample in its own group, each group average *is* that sample, so
  # the result must match ordinary per-sample ptw warping
  singletons <- correct_rt_group(Sa_pr, lambdas = 210, groups = names(Sa_pr),
                                 reference = 1, show_progress = FALSE, cl = 1)
  direct <- correct_rt(Sa_pr, lambdas = 210, alg = "ptw", reference = 1,
                       show_progress = FALSE)
  expect_equal(names(singletons), names(direct))
  for (i in seq_along(direct))
    expect_equal(singletons[[i]], direct[[i]], ignore_attr = TRUE)
})

test_that("the lptw regularization is dimensionless and portable", {
  skip_on_cran()
  # `ptw::warp.time` uses a monomial basis in raw time, so the coefficients have
  # incommensurable units. Each deviation is scaled by the weighted RMS of its
  # monomial, which puts the penalty in the same units as the fit term. The
  # value that halves the applied warp must then be the same regardless of the
  # time resolution of the data or the number of matched peaks.
  fpw <- fit_penalized_warp
  identity_coef <- c(0, 1, 0)
  # a warp with real curvature, so the penalty has something to bite on
  case <- function(scale, npk){
    t0 <- seq(100, 900, length.out = npk) * scale
    tc <- (t0 - mean(t0)) / (diff(range(t0)) / 2)
    list(samp_i = cbind(rt = t0, I = rep(100, npk)),
         ref_i  = cbind(rt = t0 + (3 + 6 * tc + 9 * tc^2) * scale,
                        I = rep(100, npk)))
  }
  displacement <- function(co, samp_rt)
    sqrt(mean((ptw::warp.time(samp_rt, co) - samp_rt)^2))
  half_point <- function(d){
    samp_rt <- d$samp_i[, "rt"]
    d0 <- displacement(fpw(d$ref_i, d$samp_i, penalty = 0), samp_rt)
    lo <- 1e-6; hi <- 1e6
    for (i in seq_len(50)){
      mid <- sqrt(lo * hi)
      if (displacement(fpw(d$ref_i, d$samp_i, penalty = mid), samp_rt) > 0.5 * d0)
        lo <- mid else hi <- mid
    }
    sqrt(lo * hi)
  }
  # invariant to the index scale, and in a human-usable range rather than ~1e12
  by_scale <- vapply(c(0.5, 1, 2), function(s) half_point(case(s, 20)), numeric(1))
  expect_equal(by_scale[1], by_scale[2], tolerance = 0.01)
  expect_equal(by_scale[2], by_scale[3], tolerance = 0.01)
  expect_gt(by_scale[2], 0.1)
  expect_lt(by_scale[2], 10)
  # and invariant to the number of matched peaks
  by_npk <- vapply(c(5, 20, 80), function(n) half_point(case(1, n)), numeric(1))
  expect_equal(by_npk[1], by_npk[2], tolerance = 0.05)
  expect_equal(by_npk[2], by_npk[3], tolerance = 0.05)
})

test_that("the lptw regularization shrinks the warp monotonically", {
  skip_on_cran()
  fpw <- fit_penalized_warp
  npk <- 20
  t0 <- seq(100, 900, length.out = npk)
  tc <- (t0 - mean(t0)) / (diff(range(t0)) / 2)
  samp_i <- cbind(rt = t0, I = rep(100, npk))
  ref_i  <- cbind(rt = t0 + 3 + 6 * tc + 9 * tc^2, I = rep(100, npk))
  disp <- function(p) sqrt(mean((ptw::warp.time(
    t0, fpw(ref_i, samp_i, penalty = p)) - t0)^2))
  d <- vapply(c(0, 0.5, 2.4, 10, 1000), disp, numeric(1))
  expect_true(all(diff(d) < 0))
  # ~2.4 is documented as halving the warp
  expect_equal(d[3] / d[1], 0.5, tolerance = 0.1)
  # a large value approaches the identity warp
  expect_lt(d[5] / d[1], 0.02)
  # and `0` is exactly unpenalized
  expect_equal(fpw(ref_i, samp_i, penalty = 0),
               fpw(ref_i, samp_i))
})

test_that("penalty_scale selects which warp coefficients are penalized", {
  skip_on_cran()
  fpw <- fit_penalized_warp
  npk <- 20
  t0 <- seq(100, 900, length.out = npk)
  tc <- (t0 - mean(t0)) / (diff(range(t0)) / 2)
  samp_i <- cbind(rt = t0, I = rep(100, npk))
  ref_i  <- cbind(rt = t0 + 3 + 6 * tc + 9 * tc^2, I = rep(100, npk))
  free <- fpw(ref_i, samp_i, penalty = 0)
  # penalizing only the quadratic term must shrink it while leaving the
  # intercept and slope comparatively free
  quad_only <- fpw(ref_i, samp_i, penalty = 10, penalty_scale = c(0, 0, 1))
  expect_lt(abs(quad_only[3]), abs(free[3]))
  all_terms <- fpw(ref_i, samp_i, penalty = 10)
  expect_lt(abs(all_terms[1] - 0), abs(quad_only[1] - 0))
  expect_error(fpw(ref_i, samp_i, penalty = 1, penalty_scale = c(1, 1)),
               regexp = "must have the same length as 'init.coef'")
})

test_that("the deprecated `n.zeros` and `n.traces` arguments still work", {
  skip_on_cran()
  data(Sa_pr)
  expect_warning(old <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw",
                                   n.zeros = 50, what = "models"),
                 regexp = "`n.zeros` argument is deprecated")
  new <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw", n_zeros = 50,
                    what = "models")
  expect_equal(attr(old, "n_zeros"), 50)
  expect_equal(lapply(old, function(x) x$warp.coef),
               lapply(new, function(x) x$warp.coef))
})

test_that("shift_warp_coef re-expresses a padded warp in unpadded indices", {
  coef <- c(-6.25, 1.01, -2e-6)
  tp <- 1:500
  expect_equal(ptw::warp.time(tp, shift_warp_coef(coef, 100)),
               ptw::warp.time(tp + 100, coef) - 100)
})

test_that("padded ptw models are stored on the unpadded axis", {
  skip_on_cran()
  data(Sa_pr)
  models <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw", n_zeros = 20,
                       what = "models", show_progress = FALSE)
  m <- models[[2]]
  expect_equal(dim(m$warp.fun), c(1, nrow(Sa_pr[[2]])))
  expect_equal(ncol(m$warped.sample), nrow(Sa_pr[[2]]))
  expect_equal(c(m$warp.fun),
               ptw::warp.time(seq_len(nrow(Sa_pr[[2]])), m$warp.coef[1, ]))
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  expect_s3_class(correct_peaks(pks, models, Sa_pr), "peak_list")
})

test_that("correct_rt records the chrom_list name when rescaling", {
  data(Sa_pr)
  m <- correct_rt(Sa_pr[1:2], lambdas = 210, scale = TRUE, what = "models",
                  cl = 1)
  expect_equal(attr(m, "chrom_list"), "Sa_pr[1:2]")
})

test_that("lptw rejects multiple wavelengths", {
  data(Sa_pr)
  expect_error(correct_rt(Sa_pr, lambdas = c(210, 220), alg = "lptw", cl = 1),
               "single wavelength")
})

test_that("a matched_peak_list is applied by sample name", {
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1)
  mt <- match_peaks(pks, chrom_list = Sa_pr, reference = 1, lambda = 210)
  rev_mt <- structure(rev(unclass(mt)), class = class(mt),
                      reference = attr(mt, "reference"))
  expect_equal(correct_rt(Sa_pr, alg = "lptw", peaks = rev_mt, cl = 1),
               correct_rt(Sa_pr, alg = "lptw", peaks = mt, cl = 1))
  short <- mt
  short[[4]] <- NULL
  expect_error(correct_rt(Sa_pr, alg = "lptw", peaks = short, cl = 1),
               "no entry for sample")
})

test_that("lptw fit warnings survive parallel workers", {
  skip_on_os("windows")
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1)
  mt <- match_peaks(pks, chrom_list = Sa_pr, reference = 1, lambda = 210)
  mt[[2]]$ref_i <- mt[[2]]$ref_i[1, , drop = FALSE]
  mt[[2]]$samp_i <- mt[[2]]$samp_i[1, , drop = FALSE]
  mt[[2]]$cor_scores <- mt[[2]]$cor_scores[1]
  expect_warning(correct_rt(Sa_pr, alg = "lptw", peaks = mt, cl = 2,
                            what = "models"),
                 "Sample 121.*insufficient")
})

test_that("a stored reference can be named by sample", {
  data(Sa_pr)
  m <- correct_rt(Sa_pr, lambdas = 210, reference = 1, what = "models", cl = 1)
  expect_no_error(correct_rt(Sa_pr, models = m, reference = names(Sa_pr)[1],
                             cl = 1))
})

test_that("fit_penalized_warp solves exactly at higher degrees", {
  t0 <- seq(40, 430, length.out = 20)
  samp_i <- cbind(rt = t0, I = rep(100, 20))
  truth <- c(-5, 1.02, -1e-4, 2e-7)
  ref_i <- cbind(rt = ptw::warp.time(t0, truth), I = rep(100, 20))
  expect_equal(fit_penalized_warp(ref_i, samp_i, init.coef = c(0, 1, 0, 0)),
               truth, tolerance = 1e-6)
})

test_that("fit_penalized_warp returns the identity when peaks share times", {
  samp_i <- cbind(rt = c(100, 100, 300), I = rep(100, 3))
  ref_i <- cbind(rt = c(105, 106, 310), I = rep(100, 3))
  expect_warning(co <- fit_penalized_warp(ref_i, samp_i), "do not determine")
  expect_equal(co, c(0, 1, 0))
})

test_that("warp_sample_ptw rejects a warp that reverses direction", {
  data(Sa_pr)
  expect_error(warp_sample_ptw(Sa_pr[[1]], c(0, 1.5, -0.003)),
               "reverses direction")
})

test_that("lptw falls back to the identity warp when a fit reverses direction", {
  skip_on_cran()
  data(Sa_pr)
  local_mocked_bindings(fit_penalized_warp = function(...) c(0, 1.5, -0.003))
  w <- capture_warnings(
    warped <- correct_rt(Sa_pr, lambdas = "210", alg = "lptw", reference = 1,
                         show_progress = FALSE, cl = 1))
  expect_match(w, "The fitted warp reverses direction")
  expect_setequal(sub(":.*", "", w), paste("Sample", names(Sa_pr)))
  expect_equal(unname(warped[[2]]), unname(Sa_pr[[2]]))
})

test_that("correct_rt applies VPdtw's own warp to every sample", {
  data(Sa_pr)
  mods <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", reference = 1,
                     what = "models", show_progress = FALSE, cl = 1)
  warped <- correct_rt(Sa_pr, models = mods, alg = "vpdtw", lambdas = "210",
                       reference = 1, show_progress = FALSE, cl = 1)
  for (s in seq_along(Sa_pr)){
    expect_equal(unname(warped[[s]][, "210"]), unname(mods$warpedQuery[, s]))
  }
  rows <- match(format_axis(get_times(Sa_pr, idx = 1)), rownames(warped[[1]]))
  expect_false(anyNA(rows))
  expect_equal(c(warped[[1]][rows, ]), c(Sa_pr[[1]]))
})

test_that("correct_rt raises a worker's error when run in parallel", {
  skip_on_cran()
  data(Sa_pr)
  mods <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw", reference = 1,
                     what = "models", show_progress = FALSE, cl = 1)
  mods[[2]]$warp.coef[1, ] <- c(0, 1, -0.01)
  expect_error(suppressWarnings(
    correct_rt(Sa_pr, models = mods, alg = "ptw", show_progress = FALSE,
               cl = 2)),
    "reverses direction")
})

test_that("fitted models are applied to samples by name", {
  skip_on_cran()
  data(Sa_pr)
  for (alg in c("ptw", "lptw", "vpdtw")){
    mods <- suppressWarnings(correct_rt(Sa_pr, lambdas = "210", alg = alg,
                                        reference = 2, what = "models", cl = 1,
                                        show_progress = FALSE))
    fwd <- correct_rt(Sa_pr, models = mods, alg = alg, lambdas = "210", cl = 1,
                      show_progress = FALSE)
    rev_out <- correct_rt(rev(Sa_pr), models = mods, alg = alg,
                          lambdas = "210", cl = 1, show_progress = FALSE)
    for (s in names(Sa_pr))
      expect_equal(unclass(rev_out[[s]]), unclass(fwd[[s]]),
                   ignore_attr = TRUE, info = paste(alg, s))
    sub <- correct_rt(Sa_pr[3:4], models = mods, alg = alg, lambdas = "210",
                      cl = 1, show_progress = FALSE)
    expect_equal(unclass(sub[[names(Sa_pr)[4]]]), unclass(fwd[[4]]),
                 ignore_attr = TRUE, info = alg)
  }
  mods <- correct_rt(Sa_pr[1:3], lambdas = "210", reference = 2,
                     what = "models", cl = 1, show_progress = FALSE)
  expect_error(correct_rt(Sa_pr, models = mods, lambdas = "210", cl = 1,
                          show_progress = FALSE),
               regexp = "no warp for sample")
})

test_that("warping rejects duplicate sample names", {
  skip_on_cran()
  data(Sa_pr)
  dn <- Sa_pr
  names(dn) <- c("a", "a", "b", "c")
  expect_error(correct_rt(dn, lambdas = "210", reference = 3, cl = 1,
                          show_progress = FALSE),
               regexp = "Duplicate sample names found in `chrom_list`")
  expect_error(correct_rt_group(dn, lambdas = "210",
                                groups = c("x", "x", "y", "y"), cl = 1,
                                show_progress = FALSE),
               regexp = "Duplicate sample names found in `chrom_list`")
})

test_that("fitted VPdtw models re-apply to multi-wavelength data without `lambdas`", {
  skip_on_cran()
  data(Sa_pr)
  m <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", what = "models",
                  show_progress = FALSE)
  warped <- correct_rt(Sa_pr, alg = "vpdtw", models = m, show_progress = FALSE)
  expect_equal(colnames(warped[[1]]), colnames(Sa_pr[[1]]))
})
