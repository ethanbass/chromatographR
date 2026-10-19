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

test_that("match_peaks picks a reference containing every target", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)
  mp <- match_peaks(pks, chrom_list = Sa_pr, targets = tg)
  ref <- attr(mp, "reference")
  expect_true(ref %in% names(pks))
  rt_ref <- attr(mp, "reference_targets")$rt
  expect_equal(attr(mp, "reference_targets")$name, c("a", "b"))
  expect_true(all(abs(rt_ref - tg$rt) <= tg$tol))
  expect_true(all(rt_ref %in% pks[[ref]][["210"]]$rt))

  only_458 <- pks
  for (s in setdiff(names(pks), "458")) {
    p <- only_458[[s]][["210"]]
    only_458[[s]][["210"]] <- p[abs(p$rt - 12.2) > 0.2, ]
  }
  mp <- match_peaks(only_458, chrom_list = Sa_pr, targets = tg)
  expect_equal(attr(mp, "reference"), "458")
})

test_that("match_peaks keeps a given reference alongside `targets`", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)
  mp <- match_peaks(pks, chrom_list = Sa_pr, reference = "121", targets = tg)
  expect_equal(attr(mp, "reference"), "121")
  expect_true(all(attr(mp, "reference_targets")$rt %in% pks[["121"]][["210"]]$rt))
  expect_equal(match_peaks(pks, chrom_list = Sa_pr, targets = c(12.2, 14.45)),
               match_peaks(pks, chrom_list = Sa_pr,
                           targets = data.frame(rt = c(12.2, 14.45))))
  no_a <- pks
  p <- no_a[["121"]][["210"]]
  no_a[["121"]][["210"]] <- p[abs(p$rt - 12.2) > 0.2, ]
  expect_error(match_peaks(no_a, chrom_list = Sa_pr, reference = "121", targets = tg),
               "reference sample .*121.* has no peak .* of a\\.")
  expect_error(match_peaks(pks, chrom_list = Sa_pr, reference = "best", targets = tg),
               "cannot be combined with `targets`")
})

test_that("match_peaks reports target coverage when no sample has every target", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "absent"), rt = c(12.2, 19.5), tol = 0.2)
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = tg),
               "No sample contains a peak for every target.*a 4 of 4, absent 0 of 4")
  tg_260 <- data.frame(rt = 12.2, lambda = 260)
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = tg_260),
               "No peaks at 260 nm: `peak_list` has peaks at 210 only")
})

test_that("match_peaks compares corrected retention times when present", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  ref <- match_peaks(pks, chrom_list = Sa_pr, reference = 1)[[2]]$matches

  late <- Sa_pr
  late[[2]] <- transfer_metadata(late[[2]], Sa_pr[[2]])
  rownames(late[[2]]) <- as.numeric(rownames(Sa_pr[[2]])) + 0.5
  shifted <- pks
  shifted[[2]][["210"]][c("rt", "start", "end")] <-
    pks[[2]][["210"]][c("rt", "start", "end")] + 0.5
  corrected <- shifted
  for (i in seq_along(corrected))
    corrected[[i]][["210"]]$rt.cor <- pks[[i]][["210"]]$rt

  n_raw <- NROW(suppressWarnings(
    match_peaks(shifted, chrom_list = late, reference = 1))[[2]]$matches)
  expect_lt(n_raw, nrow(ref))

  m <- match_peaks(corrected, chrom_list = late, reference = 1)[[2]]$matches
  expect_equal(m[, "rt_ref"], ref[, "rt_ref"])
  expect_equal(m[, "rt_samp"], ref[, "rt_samp"] + 0.5)

  n_off <- NROW(suppressWarnings(
    match_peaks(corrected, chrom_list = late, reference = 1,
                use_corrected_rt = FALSE))[[2]]$matches)
  expect_equal(n_off, n_raw)
  expect_error(match_peaks(shifted, chrom_list = late, reference = 1,
                           use_corrected_rt = TRUE), "requires corrected")

  corrected[[3]][["210"]]$rt.cor <- NULL
  expect_error(match_peaks(corrected, chrom_list = late, reference = 1),
               "others do not")
})

test_that("targets are looked up in corrected time and reported in raw time", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)
  unshifted <- match_peaks(pks, chrom_list = Sa_pr, targets = tg)

  late <- Sa_pr
  corrected <- pks
  for (i in seq_along(pks)){
    rownames(late[[i]]) <- as.numeric(rownames(Sa_pr[[i]])) + 0.5
    corrected[[i]][["210"]][c("rt", "start", "end")] <-
      pks[[i]][["210"]][c("rt", "start", "end")] + 0.5
    corrected[[i]][["210"]]$rt.cor <- pks[[i]][["210"]]$rt
  }
  mp <- match_peaks(corrected, chrom_list = late, targets = tg)
  expect_equal(attr(mp, "reference"), attr(unshifted, "reference"))
  expect_equal(attr(mp, "reference_targets")$rt,
               attr(unshifted, "reference_targets")$rt + 0.5)
})

test_that("match_peaks ranks candidate references and applies the presence floors", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)

  boosted <- pks
  boosted[["122"]][["210"]]$height <- 3 * boosted[["122"]][["210"]]$height
  expect_equal(attr(match_peaks(boosted, chrom_list = Sa_pr, targets = tg), "reference"), "122")

  near_a <- abs(boosted[["122"]][["210"]]$rt - 12.2) <= 0.2
  poor_fit <- boosted
  poor_fit[["122"]][["210"]]$r.squared[near_a] <- 0.3
  expect_false(attr(match_peaks(poor_fit, chrom_list = Sa_pr, targets = tg), "reference") == "122")

  tiny <- boosted
  tiny[["122"]][["210"]]$height[near_a] <- 0.01 * tiny[["122"]][["210"]]$height[near_a]
  expect_false(attr(match_peaks(tiny, chrom_list = Sa_pr, targets = tg), "reference") == "122")

  noisy <- boosted
  p <- noisy[["122"]][["210"]]
  blip <- p[which(near_a)[1], ]
  blip$rt <- 12.2
  blip$height <- 10 * max(p$height)
  blip$r.squared <- 0.1
  noisy[["122"]][["210"]] <- rbind(p, blip)
  expect_equal(attr(match_peaks(noisy, chrom_list = Sa_pr, targets = tg), "reference"), "122")
})

test_that("a peak within tol of two targets counts only for the nearer one", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE,
                   smooth_window = 1)
  tg <- data.frame(name = c("x", "y"), rt = c(11.36, 11.39), tol = 0.1)
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = tg), "y 0 of 4")
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = data.frame(rt = 11.37), lambda = 254),
               "No peaks at 254 nm")
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = data.frame(time = 11.37)),
               "supplied to `targets`")
})

test_that("a matched peak table names the target columns", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)
  mp <- match_peaks(pks, chrom_list = Sa_pr, targets = tg)
  tab <- get_peaktable(mp, pks = pks)
  expect_true(all(c("a", "b") %in% colnames(tab$tab)))
  expect_identical(colnames(tab$tab), colnames(tab$pk_meta))
  expect_equal(unlist(tab$pk_meta["rt", c("a", "b")]), attr(mp, "reference_targets")$rt,
               ignore_attr = TRUE)
  expect_equal(colnames(tab[, c("a", "b")]$tab), c("a", "b"))
  expect_true(all(grepl("^V[0-9]+$", setdiff(colnames(tab$tab), c("a", "b")))))

  plain <- get_peaktable(match_peaks(pks, chrom_list = Sa_pr, reference = 1), pks = pks)
  expect_true(all(grepl("^V[0-9]+$", colnames(plain$tab))))

  ref <- attr(mp, "reference")
  missing_a <- pks
  p <- missing_a[[ref]][["210"]]
  missing_a[[ref]][["210"]] <- p[abs(p$rt - attr(mp, "reference_targets")$rt[1]) > 1e-6, ]
  expect_error(get_peaktable(mp, pks = missing_a), "No reference peak found for target 'a'")
})

test_that("a matched peak table on corrected peaks looks up aligned chromatograms", {
  skip_on_cran()
  data(Sa_pr)
  mods <- correct_rt(Sa_pr, lambdas = "210", reference = 1, what = "models",
                     show_progress = FALSE)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  pks_cor <- correct_peaks(pks, mods, chrom_list = Sa_pr)
  pt <- get_peaktable(match_peaks(pks_cor, chrom_list = Sa_pr, reference = 1),
                      pks = pks_cor)
  expect_true(pt$args$use_corrected_rt)
  ref_pks <- pks_cor[[1]][["210"]]
  expect_equal(unlist(pt$pk_meta["rt", ]),
               ref_pks$rt.cor[unlist(pt$pk_meta["peak", ])], ignore_attr = TRUE)
  warped <- correct_rt(Sa_pr, models = mods, show_progress = FALSE)
  expect_false(isTRUE(all.equal(warped[[2]], Sa_pr[[2]], check.attributes = FALSE)))
  expect_equal(get_chrom_list(pt), warped, ignore_attr = TRUE)

  pks_cor <- correct_peaks(pks, mods, chrom_list = Sa_pr, chrom_list_aligned = warped)
  pt <- get_peaktable(match_peaks(pks_cor, chrom_list = Sa_pr, reference = 1),
                      pks = pks_cor)
  expect_equal(pt$args$chrom_list_aligned, "warped")
})

test_that("background-subtracted spectra ignore a sloping background", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  sloped <- Sa_pr
  t <- get_times(Sa_pr, 2)
  lam <- as.numeric(colnames(Sa_pr[[2]]))
  sloped[[2]] <- transfer_metadata(
    Sa_pr[[2]] + outer(2 * (t - min(t)) / diff(range(t)), (lam - min(lam)) / diff(range(lam))),
    Sa_pr[[2]])
  n <- function(chroms, background_correction)
    NROW(suppressWarnings(match_peaks(pks, chrom_list = chroms, reference = 1,
                                      background_correction = background_correction))[[2]]$matches)
  expect_lt(n(sloped, "none"), n(Sa_pr, "none"))
  clean <- match_peaks(pks, chrom_list = Sa_pr, reference = 1,
                       background_correction = "linear")
  expect_equal(match_peaks(pks, chrom_list = sloped, reference = 1,
                           background_correction = "linear")[[2]]$matches,
               clean[[2]]$matches)
  expect_equal(attr(clean, "background_correction"), "linear")

  no_bounds <- pks
  for (i in seq_along(no_bounds)) no_bounds[[i]][["210"]][c("start", "end")] <- NULL
  expect_error(match_peaks(no_bounds, chrom_list = Sa_pr, reference = 1,
                           background_correction = "linear"), "`start` and `end`")
})

test_that("background spectra fall back to retention time on one-dimensional chromatograms", {
  skip_on_cran()
  data(Sa_pr)
  x1 <- structure(lapply(Sa_pr, function(m) transfer_metadata(m[, "210", drop = FALSE], m)),
                  names = names(Sa_pr))
  pks <- get_peaks(x1, lambdas = 210, cl = 1, show_progress = FALSE)
  expect_warning(bg <- match_peaks(pks, chrom_list = x1, reference = 1,
                                   background_correction = "linear"),
                 "Spectral matching is unavailable")
  apex <- suppressWarnings(match_peaks(pks, chrom_list = x1, reference = 1))
  expect_equal(bg[[2]]$matches, apex[[2]]$matches)
})

test_that("a targets table can set a spectral cutoff per target", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b"), rt = c(12.2, 14.45), tol = 0.2)
  base <- match_peaks(pks, chrom_list = Sa_pr, targets = tg)
  rt_a <- attr(base, "reference_targets")$rt[1]
  n_a <- function(m) sum(vapply(m[-reference_position(m)], function(s)
    sum(s$matches[, "rt_ref"] == rt_a), numeric(1)))
  expect_gt(n_a(base), 0)

  strict <- match_peaks(pks, chrom_list = Sa_pr,
                        targets = transform(tg, spectral_cutoff = c(1, NA)))
  expect_equal(n_a(strict), 0)
  other <- function(m) lapply(m, function(s) s$matches[s$matches[, "rt_ref"] != rt_a, , drop = FALSE])
  expect_equal(other(strict), other(base))

  expect_error(match_peaks(pks, chrom_list = Sa_pr,
                           targets = transform(tg, spectral_cutoff = c(1.5, NA))),
               "between -1 and 1")
})

test_that("a targets table can set a retention time window per target", {
  times <- seq(0, 5, by = 0.005)
  run <- function(rt, h){
    y <- h * exp(-(times - rt)^2 / (2 * 0.02^2))
    matrix(c(y, 0.5 * y, 0.2 * y), ncol = 3, dimnames = list(times, c("210", "220", "230")))
  }
  x <- list(r1 = run(2.5, 100), r2 = run(2.53, 60))
  pks <- get_peaks(x, lambdas = "210", fit = "gaussian", show_progress = FALSE, cl = 1)

  wide <- match_peaks(pks, chrom_list = x, rt_window = 0.01,
                      targets = data.frame(name = "a", rt = 2.5, rt_window = 0.05))
  expect_equal(attr(wide, "reference"), "r1")
  expect_equal(nrow(wide[["r2"]]$matches), 1)

  expect_warning(narrow <- match_peaks(pks, chrom_list = x,
                                       targets = data.frame(name = "a", rt = 2.5, rt_window = 0.02)),
                 "No peaks matched for sample 2")
  expect_null(narrow[["r2"]])
  expect_warning(match_peaks(pks, chrom_list = x, rt_window = 0.01,
                             targets = data.frame(name = "a", rt = 2.5, tol = 0.05)),
                 "No peaks matched for sample 2")
  expect_equal(match_peaks(pks, chrom_list = x,
                           targets = data.frame(name = "a", rt = 2.5, rt_window = NA)),
               match_peaks(pks, chrom_list = x, targets = data.frame(name = "a", rt = 2.5)))
  expect_error(match_peaks(pks, chrom_list = x,
                           targets = data.frame(name = "a", rt = 2.5, rt_window = -1)),
               "must be numeric and non-negative")
})

test_that("a targets table can set the wavelength per target", {
  times <- seq(0, 5, by = 0.005)
  run <- function(shift, ha, hb){
    a <- ha * exp(-(times - 2.5 - shift)^2 / (2 * 0.02^2))
    b <- hb * exp(-(times - 3.5 - shift)^2 / (2 * 0.02^2))
    matrix(c(a + 0.1 * b, 0.5 * a + 0.4 * b, 0.2 * a + b), ncol = 3,
           dimnames = list(times, c("210", "220", "230")))
  }
  x <- list(r1 = run(0, 100, 80), r2 = run(0.01, 60, 50), r3 = run(-0.01, 90, 40))
  pks <- get_peaks(x, lambdas = c("210", "230"), fit = "gaussian", show_progress = FALSE, cl = 1)
  targets <- data.frame(name = c("a", "b"), rt = c(2.5, 3.5), lambda = c(210, 230), tol = 0.05)

  m <- match_peaks(pks, chrom_list = x, targets = targets)
  expect_s3_class(m, "matched_peak_list_set")
  expect_named(m, c("210", "230"))
  pt <- get_peaktable(m, pks = pks)
  expect_equal(unlist(pt$pk_meta["lambda", c("a", "b")], use.names = FALSE), c(210, 230))

  single_a <- get_peaktable(match_peaks(pks, chrom_list = x, targets = targets[1, ]), pks = pks)
  single_b <- get_peaktable(match_peaks(pks, chrom_list = x, targets = targets[2, ]), pks = pks)
  expect_equal(pt$tab$a, single_a$tab$a)
  expect_equal(pt$tab$b, single_b$tab$b)
  expect_true(all(pt$tab$a > 0) && all(pt$tab$b > 0))
  expect_false(anyDuplicated(colnames(pt$tab)) > 0)

  expect_error(match_peaks(pks, chrom_list = x,
                           targets = transform(targets, lambda = c(210, NA))),
               "Give every target a wavelength")
  expect_error(correct_rt(x, alg = "lptw", peaks = m), "takes the matches at one wavelength")
})

test_that("a spectral_cutoff column of NA uses the default cutoff", {
  times <- seq(0, 5, by = 0.005)
  run <- function(rt, h){
    y <- h * exp(-(times - rt)^2 / (2 * 0.02^2))
    matrix(c(y, 0.5 * y, 0.2 * y), ncol = 3, dimnames = list(times, c("210", "220", "230")))
  }
  x <- list(r1 = run(2.5, 100), r2 = run(2.51, 60))
  pks <- get_peaks(x, lambdas = "210", fit = "gaussian", show_progress = FALSE, cl = 1)
  tg <- data.frame(name = "a", rt = 2.5, tol = 0.05)
  expect_equal(match_peaks(pks, chrom_list = x, targets = transform(tg, spectral_cutoff = NA)),
               match_peaks(pks, chrom_list = x, targets = tg), ignore_attr = "reference_targets")
  expect_error(match_peaks(pks, chrom_list = x, targets = transform(tg, spectral_cutoff = "high")),
               "must be numeric")
})

test_that("a consensus spectrum replaces an overlapped reference spectrum", {
  times <- seq(0, 5, by = 0.005)
  wl <- seq(200, 300, by = 5)
  target <- exp(-(wl - 230)^2 / (2 * 15^2))
  interferent <- exp(-(wl - 280)^2 / (2 * 10^2))
  run <- function(h, contaminate){
    y <- h * exp(-(times - 2.5)^2 / (2 * 0.02^2))
    x <- outer(y, target) + if (contaminate) outer(0.8 * y, interferent) else 0
    withr::with_seed(h, x <- x + matrix(rnorm(length(x), sd = 0.2), nrow(x)))
    dimnames(x) <- list(times, wl)
    x
  }
  x <- list(r1 = run(100, TRUE), r2 = run(90, FALSE), r3 = run(80, FALSE), r4 = run(70, FALSE))
  pks <- get_peaks(x, lambdas = "230", fit = "gaussian", show_progress = FALSE, cl = 1)
  tg <- data.frame(name = "a", rt = 2.5, tol = 0.05)
  own <- match_peaks(pks, chrom_list = x, reference = "r1", spectral_cutoff = -1)
  con <- match_peaks(pks, chrom_list = x, reference = "r1", spectral_cutoff = -1, consensus_spectra = TRUE)
  cor_at <- function(m, s) m[[s]]$cor_scores[which.min(abs(m[[s]]$matches[, "rt_ref"] - 2.5))]
  k <- which.min(abs(attr(con, "rt_ref") - 2.5))
  expect_null(attr(own, "consensus_spectra"))
  expect_true(all(vapply(names(x)[-1], function(s) cor_at(con, s) > cor_at(own, s), logical(1))))
  own_spectrum <- x$r1[which.min(abs(times - 2.5)), ]
  expect_gt(cor(attr(con, "consensus_spectra")[, k], target), cor(own_spectrum, target) + 0.2)
  expect_equal(attr(con, "consensus_n")[k], 4)
  expect_equal(match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = -1, consensus_spectra = FALSE),
               match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = -1))
})

noise_test_runs <- function(heights, spectra, seed){
  times <- seq(0, 4, by = 0.01)
  wl <- seq(200, 300, by = 5)
  sd_wl <- 1 + 2 * exp(-(wl - 200) / 20)
  withr::with_seed(seed, lapply(seq_along(heights), function(i){
    y <- heights[i] * exp(-(times - 2)^2 / (2 * 0.03^2))
    x <- outer(y, spectra[[i]]) + sweep(matrix(rnorm(length(times) * length(wl)), length(times)), 2, sd_wl, "*") +
      rnorm(length(times), sd = 0.5)
    dimnames(x) <- list(times, wl)
    x
  }))
}
noise_test_spectra <- function(){
  wl <- seq(200, 300, by = 5)
  a <- exp(-(wl - 230)^2 / (2 * 12^2))
  list(a = a, b = a + 0.45 * exp(-(wl - 265)^2 / (2 * 10^2)))
}
matched_at <- function(m, s) {
  x <- m[[s]]
  if (is.null(x)) return(NA_real_)
  k <- which(abs(x$matches[, "rt_ref"] - 2) < 0.05)
  if (length(k)) x$matches[k[1], "rt_samp"] else NA_real_
}

test_that("the noise spectral test passes small peaks and rejects wrong compounds", {
  sp <- noise_test_spectra()
  expect_lt(cor(sp$a, sp$b), 0.95)
  expect_gt(cor(sp$a, sp$b), 0.8)
  x <- noise_test_runs(c(1000, 6, 500), list(sp$a, sp$a, sp$b), seed = 1)
  names(x) <- c("ref", "small", "wrong")
  pks <- get_peaks(x, lambdas = "230", fit = "gaussian", show_progress = FALSE, cl = 1)
  tg <- data.frame(name = "a", rt = 2, tol = 0.05)
  cut <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = 0.95))
  noi <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = 0.95,
                                      spectral_test = "noise", rt_window = 0.05))
  expect_true(is.na(matched_at(cut, "small")))
  expect_false(is.na(matched_at(noi, "small")))
  expect_true(is.na(matched_at(noi, "wrong")))
  expect_equal(attr(noi, "spectral_test"), "noise")
  expect_true(all(noi[["small"]]$spectral_p >= 0.01))

  loose <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = transform(tg, spectral_cutoff = 0.8),
                                        spectral_test = "noise"))
  expect_false(is.na(matched_at(loose, "wrong")))
  off <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = transform(tg, spectral_cutoff = -1),
                                      spectral_test = "noise"))
  expect_false(is.na(matched_at(off, "wrong")))

  expect_error(match_peaks(pks, chrom_list = x, targets = tg, spectral_test = "noise", alpha = 2),
               "`alpha` must be")
  expect_error(match_peaks(pks, chrom_list = x, targets = tg, spectral_test = "noise",
                           background_correction = "linear"),
               "needs `background_correction = \"none\"`")
  expect_equal(suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg, spectral_test = "cutoff")),
               suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg)))
})

test_that("the noise spectral test rejects real peaks at about the nominal rate", {
  skip_on_cran()
  sp <- noise_test_spectra()
  n <- 60
  x <- noise_test_runs(c(1000, rep(5, n)), rep(list(sp$a), n + 1), seed = 2)
  names(x) <- c("ref", paste0("s", seq_len(n)))
  pks <- get_peaks(x, lambdas = "230", fit = "gaussian", show_progress = FALSE, cl = 1)
  tg <- data.frame(name = "a", rt = 2, tol = 0.05)
  free <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = -1))
  noi <- suppressWarnings(match_peaks(pks, chrom_list = x, targets = tg, spectral_cutoff = 0.95,
                                      spectral_test = "noise", alpha = 0.05, rt_window = 0.05))
  found <- names(x)[-1][!is.na(vapply(names(x)[-1], matched_at, numeric(1), m = free))]
  expect_gt(length(found), 40)
  rejected <- mean(is.na(vapply(found, matched_at, numeric(1), m = noi)))
  expect_lt(rejected, 0.15)
})

test_that("trace noise is measured, not inferred from point-to-point differences, when scans are correlated", {
  set.seed(1)
  n <- 3000
  times <- seq(0.01, by = 0.01, length.out = n)
  kern <- dnorm(-8:8, sd = 2)
  noise <- as.numeric(stats::filter(rnorm(n + 16), kern / sum(kern), sides = 2))[9:(n + 8)]
  noise <- noise / sd(noise)
  rt <- seq(1, 29, by = 0.7)
  height <- exp(seq(log(2), log(500), length.out = length(rt)))
  width <- 0.03
  y <- noise + rowSums(sapply(seq_along(rt), function(i)
    height[i] * exp(-0.5 * ((times - rt[i]) / width)^2)))
  peaks <- data.frame(start = rt - 4 * width, end = rt + 4 * width, height = height)
  expect_lt(diff_noise(y), 0.5)
  expect_equal(trace_noise(y, times, peaks), 1, tolerance = 0.15)
})

test_that("trace noise is not inflated by a crowded trace", {
  set.seed(1)
  n <- 600
  times <- seq(1 / 150, by = 1 / 150, length.out = n)
  noise <- rnorm(n)
  rt <- sort(runif(60, 0, 4))
  height <- exp(runif(60, log(0.5), log(1000)))
  width <- 5 / 150 / 2.355
  y <- noise + rowSums(sapply(seq_along(rt), function(i)
    height[i] * exp(-0.5 * ((times - rt[i]) / width)^2)))
  peaks <- data.frame(start = rt - 4 * width, end = rt + 4 * width, height = height)
  expect_gt(diff_noise(y) / sd(noise), 2)
  expect_equal(trace_noise(y, times, peaks) / sd(noise), 1, tolerance = 0.3)
})

test_that("trace noise falls back to the difference estimate when peaks cover the trace", {
  set.seed(1)
  times <- seq(0.01, by = 0.01, length.out = 500)
  y <- rnorm(500)
  peaks <- data.frame(start = 0, end = 10, height = 1e6)
  expect_equal(trace_noise(y, times, peaks), diff_noise(y))
})

test_that("match_peaks rejects aligned chromatograms for peaks found in unaligned ones", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  warped <- correct_rt(Sa_pr, lambdas = "210", verbose = FALSE)
  expect_error(match_peaks(pks, chrom_list = warped, lambda = 210),
               "unaligned chromatograms")
  pks_cor <- correct_peaks(pks, correct_rt(Sa_pr, lambdas = "210",
                                           what = "models", verbose = FALSE))
  attr(pks_cor, "chrom_list") <- "not_a_real_object"
  expect_error(match_peaks(pks_cor, chrom_list = warped, lambda = 210),
               "unaligned chromatograms")
  expect_no_error(suppressWarnings(match_peaks(pks, chrom_list = Sa_pr, lambda = 210)))
})

test_that("an all-NA `spectral_cutoff` column of `targets` takes the argument", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(rt = c(12.07, 13.65, 16.72), name = c("a", "b", "c"))
  ref <- match_peaks(pks, Sa_pr, targets = tg, lambda = 210)
  tg$spectral_cutoff <- NA
  expect_equal(match_peaks(pks, Sa_pr, targets = tg, lambda = 210), ref)
})

test_that("targets at several wavelengths share one reference sample when they can", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = c(210, 254), cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("a", "b", "c"), rt = c(12.07, 13.65, 16.72),
                   lambda = c(210, 210, 254))
  m <- match_peaks(pks, chrom_list = Sa_pr, targets = tg)
  refs <- vapply(m, attr, character(1), "reference")
  expect_length(unique(refs), 1)

  near <- function(p, r) which.min(abs(p$rt - r))
  split_pks <- pks
  a1 <- split_pks[[1]][["210"]]
  split_pks[[1]][["210"]] <- a1[-near(a1, 12.07), ]
  for (s in 2:4) {
    c254 <- split_pks[[s]][["254"]]
    split_pks[[s]][["254"]] <- c254[-near(c254, 16.72), ]
  }
  expect_warning(m2 <- match_peaks(split_pks, chrom_list = Sa_pr, targets = tg),
                 "each wavelength has its own reference sample")
  expect_equal(attr(m2[["254"]], "reference"), names(Sa_pr)[1])
  expect_false(attr(m2[["210"]], "reference") == names(Sa_pr)[1])
})

test_that("matched peak lists and sets print a short summary", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = c(210, 254), cl = 1, show_progress = FALSE)
  m <- match_peaks(pks, chrom_list = Sa_pr, lambda = 210, reference = 1)
  out <- capture.output(print(m))
  expect_lte(length(out), 5)
  expect_match(out[1], "4 samples, matched at 210 nm against reference sample 119",
               fixed = TRUE)
  expect_match(out[3], "^Matched peaks per sample: 121 [0-9]+, 122 [0-9]+, 458 [0-9]+$")
  expect_output(expect_invisible(print(m)))

  tg <- data.frame(name = c("a", "b", "c"), rt = c(12.07, 13.65, 16.72),
                   lambda = c(210, 210, 254))
  out <- capture.output(print(match_peaks(pks, chrom_list = Sa_pr, targets = tg)))
  expect_equal(length(out), 4)
  expect_match(out[1], "at 2 wavelengths, reference sample 119", fixed = TRUE)
  expect_match(out[2], "210 nm: [0-9]+ reference peaks, targets a, b")
  expect_match(out[3], "254 nm: [0-9]+ reference peaks, target c")
})

test_that("lambda fills in targets without a wavelength", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = c(210, 254), cl = 1, show_progress = FALSE)
  m <- match_peaks(pks, chrom_list = Sa_pr, lambda = 254,
                   targets = data.frame(name = c("a", "c"), rt = c(12.07, 16.72),
                                        lambda = c(210, NA)))
  expect_s3_class(m, "matched_peak_list_set")
  expect_equal(names(m), c("210", "254"))
  m1 <- match_peaks(pks, chrom_list = Sa_pr, lambda = 210,
                    targets = data.frame(name = "a", rt = 12.07, lambda = NA))
  expect_s3_class(m1, "matched_peak_list")
})

test_that("targets at several wavelengths are checked as one table", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = c(210, 254), cl = 1, show_progress = FALSE)
  expect_error(match_peaks(pks, chrom_list = Sa_pr,
                           targets = data.frame(name = c("a", "b", "a"),
                                                rt = c(12.07, 13.65, 16.72),
                                                lambda = c(210, 210, 254))),
               "Target names must be unique. Duplicated: a")
  m <- match_peaks(pks, chrom_list = Sa_pr,
                   targets = data.frame(rt = c(12.07, 16.72),
                                        lambda = c("210.0", "254")))
  expect_equal(names(m), c("210", "254"))
  expect_error(expect_no_warning(
    match_peaks(pks, chrom_list = Sa_pr,
                targets = data.frame(rt = c(12.07, 16.72), lambda = c(210, 280)))),
    "No peaks at 280 nm")
})

test_that("the noise test keeps its error rate when sample and reference are equally tall", {
  withr::local_seed(2)
  wl <- seq(200, 300, by = 5)
  sd_wl <- 1 + 2 * exp(-(wl - 200) / 20)
  u <- exp(-(wl - 230)^2 / (2 * 12^2))
  noise <- function(n) sweep(matrix(rnorm(n * length(wl)), n), 2, sd_wl, "*")
  w <- 1 / sd_wl
  p <- replicate(400, {
    Ys <- noise(300); Ys <- Ys - rowMeans(Ys)
    Yr <- noise(300); Yr <- Yr - rowMeans(Yr)
    r <- 50 * u + noise(1)[1, ]; s <- 50 * u + noise(1)[1, ]
    spectral_p(s, spectral_null(r, Ys, Yr, w), w, 1)
  })
  expect_lt(mean(p < 0.05), 0.08)
})

test_that("the noise test needs at least 1 / alpha peak-free timepoints", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  expect_error(match_peaks(pks, chrom_list = Sa_pr, reference = 1,
                           spectral_test = "noise", alpha = 0.001),
               "fewer than 1000")
  null <- list(qx = qr(cbind(1, 1:5)), rc2 = 1, q_samp = rep(0, 9),
               q_ref = rep(0, 9), q_cross = rep(0, 9))
  expect_equal(spectral_p(c(0, 3, 1, 5, 2), null, rep(1, 5), 1), 1 / 10)
})

test_that("spectral_nulls matches a direct projection of the noise spectra", {
  withr::local_seed(3)
  p <- 15
  w <- runif(p, 0.5, 2)
  Ys <- matrix(rnorm(40 * p), 40)
  Yr <- matrix(rnorm(31 * p), 31)
  R <- cbind(exp(-(1:p - 6)^2 / 8), exp(-(1:p - 10)^2 / 4), NA)
  direct <- function(r, scale){
    qx <- qr(cbind(w, r * w))
    rs <- qr.resid(qx, t(Ys) * w)
    rr <- qr.resid(qx, t(Yr) * w) * scale
    n <- max(ncol(rs), ncol(rr))
    i_s <- rep_len(seq_len(ncol(rs)), n)
    i_r <- (seq_len(n) - 1 + ncol(rr) %/% 2) %% ncol(rr) + 1
    c(colSums(rs^2)[i_s], colSums(rr^2)[i_r],
      colSums(rs[, i_s] * rr[, i_r]))
  }
  out <- spectral_nulls(R, Ys, Yr, w, scale = c(1, 0.3, 1))
  for (k in 1:2){
    expect_equal(c(out[[k]]$q_samp, out[[k]]$q_ref, out[[k]]$q_cross),
                 direct(R[, k], c(1, 0.3)[k]), tolerance = 1e-10)
  }
  expect_null(out[[3]])
})
