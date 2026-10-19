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
               "matching wavelength \\(210\\)")
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
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(name = c("x", "y"), rt = c(11.36, 11.39), tol = 0.1)
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = tg), "y 0 of 4")
  expect_error(match_peaks(pks, chrom_list = Sa_pr, targets = data.frame(rt = 11.37), lambda = 254),
               "No peaks were found at wavelength")
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
  n <- function(chroms, spectrum)
    NROW(suppressWarnings(match_peaks(pks, chrom_list = chroms, reference = 1,
                                      spectrum = spectrum))[[2]]$matches)
  expect_lt(n(sloped, "apex"), n(Sa_pr, "apex"))
  clean <- match_peaks(pks, chrom_list = Sa_pr, reference = 1, spectrum = "background")
  expect_equal(match_peaks(pks, chrom_list = sloped, reference = 1,
                           spectrum = "background")[[2]]$matches,
               clean[[2]]$matches)
  expect_equal(attr(clean, "spectrum"), "background")

  no_bounds <- pks
  for (i in seq_along(no_bounds)) no_bounds[[i]][["210"]][c("start", "end")] <- NULL
  expect_error(match_peaks(no_bounds, chrom_list = Sa_pr, reference = 1,
                           spectrum = "background"), "`start` and `end`")
})

test_that("background spectra fall back to retention time on one-dimensional chromatograms", {
  skip_on_cran()
  data(Sa_pr)
  x1 <- structure(lapply(Sa_pr, function(m) transfer_metadata(m[, "210", drop = FALSE], m)),
                  names = names(Sa_pr))
  pks <- get_peaks(x1, lambdas = 210, cl = 1, show_progress = FALSE)
  expect_warning(bg <- match_peaks(pks, chrom_list = x1, reference = 1,
                                   spectrum = "background"),
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

test_that("an all-NA `spectral_cutoff` column of `targets` takes the argument", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  tg <- data.frame(rt = c(12.07, 13.65, 16.72), name = c("a", "b", "c"))
  ref <- match_peaks(pks, Sa_pr, targets = tg, lambda = 210)
  tg$spectral_cutoff <- NA
  expect_equal(match_peaks(pks, Sa_pr, targets = tg, lambda = 210), ref)
})
