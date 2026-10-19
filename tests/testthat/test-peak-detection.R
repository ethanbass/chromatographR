# Tests for find_peaks correctness using synthetic signals with known peak positions.

make_gaussian_signal <- function(n, centers, widths, heights) {
  x <- seq_len(n)
  y <- numeric(n)
  for (i in seq_along(centers))
    y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * widths[i]^2))
  y
}

# ── Single peak ───────────────────────────────────────────────────────────────

test_that("find_peaks detects a single Gaussian peak at the correct position", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = 151, widths = 10, heights = 1000)
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 1L)
  expect_equal(p$pos, 151, tolerance = 1)
  expect_true(p$lower < p$pos)
  expect_true(p$upper > p$pos)
})

# ── Two peaks ─────────────────────────────────────────────────────────────────

test_that("find_peaks detects two well-separated Gaussian peaks in the right order", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(101, 201), widths = c(10, 10),
                    heights = c(1000, 800))
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 2L)
  expect_true(p$pos[1] < p$pos[2])
  expect_equal(p$pos[1], 101, tolerance = 1)
  expect_equal(p$pos[2], 201, tolerance = 1)
})

# ── min_amplitude filtering ──────────────────────────────────────────────────────

test_that("find_peaks min_amplitude removes peaks below the amplitude threshold", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(101, 201), widths = c(10, 10),
                    heights = c(1000, 10))
  p_all  <- find_peaks(y, smooth_type = "none", min_amplitude = 0)
  p_tall <- find_peaks(y, smooth_type = "none", min_amplitude = 50)
  expect_equal(nrow(p_all),  2L)
  expect_equal(nrow(p_tall), 1L)
  expect_equal(p_tall$pos, 101, tolerance = 1)
})

# ── bounds = FALSE ────────────────────────────────────────────────────────────

test_that("find_peaks with bounds = FALSE returns a numeric vector", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = 151, widths = 10, heights = 1000)
  p <- find_peaks(y, smooth_type = "none", bounds = FALSE)
  expect_true(is.numeric(p))
  expect_equal(length(p), 1L)
  expect_equal(p, 151, tolerance = 1)
})

test_that("find_peaks keeps peak windows clear of missing values at the ends", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 301, centers = c(40, 270), widths = c(8, 8),
                            heights = c(1000, 1000))
  y[1:20] <- NA
  y[281:301] <- NA
  p <- find_peaks(y, smooth_type = "none")
  expect_equal(nrow(p), 2L)
  expect_equal(p$lower[1], 21)
  expect_equal(p$upper[2], 280)
  x <- matrix(y, ncol = 1, dimnames = list(seq_along(y) / 10, "210"))
  expect_no_warning(pks <- fit_peaks(x, lambda = "210"))
  expect_false(anyNA(pks$area))
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("find_peaks errors when y is not a vector", {
  expect_error(find_peaks(matrix(1:9, 3)), regexp = "Please provide a vector")
})

# ── Restricting peak fitting to target retention times ────────────────────────

data(Sa_pr)
lam_rt   <- "210"
pks_full <- get_peaks(Sa_pr, lambdas = lam_rt, show_progress = FALSE)
# Two apex retention times known to be present in sample 1.
rts_2 <- pks_full[[1]][[lam_rt]]$rt[c(3, 8)]

test_that("`targets = NULL` leaves get_peaks unchanged", {
  skip_on_cran()
  expect_equal(pks_full,
               get_peaks(Sa_pr, lambdas = lam_rt, targets = NULL, tol = 0.1,
                         select_by = "height", show_progress = FALSE))
})

test_that("get_peaks restricts fitting to peaks near `targets`", {
  skip_on_cran()
  tgt <- get_peaks(Sa_pr, lambdas = lam_rt, targets = rts_2, tol = 0.05,
                   show_progress = FALSE)
  expect_s3_class(tgt, "peak_list")
  expect_equal(names(tgt), names(Sa_pr))
  # At most one peak per target. Note `tol` is applied to the DETECTED apex,
  # while `rt` reports the MODEL-FITTED apex, so the reported `rt` may fall
  # slightly outside the window; the strict invariant is tested on
  # `select_peaks_by_rt` below.
  for (s in tgt){
    expect_lte(nrow(s[[lam_rt]]), length(rts_2))
  }
  # attributes are preserved
  expect_equal(attr(tgt, "fit"), attr(pks_full, "fit"))
  expect_equal(attr(tgt, "lambdas"), attr(pks_full, "lambdas"))
})

test_that("restricting by `targets` does not alter peak bounds", {
  # Guards the design decision to filter find_peaks' OUTPUT rather than
  # restrict its input. Bounds come from derivative sign changes over the whole
  # trace, so truncating the input would clip them.
  skip_on_cran()
  tgt <- get_peaks(Sa_pr, lambdas = lam_rt, targets = rts_2, tol = 0.05,
                   show_progress = FALSE)
  a <- pks_full[[1]][[lam_rt]]
  b <- tgt[[1]][[lam_rt]]
  expect_equal(a[a$rt %in% b$rt, c("rt", "start", "end")],
               b[, c("rt", "start", "end")], ignore_attr = TRUE)
})

test_that("`select_by` resolves multiple peaks in one window", {
  skip_on_cran()
  # In sample 1 the peaks at 11.20 and 11.36 both fall within 0.2 of 11.20.
  # The nearer peak (11.20) is the shorter of the two, so `height` and
  # `proximity` must disagree here or this test proves nothing.
  pick <- function(sb){
    get_peaks(Sa_pr[1], lambdas = lam_rt, targets = 11.20, tol = 0.2,
              select_by = sb, show_progress = FALSE, cl = 1,
              smooth_window = 1)[[1]][[lam_rt]]$rt
  }
  all_rt <- get_peaks(Sa_pr[1], lambdas = lam_rt, show_progress = FALSE,
                      cl = 1, smooth_window = 1)[[1]][[lam_rt]]$rt
  near <- function(t) all_rt[which.min(abs(all_rt - t))]
  expect_equal(pick("height"), near(11.36))
  expect_equal(pick("proximity"), near(11.20))
  expect_equal(sort(pick("none")), c(near(11.20), near(11.36)))
})

test_that("fitted retention times are not snapped to the time grid", {
  skip_on_cran()
  ts <- seq(10, 12, by = 0.02)
  y <- 100 * exp(-(ts - 11.013)^2 / (2 * 0.05^2))
  chrom <- list(s = matrix(y, ncol = 1,
                           dimnames = list(format(ts, nsmall = 2), "210")))
  for (f in c("gaussian", "egh", "bemg")){
    rt <- suppressWarnings(get_peaks(chrom, lambdas = 210, fit = f, cl = 1,
                                     show_progress = FALSE))[[1]][["210"]]$rt
    expect_equal(rt, 11.013, tolerance = 1e-5, info = f)
  }
})

test_that("a peak returned for a target is not reported missing", {
  skip_on_cran()
  data(Sa_warp)
  expect_no_warning(
    pks <- get_peaks(Sa_warp[1], lambdas = 210, fit = "gaussian", targets = 18.44,
                     tol = 0.01, cl = 1, show_progress = FALSE))
  expect_equal(nrow(pks[[1]][["210"]]), 1L)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(pks, targets = 18.44, tol = 0.01))
  expect_warning(get_peaks(Sa_pr, lambdas = 210, targets = c(11.36, 11.37),
                           tol = c(0.05, 0.001), cl = 1,
                           show_progress = FALSE),
                 regexp = "No peak was assigned to the following targets in any sample: 11.37\\.")
})

test_that("`tol` may be a vector with one value per target", {
  skip_on_cran()
  n <- function(...) sapply(get_peaks(Sa_pr, lambdas = lam_rt,
                                      targets = c(11.18, 12.36), ...,
                                      show_progress = FALSE),
                            function(s) nrow(s[[lam_rt]]))
  expect_equal(n(tol = 0.1), n(tol = c(0.1, 0.1)))
  # a looser tolerance on the second target recovers peaks a tight one misses
  expect_gt(sum(n(tol = c(0.02, 0.5))), sum(n(tol = 0.02)))
})

test_that("overlapping windows do not duplicate a shared peak", {
  skip_on_cran()
  expect_warning(
    d <- get_peaks(Sa_pr[1], lambdas = lam_rt, targets = c(11.36, 11.37), tol = 0.1,
                   select_by = "none", show_progress = FALSE,
                   cl = 1)[[1]][[lam_rt]],
    regexp = "No peak was assigned to the following targets in any sample: 11.36\\.")
  expect_false(any(duplicated(d$rt)))
})

test_that("targets matching nothing give empty tables and a warning", {
  skip_on_cran()
  expect_warning(
    tgt <- get_peaks(Sa_pr[1], lambdas = lam_rt, targets = 11.181, tol = 0.0001,
                     show_progress = FALSE, cl = 1),
    regexp = "No peak was assigned to the following targets in any sample")
  d <- tgt[[1]][[lam_rt]]
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 0L)
  # empty tables keep the full schema so downstream rbinds still work
  expect_true(all(c("sample", "lambda", "rt", "start", "end") %in% colnames(d)))
})

test_that("restricted peak lists compose with get_peaktable", {
  skip_on_cran()
  tgt <- get_peaks(Sa_pr, lambdas = lam_rt, targets = c(11.18, 12.36, 15.94),
                   tol = 0.1, show_progress = FALSE)
  pt <- get_peaktable(tgt, chrom_list = Sa_pr)
  expect_s3_class(pt, "peak_table")
  expect_equal(nrow(pt$tab), length(Sa_pr))
  expect_equal(ncol(pt$tab), 3L)
  expect_equal(colnames(pt$tab), colnames(pt$pk_meta))
})

test_that("targets closer than `hmax` get separate peak table columns", {
  skip_on_cran()
  pk <- get_peaks(Sa_pr, lambdas = "220", targets = c(11.63, 11.77), tol = 0.03,
                  show_progress = FALSE, cl = 1)
  expect_equal(nrow(pk[["119"]][["220"]]), 2L)
  pt <- get_peaktable(pk)
  expect_equal(ncol(pt$tab), 2L)
  expect_true(all(pt$tab["119", ] > 0))
})

test_that("a targets data.frame sets names, wavelengths and column order", {
  skip_on_cran()
  tg <- data.frame(name = c("b", "a"), rt = c(12.36, 11.63),
                   lambda = c(280, 220))
  pk <- get_peaks(Sa_pr, targets = tg, tol = 0.03, show_progress = FALSE, cl = 1)
  expect_equal(names(pk[[1]]), c("280", "220"))
  pt <- get_peaktable(pk)
  expect_equal(colnames(pt$tab), c("b", "a"))
  expect_equal(colnames(pt$pk_meta), c("b", "a"))
  expect_equal(unname(unlist(pt$pk_meta["lambda", ])), c(280, 220))

  tg3 <- rbind(tg, data.frame(name = "c", rt = 11.77, lambda = 220))
  pt3 <- get_peaktable(get_peaks(Sa_pr, targets = tg3, tol = 0.03,
                                 show_progress = FALSE, cl = 1))
  expect_equal(colnames(pt3$tab), c("b", "a", "c"))
  expect_equal(pt3$tab[, c("b", "a")], pt$tab)
})

test_that("a target detected in no sample gets no column", {
  skip_on_cran()
  tg <- data.frame(name = c("a", "ghost", "b"), rt = c(11.63, 11.181, 12.36),
                   tol = c(0.03, 0.0001, 0.03), lambda = c(220, 220, 280))
  expect_warning(pk <- get_peaks(Sa_pr, targets = tg, show_progress = FALSE, cl = 1),
                 regexp = "No peak was assigned to the following targets in any sample: ghost\\.")
  pt <- get_peaktable(pk)
  expect_equal(colnames(pt$tab), c("a", "b"))
  expect_no_error(attach_ref_spectra(pt, chrom_list = Sa_pr))
})

test_that("unnamed targets with wavelengths keep the order of `targets`", {
  skip_on_cran()
  tg <- data.frame(rt = c(11.63, 12.36, 11.77), lambda = c(220, 280, 220))
  pt <- get_peaktable(get_peaks(Sa_pr, targets = tg, tol = 0.03,
                                show_progress = FALSE, cl = 1))
  expect_equal(colnames(pt$tab), c("V1", "V2", "V3"))
  expect_equal(unname(unlist(pt$pk_meta["lambda", ])), c(220, 280, 220))
})

test_that("targets are assigned when peak bounds reach 1e5 timepoints", {
  skip_on_cran()
  i <- seq_len(200000)
  times <- (i - 1) * 0.0005
  y <- 100 * exp(-(i - 90000)^2 / (2 * 3000^2)) +
    100 * exp(-(i - 110000)^2 / (2 * 3000^2))
  m <- matrix(y, ncol = 1, dimnames = list(times, "220"))
  pk <- get_peaks(list(s1 = m, s2 = m), lambdas = "220",
                  targets = times[c(90000, 110000)], tol = 0.1, fit = "raw",
                  smooth_type = "none", show_progress = FALSE, cl = 1)
  expect_equal(pk$s1[["220"]]$target, c(1L, 2L))
  expect_equal(ncol(get_peaktable(pk)$tab), 2L)
})

test_that("named target columns work downstream", {
  skip_on_cran()
  tg <- data.frame(name = c("a", "b"), rt = c(11.63, 11.77), lambda = 220)
  pk <- get_peaks(Sa_pr, targets = tg, tol = 0.03, show_progress = FALSE, cl = 1)
  pt <- get_peaktable(filter_peaks(pk, min_area = 0.1))
  expect_equal(colnames(pt$tab), c("a", "b"))
  wide <- reshape_peaktable(pt, peaks = c("a", "b"), format = "wide")
  expect_true(all(c("a", "b") %in% colnames(wide)))
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot_spectrum("a", pt, chrom_list = Sa_pr, verbose = FALSE))
})

test_that("`mod_list` reads targets on the aligned axis", {
  skip_on_cran()
  mods <- correct_rt(Sa_pr, lambdas = "220", what = "models", alg = "ptw")
  f <- function(mod_list){
    get_peaks(Sa_pr, lambdas = "220", targets = 14.41, tol = 0.05,
              select_by = "proximity", mod_list = mod_list,
              show_progress = FALSE, cl = 1)
  }
  al <- f(mods)
  expect_true(all(vapply(al, function(s) nrow(s[["220"]]), numeric(1)) == 1))
  expect_false("rt.cor" %in% colnames(al[[1]][["220"]]))
  maps <- warp_maps(mods, Sa_pr)
  aligned <- vapply(names(al), function(s) maps[[s]](al[[s]][["220"]]$rt),
                    numeric(1))
  expect_true(all(abs(aligned - 14.41) <= 0.05))
  raw <- f(NULL)
  expect_lt(sum(vapply(raw, function(s) nrow(s[["220"]]), numeric(1))), 4)
})

test_that("with `mod_list`, `targets` may lie anywhere on the aligned axis", {
  skip_on_cran()
  mods <- correct_rt(Sa_pr, lambdas = "220", what = "models", alg = "vpdtw")
  pk <- get_peaks(Sa_pr, lambdas = "220", targets = 18.9, tol = 0.05,
                  mod_list = mods, show_progress = FALSE, cl = 1)
  expect_equal(nrow(pk[["458"]][["220"]]), 1L)
  expect_lt(pk[["458"]][["220"]]$rt, max(get_times(Sa_pr)))
})

test_that("`targets` is checked against every sample's time range", {
  skip_on_cran()
  s <- Sa_pr
  s[[1]] <- s[[1]][get_times(Sa_pr) <= 15, ]
  s <- transfer_metadata(s, Sa_pr)
  pk <- get_peaks(s, lambdas = "210", targets = 16.1, tol = 0.1,
                  show_progress = FALSE, cl = 1)
  expect_equal(nrow(pk[[1]][["210"]]), 0L)
  expect_error(get_peaks(s, lambdas = "210", targets = 19, show_progress = FALSE),
               regexp = "outside the bounds")
})

test_that("targets left without a peak are reported", {
  skip_on_cran()
  expect_warning(expect_warning(
    get_peaks(Sa_pr, lambdas = "210", max_fwhm = 0.105,
              targets = data.frame(rt = c(10.654, 11.342, 12.123),
                               name = c("a", "b", "c")),
              show_progress = FALSE, cl = 1),
    regexp = "No peak was assigned to the following targets in any sample: a, b\\."),
    regexp = "discarded for exceeding `max_fwhm`")
  expect_warning(
    get_peaks(Sa_pr, lambdas = "210",
              targets = data.frame(rt = c(12.123, 12.123), name = c("x", "y")),
              show_progress = FALSE, cl = 1),
    regexp = "No peak was assigned to the following targets in any sample: y\\.")
})

test_that("a peak detected as two apexes is fitted once", {
  times <- seq(0, 10, length.out = 1001)
  m <- matrix(pmin(100 * exp(-(times - 5)^2 / (2 * 0.2^2)), 60), ncol = 1,
              dimnames = list(times, "220"))
  expect_warning(
    pk <- get_peaks(list(s1 = m), lambdas = "220", fit = "raw",
                    targets = data.frame(rt = c(4.85, 5.15), name = c("x", "y")),
                    show_progress = FALSE, cl = 1),
    regexp = "No peak was assigned to the following targets in any sample: y\\.")
  expect_equal(nrow(pk$s1[["220"]]), 1L)
  expect_equal(pk$s1[["220"]]$target, 1L)
})

test_that("unnamed targets at one wavelength are named by position", {
  skip_on_cran()
  pk <- suppressWarnings(get_peaks(Sa_pr, lambdas = "220",
                                   targets = c(11.63, 11.181, 11.77),
                                   tol = c(0.03, 0.0001, 0.03),
                                   show_progress = FALSE, cl = 1))
  expect_equal(colnames(get_peaktable(pk)$tab), c("V1", "V3"))
})

test_that("check_targets reads factors by value and normalizes wavelengths", {
  tg <- check_targets(data.frame(rt = factor(c(12.5, 11.2)),
                                 tol = factor(c(0.1, 0.2)),
                                 lambda = c("220.0", "280")), tol = 0.1)
  expect_equal(tg$rt, c(12.5, 11.2))
  expect_equal(tg$tol, c(0.1, 0.2))
  expect_equal(tg$lambda, c("220", "280"))
  expect_equal(check_targets(factor(c(12.5, 11.2)), tol = 0.1)$rt, c(12.5, 11.2))
})

test_that("get_peaks validates a targets data.frame and `mod_list`", {
  expect_error(get_peaks(Sa_pr, lambdas = "220", targets = 11.6, bounds = FALSE),
               regexp = "cannot be used with `fit_peaks`")
  expect_error(get_peaks(Sa_pr, targets = data.frame(time = 11.6)),
               regexp = "must have an `rt` column")
  expect_error(get_peaks(Sa_pr, targets = data.frame(name = c("a", "a"),
                                                 rt = c(11.6, 12.3),
                                                 lambda = 220)),
               regexp = "must be unique. Duplicated: a")
  expect_error(get_peaks(Sa_pr, lambdas = "220",
                         targets = data.frame(rt = 11.6, lambda = 220)),
               regexp = "not both")
  expect_error(get_peaks(Sa_pr, lambdas = c("210", "220"),
                         targets = data.frame(name = "a", rt = 11.6)),
               regexp = "duplicated across wavelengths")
  expect_error(get_peaks(Sa_pr, targets = data.frame(rt = 11.6, lambda = NA)),
               regexp = "missing or empty")
  expect_error(get_peaks(Sa_pr, lambdas = "220", mod_list = list()),
               regexp = "only used to place")
})

test_that("get_peaks validates `targets`, `tol` and `select_by`", {
  expect_error(get_peaks(Sa_pr, lambdas = lam_rt, targets = 99,
                         show_progress = FALSE),
               regexp = "outside the bounds")
  expect_error(get_peaks(Sa_pr, lambdas = lam_rt, targets = c(11.18, 12.36),
                         tol = c(1, 2, 3), show_progress = FALSE),
               regexp = "single value or one value")
  expect_error(get_peaks(Sa_pr, lambdas = lam_rt, targets = 11.18,
                         select_by = "nonsense", show_progress = FALSE),
               regexp = "'arg' should be one of")
  expect_error(get_peaks(Sa_pr, lambdas = lam_rt, targets = "abc",
                         show_progress = FALSE),
               regexp = "numeric vector of retention times")
  expect_error(get_peaks(Sa_pr, lambdas = lam_rt,
                         targets = data.frame(rt = 11.18, tol = 0.1), tol = 0.2,
                         show_progress = FALSE),
               regexp = "either in `tol` or in the `tol` column")
})

test_that("a failed target fit is NA in the peak table, and a missing one 0", {
  skip_on_cran()
  pks <- get_peaks(Sa_pr[1:2], lambdas = lam_rt, targets = c(11.18, 12.36),
                   cl = 1, show_progress = FALSE, smooth_window = 1)
  pks[[1]][[lam_rt]][pks[[1]][[lam_rt]]$target == 2, c("height", "area")] <- NA
  pks[[2]][[lam_rt]] <- pks[[2]][[lam_rt]][pks[[2]][[lam_rt]]$target == 2, ]
  tab <- get_peaktable(pks)$tab
  expect_true(is.na(tab[1, 2]))
  expect_gt(tab[2, 2], 0)
  expect_equal(tab[2, 1], 0)
  expect_gt(tab[1, 1], 0)
})

test_that("select_peaks_by_rt keeps only apexes within `tol`", {
  # The helper's contract holds exactly, since it works on detected apex
  # positions rather than model-fitted retention times.
  times <- seq(0, 10, by = 0.1)
  y <- numeric(length(times))
  apex <- c(11L, 31L, 51L)            # rt 1.0, 3.0, 5.0
  y[apex] <- c(10, 100, 50)
  pos <- data.frame(pos = apex, lower = apex - 3L, upper = apex + 3L)
  f <- function(...) select_peaks_by_rt(pos, y = y,
                                                         times = times, ...)
  expect_equal(f(rts = 3.0, tol = 0.05)$pos, 31L)
  expect_equal(f(rts = c(1.0, 5.0), tol = 0.05)$pos, c(11L, 51L))
  expect_equal(nrow(f(rts = 8.0, tol = 0.05)), 0L)
  # both apexes in one window: height wins, proximity wins, none keeps both
  expect_equal(f(rts = 1.4, tol = 1.7, select_by = "height")$pos, 31L)
  expect_equal(f(rts = 1.4, tol = 1.7, select_by = "proximity")$pos, 11L)
  expect_equal(f(rts = 1.4, tol = 1.7, select_by = "none")$pos, c(11L, 31L))
  # `tol` is validated
  expect_error(f(rts = c(1, 3), tol = c(1, 2, 3)),
               regexp = "single value or one value")
  expect_error(f(rts = 1, tol = -1), regexp = "non-negative")
})

test_that("match_rts records the target that selected each peak", {
  m <- match_rts(c(0.9, 1.15), rts = c(1.0, 1.35), tol = 0.5,
                 select_by = "proximity")
  expect_equal(as.vector(m), 1L)
  expect_equal(attr(m, "target"), 1L)
  m <- match_rts(1.1, rts = c(1.0, 1.3), tol = 0.5, select_by = "proximity")
  expect_equal(attr(m, "target"), 1L)
})

test_that("overlapping windows are split between the nearest targets", {
  m <- match_rts(c(0.95, 1.06, 1.15), rts = c(1.0, 1.1), tol = 0.1,
                 weight = c(1, 10, 1), select_by = "height")
  expect_equal(as.vector(m), c(1L, 2L))
  expect_equal(attr(m, "target"), c(1L, 2L))
  for (w in list(c(100, 200), c(200, 100))){
    m <- match_rts(c(13.05, 13.15), rts = c(12.9, 13.1), tol = 0.3,
                   weight = w, select_by = "height")
    expect_equal(attr(m, "target"), 2L)
    expect_equal(as.vector(m), which.max(w))
  }
  m <- match_rts(c(1.04, 1.08), rts = c(1.0, 1.2), tol = c(0.05, 0.2),
                 weight = c(1, 1), select_by = "none")
  expect_equal(attr(m, "target"), c(1L, 2L))
  m <- match_rts(c(NA, 1.0), rts = 1.0, tol = 0.1, weight = c(1, 1))
  expect_equal(as.vector(m), 2L)
  expect_length(match_rts(numeric(0), rts = 1.0), 0L)
})

test_that("a single cluster keeps a well-formed pk_meta", {
  skip_on_cran()
  pk <- filter_peaks(get_peaks(Sa_pr, lambdas = "220", show_progress = FALSE,
                               cl = 1), min_rt = 11.55, max_rt = 11.7)
  pt <- get_peaktable(pk)
  expect_equal(ncol(pt$tab), 1L)
  expect_equal(ncol(pt$pk_meta), 1L)
  expect_true(all(c("lambda", "rt", "start", "end") %in% rownames(pt$pk_meta)))
})

# ── `bounds = FALSE` is not usable for integration ────────────────────────────

test_that("`bounds = FALSE` is rejected by fit_peaks and get_peaks", {
  data(Sa_pr)
  # `find_peaks(bounds = FALSE)` returns a bare vector of apex positions, which
  # cannot be integrated because every peak model needs a window. Previously
  # this reached `fit_peaks` through `...` and failed with the opaque message
  # "argument is of length zero".
  expect_error(fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw", bounds = FALSE),
               "cannot be used with `fit_peaks`")
  expect_error(get_peaks(Sa_pr[1], lambdas = 210, fit = "raw", bounds = FALSE,
                         cl = 1, show_progress = FALSE),
               "cannot be used with `fit_peaks`")
})

test_that("`bounds = TRUE` remains a harmless no-op and find_peaks is unaffected", {
  skip_on_cran()
  data(Sa_pr)
  expect_equal(fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw", bounds = TRUE),
               fit_peaks(Sa_pr[[1]], lambda = 210, fit = "raw"))
  # `find_peaks` still honors it when called directly, as `plot_spectrum` does
  expect_null(dim(find_peaks(Sa_pr[[1]][, "210"], bounds = FALSE)))
  expect_s3_class(find_peaks(Sa_pr[[1]][, "210"], bounds = TRUE), "data.frame")
})

# ── Apex alignment (C0) ───────────────────────────────────────────────────────

test_that("find_peaks reports the true apex, not the timepoint before it", {
  # `d[i] = y[i+1] - y[i]`, so a downward zero-crossing between `d[i]` and
  # `d[i+1]` puts the apex at `y[i+1]`. Reporting `i` made every `pos` one
  # sampling interval early.
  expect_equal(find_peaks(c(0, 1, 2, 5, 20, 100, 20, 5, 2, 1, 0),
                          smooth_type = "none")$pos, 6)
})

test_that("find_peaks pos equals which.max for every smoothing method", {
  skip_on_cran()
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  for (st in c("none", "gaussian", "box", "mva", "tmva")){
    p <- find_peaks(y, smooth_type = st)
    expect_equal(p$pos[which.max(p$prominence)], which.max(y),
                 info = st)
  }
})

test_that("min_amplitude is tested at the apex, not on the rising flank", {
  # `preprocess` used to clip negatives to exactly zero, and `y > min_amplitude`
  # is false at zero, so a peak whose preceding timepoint had been clipped was
  # dropped even at the inert `min_amplitude = 0` default.
  mk <- function(preceding){
    y <- rep(0, 50); y[20] <- preceding; y[21] <- 8; y[22] <- 3; y
  }
  expect_equal(nrow(find_peaks(mk(3), smooth_type = "none")), 1L)
  expect_equal(nrow(find_peaks(mk(0), smooth_type = "none")), 1L)
})

test_that("the deprecated `slope_thresh` and `amp_thresh` still work", {
  y <- make_gaussian_signal(n = 601, centers = c(150, 301, 450),
                            widths = c(4, 20, 40), heights = c(20, 1000, 500))
  expect_warning(old <- find_peaks(y, smooth_type = "none", amp_thresh = 50),
                 regexp = "`amp_thresh` argument is deprecated")
  expect_equal(old, find_peaks(y, smooth_type = "none", min_amplitude = 50))
  expect_warning(old <- find_peaks(y, smooth_type = "none", slope_thresh = 1),
                 regexp = "`slope_thresh` argument is deprecated")
  expect_equal(old, find_peaks(y, smooth_type = "none", min_curvature = 1))
})

test_that("min_curvature filters on curvature, excluding broad peaks first", {
  y <- make_gaussian_signal(n = 601, centers = c(150, 450),
                            widths = c(4, 40), heights = c(500, 500))
  expect_equal(nrow(find_peaks(y, smooth_type = "none", min_curvature = 0)), 2L)
  broad_gone <- find_peaks(y, smooth_type = "none", min_curvature = 1)
  expect_equal(broad_gone$pos, 150L)
})

test_that("savgol rejects a window too narrow to smooth", {
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  expect_error(find_peaks(y, smooth_type = "savgol"), regexp = "too narrow")
  expect_error(find_peaks(y, smooth_type = "savgol", smooth_window = 5),
               regexp = "returns the signal unchanged")
  p <- find_peaks(y, smooth_type = "savgol", smooth_window = 21)
  expect_equal(p$pos[which.max(p$prominence)], 301L)
  expect_equal(nrow(find_peaks(y, smooth_type = "savgol", smooth_window = 20)),
               nrow(find_peaks(y, smooth_type = "savgol", smooth_window = 21)))
})

test_that("even-width moving averages do not shift the reported apex", {
  # An even-width running mean is not centered on a sample: it creates a
  # plateau at the extremum, giving two consecutive sign changes.
  skip_on_cran()
  y <- make_gaussian_signal(n = 601, centers = 301, widths = 20, heights = 1000)
  for (st in c("mva", "tmva")){
    p <- find_peaks(y, smooth_type = st, smooth_window = 4)
    expect_equal(p$pos[which.max(p$prominence)], 301L, info = st)
  }
})

# ── Noise-referenced prominence filter (C) ────────────────────────────────────

test_that("min_prominence_snr = 0 leaves find_peaks output unchanged", {
  skip_on_cran()
  data(Sa_pr)
  for (nm in names(Sa_pr)) for (lam in c("210", "260", "310")){
    y <- Sa_pr[[nm]][, lam]
    expect_identical(find_peaks(y), find_peaks(y, min_prominence_snr = 0))
  }
})

test_that("find_peaks returns a prominence column", {
  data(Sa_pr)
  p <- find_peaks(Sa_pr[[1]][, "210"])
  expect_true("prominence" %in% names(p))
  expect_equal(nrow(p), length(p$prominence))
  expect_true(all(p$prominence >= 0))
})

test_that("raising min_prominence_snr drops the lowest-prominence apices first", {
  skip_on_cran()
  data(Sa_pr)
  y <- Sa_pr[[1]][, "210"]
  levels <- c(0, 1, 2, 4, 8, 16)
  sets <- lapply(levels, function(k) find_peaks(y, smooth_window = 5, min_prominence_snr = k))
  # monotonically fewer peaks, and each set nested in the looser one
  expect_true(all(diff(vapply(sets, nrow, integer(1))) <= 0))
  for (i in seq_len(length(sets) - 1))
    expect_true(all(sets[[i + 1]]$pos %in% sets[[i]]$pos))
  # and the retained minimum prominence rises
  mins <- vapply(sets, function(p) min(p$prominence), numeric(1))
  expect_true(all(diff(mins) > 0))
})

test_that("a noise-referenced threshold survives a dominant solvent front", {
  # A single large peak turns any threshold expressed as a fraction of the
  # largest feature into a hundreds-of-sigma filter. Here 1% of full scale is
  # 200 units against a noise sd of ~1.
  skip_on_cran()
  set.seed(42)
  tt <- seq(0, 30, by = 0.01)
  g <- function(h, m, s) h * exp(-(tt - m)^2 / (2 * s^2))
  y <- g(20000, 2, 0.05) + g(30, 12, 0.15) + g(8, 20, 0.15) +
    rnorm(length(tt))
  sw <- 0.01  # ~30 timepoints, comparable to the analyte peak width
  near <- function(p, target)
    any(abs(p$pos - which.min(abs(tt - target))) <= 5)

  all_pks <- find_peaks(y, smooth_window = sw)
  rel <- all_pks[all_pks$prominence >
                   0.01 * max(all_pks$prominence, na.rm = TRUE), ]
  snr <- find_peaks(y, smooth_window = sw, min_prominence_snr = 4)

  # a max-referenced threshold loses both analytes ...
  expect_false(near(rel, 12))
  expect_false(near(rel, 20))
  # ... while a noise-referenced one keeps them
  expect_true(near(snr, 12))
  expect_true(near(snr, 20))
})

test_that("min_prominence_snr requires bounds = TRUE", {
  expect_error(find_peaks(c(0, 1, 5, 1, 0), min_prominence_snr = 4,
                          bounds = FALSE),
               regexp = "requires `bounds = TRUE`")
  expect_type(find_peaks(c(0, 1, 5, 1, 0), bounds = FALSE), "integer")
})


test_that("the prominence filter is calibrated on correlated noise", {
  set.seed(1)
  n <- 2000
  kern <- dnorm(-8:8, sd = 2)
  noise <- as.numeric(stats::filter(rnorm(n + 16), kern / sum(kern), sides = 2))[9:(n + 8)]
  noise <- noise / sd(noise)
  rt <- seq(100, 1900, length.out = 10)
  height <- exp(seq(log(5), log(300), length.out = 10))
  sd_peak <- 20 / 2.355
  y <- noise + rowSums(sapply(seq_along(rt), function(i)
    height[i] * exp(-0.5 * ((seq_len(n) - rt[i]) / sd_peak)^2)))
  p <- find_peaks(y, smooth_window = 10, min_prominence_snr = 4)
  spurious <- sum(vapply(p$pos, function(x) min(abs(x - rt)) > 10, logical(1)))
  expect_lte(spurious, 2)
  expect_gte(nrow(p) - spurious, 9)
})

test_that("find_peaks warns when smooth_window is too narrow for the prominence filter", {
  set.seed(1)
  n <- 2000
  main <- seq(100, 1850, length.out = 12)
  rt <- c(main, main + 16)
  height <- c(exp(seq(log(20), log(300), length.out = 12)), 0.35 * exp(seq(log(20), log(300), length.out = 12)))
  sd_peak <- 20 / 2.355
  y <- rnorm(n) + rowSums(sapply(seq_along(rt), function(i)
    height[i] * exp(-0.5 * ((seq_len(n) - rt[i]) / sd_peak)^2)))
  expect_warning(find_peaks(y, smooth_window = 2, min_prominence_snr = 3),
                 regexp = "under a third of the typical peak width")
  expect_no_warning(find_peaks(y, smooth_window = 10, min_prominence_snr = 3))
  expect_no_warning(find_peaks(y, smooth_window = 2))
})

# ── exponentially modified Gaussians must not overflow ───────────────────────

test_that("log_erfc agrees with log(erfc()) and survives the far upper tail", {
  x <- seq(-3, 5, by = 0.25)
  expect_equal(log_erfc(x), log(erfc(x)), tolerance = 1e-12)
  expect_true(all(is.finite(log_erfc(seq(10, 200, by = 10)))))
  expect_true(all(is.infinite(log(erfc(seq(30, 200, by = 10))))))
})

test_that("bemg stays finite where the naive exp(A) * erfc(z) form overflows", {
  x <- seq(2.79, 3.185, by = 0.005)
  naive <- function(center, width, height, a, b){
    tm <- x - center; as2 <- a * width^2; bs2 <- b * width^2
    s <- sqrt(2 * width^2)
    height * 0.5 * (exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / s) +
                    exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / s))
  }
  expect_true(any(!is.finite(naive(2.845, 40, 960, 1, 1))))
  expect_true(all(is.finite(bemg(x, 2.845, 40, 960, 1, 1))))
  expect_true(all(is.finite(bemg(x, 2.9, 25, 500, 2, 2))))
  expect_true(all(bemg(x, 2.845, 40, 960, 1, 1) >= 0))
})

test_that("bemg is unchanged where the naive form was already finite", {
  x <- seq(-1, 1, length.out = 101)
  naive <- function(center, width, height, a, b){
    tm <- x - center; as2 <- a * width^2; bs2 <- b * width^2
    s <- sqrt(2 * width^2)
    height * 0.5 * (exp(a * as2 / 2 + a * tm) * erfc((as2 + tm) / s) +
                    exp(b * bs2 / 2 - b * tm) * erfc((bs2 - tm) / s))
  }
  for (p in list(c(0.5, 0.05, 500, 2, 3), c(-0.3, 0.1, 50, 0.5, 0.8),
                 c(0, 0.2, 1000, 0.3, 4))){
    ref <- naive(p[1], p[2], p[3], p[4], p[5])
    skip_if(any(!is.finite(ref)))
    expect_equal(bemg(x, p[1], p[2], p[3], p[4], p[5]), ref, tolerance = 1e-10)
  }
})

test_that("a peak whose model yields a non-finite intensity is dropped, or kept as NA", {
  d <- data.frame(rt = c(10, 20, 30), sd = c(1, 2, 3),
                  height = c(100, NaN, 300), area = c(5, 6, Inf),
                  converged = c(1, 1, 1))
  expect_silent(out <- drop_unfittable_peaks(d))
  expect_equal(out$rt, 10)
  kept <- drop_unfittable_peaks(d, keep = TRUE)
  expect_equal(kept$rt, d$rt)
  expect_equal(kept$height, c(100, NA, NA))
  expect_equal(kept$area, c(5, NA, NA))
  expect_equal(kept$converged, c(1, 0, 0))
})

test_that("failed fits of requested peaks are kept, and survive unrelated filters", {
  skip_on_cran()
  data(Sa_pr)
  egh <- fitpk_egh
  local_mocked_bindings(fitpk_egh = function(x, pos, ...){
    r <- egh(x, pos, ...)
    if (abs(get_times(x)[pos[[1]]] - 11.37) < 0.05) r["area"] <- NaN
    r
  })
  pks <- get_peaks(Sa_pr, lambdas = "210", targets = c(10.7, 11.37), tol = 0.1,
                   show_progress = FALSE, cl = 1)
  failed <- do.call(rbind, lapply(pks, `[[`, "210"))
  failed <- failed[failed$target == 2, ]
  expect_gt(nrow(failed), 0)
  expect_true(all(is.na(failed$area) == (failed$converged == 0)))
  n_na <- function(pl) sum(vapply(pl, function(s) sum(is.na(s[["210"]]$area)),
                                  numeric(1)))
  expect_gt(n_na(pks), 0)
  expect_equal(n_na(filter_peaks(pks, min_rt = 10)), n_na(pks))
  expect_equal(n_na(filter_peaks(pks, min_area = 0)), 0)
  expect_no_error(get_peaktable(pks, summarize_by = "max"))
  detected <- get_peaks(Sa_pr, lambdas = "210", show_progress = FALSE, cl = 1)
  expect_equal(n_na(detected), 0)
})

test_that("an out-of-range index becomes NA instead of aborting the run", {
  cl <- list(a = matrix(1, nrow = 5, ncol = 1,
                        dimnames = list(as.character(seq(0.1, 0.5, by = 0.1)), "220")))
  x <- data.frame(rt = 2, start = 1, end = 5, center = -1726,
                  sd = 1, FWHM = 2, area = 3)
  out <- convert_indices_to_times(x, timepoints = get_times(cl),
                                  tdiff = get_time_resolution(cl), tfac = 1)
  expect_true(is.na(out$center))
  expect_equal(out$rt, 0.2)
})

# ── fitted width is bounded by the integration window ───────────────────────

test_that("fit_window_bounds scales with the window and allows a margin", {
  b <- fit_window_bounds(seq(0, 10, by = 1), mult = 2)
  expect_equal(b$span, 10)
  expect_equal(b$width_max, 2 * 10 / FWHM_CONST)
  expect_equal(b$center_min, -5)
  expect_equal(b$center_max, 15)
  expect_equal(fit_window_bounds(c(4, 4))$span, 1)
  expect_equal(fit_window_bounds(seq(0, 10, by = 1), mult = 1)$center_min, 0)
})

test_that("clamp_start moves a start value inside its bounds", {
  expect_equal(clamp_start(40, 0, 4.2), 4.2)
  expect_equal(clamp_start(-3, 0, 4.2), 0)
  expect_equal(clamp_start(2, 0, 4.2), 2)
  expect_equal(clamp_start(NaN, 0, 4), 2)
})

test_that("no peak model fits wider than twice its own window", {
  set.seed(42)
  x <- seq(0, 2, length.out = 401)
  # a broad feature with a narrow window carved out of its flank: the case that
  # previously produced a fit many times wider than the data supporting it
  y <- 500 * exp(-(x - 1)^2 / (2 * 0.4^2)) + rnorm(length(x), sd = 1)
  for (w in list(c(1.30, 1.36), c(0.50, 0.58), c(1.00, 1.30))){
    i <- which(x >= w[1] & x <= w[2])
    span <- x[max(i)] - x[min(i)]
    for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
      m <- suppressWarnings(get(f)(x[i], y[i], baseline = "sloped"))
      expect_lte(FWHM_CONST * abs(unname(m$width)), FIT_WIDTH_MAX_MULT * span * 1.001)
      expect_gte(unname(m$width), 0)
      expect_true(all(is.finite(m$y)))
    }
  }
})


test_that("window_limited marks fits that came to rest on a window limit", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_false(m$window_limited, info = f)
  }

  set.seed(42)
  xb <- seq(0, 2, length.out = 401)
  broad <- 500 * exp(-(xb - 1)^2 / (2 * 0.4^2)) + rnorm(length(xb), sd = 1)
  n_pinned <- 0L
  for (w in list(c(1.30, 1.36), c(0.50, 0.58), c(1.00, 1.30), c(1.35, 1.45))){
    i <- which(xb >= w[1] & xb <= w[2])
    width_max <- FIT_WIDTH_MAX_MULT * (xb[max(i)] - xb[min(i)]) / FWHM_CONST
    for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
      m <- suppressWarnings(get(f)(xb[i], broad[i], baseline = "sloped"))
      info <- paste(f, "window", w[1], w[2])
      if (abs(unname(m$width)[1] - width_max) < 1e-6 * width_max){
        n_pinned <- n_pinned + 1L
        expect_true(m$window_limited, info = info)
      }
    }
  }
  expect_gt(n_pinned, 0)
})

test_that("converged reports the optimizer's own convergence status", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_true(m$converged, info = f)
  }

  flat <- rep(0, 20)
  for (f in c("fit_egh", "fit_bemg", "fit_gaussian")){
    m <- suppressWarnings(get(f)(seq_along(flat), flat, baseline = "none"))
    expect_false(m$converged, info = f)
  }
})

test_that("the bound does not disturb a peak that fits comfortably inside it", {
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  for (f in c("fit_egh", "fit_gaussian")){
    m <- suppressWarnings(get(f)(x, y, baseline = "none"))
    expect_equal(unname(m$width), 0.02, tolerance = 0.02)
    expect_equal(unname(m$center), 0.5, tolerance = 0.01)
    expect_lt(FWHM_CONST * abs(unname(m$width)), FIT_WIDTH_MAX_MULT * 1)
  }
})

# ── max_fwhm is in time units and agrees with filter_peaks ──────────────────

.mk_chrom <- function(dt){
  ts <- seq(0, 3, by = dt)
  y <- 1000 * exp(-(ts - 1)^2/(2*0.02^2)) + 600 * exp(-(ts - 2)^2/(2*0.06^2))
  list(s1 = matrix(y, ncol = 1, dimnames = list(as.character(ts), "220")))
}

test_that("get_peaks(max_fwhm) and filter_peaks(max_fwhm) select the same peaks", {
  skip_on_cran()
  for (dt in c(0.005, 0.01)){
    cl <- .mk_chrom(dt)
    all <- get_peaks(cl, lambdas = "220", max_fwhm = Inf, estimate_purity = FALSE,
                     show_progress = FALSE)
    for (v in c(0.03, 0.1, 1)){
      a <- suppressWarnings(get_peaks(cl, lambdas = "220", max_fwhm = v,
                                      estimate_purity = FALSE,
                                      show_progress = FALSE))[[1]][[1]]
      b <- suppressWarnings(filter_peaks(all, max_fwhm = v))[[1]][[1]]
      expect_equal(nrow(a), nrow(b))
      expect_equal(a$rt, b$rt)
    }
  }
})

test_that("the max_fwhm default applies no width filtering", {
  skip_on_cran()
  dt <- 0.005
  cl <- .mk_chrom(dt)
  d <- get_peaks(cl, lambdas = "220", estimate_purity = FALSE, show_progress = FALSE)[[1]][[1]]
  e <- get_peaks(cl, lambdas = "220", max_fwhm = Inf,
                 estimate_purity = FALSE, show_progress = FALSE)[[1]][[1]]
  expect_equal(d$rt, e$rt)
  expect_equal(d$sd, e$sd)
})

test_that("sd_max is deprecated but keeps its timepoint semantics", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  expect_warning(
    expect_warning(old <- get_peaks(cl, lambdas = "220", sd_max = 8,
                                    estimate_purity = FALSE, show_progress = FALSE),
                   "deprecated"),
    "discarded")
  new <- suppressWarnings(get_peaks(cl, lambdas = "220", max_fwhm = FWHM_CONST * 8 * 0.005,
                   estimate_purity = FALSE, show_progress = FALSE))
  expect_length(old[[1]][[1]]$rt, 1)
  expect_equal(old[[1]][[1]]$rt, new[[1]][[1]]$rt)
  expect_error(get_peaks(cl, lambdas = "220", max_fwhm = 1, sd_max = 50),
               "not both")
})

# ── FWHM is a measured width, not a rescaled sigma ──────────────────────────

.dense_fwhm <- function(f, half_span){
  g <- seq(-half_span, half_span, length.out = 1000001)
  v <- f(g)
  apex <- which.max(v)
  half <- v[[apex]] / 2
  g[apex - 1 + min(which(v[apex:length(v)] < half))] -
    g[max(which(v[seq_len(apex)] < half))]
}

test_that("egh_fwhm is the half-height width of the egh shape", {
  skip_on_cran()
  for (tau in c(0, 0.1, 1, 5, 20, 50, -1, -5, -20)){
    shape <- function(t) egh(t, 0, 1, 1, tau)
    expect_equal(egh_fwhm(1, tau), .dense_fwhm(shape, 300 + 50 * abs(tau)),
                 tolerance = 1e-3, info = paste("tau =", tau))
  }
  # and it reduces to the gaussian identity when the peak is symmetric
  expect_equal(egh_fwhm(2.5, 0), FWHM_CONST * 2.5)
})

test_that("bemg_fwhm is the half-height width of the bemg shape", {
  skip_on_cran()
  for (ab in list(c(1, 1), c(0.3, 2), c(5, 0.2), c(0.05, 0.05))){
    m <- list(center = 0, width = 1, height = 1, a = ab[1], b = ab[2])
    shape <- function(t) bemg(t, 0, 1, 1, ab[1], ab[2])
    expect_equal(bemg_fwhm(m), .dense_fwhm(shape, 500 + 100 / min(ab)),
                 tolerance = 1e-3, info = paste(ab, collapse = ", "))
  }
})

test_that("bemg_fwhm measures a peak whose gaussian core is degenerate", {
  skip_on_cran()
  m <- list(center = 0, width = 0, height = 10, a = 38.9, b = 3.125)
  expect_equal(bemg_fwhm(m),
               .dense_fwhm(function(t) bemg(t, 0, 0, 10, 38.9, 3.125), 20),
               tolerance = 1e-3)
})

test_that("FWHM is comparable across fits, unlike sd", {
  skip_on_cran()
  # each model divides a peak between core and tail differently, so sigma is a
  # property of the parameterisation while the half-height width is a property
  # of the peak
  data(Sa_pr, envir = environment())
  fits <- c("gaussian", "egh", "bemg")
  pks <- lapply(fits, function(f)
    suppressWarnings(get_peaks(Sa_pr[1], lambdas = "210", fit = f,
                               estimate_purity = FALSE, show_progress = FALSE,
                               cl = 1))[[1]][[1]])
  names(pks) <- fits
  # the same peaks, so the comparison is peak by peak
  for (f in fits[-1]) expect_equal(pks[[f]]$start, pks[[1]]$start)
  spread <- function(col){
    m <- sapply(pks, `[[`, col)
    apply(m, 1, function(v) max(v) / min(v) - 1)
  }
  expect_gt(median(spread("sd")), 0.4)
  expect_lt(median(spread("FWHM")), 0.3)
})

test_that("a raw peak has no fitted sigma and a measured FWHM", {
  skip_on_cran()
  # two peaks 0.08 apart with sigma 0.03: the valley between them sits above
  # half height, so there is no crossing to measure unless a local baseline
  # brings the window edges down
  ts <- seq(0, 3, by = 0.005)
  y <- 1000 * exp(-(ts - 1.4)^2 / (2 * 0.03^2)) +
    900 * exp(-(ts - 1.48)^2 / (2 * 0.03^2))
  cl <- list(s1 = matrix(y, ncol = 1, dimnames = list(as.character(ts), "220")))
  raw <- function(b) suppressWarnings(
    get_peaks(cl, lambdas = "220", fit = "raw", baseline = b,
              estimate_purity = FALSE, show_progress = FALSE))[[1]][[1]]
  none <- raw("none")
  expect_true(all(is.na(none$sd)))
  expect_true(all(is.na(none$FWHM)))
  sloped <- raw("sloped")
  expect_true(all(is.na(sloped$sd)))
  expect_true(all(is.finite(sloped$FWHM)))
  # an unmeasurable width is not grounds for discarding the peak
  expect_equal(nrow(none), nrow(sloped))
  none_list <- suppressWarnings(
    get_peaks(cl, lambdas = "220", fit = "raw", baseline = "none",
              estimate_purity = FALSE, show_progress = FALSE))
  expect_equal(nrow(filter_peaks(none_list, max_fwhm = 1e-9)[[1]][[1]]),
               nrow(none))
})

test_that("one max_fwhm transfers across every fit", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  n <- vapply(c("gaussian", "egh", "bemg", "raw"), function(f)
    nrow(suppressWarnings(
      get_peaks(cl, lambdas = "220", fit = f, max_fwhm = 0.2,
                estimate_purity = FALSE, show_progress = FALSE))[[1]][[1]]),
    numeric(1))
  expect_true(all(n == n[["gaussian"]]))
})

test_that("max_fwhm is ignored when windows are supplied", {
  skip_on_cran()
  cl <- .mk_chrom(0.005)
  k <- suppressMessages(get_peaks(cl, lambdas = "220",
                                  windows = matrix(c(0.9, 1.1), nrow = 1),
                                  fit = "raw", max_fwhm = 1e-9,
                                  estimate_purity = FALSE, show_progress = FALSE))
  expect_equal(nrow(k[[1]][[1]]), 1)
})

test_that("filter_peaks can filter on converged and window_limited", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa_pr, lambdas = 210, fit = "egh",
                                    baseline = "sloped", sd_max = NULL,
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  n <- function(p) sum(vapply(p, function(s) nrow(s[["210"]]), numeric(1)))
  for (col in c("converged", "window_limited")){
    kept <- do.call(filter_peaks, c(list(pks), setNames(list(TRUE), col)))
    drop <- do.call(filter_peaks, c(list(pks), setNames(list(FALSE), col)))
    expect_equal(n(kept) + n(drop), n(pks), info = col)
    expect_true(all(unlist(lapply(kept, function(s) s[["210"]][[col]])) == 1),
                info = col)
    expect_true(all(unlist(lapply(drop, function(s) s[["210"]][[col]])) == 0),
                info = col)
  }
  both <- filter_peaks(pks, converged = TRUE, window_limited = FALSE)
  expect_lte(n(both), n(filter_peaks(pks, converged = TRUE)))
  expect_s3_class(both, "peak_list")

  expect_error(filter_peaks(pks, converged = "yes"), "TRUE or FALSE")
  expect_error(filter_peaks(pks, window_limited = NA), "TRUE or FALSE")
  expect_error(filter_peaks(pks, converged = c(TRUE, FALSE)), "TRUE or FALSE")
})

test_that("filter_peaks rejects fit diagnostics for a raw peak list", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  raw <- suppressWarnings(get_peaks(Sa_pr[1], lambdas = 210, fit = "raw",
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  expect_error(filter_peaks(raw, converged = TRUE), "no `converged` column")
  expect_error(filter_peaks(raw, window_limited = FALSE),
               "no `window_limited` column")
})

test_that("bemg rejects a fit whose tail grows away from the apex", {
  skip_on_cran()
  # a rising flank with no apex of its own: the bi-exponential can fit it with
  # a negative rate, which is a wall rather than a peak
  set.seed(1)
  x <- seq(0, 2, length.out = 60)
  y <- 200 * exp(-(x - 2.4)^2 / (2 * 0.5^2)) + rnorm(60, sd = 0.5)
  first <- suppressWarnings(fit_bemg_once(x, y, NULL, NULL, NULL, NULL, NULL,
                                          NULL, "none", 1000, FALSE))
  expect_false(first$converged)
  m <- suppressWarnings(fit_bemg(x, y, baseline = "none"))
  expect_gte(unname(m$a), 0)
  expect_gte(unname(m$b), 0)

  # the rejection is what makes the width measurable
  expect_true(is.finite(bemg_fwhm(m)))
})

test_that("reject_growing_tail passes through a well behaved fit", {
  skip_on_cran()
  x <- seq(0, 1, length.out = 201)
  y <- 1000 * exp(-(x - 0.5)^2 / (2 * 0.02^2))
  m <- suppressWarnings(fit_bemg(x, y, baseline = "none"))
  expect_gte(unname(m$a), 0)
  expect_gte(unname(m$b), 0)
  expect_true(m$converged)
})

test_that("bemg recovers a peak whose first fit cannot start", {
  skip_on_cran()
  y <- Sa_pr[["119"]][, "210"]
  ts <- get_times(Sa_pr, idx = 1)
  x <- which(ts >= 15.28 - 1e-9 & ts <= 15.58 + 1e-9)
  apex <- x[which.max(y[x])]
  first <- suppressWarnings(fit_bemg_once(x, y[x], apex, NULL, NULL, NULL,
                                          y[apex], NULL, "none", 1000, FALSE))
  expect_false(first$converged)
  expect_true(is.finite(sum(first$residual^2)))
  m <- suppressWarnings(fit_bemg(x, y[x], start.center = apex,
                                 start.height = y[apex]))
  expect_true(m$converged)
  expect_equal(max(m$y), unname(y[apex]), tolerance = 0.05)
})

test_that("find_peaks keeps its column schema when nothing is detected", {
  # `sapply` over an empty vector yields `list()`, which `data.frame` drops
  p <- find_peaks(rep(0, 100), bounds = TRUE)
  expect_equal(nrow(p), 0L)
  expect_named(p, c("pos", "lower", "upper", "prominence"))
})

test_that("fit_peaks rejects a `pos` that is not a find_peaks table", {
  x <- Sa_pr[[1]]
  expect_error(fit_peaks(x, lambda = "210", pos = c(50, 100)),
               regexp = "must be a data.frame with columns")
})

test_that("fit_peaks refuses `lambda = \"max\"` with a clear message", {
  expect_error(fit_peaks(Sa_pr[[1]], lambda = "max"),
               regexp = "must be a specific wavelength")
})

test_that("get_lambda_idx checks the wavelength of a single-column chromatogram", {
  expect_error(get_lambda_idx("999", lambdas = 210, allow_max = FALSE),
               regexp = "could not be found")
  expect_equal(get_lambda_idx("210", lambdas = 210, allow_max = FALSE), 1)
  expect_equal(get_lambda_idx("999", lambdas = NA_real_, allow_max = FALSE), 1)
})

test_that("print.peak_list reports no width filter by default", {
  pks <- get_peaks(Sa_pr[1], lambdas = "210", show_progress = FALSE)
  expect_output(print(pks), "max_fwhm: none")
})

test_that("print.peak_list counts a sample with no peaks as zero", {
  pks <- get_peaks(Sa_pr[1:2], lambdas = "210", show_progress = FALSE)
  n <- nrow(pks[[1]][["210"]])
  pks[[2]] <- structure(list(), names = character(0))
  expect_output(print(pks), sprintf("Total peaks: %d", n))
})

test_that("get_peaks reports widths and areas in `time_unit` whatever the axis unit", {
  skip_on_cran()
  data(Sa_pr)
  s <- convert_time_units(Sa_pr[1], time_unit = "s")
  first <- function(x, u) get_peaks(x, lambdas = "210", time_unit = u,
                                    show_progress = FALSE, cl = 1)[[1]][["210"]][1, ]
  ref <- first(Sa_pr[1], "min")
  for (u in c("min", "s")){
    p <- first(s, u)
    f <- time_factor(u)
    expect_equal(p$area, ref$area * f, tolerance = 1e-3, info = u)
    expect_equal(p$sd, ref$sd * f, tolerance = 1e-3, info = u)
  }
  expect_equal(attr(get_peaks(s, lambdas = "210", show_progress = FALSE, cl = 1),
                    "time_unit"), "s")
})

test_that("`tol` is in minutes whatever the axis unit", {
  skip_on_cran()
  data(Sa_pr)
  s <- convert_time_units(Sa_pr[1], time_unit = "s")
  m_rt <- get_peaks(Sa_pr[1], lambdas = lam_rt, targets = rts_2, tol = 0.05,
                    show_progress = FALSE, cl = 1)[[1]][[lam_rt]]$rt
  s_rt <- get_peaks(s, lambdas = lam_rt, targets = rts_2 * 60, tol = 0.05,
                    show_progress = FALSE, cl = 1)[[1]][[lam_rt]]$rt
  expect_gt(length(m_rt), 0)
  expect_equal(s_rt, m_rt * 60, tolerance = 1e-3)
  pks <- data.frame(rt = c(60, 70), start = c(58, 68), end = c(62, 72))
  expect_equal(peak_zoom_xlim(pks, rts = 66, tol = 0.1, pad = 0,
                              new.ts = 0:100, rtf = 60), c(58, 72))
  got <- NULL
  local_mocked_bindings(peak_zoom_xlim = function(pks, rts, tol, pad,
                                                  new.ts, rtf = 1){
    got <<- rtf
    NULL
  })
  pk_s <- get_peaks(s, lambdas = lam_rt, show_progress = FALSE, cl = 1)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(pk_s, chrom_list = s, lambda = lam_rt, targets = rts_2[1] * 60)
  expect_equal(got, 60)
  expect_equal(rt_factor("hours"), 1)
})

test_that("get_peaks warns when `tol` exceeds a tenth of the run", {
  skip_on_cran()
  data(Sa_pr)
  s <- convert_time_units(Sa_pr[1], time_unit = "s")
  expect_warning(get_peaks(s, lambdas = lam_rt, targets = rts_2 * 60, tol = 6,
                           show_progress = FALSE, cl = 1),
                 regexp = "`tol` is always in minutes")
  expect_warning(get_peaks(s, lambdas = lam_rt,
                           targets = data.frame(rt = rts_2 * 60, tol = c(0.05, 6)),
                           show_progress = FALSE, cl = 1),
                 regexp = "tolerance of 6 min")
  expect_no_warning(get_peaks(s, lambdas = lam_rt, targets = rts_2 * 60, tol = 0.1,
                              show_progress = FALSE, cl = 1))
})

test_that("filter_peaks keeps each sample's time unit", {
  skip_on_cran()
  data(Sa_pr)
  s <- lapply(Sa_pr[1:2], function(m){
    rownames(m) <- as.numeric(rownames(m)) * 60
    attr(m, "time_unit") <- "s"
    m
  })
  pks <- get_peaks(s, lambdas = 210, cl = 1, show_progress = FALSE)
  expect_equal(get_time_unit(filter_peaks(pks, min_height = 1)), "s")
})

test_that("missing entries in a `tol` column of `targets` take the default", {
  out <- check_targets(data.frame(rt = c(12, 13), tol = c(0.2, NA)), tol = 0.1)
  expect_equal(out$tol, c(0.2, 0.1))
  expect_error(check_targets(data.frame(rt = 12, tol = "wide"), tol = 0.1),
               "non-negative")
})

test_that("invalid `purity_args` fail once, before fitting", {
  data(Sa_pr)
  expect_error(get_peaks(Sa_pr[1], lambdas = 210, estimate_purity = TRUE,
                         purity_args = list(reference = "nonsense"), cl = 1,
                         show_progress = FALSE),
               "purity_args$reference", fixed = TRUE)
  expect_error(get_peaks(Sa_pr[1], lambdas = 210, estimate_purity = TRUE,
                         purity_args = list(refrence = "apex"), cl = 1,
                         show_progress = FALSE),
               "Unknown `purity_args`", fixed = TRUE)
})

test_that("r_squared is 1 - SS_res / SS_tot", {
  y <- c(1, 3, 2, 5, 4)
  expect_equal(r_squared(y, y), 1)
  expect_equal(r_squared(y, y + 2), 1 - 20 / 10)
  expect_true(is.na(r_squared(rep(2, 5), 1:5)))
})
