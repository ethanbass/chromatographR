# Tests for get_peaktable clustering correctness using synthetic chromatograms.
#
# Two samples, each with two Gaussian peaks.
# Peak A: same RT in both samples (gap = 0).
# Peak B: 50-index gap between samples; on this time grid (dt ≈ 0.0167 min/index)
#         that is ~0.83 min.

.n  <- 301
.ts <- seq(0, 5, length.out = .n)   # dt ≈ 0.0167 min/index

.gauss <- function(centers, heights) {
  x <- seq_len(.n); y <- numeric(.n)
  for (i in seq_along(centers))
    y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * 10^2))
  matrix(y, ncol = 1, dimnames = list(as.character(.ts), "220"))
}

.cl <- list(
  s1 = .gauss(centers = c(101, 201), heights = c(1000, 800)),
  s2 = .gauss(centers = c(101, 251), heights = c(1000, 800))
)

.pks <- get_peaks(.cl, lambdas = "220", fit = "gaussian", smooth_type = "none",
                  estimate_purity = FALSE, show_progress = FALSE)

test_that("get_peaktable carries metadata added to the chromatograms", {
  cl <- .cl
  attr(cl$s1, "dose") <- 1
  attr(cl$s2, "dose") <- 2
  attr(cl$s1, "added_metadata") <- attr(cl$s2, "added_metadata") <- "dose"
  pks <- get_peaks(cl, lambdas = "220", fit = "gaussian", smooth_type = "none",
                   estimate_purity = FALSE, show_progress = FALSE)
  pt <- get_peaktable(pks)
  expect_equal(pt$sample_meta$name, rownames(pt$tab))
  expect_equal(pt$sample_meta$dose, c(1, 2))

  pt <- attach_metadata(pt, data.frame(id = c("s2", "s1"), site = c("b", "a")),
                        column = "id")
  expect_equal(pt$sample_meta$site, c("a", "b"))
  expect_equal(pt$sample_meta$dose, c(1, 2))

  expect_false(is.data.frame(get_peaktable(.pks)$sample_meta))
})

# ── hmax = 0.2: nearby peaks merge, distant peaks stay separate ───────────────

test_that("get_peaktable clusters nearby peaks into one column and separates distant peaks", {
  skip_on_cran()
  pt <- get_peaktable(.pks, hmax = 0.2)

  # 3 distinct peaks total → 3 columns
  expect_equal(ncol(pt$tab), 3L)

  # The shared peak A column: both samples are non-zero
  rt_A  <- as.numeric(.ts[101])
  col_A <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_A))
  expect_gt(pt$tab["s1", col_A], 0)
  expect_gt(pt$tab["s2", col_A], 0)

  # Peak B lands in two separate columns, one per sample
  rt_B_s1  <- as.numeric(.ts[201])
  rt_B_s2  <- as.numeric(.ts[251])
  col_B_s1 <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_B_s1))
  col_B_s2 <- which.min(abs(as.numeric(pt$pk_meta["rt", ]) - rt_B_s2))
  expect_true(col_B_s1 != col_B_s2)
  expect_equal(pt$tab["s2", col_B_s1], 0)
  expect_equal(pt$tab["s1", col_B_s2], 0)
})

# ── hmax = 1.0: large enough to merge even the distant peaks ─────────────────

test_that("get_peaktable merges distant peaks when hmax is large enough", {
  skip_on_cran()
  pt <- get_peaktable(.pks, hmax = 1.0)

  # hmax = 1.0 min >> 0.83 min gap → all peaks collapse to 2 columns
  expect_equal(ncol(pt$tab), 2L)
  expect_true(all(pt$tab > 0))
})

test_that("get_peaktable rejects misspelled arguments", {
  skip_on_cran()
  expect_error(get_peaktable(.pks, hmx = 1.0), regexp = "unused argument.*hmx")
})

# ── metadata reported by only some peaks must survive summarization ──────────

test_that("weighted.mean summarization keeps metadata that only some peaks report", {
  skip_on_cran()
  data(Sa, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa, lambdas = "210", fit = "raw", cl = 1,
                                    show_progress = FALSE))
  fwhm <- unlist(lapply(pks, function(s) s[["210"]]$FWHM))
  expect_true(any(is.finite(fwhm)))
  expect_true(any(!is.finite(fwhm)))

  m <- get_peaktable(pks, summarize_by = "weighted.mean")$pk_meta
  summarized <- as.numeric(m["FWHM", ])
  kept <- summarized[is.finite(summarized)]
  expect_gt(length(kept), sum(!is.finite(summarized)))
  expect_true(all(kept >= min(fwhm, na.rm = TRUE)))
  expect_true(all(kept <= max(fwhm, na.rm = TRUE)))

  expect_true(all(!is.finite(as.numeric(m["sd", ]))))
})

test_that("weighted_mean_finite ignores unusable weights and values", {
  expect_equal(weighted_mean_finite(c(1, 3), c(1, 1)), 2)
  expect_equal(weighted_mean_finite(c(1, 3), c(1, NaN)), 1)
  expect_equal(weighted_mean_finite(c(1, NaN), c(1, 5)), 1)
  expect_equal(weighted_mean_finite(c(2, 4), c(0, 0)), 3)
  expect_equal(weighted_mean_finite(c(2, 4), c(-1, -1)), 3)
  expect_identical(weighted_mean_finite(c(NA_real_, NaN), c(1, 1)), NA_real_)
})

test_that("get_peaktable carries window_limited through to pk_meta", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa_pr, lambdas = 210, fit = "egh",
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  expect_true("window_limited" %in% colnames(pks[[1]][[1]]))
  expect_true(all(pks[[1]][[1]]$window_limited %in% c(0, 1)))
  expect_true("converged" %in% colnames(pks[[1]][[1]]))
  expect_true(all(pks[[1]][[1]]$converged %in% c(0, 1)))
  pt <- suppressWarnings(get_peaktable(pks))
  expect_true("window_limited" %in% rownames(pt$pk_meta))
  expect_true("converged" %in% rownames(pt$pk_meta))
  v <- pt$pk_meta["window_limited", ]
  expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
})

test_that("filter_peaktable can filter on min_converged", {
  skip_on_cran()
  data(Sa_pr, envir = environment())
  pks <- suppressWarnings(get_peaks(Sa_pr, lambdas = 210, fit = "egh",
                                    baseline = "sloped", sd_max = NULL,
                                    estimate_purity = FALSE, cl = 1,
                                    show_progress = FALSE))
  pt <- suppressWarnings(get_peaktable(pks))
  cv <- as.numeric(pt$pk_meta["converged", ])
  for (th in c(0, 0.8, 1)){
    f <- suppressWarnings(filter_peaktable(pt, min_converged = th))
    expect_equal(ncol(f$tab), sum(cv >= th), info = paste("th =", th))
    expect_equal(ncol(f$pk_meta), ncol(f$tab), info = paste("th =", th))
    expect_s3_class(f, "peak_table")
  }
  expect_lte(ncol(suppressWarnings(
    filter_peaktable(pt, min_converged = 1))$tab), ncol(pt$tab))

  for (bad in list(1.5, -0.1, "x", NA, c(0.5, 0.6))){
    expect_error(filter_peaktable(pt, min_converged = bad),
                 "between 0 and 1")
  }
})

test_that("filter_peaktable rejects min_converged without a converged row", {
  skip_on_cran()
  data(pk_tab, envir = environment())
  if (!"converged" %in% rownames(pk_tab$pk_meta)){
    expect_error(filter_peaktable(pk_tab, min_converged = 0.5),
                 "no `converged` row")
  } else {
    skip("stored pk_tab already carries a converged row")
  }
})

test_that("filter_peaktable matches one peak per requested retention time", {
  data(pk_tab, envir = environment())
  rt <- as.numeric(pk_tab$pk_meta["rt", ])
  one <- filter_peaktable(pk_tab, rts = rt[2], tol = 5)
  expect_equal(ncol(one$tab), 1)
  expect_equal(as.numeric(one$pk_meta["rt", ]), rt[2])

  two <- filter_peaktable(pk_tab, rts = rt[c(2, 10)], tol = 5)
  expect_equal(ncol(two$tab), 2)
  expect_setequal(as.numeric(two$pk_meta["rt", ]), rt[c(2, 10)])

  expect_warning(none <- filter_peaktable(pk_tab, rts = 99, tol = 0.01),
                 regexp = "not identified in the peak table")
  expect_equal(ncol(none$tab), 0)
})

test_that("filter_peaktable rejects `rts` combined with `min_rt`/`max_rt`", {
  data(pk_tab, envir = environment())
  rt <- as.numeric(pk_tab$pk_meta["rt", ])
  expect_error(filter_peaktable(pk_tab, rts = rt[2], min_rt = 10),
               "not both")
  expect_error(filter_peaktable(pk_tab, rts = rt[2], max_rt = 10),
               "not both")
})

test_that("filter_peaktable keeps ref_spectra aligned with tab", {
  data(pk_tab, envir = environment())
  pk_tab$ref_spectra <- as.data.frame(matrix(
    seq_len(3 * ncol(pk_tab$tab)), 3, dimnames = list(NULL, colnames(pk_tab$tab))))
  f <- filter_peaktable(pk_tab, min_rt = 12)
  expect_lt(ncol(f$tab), ncol(pk_tab$tab))
  expect_equal(colnames(f$ref_spectra), colnames(f$tab))
})

test_that("filter_peaks retains peaks whose flag is NA", {
  pks <- structure(list(s = list("220" = data.frame(
    rt = 1:4, FWHM = 0.1, height = 10, area = 10,
    converged = c(1, 0, NA, 1), window_limited = c(0, 1, NA, 0)))),
    class = c("peak_list", "list"))
  expect_equal(filter_peaks(pks, converged = TRUE)$s[["220"]]$rt, c(1, 3, 4))
  expect_equal(filter_peaks(pks, converged = FALSE)$s[["220"]]$rt, c(2, 3))
  expect_equal(filter_peaks(pks, window_limited = FALSE)$s[["220"]]$rt,
               c(1, 3, 4))
})

test_that("filter calls are recorded as code that repeats them", {
  f <- suppressWarnings(filter_peaks(.pks, min_fwhm = 0.05))
  f <- suppressWarnings(filter_peaks(f, converged = TRUE, max_rt = 2))
  expect_equal(attr(f, "filters"), c("filter_peaks(min_fwhm = 0.05)",
                                     "filter_peaks(max_rt = 2, converged = TRUE)"))
  expect_output(print(f), "Filtered with: filter_peaks\\(min_fwhm = 0.05\\)")
  pt <- get_peaktable(f)
  expect_equal(pt$args$filters, paste(attr(f, "filters"), collapse = "; "))
  pt <- filter_peaktable(pt, min_rt = 0)
  expect_match(pt$args$filters, "; filter_peaktable\\(min_rt = 0\\)$")
  expect_true(is.na(get_peaktable(.pks)$args$filters))
  expect_equal(filter_peaktable(pk_tab, min_rt = 12)$args$filters,
               "filter_peaktable(min_rt = 12)")
})

test_that("peak tables record how their peaks were fitted", {
  pt <- get_peaktable(.pks)
  expect_equal(pt$args[c("fit", "baseline", "max_iter")],
               list(fit = "gaussian", baseline = "none", max_iter = 1000))
  expect_null(attr(pt, "windows"))
  expect_null(attr(pt, "targets"))
  w <- rbind(c(1.2, 2.1), c(2.9, 3.8))
  pw <- get_peaktable(get_peaks(.cl, lambdas = "220", windows = w, fit = "raw",
                                estimate_purity = FALSE, show_progress = FALSE))
  again <- get_peaktable(get_peaks(.cl, lambdas = "220",
                                   windows = attr(pw, "windows"), fit = "raw",
                                   estimate_purity = FALSE,
                                   show_progress = FALSE))
  expect_equal(again$tab, pw$tab)
  expect_equal(attr(filter_peaktable(pw, max_rt = 2), "windows"),
               attr(pw, "windows"))
  tg <- data.frame(rt = c(1.67, 3.33), name = c("a", "b"))
  pt <- get_peaktable(get_peaks(.cl, lambdas = "220", targets = tg,
                                fit = "gaussian", estimate_purity = FALSE,
                                show_progress = FALSE))
  expect_equal(attr(pt, "targets")$name, c("a", "b"))
  for (tg in list(tg, c(1.67, 3.33))){
    pt <- get_peaktable(get_peaks(.cl, lambdas = "220", targets = tg,
                                  fit = "gaussian", estimate_purity = FALSE,
                                  show_progress = FALSE))
    again <- get_peaktable(get_peaks(.cl, lambdas = "220",
                                     targets = attr(pt, "targets"),
                                     fit = "gaussian", estimate_purity = FALSE,
                                     show_progress = FALSE))
    expect_equal(again$tab, pt$tab)
  }
})

test_that("write_peaktable writes args with several values as one field", {
  pt <- get_peaktable(suppressWarnings(
    filter_peaks(.pks, min_fwhm = 0.05, max_rt = 2)))
  pt$args$lambdas <- list(c("210", "260"))
  d <- withr::local_tempdir()
  write_peaktable(pt, d, what = "args")
  a <- utils::read.csv(file.path(d, "peak_table-args.csv"), row.names = 1)
  expect_equal(ncol(a), 1)
  expect_equal(a["lambdas", 1], "210, 260")
  expect_equal(a["filters", 1], pt$args$filters)
})

test_that("summary of a peak table reports the retention time unit", {
  cl <- lapply(.cl, function(m) structure(m, time_unit = "Minutes"))
  pks <- get_peaks(cl, lambdas = "220", fit = "gaussian", smooth_type = "none",
                   estimate_purity = FALSE, show_progress = FALSE)
  pt <- get_peaktable(pks)
  expect_equal(pt$args$rt_unit, "min")
  expect_output(print(summary(pt)), "\\] min")
  expect_output(print(summary(get_peaktable(.pks))), "\\]\n")
})

# ── corrected retention times resolve to aligned chromatograms ────────────────

.aligned_fixture <- function(alg = "ptw") {
  models <- suppressWarnings(correct_rt(Sa_pr, lambdas = 210, alg = alg,
                                        what = "models", show_progress = FALSE))
  warped <- suppressWarnings(correct_rt(Sa_pr, models = models, alg = alg,
                                        lambdas = 210, show_progress = FALSE))
  pks <- suppressWarnings(get_peaks(Sa_pr, lambdas = 210, show_progress = FALSE))
  list(models = models, warped = warped, pks = pks)
}

test_that("a peak table from corrected times rebuilds the aligned chromatograms", {
  skip_on_cran()
  for (alg in c("ptw", "vpdtw")) {
    f <- .aligned_fixture(alg)
    pt <- get_peaktable(correct_peaks(f$pks, f$models, Sa_pr))
    expect_true(pt$args$use_corrected_rt)
    expect_equal(unclass(get_chrom_list(pt)), unclass(f$warped),
                 ignore_attr = TRUE)
  }
})

test_that("stored ptw models keep only the warp", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  pt <- get_peaktable(correct_peaks(f$pks, f$models, Sa_pr))
  stored <- attr(pt, "warp_models")
  expect_s3_class(stored, "ptw_list")
  expect_setequal(names(stored[[1]]), c("warp.coef", "mode"))
})

test_that("stored VPdtw models keep only the warp", {
  skip_on_cran()
  f <- .aligned_fixture("vpdtw")
  pt <- get_peaktable(correct_peaks(f$pks, f$models, Sa_pr))
  stored <- attr(pt, "warp_models")
  expect_s3_class(stored, "VPdtw")
  expect_setequal(names(stored), c("xVals", "shift", "samples", "n_timepoints"))
  expect_equal(attr(stored, "reference"), attr(f$models, "reference"))
  expect_equal(attr(stored, "parameters"), attr(f$models, "parameters"))
})

test_that("a table of some samples rebuilds them from models fitted to all", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  pks <- suppressWarnings(get_peaks(Sa_pr[1:2], lambdas = 210,
                                    show_progress = FALSE))
  pt <- get_peaktable(correct_peaks(pks, f$models, Sa_pr))
  expect_equal(unclass(get_chrom_list(pt)), unclass(f$warped[1:2]),
               ignore_attr = TRUE)
})

test_that("corrected peak tables store corrected times in the rt rows", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  cp <- correct_peaks(f$pks, f$models, Sa_pr)
  pt <- get_peaktable(cp)
  raw <- get_peaktable(cp, use_corrected_rt = FALSE)
  expect_equal(rownames(pt$pk_meta), rownames(raw$pk_meta))
  expect_false(isTRUE(all.equal(pt$pk_meta["rt", ], raw$pk_meta["rt", ])))
  aligned <- get_chrom_list(pt)
  win <- get_retention_idx(unlist(pt$pk_meta[c("start", "end"), "V2"]),
                           get_times(aligned))
  coords <- get_peak_coordinates(pt, aligned, loc = "V2", idx = "max",
                                 lambda = "210")
  expect_equal(coords$coordinates[[2]],
               max(aligned[[which.max(pt$tab$V2)]][win[1]:win[2], "210"]))
  long <- reshape_peaktable(pt, peaks = "V2")
  expect_equal(colnames(long), c("peak", "sample", "lambda", "rt", "area"))
  expect_equal(unique(long$rt), pt$pk_meta["rt", "V2"])
})

test_that("use.cor still works, and old peak tables recording it stay aligned", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  cp <- correct_peaks(f$pks, f$models, Sa_pr)
  expect_warning(old <- get_peaktable(cp, use.cor = FALSE),
                 "`use.cor` argument is deprecated")
  expect_equal(old, get_peaktable(cp, use_corrected_rt = FALSE))

  pt <- get_peaktable(cp)
  names(pt$args)[names(pt$args) == "use_corrected_rt"] <- "use.cor"
  expect_equal(get_chrom_list(pt), get_chrom_list(get_peaktable(cp)))
})

test_that("a named aligned chrom_list takes precedence over the models", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  assign("aligned_test_chroms", f$warped, envir = globalenv())
  withr::defer(rm("aligned_test_chroms", envir = globalenv()))
  cp <- correct_peaks(f$pks, f$models, Sa_pr,
                      chrom_list_aligned = aligned_test_chroms)
  pt <- get_peaktable(cp)
  expect_equal(pt$args$chrom_list_aligned, "aligned_test_chroms")
  attr(pt, "warp_models") <- NULL
  expect_identical(get_chrom_list(pt), f$warped)
})

test_that("uncorrected times resolve to the unaligned chromatograms", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  pt <- get_peaktable(correct_peaks(f$pks, f$models, Sa_pr), use_corrected_rt = FALSE)
  expect_identical(get_chrom_list(pt), Sa_pr)
})

test_that("corrected times without models or aligned chromatograms warn", {
  skip_on_cran()
  f <- .aligned_fixture("ptw")
  pt <- get_peaktable(correct_peaks(f$pks, f$models, Sa_pr))
  attr(pt, "warp_models") <- NULL
  expect_warning(out <- get_chrom_list(pt), "unaligned `Sa_pr`")
  expect_identical(out, Sa_pr)
})

test_that("filtering a table of corrected peaks does not add rows", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = "210", cl = 1, show_progress = FALSE)
  mods <- correct_rt(Sa_pr, lambdas = "210", reference = 1, what = "models",
                     cl = 1, show_progress = FALSE)
  pt <- get_peaktable(correct_peaks(pks, mods, chrom_list = Sa_pr),
                      chrom_list = Sa_pr)
  filtered <- filter_peaktable(pt, min_rt = 11, max_rt = 13)
  expect_equal(rownames(filtered$pk_meta), rownames(pt$pk_meta))
  rt_cor <- unlist(filtered$pk_meta["rt", ])
  expect_true(all(rt_cor >= 11 & rt_cor <= 13))
})

test_that("hmax, sigma_t and combine_peaks' tol are minutes whatever the time unit", {
  skip_on_cran()
  data(Sa_pr)
  to_seconds <- function(pks){
    for (i in seq_along(pks)){
      for (l in names(pks[[i]]))
        pks[[i]][[l]][c("rt", "start", "end")] <- pks[[i]][[l]][c("rt", "start", "end")] * 60
      attr(pks[[i]], "time_unit") <- "s"
    }
    pks
  }
  pm <- get_peaks(Sa_pr, lambdas = "210", show_progress = FALSE, cl = 1)
  ps <- to_seconds(pm)
  withr::local_seed(1)
  jittered <- pm
  for (i in seq_along(jittered))
    jittered[[i]][["210"]]$rt <- jittered[[i]][["210"]]$rt +
      runif(nrow(jittered[[i]][["210"]]), 0, 1e-4)
  expect_equal(unname(as.matrix(get_peaktable(to_seconds(jittered))$tab)),
               unname(as.matrix(get_peaktable(jittered)$tab)))
  sp_m <- get_peaktable(pm, chrom_list = Sa_pr, clust = "sp.rt", sigma_t = 0.1)
  sp_s <- get_peaktable(ps, chrom_list = convert_time_units(Sa_pr, "s"),
                        clust = "sp.rt", sigma_t = 0.1)
  expect_equal(ncol(sp_s$tab), ncol(sp_m$tab))
  expect_warning(
    old <- get_peaktable(ps, chrom_list = convert_time_units(Sa_pr, "s"),
                         clust = "sp.rt", sigma.t = 6),
    "`sigma.t` argument is deprecated")
  expect_equal(old$tab, sp_s$tab)
  tm <- attach_ref_spectra(get_peaktable(pm), chrom_list = Sa_pr)
  ts <- tm
  ts$pk_meta["rt", ] <- as.numeric(tm$pk_meta["rt", ]) * 60
  ts$args$rt_unit <- "s"
  merged <- colnames(combine_peaks(tm, tol = 0.3, min_cor = 0)$tab)
  expect_lt(length(merged), ncol(tm$tab))
  expect_equal(colnames(combine_peaks(ts, tol = 0.3, min_cor = 0)$tab), merged)
})

test_that("get_peaktable errors clearly when no target is found", {
  skip_on_cran()
  data(Sa_pr)
  p <- suppressWarnings(get_peaks(Sa_pr, lambdas = 220, targets = 11.011, tol = 1e-4,
                                  cl = 1, show_progress = FALSE))
  expect_error(get_peaktable(p), "None of the targets was found in any sample")
})
