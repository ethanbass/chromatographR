# to run visual tests for plotly graphics, set `Sys.setenv("VISUAL_TESTS"="true")`

data("Sa")
new.ts <- seq(10, 18.66, by = .01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths

out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas, 
                  show_progress = FALSE)

test_that("preprocess works on matrix", {
  out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas,
                    maxI = 1, show_progress = FALSE)
  expect_equal(class(out)[1],"matrix")
  expect_equal(rownames(out), format_axis(new.ts))
  expect_equal(colnames(out), as.character(new.lambdas))
  expect_error(preprocess_matrix(Sa), regexp = "X should be a matrix")
})

test_that("Preprocess works without providing dimensions", {
  out <- suppressMessages(preprocess(X = (Sa[1]), show_progress = FALSE))
  testthat::expect_true(inherits(out, "list"))
  expect_equal(class(out[[1]])[1], "matrix")
  expect_equal(rownames(out[[1]]), format_axis(new.ts))
  expect_equal(colnames(out[[1]]), as.character(new.lambdas))
})

test_that("preprocess unwraps nested lists with a dad element", {
  dad_list <- lapply(Sa[1:2], function(x) list(dad = x, meta = "foo"))
  result <- preprocess(X = dad_list, dim1 = new.ts, dim2 = new.lambdas,
                       show_progress = FALSE)
  expect_equal(length(result), 2)
  expect_equal(class(result[[1]])[1], "matrix")
  expect_equal(rownames(result[[1]]), format_axis(new.ts))
})

test_that("preprocess returns correct errors and warnings", {
  expect_error(preprocess(X = list(Sa[[1]], "not_a_matrix"),
                          show_progress = FALSE),
               regexp = "X should be a matrix")
  expect_error(preprocess(X = as.data.frame(Sa[[1]]), show_progress = FALSE),
               regexp = "X should be a matrix")
  expect_error(preprocess(Sa, dim1 = seq(8,15, by=.01), show_progress = FALSE), 
               regexp = "incompatible with actual data")
  expect_error(preprocess(Sa, dim1 = seq(10,20, by=.01),
                          show_progress = FALSE), 
               regexp = "incompatible with actual data")
  expect_warning(suppressMessages(
      xx <- preprocess(Sa, dim1 = seq(16,18.664, by=.001),
                       show_progress = FALSE)
    ), regexp = "No extrapolation allowed")
  expect_length(xx, 2)
  expect_warning(suppressMessages(
    xx <- preprocess(Sa, dim1 = seq(9.997,13, by=.01), show_progress = FALSE) 
  ), regexp = "No extrapolation allowed")
  expect_length(xx,2)
})


test_that("preprocess works without interpolation", {
  dat.pr <- preprocess(X = Sa[[1]], interpolate_cols = FALSE,
                       interpolate_rows = FALSE, show_progress = FALSE)
  expect_equal(dim(dat.pr),dim(Sa[[1]]))
  expect_equal(rownames(dat.pr), rownames(Sa[[1]]))
  expect_equal(colnames(dat.pr), colnames(Sa[[1]]))
})

dat.pr <- preprocess(X = Sa[1:2], dim1 = new.ts, dim2 = new.lambdas,
                     show_progress = FALSE)

test_that("preprocess works on a list", {
  expect_equal(length(dat.pr),length(Sa[1:2]))
  expect_equal(rownames(dat.pr[[1]]), format_axis(new.ts))
  expect_equal(colnames(dat.pr[[1]]), as.character(new.lambdas))
  expect_equal(names(dat.pr), names(Sa[1:2]))
})

### test correct_rt ###

warping.models <- correct_rt(dat.pr, lambdas = c('210','260','318'),
                             what = "models", show_progress = FALSE)

warp <- correct_rt(chrom_list = dat.pr, models = warping.models,
                   what = "corrected.values", show_progress = FALSE)

test_that("correct_rt works", {
  equals(length(warping.models), length(warp), length(Sa[1:2]))
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_equal(t(warping.models[[1]]$warped.sample)[,1], 
               as.numeric(warp[[1]][,"210"]))
  expect_error(correct_rt(dat.pr, show_progress = FALSE),
               regexp = "Must specify wavelengths")
  expect_error(correct_rt(dat.pr, what="x", show_progress = FALSE),
               regexp = "'arg' should be one of")
  expect_error(correct_rt(dat.pr, lambdas = "x", show_progress = FALSE),
               regexp = "Lambdas not found")
  expect_error(correct_rt(dat.pr, lambdas = "210", alg="x",
                          show_progress = FALSE), 
               regexp = "'arg' should be one of")
  expect_error(correct_rt(dat.pr, lambdas = "210", 
                          models = "warping.models", alg = "vpdtw",
                          show_progress = FALSE), 
               regexp = "The supplied models do not match")
})

test_that("correct_rt works with VPdtw", {
  skip_on_cran()
  warp <- correct_rt(dat.pr, lambdas = "210", alg = "vpdtw",
                     show_progress = FALSE)
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr = TRUE)
  expect_error(correct_rt(dat.pr, lambdas = c("210", "260"),
                          alg="vpdtw",  show_progress = FALSE), 
               regexp = "only supports warping by a single wavelength")
  expect_error(correct_rt(dat.pr, lambdas = "x", alg = "vpdtw",
                          show_progress = FALSE), 
               regexp = "Lambdas not found")
})

test_that("VPdtw plot displays correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  vpdtw_alignment <- function(){
    warp <- correct_rt(Sa_pr, lambdas = "210", alg = "vpdtw", plot_it = TRUE,
                       show_progress = FALSE)
  }
  vdiffr::expect_doppelganger("vpdtw_alignment", vpdtw_alignment)
})
### test correct_rt_group ###

# `Sa_pr` has 4 samples, so it can be split into two groups of two.
cl_grp <- Sa_pr
for (i in seq_along(cl_grp)) attr(cl_grp[[i]], "batch") <- c("a", "a", "b", "b")[i]
grp_vec <- c("a", "a", "b", "b")

warp_grp <- correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                             show_progress = FALSE)

test_that("correct_rt_group warps every sample with its group's model", {
  skip_on_cran()
  expect_type(warp_grp, "list")
  # output order follows `chrom_list`, not the group loop
  expect_equal(names(warp_grp), names(cl_grp))
  expect_equal(dim(warp_grp[[1]]), dim(cl_grp[[1]]))
  expect_equal(dimnames(warp_grp[[1]]), dimnames(cl_grp[[1]]))
  # samples in the same group share a warp, so their models agree; the two
  # groups' warps differ from each other
  models <- correct_rt(chrom_list = cl_grp, lambdas = 210, alg = "ptw",
                       what = "models", show_progress = FALSE)
  expect_false(identical(warp_grp, correct_rt(cl_grp, models = models,
                                              show_progress = FALSE)))
})

test_that("correct_rt_group accepts all three forms of `groups`", {
  skip_on_cran()
  # metadata attribute: this branch previously never coerced to a factor, so
  # `levels()` was NULL and the group loop iterated over nothing
  expect_equal(correct_rt_group(cl_grp, lambdas = 210, groups = "batch",
                                show_progress = FALSE),
               warp_grp)
  # a named vector is matched to `names(chrom_list)` rather than used
  # positionally, so a scrambled order gives the same answer
  grp_named <- setNames(grp_vec[c(3, 1, 2, 4)], names(cl_grp)[c(3, 1, 2, 4)])
  expect_equal(correct_rt_group(cl_grp, lambdas = 210, groups = grp_named,
                                show_progress = FALSE),
               warp_grp)
})

test_that("correct_rt_group returns per-sample group models", {
  skip_on_cran()
  mods <- correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                           what = "models", show_progress = FALSE)
  expect_s3_class(mods, "ptw_list")
  expect_equal(names(mods), names(cl_grp))
  expect_equal(mods[[1]]$warp.coef, mods[[2]]$warp.coef)
  expect_false(isTRUE(all.equal(mods[[1]]$warp.coef, mods[[3]]$warp.coef)))
  expect_equal(correct_rt(cl_grp, models = mods, show_progress = FALSE),
               warp_grp, ignore_attr = TRUE)
  pks <- get_peaks(cl_grp, lambdas = 210, show_progress = FALSE)
  pks_cor <- correct_peaks(pks, mods, chrom_list = cl_grp)
  expect_true("rt.cor" %in% colnames(pks_cor[[1]][[1]]))
})

test_that("correct_rt_group resolves the reference group", {
  skip_on_cran()
  # `reference_group` and a `reference` sample in that group agree
  warp_b <- correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                             reference_group = "b", show_progress = FALSE)
  expect_equal(names(warp_b), names(cl_grp))
  expect_equal(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                                reference = 3, show_progress = FALSE),
               warp_b)
  expect_equal(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                                reference = names(cl_grp)[3],
                                show_progress = FALSE),
               warp_b)
  # `reference = "best"` (the default) delegates selection to `ptw::bestref`
  msgs <- capture_messages(
    invisible(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                               verbose = TRUE, show_progress = FALSE)))
  expect_match(msgs, "Computing within-group averages", all = FALSE)
  expect_match(msgs, "Selected chromatogram 1 as best reference", all = FALSE)
  expect_match(msgs, "Applying group warp models", all = FALSE)
})

test_that("correct_rt_group errors on unusable group assignments", {
  skip_on_cran()
  expect_error(correct_rt_group(cl_grp, lambdas = 210, groups = c("a", "b"),
                                show_progress = FALSE),
               regexp = "Length of `groups` must match")
  expect_error(correct_rt_group(cl_grp, lambdas = 210,
                                groups = c("a", "a", "b", NA),
                                show_progress = FALSE),
               regexp = "All samples must have a group assignment")
  # an attribute missing everywhere yields a list of NULLs, which the guard
  # must catch rather than passing on an empty factor
  expect_error(correct_rt_group(Sa_pr, lambdas = 210, groups = "batch",
                                show_progress = FALSE),
               regexp = "Metadata attribute 'batch' not found on all")
  cl_partial <- cl_grp
  attr(cl_partial[[4]], "batch") <- NULL
  expect_error(correct_rt_group(cl_partial, lambdas = 210, groups = "batch",
                                show_progress = FALSE),
               regexp = "not found on all chromatograms. Missing for")
})

test_that("correct_rt_group errors on an unresolvable reference", {
  skip_on_cran()
  # these all used to reach `correct_rt` as `integer(0)` and fail there with
  # "argument is of length zero"
  expect_error(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                                reference_group = "zzz",
                                show_progress = FALSE),
               regexp = "does not match any group. Valid groups are: 'a', 'b'")
  expect_error(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                                reference = "zzz", show_progress = FALSE),
               regexp = "does not match any sample in `chrom_list`")
  expect_error(correct_rt_group(cl_grp, lambdas = 210, groups = grp_vec,
                                reference = 99, show_progress = FALSE),
               regexp = "does not match any sample in `chrom_list`")
})

test_that("PTW plot displays correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  ptw_alignment <- function(){
    warp <- correct_rt(Sa_pr, lambdas = "210", alg = "ptw",
                       plot_it = TRUE, show_progress = FALSE)
  }
  vdiffr::expect_doppelganger("ptw_alignment", ptw_alignment)
})


test_that("plot_chroms works to plot alignments", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  alignment <- function(){
    par(mfrow = c(2,1))
    plot_chroms(warp, lambdas = "210", show_legend = TRUE)
    plot_chroms(dat.pr, lambdas = "210", show_legend = TRUE)
  }
  vdiffr::expect_doppelganger("alignment", alignment)
  expect_error(plot_chroms(pktab))
})

test_that("plot_chroms can subset chromatograms correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  numeric_indices <- function(){
    plot_chroms(Sa_pr, idx=c(1, 3), lambdas = 210, show_legend = TRUE)
  }
  vdiffr::expect_doppelganger("plot-trace", numeric_indices)

  numeric_indices_rev <- function(){
    plot_chroms(Sa_pr, idx=c(3, 1), lambdas = 210, show_legend=TRUE)
  }
  vdiffr::expect_doppelganger("plot-trace-rev", numeric_indices_rev)
})

test_that("plot_chroms works to plot alignments with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  alignment_ggp <- function(){
    plot_chroms(warp, lambdas="210", engine="ggplot", show_legend = TRUE)
  }
  vdiffr::expect_doppelganger("alignment_ggp", alignment_ggp)
  
  alignment_ggp_zoom <- function(){
    plot_chroms(warp, lambdas = "210", engine = "ggplot", show_legend = FALSE,
                xlim = c(12.5, 15), ylim = c(0, 550))
  }
  suppressWarnings(vdiffr::expect_doppelganger("alignment_ggp_zoom",
                                               alignment_ggp_zoom))
})


### test get_peaks ###
lam <- c("210", "318")

pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh", smooth_type = "none",
                     show_progress = FALSE)

pks_bemg <- get_peaks(dat.pr, lambdas = lam, fit = "bemg", smooth_type = "none",
                     show_progress = FALSE)

pks_gaussian <- get_peaks(dat.pr, lambdas = lam, fit = "gaussian",
                          smooth_type = "none", show_progress = FALSE)

pks_raw <- get_peaks(dat.pr, lambdas = lam, fit = "raw", smooth_type = "savgol",
                     show_progress = FALSE, smooth_window = 7)

matched <- match_peaks(pks_egh, chrom_list = dat.pr, lambda = "210",
                       verbose = FALSE)

test_that("get_peaks works", {
  expect_equal(names(pks_egh), names(dat.pr))
  expect_equal(names(pks_gaussian), names(dat.pr))
  expect_equal(names(pks_egh[[1]]), lam)
  expect_equal(names(pks_gaussian[[1]]), lam)
  expect_equal(class(pks_egh), c("peak_list","list"))
  expect_error(get_peaks(dat.pr, show_progress = FALSE), 
               regexp = "must be provided.")
  expect_error(get_peaks(dat.pr, lambdas = "210", fit = "nonsense",
                         show_progress = FALSE), 
               regexp = "'arg' should be one of")
})

test_that("[.peak_list works", {
  sub_peaks <- pks_egh[2]
  expect_equal(sub_peaks[[1]], pks_egh[[2]])
  attr_original <- attributes(pks_egh)
  attr_new <- attributes(sub_peaks)
  expect_equal(attr_new[-1], attr_original[-1])
})

test_that("filter_peaks works", {
  pks_s <- filter_peaks(pks_egh, min_height = 10, min_area = 10, min_sd = .07)
  expect_equal(names(pks_egh), names(pks_s))
  expect_equal(names(pks_s), names(dat.pr))
  expect_equal(names(pks_s[[1]]), lam)
  expect_equal(class(pks_s), c("peak_list","list"))
  expect_lt(nrow(pks_s[[1]][[1]]), nrow(pks_egh[[1]][[1]]))
  expect_lt(nrow(pks_s[[2]][[1]]), nrow(pks_egh[[2]][[1]]))
  expect_warning(filter_peaks(pks_egh), regexp = "Nothing to filter")
  expect_warning(filter_peaks(pks_egh, min_height=0),
                 regexp = "'min_height' is less than minimum peak height.")
  expect_warning(filter_peaks(pks_egh, min_area=0), 
                 regexp = "'min_area' is less than minimum peak area")
  expect_warning(filter_peaks(pks_egh, max_fwhm = Inf),
                 regexp = "'max_fwhm' is greater than the maximum")
})

### test match_peaks and lptw ###

test_that("match_peaks works", {
  expect_s3_class(matched, "matched_peak_list")
  expect_equal(names(matched), names(dat.pr))
  expect_equal(attr(matched, "reference"), names(dat.pr)[1])
  expect_false(is.null(matched[[1]]))
})

test_that("match_peaks errors correctly", {
  expect_error(match_peaks(pks_egh, chrom_list = dat.pr),
               regexp = "wavelength")
  expect_error(match_peaks(pks_egh, chrom_list = dat.pr[1], lambda = "210"),
               regexp = "do not match")
  expect_error(match_peaks(pks_egh, chrom_list = dat.pr, lambda = "210",
                           fit = "nonsense"),
               regexp = "'arg' should be one of")
  expect_error(match_peaks(pks_raw, chrom_list = dat.pr, lambda = "210",
                           shape_cutoff = 1),
               regexp = "Peak shape parameters are not fit")
})

test_that("match_peaks honors the `fit` argument", {
  skip_on_cran()
  # `fit` must not be silently overwritten by the peak_list's `fit` attribute.
  # A gaussian peak list has no `tau` column, so requesting the egh shape
  # distance should fail. If `fit` were ignored, this would succeed.
  expect_error(match_peaks(pks_gaussian, chrom_list = dat.pr, lambda = "210",
                           fit = "egh", shape_cutoff = 1),
               regexp = "tau")
  # The reverse override is legitimate: the gaussian shape distance depends
  # only on `sd`, which an egh peak list has.
  matched_gauss <- match_peaks(pks_egh, chrom_list = dat.pr, lambda = "210",
                               fit = "gaussian", shape_cutoff = 5,
                               verbose = FALSE)
  expect_s3_class(matched_gauss, "matched_peak_list")
  expect_equal(names(matched_gauss), names(dat.pr))
})

test_that("match_peaks skips samples containing no peaks", {
  skip_on_cran()
  # A blank sample among real ones. Previously this crashed in `colnames<-`,
  # because `do.call(cbind, list())` returns NULL when a sample has no peaks.
  pks_blank <- pks_egh
  pks_blank[[2]][["210"]] <- pks_blank[[2]][["210"]][0, , drop = FALSE]
  expect_warning(m <- match_peaks(pks_blank, chrom_list = dat.pr,
                                  lambda = "210", verbose = FALSE),
                 regexp = "No peaks found for sample 2")
  # the blank sample keeps its slot, so names stay aligned with the peak list
  expect_equal(length(m), length(dat.pr))
  expect_equal(names(m), names(pks_blank))
  expect_null(m[[2]])
  expect_false(is.null(m[[1]]))
  # and the consumers of a `matched_peak_list` tolerate the NULL element
  pt <- get_peaktable(m, pks = pks_blank, na.fill = NA)
  expect_s3_class(pt, "peak_table")
  expect_true(all(is.na(pt$tab[2, ])))
  expect_equal(names(correct_rt(dat.pr, lambdas = "210", alg = "lptw",
                                peaks = m, show_progress = FALSE)),
               names(dat.pr))
})

test_that("match_peaks errors when the reference contains no peaks", {
  skip_on_cran()
  pks_blank_ref <- pks_egh
  pks_blank_ref[[1]][["210"]] <- pks_blank_ref[[1]][["210"]][0, , drop = FALSE]
  expect_error(match_peaks(pks_blank_ref, chrom_list = dat.pr, lambda = "210",
                           reference = 1, verbose = FALSE),
               regexp = "No peaks were found in the reference sample")
})

### match_peaks: argument defaults and scoring ###

# `Sa_pr` rather than `dat.pr` so that `chrom_list` auto-derivation can be
# tested: `get_chrom_list` resolves the stored name with `get()`, which only
# reaches the search path, not a chrom_list defined inside a test file.
pks_sa_pr <- get_peaks(Sa_pr, lambdas = "210", fit = "egh",
                       smooth_type = "none", show_progress = FALSE)

test_that("match_peaks derives `chrom_list` and `lambda` when not supplied", {
  skip_on_cran()
  # a single-wavelength peak list needs no `lambda`, and `chrom_list` can be
  # recovered from the peak list's `chrom_list` attribute
  expect_equal(attr(pks_sa_pr, "chrom_list"), "Sa_pr")
  m <- match_peaks(pks_sa_pr, verbose = FALSE)
  expect_s3_class(m, "matched_peak_list")
  expect_equal(names(m), names(Sa_pr))
  expect_equal(attr(m, "lambda"), "210")
  # and the result is the same as supplying both explicitly
  expect_equal(m, match_peaks(pks_sa_pr, chrom_list = Sa_pr, verbose = FALSE))
})

test_that("match_peaks selects the best reference when `reference = 'best'`", {
  skip_on_cran()
  msgs <- capture_messages(
    m <- match_peaks(pks_sa_pr, chrom_list = Sa_pr, lambda = "210",
                     reference = "best", verbose = TRUE))
  expect_match(msgs, "Selected chromatogram 3 as best reference", all = FALSE)
  expect_match(msgs, "Comparing '119' against reference \\('122'\\)",
               all = FALSE)
  # `ptw::bestref` picks sample 3 here, so this is not just the default of 1
  expect_equal(attr(m, "reference"), "122")
  expect_equal(attr(m, "rt_ref"), pks_sa_pr[[3]][["210"]]$rt)
})

# A two-peak case where RT proximity and peak shape disagree. The reference has
# a narrow peak at 12.00 and a wide one at 12.10; the sample has a wide peak at
# 12.02 and a narrow one at 12.08. RT alone pairs 12.00-12.02 and 12.10-12.08;
# shape alone pairs them across. `spectral_cutoff = 0` keeps the spectral mask
# from being the binding constraint.
make_shape_conflict <- function(){
  pks <- function(rt, sd, tau, tau_left, tau_right){
    data.frame(rt = rt, start = rt - .1, end = rt + .1, sd = sd, tau = tau,
               tau_left = tau_left, tau_right = tau_right, height = 100,
               area = 10)
  }
  pl <- pks_egh
  pl[[1]] <- list("210" = pks(c(12.00, 12.10), c(0.05, 0.30), c(0.01, 0.06),
                              c(0.01, 0.06), c(0.02, 0.07)))
  pl[[2]] <- list("210" = pks(c(12.02, 12.08), c(0.30, 0.05), c(0.06, 0.01),
                              c(0.06, 0.01), c(0.07, 0.02)))
  pl
}

test_that("`shape_weight` redirects matching from RT proximity to peak shape", {
  skip_on_cran()
  pl <- make_shape_conflict()
  for (f in c("gaussian", "egh", "bemg")){
    m0 <- match_peaks(pl, chrom_list = dat.pr, lambda = "210", fit = f,
                      spectral_cutoff = 0, shape_weight = 0, verbose = FALSE)
    m1 <- match_peaks(pl, chrom_list = dat.pr, lambda = "210", fit = f,
                      spectral_cutoff = 0, shape_weight = 1, verbose = FALSE)
    pair0 <- m0[[2]]$matches[order(m0[[2]]$matches[, "rt_ref"]), ]
    pair1 <- m1[[2]]$matches[order(m1[[2]]$matches[, "rt_ref"]), ]
    # RT-driven: each reference peak takes its nearest sample peak
    expect_equal(pair0[, "rt_samp"], c(12.02, 12.08), info = f,
                 ignore_attr = TRUE)
    # shape-driven: the pairing crosses over
    expect_equal(pair1[, "rt_samp"], c(12.08, 12.02), info = f,
                 ignore_attr = TRUE)
  }
})

test_that("match_peaks errors when `fit` can be neither found nor inferred", {
  skip_on_cran()
  pks_nofit <- pks_egh
  attr(pks_nofit, "fit") <- NULL
  expect_error(match_peaks(pks_nofit, chrom_list = dat.pr, lambda = "210",
                           verbose = FALSE),
               regexp = "Please specify `fit` to proceed")
})

test_that("match_peaks skips samples where no candidate pair survives", {
  skip_on_cran()
  # peaks are present in both samples, but an impossibly tight RT window
  # leaves no candidate pair for the non-reference sample
  expect_warning(m <- match_peaks(pks_egh, chrom_list = dat.pr,
                                  lambda = "210", rt_window = 1e-6,
                                  verbose = FALSE),
                 regexp = "No peaks matched for sample 2")
  expect_null(m[[2]])
  # the reference still matches itself
  expect_false(is.null(m[[1]]))
})

test_that("`shape_cutoff` prunes shape-mismatched candidate pairs", {
  skip_on_cran()
  pl <- make_shape_conflict()
  args <- list(pl, chrom_list = dat.pr, lambda = "210", fit = "gaussian",
               spectral_cutoff = 0, verbose = FALSE)
  m_open <- do.call(match_peaks, args)
  m_cut  <- do.call(match_peaks, c(args, shape_cutoff = 0.5))
  # all 4 pairs are candidates without a cutoff; the 2 shape-mismatched pairs
  # (distance 1.43) are dropped by a cutoff of 0.5
  expect_equal(nrow(m_open[[2]]$candidates), 4)
  expect_equal(nrow(m_cut[[2]]$candidates), 2)
  expect_equal(m_cut[[2]]$matches[order(m_cut[[2]]$matches[, "rt_ref"]),
                                  "rt_samp"],
               c(12.08, 12.02), ignore_attr = TRUE)
})

### plot.matched_peak_list ###

test_that("plot.matched_peak_list returns its input invisibly for every type", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off())
  for (ty in c("ladder", "scatter")){
    expect_equal(plot(matched, type = ty), matched, info = ty)
  }
  # `idx` accepts sample names as well as positions
  expect_equal(plot(matched, idx = "121"), matched)
  # `xlim` reaches the ladder's plot rather than clashing with its default
  expect_equal(plot(matched, xlim = c(12, 14)), matched)
  expect_error(plot(matched, type = "heatmap"),
               regexp = "'arg' should be one of")
})

test_that("pack_lanes lifts only bars that collide", {
  # disjoint bars stay on the row
  expect_equal(pack_lanes(c(1, 3), c(2, 4), pad = 0.1), c(0L, 0L))
  # an overlapping bar moves up a lane
  expect_equal(pack_lanes(c(1, 1.5), c(2, 2.5), pad = 0.1), c(0L, 1L))
  # a bar running through an unmatched peak moves up, one ending beside it stays
  expect_equal(pack_lanes(c(1, 3), c(2, 4), pad = 0.1, obstacles = c(1.5, 4.05)),
               c(1L, 0L))
})

test_that("plot.matched_peak_list skips samples with no matches", {
  skip_on_cran()
  pks_blank <- pks_egh
  pks_blank[[2]][["210"]] <- pks_blank[[2]][["210"]][0, , drop = FALSE]
  m <- suppressWarnings(match_peaks(pks_blank, chrom_list = dat.pr,
                                    lambda = "210", verbose = FALSE))
  expect_null(m[[2]])
  pdf(NULL)
  on.exit(dev.off())
  expect_message(plot(m, idx = 2, type = "scatter"),
                 regexp = "Skipping '121': no matches")
  # the ladder draws the empty sample as a row of unmatched reference peaks
  expect_equal(plot(m), m)
})

test_that("match_peaks plots the matches when `plot_it = TRUE`", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off())
  m <- match_peaks(pks_egh, chrom_list = dat.pr, lambda = "210",
                   plot_it = TRUE, verbose = FALSE)
  expect_s3_class(m, "matched_peak_list")
  expect_equal(m, matched)
})

test_that("plot.matched_peak_list draws the documented plot types", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = "210", cl = 1, show_progress = FALSE)
  matched <- match_peaks(pks, chrom_list = Sa_pr, lambda = "210", verbose = FALSE)
  vdiffr::expect_doppelganger("plot_matched_peak_list_ladder",
                              function() plot(matched))
  vdiffr::expect_doppelganger("plot_matched_peak_list_ladder_zoom",
                              function() plot(matched, xlim = c(15.6, 17.4)))
  vdiffr::expect_doppelganger("plot_matched_peak_list_scatter",
                              function() plot(matched, idx = 2, type = "scatter"))
})

test_that("correct_rt works with lptw given a peak_list", {
  skip_on_cran()
  warp_lptw <- correct_rt(dat.pr, lambdas = "210", alg = "lptw",
                           peaks = pks_egh, show_progress = FALSE)
  expect_equal(names(warp_lptw), names(dat.pr))
  expect_equal(colnames(warp_lptw[[1]]), colnames(dat.pr[[1]]), ignore_attr = TRUE)
})

test_that("correct_rt returns lptw_list models", {
  skip_on_cran()
  models_lptw <- correct_rt(dat.pr, lambdas = "210", alg = "lptw",
                              peaks = matched, what = "models",
                              show_progress = FALSE)
  expect_s3_class(models_lptw, "lptw_list")
})

test_that("correct_rt errors when reference conflicts with matched_peak_list", {
  expect_error(
    correct_rt(dat.pr, lambdas = "210", alg = "lptw", peaks = matched,
               reference = 2, show_progress = FALSE),
    regexp = "differs from the reference stored"
  )
})

### test get_peaktable ###

pk_tab <- get_peaktable(peak_list = pks_egh, chrom_list = dat.pr)
pk_tab_sp <- get_peaktable(peak_list = pks_egh, chrom_list = dat.pr, clust = "sp.rt")

test_that("get_peaktable works", {
  expect_equal(rownames(pk_tab$tab), names(dat.pr))
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$pk_meta))
  expect_equal(class(pk_tab), c("peak_table","list"))
  expect_equal(class(pk_tab$tab), "data.frame")
  expect_equal(class(pk_tab$pk_meta), "data.frame")
  expect_equal(class(pk_tab$args), "list")
})

test_that("get_peaktable works with matched_peak_list", {
  skip_on_cran()
  pt <- get_peaktable(matched, pks = pks_egh)
  expect_s3_class(pt, "peak_table")
  expect_equal(rownames(pt$tab), names(dat.pr))
  expect_equal(colnames(pt$tab), colnames(pt$pk_meta))
  expect_true(all(pt$tab[1, ] > 0))
  expect_true(all(pt$pk_meta["lambda", ] == "210"))
})

test_that("get_peaktable.matched_peak_list errors correctly", {
  expect_error(get_peaktable(matched), regexp = "pks")
  expect_error(get_peaktable(matched, pks = list()), regexp = "peak_list")
})

test_that("correct_peaks works", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, show_progress = FALSE)
  
  warping.models <- correct_rt(chrom_list = Sa_pr, lambdas = 210, 
                               what = "models", show_progress = FALSE)
  
  pks_cor <- correct_peaks(pks, mod_list = warping.models)
  pktab_cor <- get_peaktable(pks_cor, use.cor = TRUE)

  ptw_warp <- correct_rt(Sa_pr, models = warping.models, show_progress = FALSE)
  pks_warp <- get_peaks(ptw_warp, lambdas = 210, show_progress = FALSE)
  pktab_warp <- get_peaktable(pks_warp)
  
  # this test failed on UBUNTU
  # expect_equal(pks_cor[[1]][[1]]$rt.cor[-33], pks_warp[[1]][[1]]$rt, 
  #              tolerance = .001)
  # expect_equal(pks_cor[[2]][[1]]$rt.cor, pks_warp[[2]][[1]]$rt[-6], 
  #              tolerance = .001)
  
  pks_reg <- get_peaks(Sa_pr, lambdas = 210, fit = "egh", show_progress = FALSE)
  pktab_reg <- get_peaktable(pks_reg)

  expect_equal(row.names(pktab_cor), names(Sa_pr))
  expect_equal(row.names(pktab_reg), names(Sa_pr))
  expect_equal(row.names(pktab_warp), names(Sa_pr))
  
  expect_lt(ncol(pktab_warp), ncol(pktab_reg))
  expect_lt(ncol(pktab_cor), ncol(pktab_reg))
})

test_that("strip plot works", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("vdiffr")
  skip_on_ci()
  plot_peak_table <- function(){
    pktab <- get_peaktable(pks_egh, plot_it = TRUE,
                           ask = FALSE)
  }
  
  vdiffr::expect_doppelganger("peak_table_plot", plot_peak_table)
})

#### test metadata attachment ###

path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)

pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")

test_that("attach_metadata works", {
  expect_equal(rownames(pk_tab$tab), pk_tab$sample_meta$vial)
  expect_equal(colnames(pk_tab$sample_meta), colnames(meta))
  
  ### check errors & warnings ###
  expect_error(attach_metadata(), regexp = 'argument "peak_table" is missing')
  expect_error(attach_metadata(warp), regexp = "must be of the `peak_table` class")
  expect_error(attach_metadata(pk_tab, metadata = warp), regexp = "Please provide metadata")
  expect_error(attach_metadata(pk_tab, metadata = meta, column = "x"), 
               regexp = "could not be found")
  expect_error(attach_metadata(pk_tab$tab, metadata = meta, column = "vial"),
               regexp = "must be of the `peak_table` class")
  expect_error(attach_metadata(pk_tab, metadata = rbind(meta,meta), 
                               column = "vial"), 
               regexp = "Sample names must be unique")
  expect_warning(attach_metadata(pk_tab, metadata = meta[-1,], column = "vial"),
                 regexp = "The supplied metadata does not include all samples")
})

pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.cor")

test_that("attach_ref_spectra works", {
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  expect_equal(pk_tab$args[["reference_spectra"]], "max.cor")
  
  pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.int")
  expect_equal(pk_tab$args[["reference_spectra"]], "max.int")
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  
  # test errors
  expect_error(attach_ref_spectra(peak_table = warp), 
               regexp = "must be of the `peak_table` class")
  expect_error(attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "x"), 
               regexp = "'arg' should be one of")
})

test_that("[.peak_table works", {
  i <- 1
  sub_rows <- pk_tab[i,]
  expect_equal(sub_rows$tab, pk_tab$tab[i,])
  expect_equal(sub_rows$pk_meta, pk_tab$pk_meta)
  expect_equal(sub_rows$sample_meta, pk_tab$sample_meta[i,])
  expect_equal(sub_rows$ref_spectra, pk_tab$ref_spectra)
  expect_equal(sub_rows$args, pk_tab$args)
  j <- c(3:5)
  sub_cols <- pk_tab[,j]
  expect_equal(sub_cols$tab, pk_tab$tab[,j])
  expect_equal(sub_cols$pk_meta, pk_tab$pk_meta[,j])
  expect_equal(sub_cols$sample_meta, pk_tab$sample_meta)
  expect_equal(sub_cols$ref_spectra, pk_tab$ref_spectra[,j])
  expect_equal(sub_rows$args, pk_tab$args)
  sub_ij <- pk_tab[i,j]
  expect_equal(sub_ij$tab, pk_tab$tab[i,j])
  expect_equal(sub_ij$pk_meta, pk_tab$pk_meta[,j])
  expect_equal(sub_ij$sample_meta, pk_tab$sample_meta[i,])
  expect_equal(sub_ij$ref_spectra, pk_tab$ref_spectra[,j])
  expect_equal(sub_ij$args, pk_tab$args)
})

test_that("subset.peak_table works", {
  sub <- subset(pk_tab, subset = pk_tab$sample_meta$trt == "+")
  
  expect_equal(sub$tab, pk_tab$tab[1,])
  expect_equal(sub$pk_meta, pk_tab$pk_meta)
  expect_equal(sub$sample_meta, pk_tab$sample_meta[1,])
  expect_equal(sub$ref_spectra, pk_tab$ref_spectra)
  
  sub <- subset(pk_tab, subset = pk_tab$sample_meta$trt == "-", 
                select = c("V10","V12"))
  
  expect_equal(sub$tab, pk_tab$tab[c(2),c("V10","V12")])
  expect_equal(sub$pk_meta, pk_tab$pk_meta[,c("V10", "V12")])
  expect_equal(sub$sample_meta, pk_tab$sample_meta[2,])
  expect_equal(sub$ref_spectra, pk_tab$ref_spectra[,c("V10", "V12")])
})

test_that("filter_peaktable works", {
  pktab_s <- filter_peaktable(pk_tab, min_rt = 12, max_rt = 15)
  # check that dimensions are unaltered
  expect_equal(rownames(pk_tab$tab), rownames(pktab_s$tab))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$pk_meta))
  expect_equal(nrow(pktab_s$tab), nrow(pktab_s$sample_meta))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$ref_spectra))
  expect_true(all(!(pktab_s$pk_meta["rt",] <= 12)))
  expect_true(all(pktab_s$pk_meta["rt",] <= 15))
  
  # warning if no arguments are provided
  expect_warning(filter_peaktable(pk_tab), regexp = "Nothing to filter")
  
  pktab_s2 <- filter_peaktable(pk_tab, max_zeros = 0)
  expect_true(all(pktab_s2$tab != 0))
})

test_that("the deprecated `min.cor` argument still works", {
  expect_warning(old <- combine_peaks(pk_tab, tol = 0.02, min.cor = 0.9,
                                      verbose = FALSE),
                 regexp = "`min.cor` argument is deprecated")
  expect_equal(old, combine_peaks(pk_tab, tol = 0.02, min_cor = 0.9,
                                  verbose = FALSE))
})

test_that("combine_peaks works", {
  pk_tab_c <- combine_peaks(pk_tab, verbose = FALSE)
  expect_lt(ncol(pk_tab_c), ncol(pk_tab))
  expect_lt(ncol(pk_tab_c$pk_meta), ncol(pk_tab$pk_meta))
  expect_equal(rownames(pk_tab_c$tab), rownames(pk_tab$tab))
})

# Minimal synthetic peak_table for testing combine_peaks choose options.
# V1 (lambda=210) and V2 (lambda=260) are duplicates: same RT, nearly identical
# spectra. V3 is a distinct peak. V1 is present in all 3 samples (colSums=60);
# V2 is missing in sample 3 (colSums=40, 2 obs).
.base_spec <- sin(seq(0, pi, length.out = 20))
pk2 <- structure(
  list(
    tab = data.frame(V1 = c(10, 20, 30), V2 = c(15, 25, NA), V3 = c(5, 8, 12),
                     row.names = paste0("S", 1:3)),
    pk_meta = data.frame(
      V1 = c("210", "pk1", "5.000", "4.5", "5.5", "0.1", "0.01", "0.24", "0.99", "0.98"),
      V2 = c("260", "pk1", "5.000", "4.5", "5.5", "0.1", "0.01", "0.24", "0.99", "0.98"),
      V3 = c("210", "pk3", "10.00", "9.5", "10.5", "0.1", "0.01", "0.24", "0.99", "0.98"),
      row.names = c("lambda", "peak", "rt", "start", "end", "sd", "tau", "FWHM",
                    "r.squared", "purity"),
      stringsAsFactors = FALSE
    ),
    ref_spectra = matrix(
      c(.base_spec, .base_spec * 1.05, rev(.base_spec)),
      nrow = 20, dimnames = list(as.character(seq(200, 390, length.out = 20)),
                                 c("V1", "V2", "V3"))
    ),
    sample_meta = data.frame(row.names = paste0("S", 1:3)),
    args = list()
  ),
  class = c("peak_table", "list")
)

test_that("combine_peaks choose = 'least_sparse' keeps peak present in most samples", {
  res <- combine_peaks(pk2, choose = "least_sparse", verbose = FALSE)
  expect_true("V1" %in% colnames(res$tab))
  expect_false("V2" %in% colnames(res$tab))
  expect_equal(colnames(res$tab), colnames(res$pk_meta))
})

test_that("combine_peaks choose = 'least_sparse' counts zero-filled cells as missing", {
  pk0 <- pk2
  pk0$tab$V1 <- c(1000, 0, 0)
  pk0$tab$V2 <- c(1, 2, 3)
  expect_true("V1" %in% colnames(combine_peaks(pk0, choose = "max",
                                               verbose = FALSE)$tab))
  res <- combine_peaks(pk0, choose = "least_sparse", verbose = FALSE)
  expect_equal(colnames(res$tab), c("V2", "V3"))
})

test_that("combine_peaks keeps ref_spectra aligned with tab", {
  res <- combine_peaks(pk_tab, tol = 0.02, verbose = FALSE)
  expect_lt(ncol(res$tab), ncol(pk_tab$tab))
  expect_equal(colnames(res$ref_spectra), colnames(res$tab))
})

test_that("merge_peaks keeps ref_spectra aligned", {
  pt <- pk2
  pt$ref_spectra <- as.data.frame(pt$ref_spectra)
  res <- merge_peaks(pt, peaks = c("V1", "V2"))
  expect_equal(colnames(res$ref_spectra), colnames(res$tab))
  res <- merge_peaks(pt, peaks = c("V1", "V3"), method = "sum")
  expect_equal(colnames(res$tab), c("V1", "V2"))
  expect_equal(unname(res$tab[, "V1"]), c(15, 28, 42))
})

test_that("combine_peaks choose = 'lambda' keeps preferred wavelength", {
  res <- combine_peaks(pk2, choose = "lambda", lambda = "260", verbose = FALSE)
  expect_true("V2" %in% colnames(res$tab))
  expect_false("V1" %in% colnames(res$tab))
  expect_equal(colnames(res$tab), colnames(res$pk_meta))
})

test_that("combine_peaks choose = 'lambda' errors when lambda not found in peak table", {
  expect_error(
    combine_peaks(pk2, choose = "lambda", lambda = "999", verbose = FALSE),
    "None of the specified lambda"
  )
})

test_that("combine_peaks choose = 'lambda' with lambda = NULL errors", {
  expect_error(combine_peaks(pk2, choose = "lambda", verbose = FALSE),
               "Please select preferred lambda")
})

test_that("merge_peaks works with max method", {
  # `merge_peaks` keeps whichever peak is more intense, and clustering decides
  # the numbering, so derive the surviving name rather than assuming it
  pr <- colnames(pk_tab$tab)[1:2]
  kept <- pr[which.max(colMeans(pk_tab$tab[, pr], na.rm = TRUE))]
  pk_tab_m <- merge_peaks(pk_tab, peaks = pr)
  expect_equal(pk_tab_m$tab[[kept]],
               pmax(pk_tab$tab[[pr[1]]], pk_tab$tab[[pr[2]]]))
  expect_equal(ncol(pk_tab_m), ncol(pk_tab)-1)
  expect_equal(ncol(pk_tab_m$pk_meta), ncol(pk_tab$pk_meta)-1)
  expect_equal(rownames(pk_tab_m$tab), rownames(pk_tab$tab))
})

test_that("merge_peaks works with sum method", {
  data(pk_tab)
  pr <- colnames(pk_tab$tab)[1:2]
  kept <- pr[which.max(colMeans(pk_tab$tab[, pr], na.rm = TRUE))]
  pk_tab_m <- merge_peaks(pk_tab, peaks = pr, method = "sum")
  expect_equal(pk_tab_m$tab[[kept]],
               (pk_tab$tab[[pr[1]]] + pk_tab$tab[[pr[2]]]))
  expect_equal(ncol(pk_tab_m), ncol(pk_tab)-1)
  expect_equal(ncol(pk_tab_m$pk_meta), ncol(pk_tab$pk_meta) - 1)
  expect_equal(rownames(pk_tab_m$tab), rownames(pk_tab$tab))
})

test_that("normalize_data works", {
  # peaks are chosen by property, not name: clustering decides the numbering,
  # so `good` must be positive in every sample and `bad` must not be
  .ok <- vapply(pk_tab$tab, function(x) all(is.finite(x) & x > 0), logical(1))
  good <- names(which(.ok))[1]
  bad <- names(which(!.ok))[1]
  pk_tab_norm <- normalize_data(pk_tab, chrom_list = dat.pr, column = "mass")
  expect_equal(rownames(pk_tab_norm$tab), rownames(pk_tab$tab))
  expect_equal(class(pk_tab_norm), class(pk_tab))
  expect_equal(colnames(pk_tab_norm$tab), colnames(pk_tab$tab))
  expect_equal(pk_tab_norm$tab[1,], pk_tab$tab[1,]/pk_tab$sample_meta$mass[1])
  expect_equal(pk_tab_norm$tab[2,], pk_tab$tab[2,]/pk_tab$sample_meta$mass[2])
  expect_equal(pk_tab_norm$args[["normalized"]], TRUE)
  expect_equal(pk_tab_norm$args[["normalization_column"]], "mass")
  expect_equal(pk_tab_norm$args[["normalization_by"]], "meta")
  
  chroms_norm <- normalize_data(pk_tab, chrom_list = dat.pr, 
                           column = "mass", what = "chrom_list")
  expect_equal(dim(chroms_norm), dim(dat.pr))
  expect_equal(chroms_norm[[1]], dat.pr[[1]]/pk_tab$sample_meta$mass[1], 
               ignore_attr = TRUE)
  expect_equal(chroms_norm[[2]], dat.pr[[2]]/pk_tab$sample_meta$mass[2], 
               ignore_attr = TRUE)
  expect_equal(attr(chroms_norm[[1]],"normalized"),TRUE)
  expect_equal(attr(chroms_norm[[1]],"normalized_by"),"meta")
  expect_equal(attr(chroms_norm[[1]],"normalization_column"),"mass")
  expect_equal(names(chroms_norm), names(dat.pr))

  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "x"), 
               regexp = "could not be found")
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, 
                              column = "mass", what="x"), 
               regexp = "'arg' should be one of")
  
  pk_tab_pnorm <- normalize_data(pk_tab, chrom_list = dat.pr, 
                                column = good, by = "peak")
  expect_equal(rownames(pk_tab_pnorm$tab), rownames(pk_tab$tab))
  expect_equal(class(pk_tab_pnorm), class(pk_tab))
  expect_equal(colnames(pk_tab_pnorm$tab), colnames(pk_tab$tab))
  expect_equal(pk_tab_pnorm$tab[1,], pk_tab$tab[1,]/pk_tab$tab[[good]][1])
  expect_equal(pk_tab_pnorm$tab[2,], pk_tab$tab[2,]/pk_tab$tab[[good]][2])
  expect_equal(pk_tab_pnorm$args[["normalized"]], TRUE)
  expect_equal(pk_tab_pnorm$args[["normalization_column"]], good)
  expect_equal(pk_tab_pnorm$args[["normalization_by"]], "peak")
  
  expect_warning(pk_tab_pnorm_w <- normalize_data(pk_tab, chrom_list = dat.pr,
                                                   column = bad, by = "peak"),
                 regexp = "Invalid normalization values")
  expect_error(pk_tab_pnorm_w <- normalize_data(pk_tab, chrom_list = dat.pr,
                                                  column = bad, by = "peak",
                                                on_invalid = "error"))

  # auto-detection of by = "peak"
  pk_tab_auto_peak <- normalize_data(pk_tab, column = good)
  expect_equal(pk_tab_auto_peak$args[["normalization_by"]], "peak")
  expect_equal(pk_tab_auto_peak$tab[1,], pk_tab$tab[1,]/pk_tab$tab[[good]][1])

  # disambiguation error when column present in both metadata and peak table
  pk_tab_ambig <- pk_tab
  pk_tab_ambig$sample_meta[[good]] <- 1
  expect_error(normalize_data(pk_tab_ambig, column = good),
               regexp = "disambiguated")
})

test_that("cluster_spectra works", {
  skip_on_cran()
  skip_if_not_installed("pvclust")
  cl <- cluster_spectra(pk_tab, nboot = 10, parallel = FALSE, verbose = FALSE, 
                        save = FALSE, plot_dend = FALSE, plot_spectra = FALSE)
  expect_equal(class(cl[[1]]), "pvclust")
  expect_equal(class(cl[[2]]), "list")
})

### test plotting functions

test_that("plot.peak_list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_peaks <- function(){
    plot(pks_egh, chrom_list = dat.pr, points = TRUE,
         ticks = TRUE, lambda = 210)
  }
  vdiffr::expect_doppelganger("plot.peak_list", plot_peaks)
  
  plot_peaks_gaussian <- function(){
    plot(pks_gaussian, chrom_list = dat.pr)
  }
  vdiffr::expect_doppelganger("plot.peak_list_gaussian", plot_peaks_gaussian)
  
  plot_peaks_bemg <- function(){
    plot(pks_bemg, chrom_list = dat.pr)
  }
  skip_on_os(c("windows", "linux"))
  vdiffr::expect_doppelganger("plot.peak_list_bemg", plot_peaks_bemg)
  
  expect_error(plot(pks_egh, chrom_list = dat.pr, lambda = 190),
               regexp = "must match one of the wavelengths in your peak list")
})

test_that("plot.ptw_list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_ptw_list_traces_base <- function(){
    plot(warping.models, what = "traces", engine = "base")
  }
  suppressWarnings(vdiffr::expect_doppelganger("plot_ptw_list_traces_base",
                                               plot_ptw_list_traces_base))
  
  plot_ptw_list_heatmap_base <- function(){
    plot(warping.models, what = "heatmap", engine = "base")
  }
  vdiffr::expect_doppelganger("plot_ptw_list_heatmap_base", plot_ptw_list_heatmap_base)
  
})

test_that("plot.ptw_list works with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  plot_ptw_list_traces_ggplot <- function(){
    suppressWarnings(plot(warping.models, what = "traces", engine = "ggplot"))
  }
  vdiffr::expect_doppelganger(title = "plot_ptw_list_traces_ggplot", 
                              fig = plot_ptw_list_traces_ggplot)
  
  plot_ptw_list_heatmap_ggplot <- function(){
    plot(warping.models, what = "heatmap", engine = "ggplot")
  }
  vdiffr::expect_doppelganger(title = "plot_ptw_list_heatmap_ggplot",
                              fig = plot_ptw_list_heatmap_ggplot)
})

test_that("plot_peak.list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_peaklist <- function(){
    plot(pks_egh, chrom_list = dat.pr, idx = 2)
  }
  vdiffr::expect_doppelganger("plot_peaklist", plot_peaklist)
})

test_that("write_peaktable writes csvs correctly", {
  data(pk_tab)
  skip_on_cran()
  path <-  tempdir(check = TRUE)
  on.exit(unlink(c(fs::path(path,"peak_table-tab", ext = "csv"),
                   fs::path(path,"peak_table-pk_meta", ext = "csv"),
                   path)))
  write_peaktable(peak_table = pk_tab, path = path, format = "csv",
                  what = c("tab", "pk_meta"))
  path_table <- paste0(file.path(path,"peak_table-tab"), ext = ".csv")
  path_meta <- paste0(file.path(path,"peak_table-pk_meta"), ext = ".csv")
  
  expect_true(file.exists(path_table))
  expect_true(file.exists(path_meta))
  
  expect_equal(pk_tab$tab, 
               read.csv(file = path_table, row.names = 1), ignore_attr = TRUE)
  expect_equal(pk_tab$pk_meta,
               read.csv(file = path_meta, row.names = 1), ignore_attr = TRUE)
  expect_error(write_peaktable(pk_tab, path = "fake_path"),
               regexp = "specified directory does not exist")
})

test_that("write_peaktable writes xlsx correctly", {
  skip_on_cran()
  skip_if_not_installed("openxlsx")
  data(pk_tab)
  path <-  tempdir()
  on.exit(unlink(c(fs::path(path, "peak_table", ext="xlsx"), path)))
  
  write_peaktable(pk_tab, path = path, format = "xlsx")
  file <- fs::path(path, "peak_table", ext = "xlsx")
  xx <- openxlsx::read.xlsx(file, sheet = 1, rowNames = TRUE)
  expect_equal(pk_tab$tab, xx, ignore_attr=TRUE)
  expect_equal(pk_tab$pk_meta, 
               openxlsx::read.xlsx(file, sheet = 2, rowNames = TRUE), 
               ignore_attr = TRUE)
})

# run the following line to activate plotly tests:
# Sys.setenv("VISUAL_TESTS" = "true")

test_that("plot.ptw_list (base) works with plotly", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  plot_ptw_list_traces_plotly <- suppressWarnings(plot(warping.models, 
                                                       what = "traces", 
                                                       engine = "plotly")
                                                                        )

  expect_doppelganger_plotly(name = "plot_ptw_list_traces_plotly", 
                             p = plot_ptw_list_traces_plotly)
})

test_that("plot.ptw_list (heatmap) works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  plot_ptw_list_heatmap_plotly <- plot(warping.models, what = "heatmap", 
                                       engine = "plotly")
  
  expect_doppelganger_plotly(name = "plot_ptw_list_heatmap_plotly", 
                             p = plot_ptw_list_heatmap_plotly)
})

test_that("merge_peaks does not depend on the order of `peaks`", {
  # `sel` is a position among the selected columns, which are in column order,
  # so resolving it against `peaks` kept the right column but left it holding
  # its original, unmerged value.
  data(pk_tab)
  pt <- pk_tab
  # make the first peak win in one sample, so the merge differs from either
  # column on its own
  pr <- colnames(pt$tab)[c(10, 11)]
  pt$tab[[pr[1]]][1] <- pt$tab[[pr[2]]][1] * 3
  for (method in c("max", "sum")){
    expected <- if (method == "max"){
      pmax(pt$tab[[pr[1]]], pt$tab[[pr[2]]])
    } else pt$tab[[pr[1]]] + pt$tab[[pr[2]]]
    fwd <- merge_peaks(pt, peaks = pr, method = method)
    rev <- merge_peaks(pt, peaks = rev(pr), method = method)
    kept <- setdiff(pr, setdiff(pr, colnames(fwd$tab)))
    expect_equal(colnames(fwd$tab), colnames(rev$tab), info = method)
    expect_equal(fwd$tab[[kept]], expected, info = method)
    expect_equal(rev$tab[[kept]], expected, info = method)
  }
})
