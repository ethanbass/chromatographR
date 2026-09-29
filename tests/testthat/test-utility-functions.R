test_that("get_times works as expected", {
  data(Sa_pr)
  ts <- get_times(Sa_pr)
  expect_equal(ts, as.numeric(rownames(Sa_pr[[1]])))
})

test_that("get_lambdas works as expected", {
  data(Sa_pr)
  lambdas <- get_lambdas(Sa_pr)
  expect_equal(lambdas, as.numeric(colnames(Sa_pr[[1]])))
})

test_that("choose_apply_fnc works as expected", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  fn <- choose_apply_fnc(show_progress = TRUE)
  expect_equal(class(fn), c("purrr_function_partial", "function"))
  
  for (fn in list(choose_apply_fnc(show_progress = FALSE, parallel = FALSE),
                  choose_apply_fnc(show_progress = FALSE, cl = 1),
                  choose_apply_fnc(show_progress = FALSE, cl = NULL))){
    expect_warning(expect_equal(fn(c(a = 1, b = 4), function(x){
      warning("w"); sqrt(x)
    }), list(a = 1, b = 2)), "Samples a, b: w")
  }
})

test_that("choose apply_fnc works as expected on unix/linux", {
  skip_on_cran()
  skip_on_os("windows")
  
  fn <- choose_apply_fnc(show_progress = TRUE, parallel = TRUE)
  expect_equal(class(fn), c("purrr_function_partial","function"))
  
  fn <- choose_apply_fnc(show_progress = FALSE, parallel = TRUE)
  expect_equal(class(fn), c("purrr_function_partial","function"))
  
  fn<-choose_apply_fnc(show_progress=NULL, parallel = TRUE)
  pbapply_exists <- check_for_pkg("pbapply", return_boolean = TRUE)
  expect_equal(fn, choose_apply_fnc(show_progress = pbapply_exists, parallel = TRUE))
  
  cl <- parallel::makeCluster(2)
  fn <- choose_apply_fnc(show_progress = NULL, cl = cl)
  expect_equal(class(fn), c("purrr_function_partial", "function"))
  parallel::stopCluster(cl)
  
  fn <- choose_apply_fnc(show_progress = NULL, cl = 2)
  expect_equal(class(fn), c("purrr_function_partial", "function"))
})

test_that("choose apply_fnc works as expected on windows", {
  skip_on_cran()
  skip_on_os(c("mac","linux","solaris"))
  skip_if_not_installed("pbapply")

  cl <- parallel::makeCluster(2)
  fn <- choose_apply_fnc(show_progress = NULL, cl = cl)
  parallel::stopCluster(cl)
  expect_equal(class(fn), c("purrr_function_partial", "function"))
})

test_that("choose_apply_fnc works as expected with NULL value", {
  skip_on_cran()
  pbapply_exists <- check_for_pkg("pbapply", return_boolean = TRUE)
  
  fn<-choose_apply_fnc(show_progress=NULL, parallel = FALSE)
  expect_equal(fn, choose_apply_fnc(show_progress = pbapply_exists, parallel = FALSE))
  
  fn<-choose_apply_fnc(show_progress = NULL, cl = 1)
  expect_equal(fn, choose_apply_fnc(show_progress = pbapply_exists, parallel = FALSE))
  
  fn<-choose_apply_fnc(show_progress = NULL, cl = NULL)
  expect_equal(fn, choose_apply_fnc(show_progress = pbapply_exists, parallel = FALSE))
})

test_that("check_peaktable works as expected", {
  data(Sa_pr)
  data(pk_tab)
  expect_null(check_peaktable(pk_tab))
  expect_error(check_peaktable("xxx"), 
               regexp = "must be of the `peak_table` class")
  expect_error(check_peaktable("Sa_pr"), 
               regexp = "must be of the `peak_table` class")
})

test_that("summary.peak_table works as expected", {
  data(pk_tab)
  expect_snapshot(print(summary(pk_tab)))
})

# test_that("get_chromlist works as expected", {
#   expect_equal(get_chrom_list(x = pk_tab)[[1]], dat.pr[[1]])
# })

test_that("get_retention_idx works as expected", {
  data(Sa_pr)
  idx <- get_retention_idx(RT = 18.49, times = get_times(Sa_pr))
  expect_equal(idx, 425)
  expect_error(get_retention_idx(RT = 9, times = get_times(Sa_pr)), 
               regexp = "falls outside the bounds of the chromatogram")
  expect_error(get_retention_idx(RT = 20, times = get_times(Sa_pr)), 
               regexp = "falls outside the bounds of the chromatogram")
})

test_that("check_idx works as expected", {
  data(Sa_pr)
  expect_error(check_idx(idx = 1000, chrom_list = Sa_pr))
})

test_that("get_lambda_idx works as expected", {
  data(Sa_pr)
  lam <- get_lambda_idx(lambda = 200, lambdas = get_lambdas(x = Sa_pr))
  expect_equal(lam, 1)
  expect_error(get_lambda_idx(lambda=190, lambdas = get_lambdas(x = Sa_pr)), 
               regexp = "specified wavelength \\(`190`\\) could not be found")
  expect_error(get_lambda_idx(lambda=400, lambdas = get_lambdas(x = Sa_pr)), 
               regexp = "specified wavelength \\(`400`\\) could not be found")
  expect_error(get_lambda_idx(lambda="max", lambdas = get_lambdas(x = Sa_pr),
                              allow_max = FALSE), 
               regexp = "must be a specific wavelength")
})

test_that("check_for_pkg functions as expected",{
  expect_error(check_for_pkg(pkg = "nonsense-package-that-doesn't-exist"),
               regexp = "must be installed to perform this action")
})

test_that("reshape_peaktable works as expected",{
  data(pk_tab)
  path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
  meta <- read.csv(path)
  pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")
  
  pktab_long <- reshape_peaktable(pk_tab)
  expect_equal(ncol(pktab_long), 5 + ncol(pk_tab$sample_meta))
  expect_equal(nrow(pktab_long), nrow(pk_tab)*ncol(pk_tab))
  
  pktab_long <- reshape_peaktable(pk_tab, peaks = c("V51","V60"), metadata = "trt")
  expect_equal(ncol(pktab_long), 5 + 1)
  expect_equal(nrow(pktab_long), nrow(pk_tab)*2)

  pktab_named <- reshape_peaktable(pk_tab, peaks = c(a = "V51", b = "V60"))
  expect_equal(levels(pktab_named$peak), c("a", "b"))
  expect_equal(unique(pktab_named$rt[pktab_named$peak == "b"]),
               pk_tab$pk_meta["rt", "V60"])

  pktab_wide <- reshape_peaktable(pk_tab, peaks = c(a = "V51", b = "V60"),
                                  metadata = "trt", format = "wide")
  expect_equal(colnames(pktab_wide), c("sample", "trt", "a", "b"))
  expect_equal(pktab_wide$sample, rownames(pk_tab$tab))
  expect_equal(pktab_wide$b, pk_tab$tab$V60)
  expect_equal(pktab_wide$trt, pk_tab$sample_meta$trt)
})

test_that("reshape_chroms works as expected", {
  data(Sa_pr)
  chrom_list_long <- reshape_chroms(Sa_pr, lambdas = "210")
  expect_equal(ncol(chrom_list_long), 4)
  expect_equal(nrow(chrom_list_long), length(Sa_pr)*nrow(Sa_pr[[1]]))
  expect_equal(unique(chrom_list_long$lambda), 210)
  expect_equal(levels(chrom_list_long$sample), names(Sa_pr))
  chrom_list_long_subset <- reshape_chroms(Sa_pr, lambdas = 210, idx = c(1:2))
  expect_equal(ncol(chrom_list_long_subset), 4)
  expect_equal(nrow(chrom_list_long_subset), 2*nrow(Sa_pr[[1]]))
  expect_equal(levels(chrom_list_long_subset$sample), names(Sa_pr)[1:2])
})

test_that("extract_idx function works as expected", {
  expect_equal(extract_idx("dat.pr[[1]]"), 1)
  expect_equal(extract_idx("dat.pr[1:3]"), 1:3)
  expect_equal(extract_idx("dat.pr[c(1,5,7)]"), c(1,5,7))
})

test_that("head.peak_table works as expected", {
  data(pk_tab)
  expect_snapshot(head(pk_tab, 1))
  expect_snapshot(tail(pk_tab, 1))
})
  
test_that("resolve_deprecated works as expected", {
  max_iter <- 100
  max.iter <- 50
  expect_warning(resolve_deprecated(max.iter, max_iter), 
                 regexp = "The `max.iter` argument is deprecated")
  expect_equal(suppressWarnings(resolve_deprecated(max.iter, max_iter)),50)
  
  max_iter <- 100
  max.iter <- NULL
  expect_equal(resolve_deprecated(max.iter, max_iter),100)
})

test_that("the deprecated `spec.smooth` argument still works", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 18.66, by = .01), dim2 = seq(200, 318, by = 2),
               show_progress = FALSE, cl = 1)
  expect_warning(old <- do.call(preprocess, c(list(Sa[[1]]), args,
                                              spec.smooth = FALSE)),
                 regexp = "`spec.smooth` argument is deprecated")
  expect_equal(as.vector(old),
               as.vector(do.call(preprocess, c(list(Sa[[1]]), args,
                                               smooth_spectra = FALSE))))
})

test_that("check_norm_values behaves correctly for all on_invalid options", {
  samples <- c("s1", "s2", "s3")
  valid   <- c(1.5, 2.0, 3.0)
  
  # valid values returned unchanged
  expect_identical(check_norm_values(valid, samples, "error"), valid)
  expect_identical(check_norm_values(valid, samples, "warn"),  valid)
  expect_identical(check_norm_values(valid, samples, "silent"), valid)
  
  # error throws and returns nothing
  expect_error(check_norm_values(c(NA,  2, 3), samples, "error"), "NA")
  expect_error(check_norm_values(c(0,   2, 3), samples, "error"), "0")
  expect_error(check_norm_values(c(-1,  2, 3), samples, "error"), "negative")
  
  # warn emits message and replaces invalid with NA
  expect_warning(w <- check_norm_values(c(NA,  2, 3), samples, "warn"), "NA")
  expect_identical(w, c(NA_real_, 2, 3))
  
  expect_warning(w <- check_norm_values(c(0,   2, 3), samples, "warn"), "0")
  expect_identical(w, c(NA_real_, 2, 3))
  
  expect_warning(w <- check_norm_values(c(-1,  2, 3), samples, "warn"), "negative")
  expect_identical(w, c(NA_real_, 2, 3))
  
  # silent replaces invalid with NA, no condition
  expect_silent(w <- check_norm_values(c(NA, 2, 3), samples, "silent"))
  expect_identical(w, c(NA_real_, 2, 3))
  
  expect_silent(w <- check_norm_values(c(0,  2, 3), samples, "silent"))
  expect_identical(w, c(NA_real_, 2, 3))
  
  expect_silent(w <- check_norm_values(c(-1, 2, 3), samples, "silent"))
  expect_identical(w, c(NA_real_, 2, 3))
  
  # sample names appear in messages
  expect_error(check_norm_values(c(NA, 2, 3), samples, "error"), "s1")
  
  # multiple invalid values - all reasons appear in message
  w <- testthat::capture_warnings(
    v <- check_norm_values(c(NA, 0, -1), samples, "warn")
  )
  expect_match(w, "NA",       all = FALSE)
  expect_match(w, "0",     all = FALSE)
  expect_match(w, "negative", all = FALSE)
  expect_identical(v, rep(NA_real_, 3))
})

test_that("format_list shows ten items and counts the rest", {
  expect_equal(format_list(c("a", "b"), quote = FALSE), "a, b")
  expect_equal(format_list(letters[1:10], quote = FALSE),
               paste(letters[1:10], collapse = ", "))
  expect_equal(format_list(letters[1:13], quote = FALSE),
               paste0(paste(letters[1:10], collapse = ", "), " and 3 more"))
})

test_that("check_duplicated_names works as expected", {
  expect_error(check_duplicated_names(c("X","X")),
               regexp = "Duplicate sample names found")
  expect_no_error(check_duplicated_names(c("X","Y")))
  expect_warning(check_duplicated_names(c("X","X"), warn = TRUE),
                 regexp = "Duplicate sample names found")
  expect_error(check_duplicated_names(c("BLANK", "A", "BLANK", "blank",
                                        "BLANK", "blank")),
               regexp = "BLANK (3), blank (2).", fixed = TRUE)
  expect_error(check_duplicated_names(c("X","X")), regexp = "make.unique",
               fixed = TRUE)
})

test_that("parallel_lapply re-raises a worker's try-error", {
  as_mclapply <- function(X, FUN, ...)
    lapply(X, function(x) try(FUN(x, ...), silent = TRUE))
  fit_sample <- function(i) if (i == 2) stop("wavelength not found") else list()
  expect_error(parallel_lapply(1:2, fit_sample, apply_fn = as_mclapply),
               "wavelength not found")
  expect_equal(parallel_lapply(1:2, function(i) i, apply_fn = as_mclapply),
               list(1L, 2L))
})

test_that("parallel_lapply re-raises each worker warning once, naming its samples", {
  skip_on_cran()
  skip_on_os("windows")
  f <- function(i){
    warning("shared")
    if (i == 2) warning("only b")
    i
  }
  w <- character(0)
  out <- withCallingHandlers(
    parallel_lapply(c(a = 1, b = 2, c = 3), f, apply_fn = parallel::mclapply,
                    mc.cores = 2),
    warning = function(cond){
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_equal(out, list(a = 1, b = 2, c = 3))
  expect_equal(w, c("Samples a, b, c: shared", "Sample b: only b"))
})

test_that("parallel_lapply stops when a worker process dies", {
  skip_on_cran()
  skip_on_os("windows")
  die <- function(i) if (i == 2) tools::pskill(Sys.getpid(), tools::SIGKILL) else i
  expect_error(suppressWarnings(
    parallel_lapply(c(a = 1, b = 2), die, apply_fn = parallel::mclapply,
                    mc.cores = 2)),
    "worker process died, so no result was returned for b\\.")
})

test_that("get_peaks reports the underlying error from a parallel worker", {
  skip_on_os("windows")
  data(Sa_pr)
  x <- lapply(Sa_pr[1:2], function(m) m[, "210", drop = FALSE])
  expect_error(get_peaks(x, lambdas = "280", cl = 2, show_progress = FALSE),
               "wavelength \\(`280`\\) could not be found")
})

test_that("get_time_resolution converts to the requested unit", {
  data(Sa_warp, envir = environment())
  sec <- Sa_warp
  for (i in seq_along(sec)){
    rownames(sec[[i]]) <- as.numeric(rownames(sec[[i]])) * 60
    attr(sec[[i]], "time_unit") <- "Seconds"
  }
  expect_equal(get_time_resolution(sec), 1.2)
  expect_equal(get_time_resolution(sec, unit = "min"), 0.02)
  expect_equal(get_time_resolution(sec, unit = "ms"), 1200)
  expect_equal(get_time_resolution(Sa_warp, unit = "min"),
               get_time_resolution(Sa_warp))
  attr(sec[[1]], "time_unit") <- "furlongs"
  expect_error(get_time_resolution(sec, unit = "min"), "unrecognized")
})

test_that("normalize_data matches `column` exactly when inferring `by`", {
  data(pk_tab)
  pk <- pk_tab
  pk$sample_meta <- data.frame(V1x = 1:4, row.names = rownames(pk$tab))
  norm <- normalize_data(pk, "V1", on_invalid = "silent")
  expect_equal(norm$args$normalization_by, "peak")
  expect_error(normalize_data(pk_tab, "zzz"), "could not be found")
})

test_that("get_peaktable rejects objects it has no method for", {
  expect_error(get_peaktable("nonsense"),
               regexp = "No get_peaktable method for object of class character")
  expect_error(get_peaktable(1:3),
               regexp = "No get_peaktable method for object of class integer")
})

test_that("get_peaktable accepts a peak list that has lost its class", {
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  bare <- lapply(pks, function(s) lapply(s, subset, height > 100))
  expect_false(inherits(bare, "peak_list"))
  expect_s3_class(get_peaktable(bare), "peak_table")
})

test_that("pad_zeros and strip_zeros round-trip", {
  x <- matrix(seq_len(12), ncol = 3)
  for (side in c("both", "left", "right")){
    padded <- pad_zeros(x, n_zeros = 2, side = side)
    n_added <- if (side == "both") 4 else 2
    expect_equal(nrow(padded), nrow(x) + n_added, info = side)
    expect_equal(strip_zeros(padded, n_zeros = 2, side = side), x, info = side)
  }
  # a vector is promoted to a one-column matrix
  padded_vec <- pad_zeros(1:4, n_zeros = 1, side = "left")
  expect_equal(dim(padded_vec), c(5L, 1L))
  expect_equal(padded_vec[1, 1], 0)
})

test_that("predict.ptw matches the method it overrides in ptw", {
  skip_on_cran()
  # `chromatographR` registers its own `predict.ptw`, overwriting `ptw`'s. The
  # two must agree, or loading the package silently changes `ptw`'s behavior.
  data(Sa_pr)
  ptw_predict <- getS3method("predict", "ptw", envir = asNamespace("ptw"))
  tr <- t(sapply(Sa_pr, function(x) x[, "210"]))
  colnames(tr) <- rownames(Sa_pr[[1]])
  pw_global <- ptw::ptw(ref = tr[3, , drop = FALSE], samp = tr[1:2, ],
                        warp.type = "global")
  pw_indiv  <- ptw::ptw(ref = tr[3, , drop = FALSE], samp = tr[1:2, ],
                        warp.type = "individual")

  for (pw in list(global = pw_global, individual = pw_indiv)){
    # response, with and without `newdata`
    expect_equal(predict(pw), ptw_predict(pw))
    expect_equal(predict(pw, newdata = tr[1:2, ]),
                 ptw_predict(pw, newdata = tr[1:2, ]))
    # a bare vector is promoted to a one-row matrix
    expect_equal(predict(pw, newdata = tr[1, ]),
                 ptw_predict(pw, newdata = tr[1, ]))
    # corrected times, with and without `newdata` and `RTref`
    expect_equal(predict(pw, what = "time"), ptw_predict(pw, what = "time"))
    expect_equal(predict(pw, what = "time", newdata = c(12, 14, 16)),
                 ptw_predict(pw, what = "time", newdata = c(12, 14, 16)))
    expect_equal(predict(pw, what = "time", RTref = as.numeric(colnames(tr))),
                 ptw_predict(pw, what = "time",
                             RTref = as.numeric(colnames(tr))))
    expect_error(predict(pw, what = "nonsense"),
                 regexp = "'arg' should be one of")
  }
  # `warp.type = "individual"` needs one warp function per row of `newdata`
  expect_error(predict(pw_indiv, newdata = tr[1:3, ]),
               regexp = "Incorrect number of rows in newdata")
  # the row-count guard only covers `warp.type = "individual"`; a global warp
  # is recycled to `nrow(object$sample)` rows, so a mismatched `newdata`
  # fails downstream instead. `ptw`'s own method does the same, so this is
  # inherited rather than introduced here.
  expect_error(predict(pw_global, newdata = tr[1:3, ]),
               regexp = "subscript out of bounds")
  expect_error(ptw_predict(pw_global, newdata = tr[1:3, ]),
               regexp = "subscript out of bounds")
})
