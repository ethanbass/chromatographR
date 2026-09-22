#' Check that peak_table is of proper class
#' 
#' This function validates that the provided object is of the `peak_table`
#' class.
#' 
#' @author Ethan Bass
#' @noRd
check_peaktable <- function(peak_table){
  if (!inherits(peak_table, "peak_table"))
    stop("The provided peak_table must be of the `peak_table` class")
}

#' Retrieve chrom_list from peak_table object and check properties
#' @param peak_table Peak table object
#' @author Ethan Bass
#' @noRd
get_chrom_list <- function(x, chrom_list = NULL, verbose = FALSE){
  if (is.null(chrom_list)){
    if (inherits(x, "peak_table")){
      string <- x$args[["chrom_list"]]
    } else if (inherits(x, "peak_list") | inherits(x, "ptw_list")){
      string <- attr(x, "chrom_list")
    }
    if (grepl("\\[*\\]", string)){
      subsetted <- TRUE
      idx <- extract_idx(string)
      string <- gsub("\\[\\[?(.*?)\\]?\\]", "", string)
    } else {
      subsetted <- FALSE
    }
    chrom_list <- try(get(string))
    if (inherits(chrom_list, "try-error")){
      stop("Chromatograms not found! Please make sure the appropriate chrom_list
             object has been loaded")
    }
    if (subsetted){
      chrom_list <- chrom_list[idx]
    }
  }
  check_chrom_list(x, chrom_list, verbose = verbose)
  chrom_list
}

#' Extracts idx from string
#' @noRd
extract_idx <- function(string, chrom_names){
  idx <- sub(".*?\\[(.*?)\\].*", "\\1", string)
  if (any(grepl(":", idx))){
    if (grepl("c(.*?)", idx)){
      matches <- regmatches(idx, gregexpr('\\(.*\\)', idx))[[1]]
      idx <- gsub('[()]', '', matches)
    }
    split <- as.numeric(strsplit(idx, ":")[[1]])
    idx <- split[1]:split[2]
  } else if (grepl("c(.*?)", idx)){
    if (grepl("[']|[\"]", idx)){
      matches <- regmatches(idx, gregexpr("'([^']*)'|\"([^\"]*)\"", idx))[[1]]
      idx <- gsub('[\'\"]', '', matches)
    } else{
      matches <- regmatches(idx, gregexpr('\\(.*\\)', idx))[[1]]
      idx <- gsub('[()]', '', matches)
      idx <- strsplit(idx,",")[[1]]
    }
  }
  if (any(!is.na(suppressWarnings(as.numeric(idx))))){
    idx <- as.numeric(idx)
  }
  idx
}

#' Check that chrom_list is matching a peak table or peak list.
#' @noRd
check_chrom_list <- function(x, chrom_list, verbose = FALSE){
  if (inherits(x, "peak_table")){
    if (length(chrom_list) != nrow(x$tab)){
      stop("Dimensions of chrom_list and peak_table do not match.")
    } else{
      if (verbose & any(names(chrom_list) != rownames(x$tab)))
        warning("Names of chromatograms do not match peak_table")
    }
  } else if (inherits(x, "peak_list")){
    if (length(chrom_list) != length(x)){
      stop("Dimensions of chrom_list and peak_list do not match.")
    } else{
      if (verbose & any(names(chrom_list) != names(x)))
        warning("Names of chromatograms do not match peak_list")
    }
  }
}

#' Get retention index
#' @author Ethan Bass
#' @noRd
get_retention_idx <- function(RT, times){
  if (!is.numeric(RT))
    stop("Retention time not found!")
  if (RT > tail(times, 1) | RT < head(times, 1))
    stop("The supplied retention time falls outside the bounds of the chromatogram.")
  idx <- which.min(abs(RT - times))
  idx
}

#' Check index
#' @author Ethan Bass
#' @noRd
check_idx <- function(idx, chrom_list){
  if (idx > nrow(chrom_list[[1]]))
    stop("The supplied index time exceeds the bounds of the chromatogram.")
}

#' Get wavelength index
#' @param lambda lambda to be found
#' @param lambdas lambdas to be matched
#' @param y Signal as numeric vector.
#' @param allow_max Logical. Whether to return lambda at maximum signal intensity
#' of vector `y`. Defaults to `TRUE.`
#' @author Ethan Bass
#' @noRd
get_lambda_idx <- function(lambda, lambdas, y, allow_max = TRUE){
  if (identical(lambda, "max")){
    if (allow_max){
      lambda.idx <- which.max(if (is.null(dim(y))) y else rowMeans(y, na.rm = TRUE))
    } else{
      stop("Wavelength (`lambda`) must be specified for interactive scanning.")
    }
  } else{
    lambda.idx <- ifelse((length(lambdas) == 1), 1,
                                which(lambdas == as.numeric(lambda))
    )
  }
  if (is.na(lambda.idx) | length(lambda.idx) == 0)
    stop(sprintf("The specified wavelength (`%s`) could not be found!", lambda))
  lambda.idx
}

#' Check chromatogram
#' @noRd
check_chr <- function(chr, loc = NULL, peak_table, chrom_list, allow_max = TRUE){
  if (identical(chr, "max")){
    if (allow_max){
      if (length(loc) > 1) {
        chr <- which.max(rowSums(apply(peak_table$tab[, loc, drop = FALSE], 2, 
                                       function(x) x / max(x, na.rm = TRUE))))
      } else{
        chr <- which.max(peak_table$tab[, loc])
      }
    } else{
      stop("Chromatogram must be specified for scan function.")
    }
  }
  if (is.character(chr)){
    if (!(chr %in% names(chrom_list))){
      stop(sprintf("Chromatogram `%s` could not be found. Please check `idx` argument and try again!", chr))
    } else{
        chr <- which(names(chrom_list) == chr)
      }
  }
  chr
}

#' Elementwise all equal function
#' @author Brian Diggs
#' @references https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
#' @noRd
elementwise.all.equal <- Vectorize(function(x, y, ...){
  isTRUE(all.equal(x, y, ...))
})

#' Get retention times
#' 
#' Get retention times from a list of chromatograms or a `peak_table` object.
#' 
#' @param x A list of chromatograms or `peak_table` object.
#' @param idx Index of chromatogram from which to extract times.
#' @return Numeric vector of retention times from the chromatogram specified by
#' `idx`.
#' @family accessors
#' @examples
#' data(Sa_warp)
#' get_times(Sa_warp)
#' @export

get_times <- function(x, idx = 1){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    as.numeric(rownames(x[[idx]]))
  } else if (inherits(x, c("matrix","data.frame"))){
    as.numeric(rownames(x))
  }
}

#' Get lambdas
#' 
#' Get wavelengths from a list of chromatograms or a `peak_table` object.
#' 
#' @param x A list of chromatograms or `peak_table` object.
#' @return A numeric vector of wavelengths.
#' @family accessors
#' @examples
#' data(Sa_warp)
#' get_lambdas(Sa_warp)
#' @export

get_lambdas <- function(x){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    as.numeric(colnames(x[[1]]))
  } else if (inherits(x, "matrix")){
    as.numeric(colnames(x))
  }
}

#' Retrieve estimated baseline
#'
#' Retrieves the baseline estimated by [`preprocess`] from a chromatogram, or a
#' list of chromatograms, preprocessed with `attach_baseline = TRUE`. This
#' function does not estimate a baseline itself; to compute one, see the
#' `baseline` argument of [`preprocess`] or [`correct_baseline`].
#'
#' The baseline is indexed by the time axis, so it is deliberately dropped by
#' any transform that resamples the signal ([`correct_rt`] and
#' [`subtract_blanks`]). Call this before aligning, or re-run [`preprocess`].
#'
#' @param x A chromatogram, or a list of chromatograms, returned by
#' [`preprocess`] with `attach_baseline = TRUE`.
#' @return A matrix of the same dimensions as `x` containing the estimated
#' baseline, or a list of such matrices if `x` is a list.
#' @family accessors
#' @examples
#' data(Sa)
#' pr <- preprocess(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
#'                  dim2 = seq(200, 318, by = 2), attach_baseline = TRUE)
#' bl <- get_baseline(pr)
#' @export

get_baseline <- function(x){
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    out <- lapply(x, get_baseline)
    class(out) <- class(x)
    return(out)
  }
  bl <- attr(x, "baseline", exact = TRUE)
  if (is.null(bl))
    stop("No baseline found. `preprocess` and `correct_baseline` retain the ",
         "estimated baseline only when called with `attach_baseline = TRUE`.")
  bl
}

#' Get time resolution
#'
#' Get the sampling interval of a chromatogram: the median gap between
#' consecutive retention times.
#'
#' @param x A list of chromatograms or `peak_table` object.
#' @param idx Index of chromatogram from which to derive the interval.
#' @param digits Number of significant digits to round to. A regular time axis
#' carries floating-point noise in the last digits of `diff`, so the median is
#' rounded by default.
#' @param unit Unit for the returned interval: `"min"`, `"s"` or `"ms"`.
#' Defaults to `NULL`, giving the interval in the units of the chromatogram's
#' own time axis, which is what [`preprocess`]'s `dim1` expects. Supply
#' `"min"` to pair with the `time_resolution` arguments of [`plot_chroms`] and
#' [`reshape_chroms`], which are in minutes whatever the units of `x`. The
#' units of `x` are read from its `time_unit` attribute, falling back to
#' minutes when it is absent, as [`plot_chroms`] does.
#' @return A length-one numeric giving the time between consecutive timepoints,
#' in `unit` if supplied and otherwise in the units of the chromatogram's time
#' axis.
#' @family accessors
#' @examples
#' data(Sa_warp)
#' get_time_resolution(Sa_warp)
#' @export

get_time_resolution <- function(x, idx = 1, digits = 6, unit = NULL){
  ts <- get_times(x = x, idx = idx)
  res <- median(diff(ts))
  if (!is.null(unit)){
    unit <- match.arg(unit, c("min", "s", "ms"))
    axis_unit <- get_time_unit(x, idx = idx, na_value = "min")
    if (!isTRUE(axis_unit %in% c("min", "s", "ms"))){
      stop(sprintf("Cannot convert the sampling interval to `%s`: the time unit of `x` (`%s`) is unrecognized.",
                   unit, axis_unit))
    }
    res <- res * time_factor(unit) / time_factor(axis_unit)
  }
  signif(res, digits = digits)
}

#' Format numbers for an error message, one at a time
#' @noRd
num_label <- function(v){
  vapply(v, function(z) format(z, trim = TRUE, scientific = FALSE,
                               drop0trailing = TRUE), character(1),
         USE.NAMES = FALSE)
}

#' Validate an interpolation grid
#'
#' `preprocess` reads the range of `dim1`/`dim2` with `head`/`tail` and
#' interpolates onto them with `approx`, both of which assume a sorted grid, so
#' the assumption is checked once here rather than surfacing as a confusing
#' extrapolation error. A named `by`, `length.out` or `along.with` element is
#' the fingerprint of `c(from, to, by = x)` written where `seq(from, to,
#' by = x)` was meant: `c()` silently returns a three-element vector.
#' @noRd
check_interp_grid <- function(v, arg, check_order = TRUE){
  nm <- names(v)
  seq_arg <- intersect(nm, c("by", "length.out", "along.with", "from", "to"))
  if (length(seq_arg) > 0){
    # formatted element-wise: `format` on a vector pads to a common width and
    # switches to scientific notation, neither of which echoes what was typed
    args <- paste0(ifelse(nzchar(nm), paste0(nm, " = "), ""), num_label(v),
                   collapse = ", ")
    stop(sprintf(paste0(
      "`%s` was given `c(%s)`, which is a %d-element vector rather than a ",
      "sequence.\n  Did you mean `seq(%s)`?"),
      arg, args, length(v), args))
  }
  if (!is.numeric(v) || length(v) == 0 || !all(is.finite(v)))
    stop(sprintf("`%s` must be a numeric vector of finite values.", arg))
  if (check_order && length(v) < 2)
    stop(sprintf(paste(
      "`%s` must contain at least two values to interpolate onto, but it has",
      "%d. To select a single %s without interpolating, set `%s = FALSE`."),
      arg, length(v), if (arg == "dim1") "timepoint" else "wavelength",
      if (arg == "dim1") "interpolate_rows" else "interpolate_cols"))
  if (check_order && length(v) > 1){
    bad <- which(diff(v) <= 0)
    if (length(bad) > 0)
      stop(sprintf(paste0(
        "`%s` must be strictly increasing, but element %d (%s) does not ",
        "precede element %d (%s)."),
        arg, bad[1], num_label(v[bad[1]]), bad[1] + 1L,
        num_label(v[bad[1] + 1L])))
  }
  invisible(v)
}

#' Sampling interval of a numeric axis
#' @noRd
axis_dt <- function(v){
  v <- suppressWarnings(as.numeric(v))
  if (length(v) > 1 && all(is.finite(v))) median(diff(v)) else NA_real_
}


#' Multiplier from minutes into a given time unit
#' @noRd
time_factor <- function(time_unit = c("min", "s", "ms")){
  time_unit <- match.arg(time_unit)
  switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
}

#' Plural suffix
#' @noRd
plural <- function(n) if (n > 1) "s" else ""

#' Get time unit
#' 
#' Get time units from a list of chromatograms or a `peak_table` object.
#' 
#' @param x A list of chromatograms or `peak_table` object.
#' @param idx Index of chromatogram from which to extract times.
#' @param na_value What to return if `time_unit` attribute is not present.
#' @return A string specifying the time units for the provided chromatograms.
#' @noRd
get_time_unit <- function(x, idx = 1, na_value = NA){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    chrom <- x[[1]]
  } else if (inherits(x, "matrix")){
    chrom <- x
  } else{
    chrom <- NA
  }
  time_unit <- attr(chrom, "time_unit")
  if (is.null(time_unit)){
    return(na_value)
  } else{
    switch(time_unit, "Minutes" = "min", 
           "Seconds" = "s", 
           "Milliseconds" = "ms",
            time_unit)
  }
}

#' Convert time units
#' @noRd
convert_time_units <- function(chrom_list, time_unit = c("min", "s", "ms")){
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
  lapply(chrom_list, function(x){
    rownames(x) <- format_axis(as.numeric(rownames(x))*tfac)
    attr(x, "time_unit") <- time_unit
    x
  })
}


#' Check for suggested package
#' 
#' This function checks for a suggested package and returns an error if the 
#' package is not installed (if `return_boolean` is `FALSE.` Otherwise, it 
#' returns a boolean value.
#' 
#' @noRd
check_for_pkg <- function(pkg, return_boolean = FALSE){
  pkg_exists <- requireNamespace(pkg, quietly = TRUE)
  if (return_boolean){
    return(pkg_exists)
  } else if (!pkg_exists) {
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
}

#' Extract variables from the left-hand-side of a formula.
#' @param formula A [`formula`] object.
#' @importFrom Formula Formula
#' @noRd
#' @note Adapted from https://github.com/adibender/pammtools/blob/master/R/formula-utils.R
#' (c) Copyright © 2017 Andreas Bender and Fabian Scheipl under MIT license:
#' Permission is hereby granted, free of charge, to any person obtaining a copy of
#' this software and associated documentation files (the “Software”), to deal in 
#' the Software without restriction, including without limitation the rights to use, 
#' copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
#' Software, and to permit persons to whom the Software is furnished to do so, 
#' subject to the following conditions:
#'   
#' The above copyright notice and this permission notice shall be included in all 
#' copies or substantial portions of the Software.
#' 
#' THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#' WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#' CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

get_lhs_vars <- function(formula) {
  if (is.character(formula) ) formula <- as.formula(formula)
  form <- formula(Formula::Formula(formula), lhs = TRUE, rhs = FALSE)
  all.vars(form)
}

#' Coerce preprocessing input to a list of matrices
#' @noRd
as_chrom_list <- function(x, arg = "X"){
  x_class <- class(x)
  if (is.matrix(x)){
    x <- list(x)
    return_matrix <- TRUE
  } else return_matrix <- FALSE
  if (inherits(x, "list") &&
      all(sapply(x, function(xx){
        inherits(xx, "list") && "dad" %in% names(xx) && is.matrix(xx$dad)
      }))){
    x <- lapply(x, function(xx) xx$dad)
  }
  if (!(inherits(x, "matrix") || (inherits(x, "list") && all(sapply(x, is.matrix)))))
    stop(sprintf("%s should be a matrix or a list of matrices", arg))
  list(chrom_list = x, return_matrix = return_matrix, class = x_class)
}

#' Re-emit baseline warnings stashed by the per-matrix workers
#' @noRd
flush_baseline_warnings <- function(x){
  bl_warnings <- unique(unlist(lapply(x, attr, "baseline_warnings")))
  x <- lapply(x, function(xx){
    attr(xx, "baseline_warnings") <- NULL
    xx
  })
  for (w in bl_warnings) warning(w, call. = FALSE, immediate. = TRUE)
  x
}

#' Transfer metadata
#' Transfers metadata attributes between objects.
#' @details
#' `"baseline"` is excluded by default because it is indexed by the time axis
#' and is therefore invalidated by any transform that resamples or rescales the
#' signal. This matters most for warping: `correct_rt` transfers metadata over
#' the warped matrices, and under `vpdtw` the baseline and the data end up with
#' different numbers of rows outright, while under `ptw` the shapes and
#' rownames still agree even though the signal has been mapped onto a new time
#' axis and the baseline has not -- so it looks valid while being silently
#' misaligned. Note that `exclude` only filters what is copied from
#' `old_object`: an attribute that arithmetic already carried onto `new_object`
#' survives, which is why `subtract_blanks` drops the baseline itself.
#' @noRd
transfer_metadata <- function(new_object, old_object, transfer_class = TRUE,
                              exclude = c('names','row.names','dim','dimnames',
                                          'baseline')){
  if (!transfer_class){
    exclude <- c(exclude, "class")
  }
  a <- attributes(old_object)
  a[exclude] <- NULL
  attributes(new_object) <- c(attributes(new_object), a)
  new_object
}

#' Choose apply function
#' This function chooses an apply function based on arguments provided by the
#' user. The options are `lapply`, `pblapply` and `mclapply`.
#' @importFrom parallel mclapply
#' @return Returns [`pbapply`][pbapply::pbapply] if `progress_bar == TRUE`,
#' otherwise returns [`lapply`].
#' @noRd
choose_apply_fnc <- function(show_progress, parallel = NULL, cl = 2){
  pbapply_installed <- check_for_pkg("pbapply", return_boolean = TRUE)
  is_not_windows <- .Platform$OS.type != "windows"
  
  if (is.null(show_progress)){
    show_progress <- pbapply_installed
  }
  
  if (is.null(parallel)){
    if (!is.null(cl) && 
        (class(cl)[1] == "SOCKcluster" || (is_not_windows && cl > 1))){
      parallel <- TRUE
    } else {
      parallel <- FALSE
    }
  }
  
  if (!parallel){
    cl <- 1
  }

  if (.Platform$OS.type == "windows" && parallel && !pbapply_installed){
    stop("Please install the pbapply package to use parallel processing on Windows.")
  }
  
  if (show_progress){
    check_for_pkg("pbapply")
    pblapply <- pbapply::pblapply
    fn <- purrr::partial(pblapply, cl = cl)
  } else if (!show_progress && is_not_windows && is.numeric(cl) && cl > 1){
      fn <- purrr::partial(mclapply, mc.cores = cl)
  } else {
    fn <- lapply
  }
  fn
}

#' Get y bounds
#' @noRd
get_y_bounds <- function(x, idx, lambdas.idx, pad = 1.1){
  mx <- get_maximum(x, idx, lambdas.idx = lambdas.idx)*pad
  c(0, mx)
}

#' Get maximum
#' @noRd
get_maximum <- function(x, idx, lambdas.idx){
  max(sapply(x[idx], function(xx){
    max(xx[, lambdas.idx], na.rm = TRUE)
  }))
}

#' Get minimum
#' @noRd
get_minimum <- function(x, idx, lambdas.idx){
  min(sapply(x[idx], function(xx){
    min(xx[, lambdas.idx], na.rm = TRUE)
  }))
}

#' Simple cap
#' Copied from `toupper` examples.
#' @noRd
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#' Deprecated function warning
#' @noRd
resolve_deprecated <- function(arg, new_arg){
  if (!is.null(arg)){
    warning(sprintf("The `%s` argument is deprecated. Please use `%s` instead.",
                    deparse(substitute(arg)), deparse(substitute(new_arg))),
            immediate. = TRUE)
    arg
  } else{
    new_arg
  }
}

#' Resolve the deprecated `remove.time.baseline` argument
#' @noRd
resolve_remove_time_baseline <- function(remove.time.baseline, baseline,
                                        baseline_supplied){
  if (is.null(remove.time.baseline)) return(baseline)
  if (!isTRUE(remove.time.baseline) && !isFALSE(remove.time.baseline))
    stop("`remove.time.baseline` must be either `TRUE` or `FALSE`.")
  warning("The `remove.time.baseline` argument is deprecated. Please use ",
          "`baseline = \"none\"` instead.", immediate. = TRUE)
  if (isFALSE(remove.time.baseline)){
    if (baseline_supplied && baseline != "none")
      stop("`remove.time.baseline = FALSE` conflicts with `baseline = \"",
           baseline, "\"`. Supply `baseline = \"none\"` alone instead.")
    return("none")
  }
  if (baseline_supplied && baseline == "none")
    stop("`remove.time.baseline = TRUE` conflicts with `baseline = \"none\"`. ",
         "Supply `baseline` alone instead.")
  baseline
}

#' Get metadata attribute
#' @noRd
get_metadata_attribute <- function(x, which, ...){
  result <- attr(x, which = which, ...)
  if (is.null(result)){
    result <- NA
  }
  result
}

#' @noRd
seq_nrow <- function(x){
  seq_len(nrow(x))
}

#' @noRd
seq_ncol <- function(x){
  seq_len(ncol(x))
}

#' Check for duplicated names
#' @noRd
check_duplicated_names <- function(sample_names, warn = FALSE){
  if (anyDuplicated(sample_names)) {
    fn <- ifelse(warn, warning, stop)
    dups <- paste(sample_names[duplicated(sample_names)], collapse = ", ")
    limit <- getOption("warning.length") - 60
    if (nchar(dups) > limit) {
      last_comma <- max(gregexpr(",", substr(dups, 1, limit))[[1]])
      dups <- paste0(substr(dups, 1, last_comma - 1), ", etc.")
    }
    fn("Duplicate sample names found in peak_list: ", dups)
  }
}

#' Resolve idx and lambda
#' @noRd
resolve_idx_lambda <- function(idx, lambda, loc, peak_table, chrom_list) {
  # resolve idx_max for spectral slice (always scalar)
  idx_max <- which.max(rowSums(apply(peak_table$tab[, loc, drop = FALSE], 2,
                                     function(x) x / max(x, na.rm = TRUE))))
  
  # resolve idx (may be vector)
  idx <- check_chr(idx, loc, peak_table, chrom_list)
  
  # resolve lambda using idx_max
  lambdas <- get_lambdas(chrom_list)
  times   <- get_times(chrom_list)
  lambda_idx <- if (lambda == "max") {
    spectra <- vapply(loc, function(l) {
      rt     <- peak_table$pk_meta["rt", l]
      rt_idx <- get_retention_idx(rt, times)
      y      <- chrom_list[[idx_max]][rt_idx, ]
      y / max(y, na.rm = TRUE)
    }, numeric(length(lambdas)))
    get_lambda_idx("max", lambdas, y = rowSums(spectra))
  } else {
    get_lambda_idx(lambda, lambdas, y = NULL)
  }
  
  list(idx = idx, lambda_idx = lambda_idx, lambda = lambdas[lambda_idx])
}

#' Get peak coordinates
#' @noRd
get_peak_coordinates <- function(peak_table, chrom_list=NULL, loc, idx, lambda){
  chrom_list <- get_chrom_list(peak_table, chrom_list = chrom_list)
  idx_max <- which.max(peak_table$tab[[loc]])
  if (identical(idx, "max")){
    idx <- idx_max
  }
  rt <- peak_table$pk_meta["rt", loc]
  times <- get_times(chrom_list)
  rt_idx <- get_retention_idx(rt, times)
  has_bounds <- all(c("start", "end") %in% rownames(peak_table$pk_meta))
  if (has_bounds){
    start_idx <- get_retention_idx(as.numeric(peak_table$pk_meta["start", loc]), times)
    end_idx   <- get_retention_idx(as.numeric(peak_table$pk_meta["end",   loc]), times)
  }

  lambda_indices <- sapply(lambda, function(l) {
    get_lambda_idx(l, get_lambdas(chrom_list),
                   y = chrom_list[[idx_max]][rt_idx, ])
  })
  lambda_idx <- lambda_indices[which.max(chrom_list[[idx_max]][rt_idx, lambda_indices])]

  int <- max(vapply(chrom_list[idx], function(chr){
    if (has_bounds)
      max(chr[start_idx:end_idx, lambda_idx])
    else
      chr[rt_idx, lambda_idx]
  }, numeric(1)))
  list(coordinates = c(rt, int), lambda_idx = lambda_idx, 
       lambda = get_lambdas(chrom_list)[lambda_idx])
}

#' Default operator for `NULL`
#'
#' `%||%` entered base R in 4.4.0, but the package supports R >= 4.1.0. The
#' right-hand side must stay lazy, since every call site pairs it with
#' `skip_underdetermined` to avoid an expensive fit.
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Format a numeric axis as fixed-width dimnames
#'
#' Gives every element of an axis the same number of decimal places: the fewest
#' that still resolve neighboring values, which on a regular grid is the
#' precision the axis was written in. Assigning a numeric axis to `dimnames`
#' goes through `as.character`, which formats each element independently and
#' drops trailing zeros, so a regular grid comes out ragged (`10.09`, `10.1`,
#' `10.11`) and lookup by name depends on which values happen to end in a
#' zero.
#'
#' @param x A numeric vector, or a character vector of numbers.
#' @return A character vector of the same length, using a common number of
#' decimal places.
#' @author Ethan Bass
#' @noRd
format_axis <- function(x){
  if (is.null(x) || length(x) == 0) return(x)
  num <- suppressWarnings(as.numeric(x))
  if (!all(is.finite(num))) return(as.character(x))
  gaps <- diff(sort(unique(num)))
  tol <- min(gaps) / 1000
  for (digits in 0:12){
    out <- formatC(num, format = "f", digits = digits)
    if (max(abs(as.numeric(out) - num)) < tol) return(out)
  }
  formatC(num, format = "f", digits = 12)
}
