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
#' @param aligned Whether to return chromatograms on the aligned time axis.
#' Defaults to `TRUE` for a `peak_table` built from corrected retention times.
#' @author Ethan Bass
#' @noRd
get_chrom_list <- function(x, chrom_list = NULL, verbose = FALSE,
                           aligned = inherits(x, "peak_table") &&
                             isTRUE(x$args$use_corrected_rt %||% x$args$use.cor)){
  if (is.null(chrom_list)){
    if (inherits(x, "peak_table")){
      string <- x$args[["chrom_list"]]
    } else if (inherits(x, c("peak_list", "ptw_list", "lptw_list", "VPdtw"))){
      string <- attr(x, "chrom_list")
    }
    chrom_list <- if (aligned) get_aligned_chrom_list(x, string)
                  else lookup_chrom_list(string)
  }
  chrom_list <- unnest_chrom_list(chrom_list)
  check_chrom_list(x, chrom_list, verbose = verbose)
  chrom_list
}

#' Look up a chrom_list by its recorded name
#' @noRd
lookup_chrom_list <- function(string){
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
  chrom_list
}

#' Get the chromatograms matching corrected retention times
#'
#' Uses the aligned chromatograms recorded by [`correct_peaks`] if there are
#' any; otherwise re-applies the stored warping models to the unaligned
#' chromatograms.
#' @noRd
get_aligned_chrom_list <- function(x, string){
  aligned_name <- x$args[["chrom_list_aligned"]]
  if (!is.null(aligned_name) && !is.na(aligned_name))
    return(lookup_chrom_list(aligned_name))
  models <- attr(x, "warp_models")
  raw <- lookup_chrom_list(string)
  if (is.null(models)){
    warning(sprintf(paste0(
      "This peak table uses corrected retention times, but records neither ",
      "warping models nor aligned chromatograms, so spectra and traces are ",
      "read from the unaligned `%s`. Rebuild the peak list with ",
      "`correct_peaks`, or pass the aligned chromatograms as `chrom_list`."),
      string), call. = FALSE)
    return(raw)
  }
  alg <- if (inherits(models, "VPdtw")) "vpdtw"
         else if (inherits(models, "lptw_list")) "lptw" else "ptw"
  correct_rt(raw, models = models, alg = alg, show_progress = FALSE, cl = 1)
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
  check_chrom_list_class(chrom_list)
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

#' @noRd
check_chrom_list_class <- function(chrom_list){
  if (inherits(chrom_list, "peak_table") ||
      !inherits(chrom_list, c("chrom_list", "list", "matrix")))
    stop(sprintf(paste("`chrom_list` must be a list of chromatograms or a",
                       "single chromatogram matrix, not a `%s`."),
                 class(chrom_list)[1]), call. = FALSE)
}

#' Select detected peaks by target retention time
#'
#' Restricts a table of detected peak positions to those whose apex falls
#' within `tol` of one of the retention times in `rts`. Called by [`get_peaks`].
#'
#' Filters the output of [`find_peaks`] rather than restricting its input,
#' because peak bounds are derived from derivative sign changes across the whole
#' trace. Truncating the input would clip the bounds of any peak whose true
#' start precedes the window, alter the smoothing bandwidth (which is a fraction
#' of `length(y)`), and introduce smoothing edge artifacts at every window
#' boundary. Filtering afterwards is also where essentially all of the speedup
#' lies, since peak detection is a small fraction of the cost of model fitting.
#' @author Ethan Bass
#' @noRd
select_peaks_by_rt <- function(pos, y, times, rts, tol = 0.1,
                               select_by = c("height", "proximity", "none")){
  select_by <- match.arg(select_by)
  tol <- check_rts(rts, tol)
  apex_idx <- pos[, "pos"]
  keep <- match_rts(times[apex_idx], rts, tol = tol, weight = y[apex_idx],
                    select_by = select_by)
  structure(pos[keep, , drop = FALSE], target = attr(keep, "target"))
}

#' Validate target retention times and tolerances
#'
#' Returns `tol` recycled to the length of `rts`.
#' @author Ethan Bass
#' @noRd
check_rts <- function(rts, tol){
  rts <- suppressWarnings(as.numeric(rts))
  if (length(rts) == 0 || anyNA(rts)){
    stop("Please provide a numeric vector of retention times to `targets`.")
  }
  if (!length(tol) %in% c(1L, length(rts))){
    stop(sprintf(
      "`tol` must be a single value or one value for each of the %d retention time%s in `targets`, not %d.",
      length(rts), plural(length(rts)), length(tol)))
  }
  if (any(!is.finite(tol)) || any(tol < 0)){
    stop("`tol` must be a non-negative number.")
  }
  rep_len(tol, length(rts))
}

#' Reject `bounds = FALSE`, which `find_peaks` accepts but integration cannot use
#' @noRd
check_bounds_arg <- function(...){
  dots <- list(...)
  if ("bounds" %in% names(dots) && !isTRUE(dots$bounds)){
    stop("`bounds = FALSE` cannot be used with `fit_peaks`, since peak boundaries are required in order to integrate. It is only meaningful when calling `find_peaks` directly.")
  }
}

#' Standardize the `targets` argument to `get_peaks` into a table of targets
#'
#' Returns a data.frame with columns `rt`, `tol`, `lambda` and `name`, where
#' `lambda` is `NA` for a target measured at every wavelength and `name` is
#' `NA` for an unnamed target.
#' @author Ethan Bass
#' @noRd
check_targets <- function(rts, tol){
  as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  if (!is.data.frame(rts)){
    rt <- as_number(rts)
    return(data.frame(rt = rt, tol = check_rts(rt, tol),
                      lambda = NA_character_, name = NA_character_))
  }
  if (!"rt" %in% colnames(rts)){
    stop("A data.frame supplied to `targets` must have an `rt` column.")
  }
  if ("tol" %in% colnames(rts)){
    fill <- is.na(rts$tol)
    tol <- replace(as_number(rts$tol), fill, rep_len(tol, nrow(rts))[fill])
  }
  rt <- as_number(rts$rt)
  out <- data.frame(rt = rt, tol = check_rts(rt, tol),
                    lambda = NA_character_, name = NA_character_)
  for (col in c("lambda", "name")){
    if (col %in% colnames(rts)){
      v <- as.character(rts[[col]])
      if (anyNA(v) || any(v == "")){
        stop(sprintf("The `%s` column of `targets` must not contain missing or empty values.", col))
      }
      out[[col]] <- v
    }
  }
  num <- as_number(out$lambda)
  out$lambda[!is.na(num)] <- as.character(num[!is.na(num)])
  if (anyDuplicated(out$name[!is.na(out$name)])){
    stop("Target names must be unique. Duplicated: ",
         paste(unique(out$name[duplicated(out$name)]), collapse = ", "), ".")
  }
  out
}

#' Match peak apices to target retention times
#'
#' Returns the indices of `apex_rt` falling within `tol` of any element of
#' `rts`, resolving ties within a window according to `select_by`. An apex
#' within `tol` of several targets is a candidate only for the nearest of them.
#' The `target` attribute gives the element of `rts` each index was matched to.
#' @author Ethan Bass
#' @noRd
match_rts <- function(apex_rt, rts, tol = 0.1, weight = NULL,
                      select_by = c("height", "proximity", "none")){
  select_by <- match.arg(select_by)
  tol <- rep_len(tol, length(rts))
  d <- abs(outer(apex_rt, rts, "-"))
  d[is.na(d) | d > rep(tol, each = length(apex_rt))] <- Inf
  owner <- max.col(-d, ties.method = "first")
  owner[!is.finite(d[cbind(seq_along(apex_rt), owner)])] <- NA
  keep <- target <- integer(0)
  for (k in seq_along(rts)){
    cand <- which(owner == k)
    if (length(cand) == 0) next
    sel <- switch(select_by,
                  none = cand,
                  proximity = cand[which.min(d[cand, k])],
                  height = cand[which.max(weight[cand])])
    keep <- c(keep, sel)
    target <- c(target, rep(k, length(sel)))
  }
  o <- order(keep)
  structure(keep[o], target = target[o])
}

#' Get retention index
#'
#' Matches retention times to the nearest positions on a sorted time axis.
#' @param RT Numeric vector of retention times to locate.
#' @param times Sorted numeric vector of time points, e.g. from [`get_times`].
#' @param clamp Logical. Whether to snap retention times falling outside
#' `times` to the ends of the axis instead of erroring. Defaults to `FALSE`.
#' @return Returns an integer vector of positions in `times`, one per element
#' of `RT`.
#' @author Ethan Bass
#' @noRd
get_retention_idx <- function(RT, times, clamp = FALSE){
  if (!is.numeric(RT))
    stop("Retention time not found!")
  if (!clamp && any(RT > tail(times, 1) | RT < head(times, 1)))
    stop("The supplied retention time falls outside the bounds of the chromatogram.")
  if (length(times) < 2) return(rep(1L, length(RT)))
  i <- findInterval(RT, times, all.inside = TRUE)
  j <- pmin(i + 1L, length(times))
  ifelse(abs(RT - times[i]) <= abs(times[j] - RT), i, j)
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
      stop("`lambda` must be a specific wavelength here; `\"max\"` is not supported.")
    }
  } else if (length(lambdas) == 1 && is.na(lambdas)){
    lambda.idx <- 1
  } else{
    lambda.idx <- which(lambdas == as.numeric(lambda))
  }
  if (length(lambda.idx) == 0 || anyNA(lambda.idx))
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
      stop(paste0("Please specify a chromatogram with `idx`. `idx = \"max\"` ",
                  "is available only in `plot_spectrum(what = \"peak\")`."))
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
#' Get wavelengths from a chromatogram, a list of chromatograms, or a
#' `peak_table` object.
#' 
#' @param x A chromatogram, a list of chromatograms, or a `peak_table` object.
#' @return A numeric vector of wavelengths, taken from the first chromatogram
#' if `x` is a list.
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
    suppressWarnings(as.numeric(colnames(x[[1]])))
  } else if (inherits(x, "matrix")){
    suppressWarnings(as.numeric(colnames(x)))
  }
}

#' Retrieve estimated baseline
#'
#' Retrieves the baseline retained by [`preprocess`] or [`correct_baseline`]
#' when called with `attach_baseline = TRUE`. This function does not estimate
#' a baseline itself.
#'
#' The baseline is dropped by [`correct_rt`], which warps the signal along the
#' time axis, and by [`subtract_blanks`], after which it no longer describes
#' the signal. Retrieve it before either step. [`crop_chroms`] keeps it,
#' cropped to match.
#'
#' @param x A chromatogram, or a list of chromatograms, returned by
#' [`preprocess`] or [`correct_baseline`] with `attach_baseline = TRUE`.
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
  bl <- attr(x, "estimated_baseline", exact = TRUE)
  if (is.null(bl)){
    if (isTRUE(attr(attr(x, "baseline_param", exact = TRUE), "attached")))
      stop("No baseline found. It was discarded by a later step that changes ",
           "the signal or its time axis, such as `correct_rt` or ",
           "`subtract_blanks`.")
    stop("No baseline found. `preprocess` and `correct_baseline` retain the ",
         "estimated baseline only when called with `attach_baseline = TRUE`.")
  }
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


#' Check manual integration windows
#'
#' Validates and standardizes the `windows` argument to [`get_peaks`],
#' returning a two-column numeric matrix of start and end retention times.
#' @author Ethan Bass
#' @noRd
check_windows <- function(windows, times){
  if (is.null(dim(windows))){
    if (length(windows) != 2L){
      stop("`windows` must be a two-column matrix or data.frame of start and end retention times, or a single pair of times.")
    }
    windows <- matrix(windows, nrow = 1L)
  }
  if (ncol(windows) != 2L){
    stop(sprintf(
      "`windows` must have exactly two columns (start and end retention times), not %d.",
      ncol(windows)))
  }
  # `data.matrix` is deliberately avoided here: on a character matrix it
  # silently returns factor codes rather than `NA`.
  values <- suppressWarnings(as.numeric(as.matrix(windows)))
  if (anyNA(values)){
    stop("Please provide numeric retention times to `windows`.")
  }
  windows <- matrix(values, ncol = 2L)
  bad <- which(windows[, 1] >= windows[, 2])
  if (length(bad) > 0){
    stop(sprintf(
      "The start of each window in `windows` must precede its end. Check window%s %s.",
      plural(length(bad)), paste(bad, collapse = ", ")))
  }
  ts_range <- range(times, na.rm = TRUE)
  oob <- which(windows < ts_range[1] | windows > ts_range[2], arr.ind = TRUE)
  if (nrow(oob) > 0){
    stop(sprintf(
      "The following windows fall outside the bounds of the chromatograms (%s-%s):\n%s",
      signif(ts_range[1], 4), signif(ts_range[2], 4),
      paste(sprintf("         %d (%s-%s)", sort(unique(oob[, "row"])),
                    sQuote(windows[sort(unique(oob[, "row"])), 1]),
                    sQuote(windows[sort(unique(oob[, "row"])), 2])),
            collapse = ", \n")))
  }
  colnames(windows) <- c("start", "end")
  windows
}

#' Convert manual integration windows to peak positions
#'
#' Converts a two-column matrix of start and end retention times into the
#' `pos`/`lower`/`upper` table consumed by [`fit_peaks`], snapping each bound to
#' the nearest time point and locating the apex within each window.
#'
#' The apex must be derived here rather than supplied by the user, since `rt`
#' and `height` are read from it. Two windows sharing a boundary time resolve to
#' the same index, which is what makes a shared boundary behave as a vertical
#' drop line: the trapezoidal sum telescopes over the shared node, so the areas
#' of the two windows sum exactly to the area of their union.
#' @author Ethan Bass
#' @noRd
windows_to_pos <- function(windows, y, times){
  windows <- check_windows(windows, times = times)
  lower <- get_retention_idx(windows[, 1], times = times)
  upper <- get_retention_idx(windows[, 2], times = times)
  apex <- mapply(function(l, u) l - 1L + which.max(y[seq.int(l, u)]), lower, upper)
  data.frame(pos = unname(apex), lower = unname(lower), upper = unname(upper))
}

#' Multiplier from minutes into a given time unit
#' @noRd
time_factor <- function(time_unit = c("min", "s", "ms")){
  time_unit <- match.arg(time_unit)
  switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
}

#' Multiplier from minutes to the units of a retention-time axis
#'
#' An unrecorded or unrecognized unit is taken to be minutes.
#' @noRd
rt_factor <- function(unit){
  if (is.null(unit) || is.na(unit) || !unit %in% c("min", "s", "ms")) return(1)
  time_factor(unit)
}

#' Indices retained when thinning a time axis
#'
#' `time_resolution` must already be expressed in the units of `times`.
#' @noRd
thin_idx <- function(times, time_resolution){
  if (is.null(time_resolution)) return(seq_along(times))
  time_diff <- mean(diff(head(times, 10)))
  seq(1, length(times), by = max(1, round(time_resolution / time_diff)))
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
  if (inherits(x, "list")) x <- unnest_chrom_list(x)
  if (!(inherits(x, "matrix") || (inherits(x, "list") && all(sapply(x, is.matrix)))))
    stop(sprintf("%s should be a matrix or a list of matrices", arg))
  list(chrom_list = x, return_matrix = return_matrix, class = x_class)
}

#' Reduce nested chromatograms to one matrix per sample
#'
#' A list of matrices is returned unchanged. Otherwise `what` (default `"dad"`)
#' selects an element of each sample with `get_element()`.
#' @noRd
unnest_chrom_list <- function(x, what = NULL){
  if (length(x) == 0 || !is.list(x[[1]]) || is.data.frame(x[[1]])) return(x)
  get_element(x, what %||% "dad")
}

#' Extract one element from each chromatogram in a nested list
#'
#' @param x A list of chromatograms, each a list of elements.
#' @param what A single name (e.g. `"dad"`), a vector of names for a nested
#' element, or one string split at the first colon (e.g.
#' `"instrument:PMP1, Pressure"`), unless that string is itself an element name.
#' @return A list of matrices, with the names and class of `x`.
#' @noRd
get_element <- function(x, what){
  colon <- length(what) == 1 && grepl(":", what, fixed = TRUE) &&
    !(what %in% names(x[[1]]))
  if (colon)
    what <- trimws(c(sub(":.*$", "", what), sub("^[^:]*:", "", what)))
  path <- function(w){
    sprintf("`what = %s`", if (colon) deparse(paste(w, collapse = ":")) else deparse(w))
  }
  out <- lapply(x, function(xx){
    el <- if (is.list(xx)) tryCatch(xx[[what]], error = function(e) NULL)
    if (is.matrix(el)) return(el)
    if (is.list(el))
      stop(sprintf("%s is a list, not a matrix. Choose one of its elements, e.g. %s. Elements: %s.",
                   path(what), path(c(what, names(el)[1])),
                   paste(sQuote(names(el), FALSE), collapse = ", ")), call. = FALSE)
    parent <- if (length(what) > 1) tryCatch(xx[[what[-length(what)]]], error = function(e) NULL) else xx
    stop(sprintf("%s was not found. Elements available%s: %s.", path(what),
                 if (length(what) > 1) sprintf(" in %s", path(what[-length(what)])) else "",
                 paste(sQuote(names(parent), FALSE), collapse = ", ")), call. = FALSE)
  })
  class(out) <- class(x)
  out
}

#' Transfer metadata
#' Transfers metadata attributes between objects.
#' @details
#' `exclude` only filters what is copied from `old_object`: an attribute that
#' arithmetic already carried onto `new_object` survives, which is why
#' `subtract_blanks` drops the baseline itself.
#' @noRd
transfer_metadata <- function(new_object, old_object, transfer_class = TRUE,
                              exclude = c('names','row.names','dim','dimnames',
                                          'estimated_baseline')){
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
#' @importFrom stats weighted.mean
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
    fn <- purrr::partial(parallel_lapply, apply_fn = pblapply, cl = cl)
  } else if (!show_progress && is_not_windows && is.numeric(cl) && cl > 1){
      fn <- purrr::partial(parallel_lapply, apply_fn = mclapply, mc.cores = cl)
  } else {
    fn <- purrr::partial(parallel_lapply, apply_fn = lapply)
  }
  fn
}

#' Apply a function in parallel, raising the first failure a worker returned
#'
#' `mclapply` returns a worker's error as a `try-error` element and a worker
#' that died as `NULL`, and every other element scheduled on that worker
#' shares its fate. Warnings are re-raised once per distinct message, naming
#' the elements of `X` that raised it.
#' @noRd
parallel_lapply <- function(X, FUN, ..., apply_fn){
  out <- withCallingHandlers(
    apply_fn(X, collect_warnings, .f = FUN, ...),
    warning = function(cond){
      if (grepl("encountered errors in user code|resulted in an error",
                conditionMessage(cond)))
        invokeRestart("muffleWarning")
    })
  errored <- which(vapply(out, inherits, logical(1), "try-error"))
  if (length(errored) > 0){
    cnd <- attr(out[[errored[1]]], "condition")
    stop(if (is.null(cnd)) as.character(out[[errored[1]]]) else conditionMessage(cnd),
         call. = FALSE)
  }
  lost <- which(vapply(out, is.null, logical(1)))
  if (length(lost) > 0){
    ids <- if (is.null(names(X))) lost else names(X)[lost]
    stop(sprintf(paste(
      "A parallel worker process died, so no result was returned for %s.",
      "Items handled by the same worker are lost with it. If memory ran out,",
      "a lower `cl` may help."), format_list(ids, quote = FALSE)),
      call. = FALSE)
  }
  replay_warnings(lapply(out, `[[`, "warnings"), names(X))
  lapply(out, `[[`, "value")
}

#' Run `.f` on one element, returning its value and the warnings it raised
#' @noRd
collect_warnings <- function(x, ..., .f){
  warnings <- character(0)
  value <- withCallingHandlers(.f(x, ...), warning = function(cond){
    warnings <<- c(warnings, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(warnings))
}

#' Raise each distinct warning once, naming the samples that raised it
#' @noRd
replay_warnings <- function(warnings, ids){
  for (msg in unique(unlist(warnings))){
    if (!is.null(ids)){
      who <- ids[vapply(warnings, function(w) msg %in% w, logical(1))]
      msg <- if (length(ids) > 1 && length(who) == length(ids))
        paste0("All samples: ", msg) else
        sprintf("Sample%s %s: %s", plural(length(who)),
                format_list(who, quote = FALSE), msg)
    }
    warning(msg, call. = FALSE, immediate. = TRUE)
  }
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

#' Format a list of names for a message
#'
#' Shows the first `n` and counts the rest, so a long list cannot push the
#' remainder of a message past `warning.length`.
#' @noRd
format_list <- function(x, n = 10, quote = TRUE){
  x <- if (quote) sQuote(x) else as.character(x)
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(n)], collapse = ", "), " and ", length(x) - n, " more")
}

#' Check for duplicated names
#' @noRd
check_duplicated_names <- function(sample_names, warn = FALSE,
                                   what = "peak_list"){
  if (anyDuplicated(sample_names)) {
    fn <- ifelse(warn, warning, stop)
    dup <- unique(sample_names[duplicated(sample_names)])
    n <- vapply(dup, function(d) sum(sample_names == d), integer(1))
    dups <- format_list(paste0(dup, " (", n, ")"), quote = FALSE)
    fn("Duplicate sample names found in `", what, "`: ", dups, ".\n",
       "Rename them with `names(x) <- make.unique(names(x))`, after deciding ",
       "what to do with any sample that was run more than once.",
       call. = FALSE)
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

#' Calculate apex prominence
#'
#' Computes the height of each apex above the higher of its two flanking
#' minima. Prominence is measured against the *local* surrounding signal, so
#' unlike raw amplitude it is not inflated by a rising baseline or deflated by
#' a neighboring peak.
#'
#' @param y A numeric vector (a chromatographic signal).
#' @param peaks Output of [`find_peaks`], including `pos`, `lower` and `upper`.
#' @return A numeric vector of prominences, one per row of `peaks`.
#' @author Ethan Bass
#' @noRd
get_prominence <- function(y, peaks){
  vapply(seq_nrow(peaks), function(i){
    left <- min(y[peaks$lower[i]:peaks$pos[i]], na.rm = TRUE)
    right <- min(y[peaks$pos[i]:peaks$upper[i]], na.rm = TRUE)
    y[peaks$pos[i]] - max(left, right)
  }, numeric(1))
}

#' Round a smoothing window up to an odd integer
#'
#' An even-width running mean is not centered on a sample: it creates plateaus
#' at extrema, which produce two consecutive sign changes in the derivative and
#' so report an apex one timepoint off. Rounding up to an odd width of at least
#' three keeps the filter centered.
#'
#' @param w A smoothing window, in timepoints.
#' @return An odd integer of at least `3`.
#' @author Ethan Bass
#' @noRd
odd_window <- function(w){
  w <- max(round(w), 3)
  if (w %% 2 == 0) w <- w + 1L
  w
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
  tol <- if (length(gaps)) min(gaps) / 1000 else max(abs(num), 1) * 1e-9
  u <- unique(num)
  for (digits in 0:12){
    if (max(abs(round(u, digits) - u)) < tol)
      return(formatC(num, format = "f", digits = digits))
  }
  formatC(num, format = "f", digits = 12)
}

#' Weighted mean tolerant of non-finite values
#'
#' Computes a weighted mean, ignoring elements with a non-finite value or a
#' non-positive or non-finite weight. Falls back to an unweighted mean of the
#' finite values when no usable weight remains, and returns `NA_real_` only
#' when every value is non-finite.
#'
#' @param x A numeric vector of values.
#' @param w A numeric vector of weights, the same length as `x`.
#' @return A length-one numeric.
#' @author Ethan Bass
#' @noRd
weighted_mean_finite <- function(x, w){
  usable <- is.finite(x) & is.finite(w) & w > 0
  if (any(usable)) return(weighted.mean(x[usable], w = w[usable]))
  finite_x <- is.finite(x)
  if (any(finite_x)) mean(x[finite_x]) else NA_real_
}
