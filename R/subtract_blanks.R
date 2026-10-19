#' Subtract blank chromatograms
#'
#' Subtracts a blank chromatogram from each sample in a `chrom_list`,
#' using either the mean of all blanks, the nearest blank in run order,
#' or the most recent preceding blank.
#'
#' @param chrom_list A `chrom_list` object.
#' @param pattern A regular expression matched against
#'   `names(chrom_list)` to identify blanks (e.g. `"^Blank"` for samples
#'   named `Blank1`, `Blank2`, etc.). Takes precedence over `blank_id`
#'   if both are supplied.
#' @param blank_id Names or indices of the elements in `chrom_list`
#'   that are blanks. Ignored if `pattern` is supplied.
#' @param zero_floor Logical. If `TRUE`, values below zero after subtraction
#' are set to zero, with a warning if this affects more than 5% of a sample's
#' values. Defaults to `FALSE`: a blank that matches the samples' baseline
#' leaves them scattered around zero, and flooring would discard that
#' signal before baseline correction.
#' @param method How to select the blank to subtract for each sample:
#'   `"mean"` subtracts the elementwise mean of all blanks from every
#'   sample; `"nearest"` subtracts the blank closest in run order
#'   (list position), taking the earlier one when two are equally close;
#'   `"preceding"` subtracts the closest blank occurring before the sample
#'   in run order, and raises an error if a sample has none.
#'
#' @details Blanks recorded on a different time grid from a sample (for
#' example, one or two points longer or shorter) are linearly interpolated
#' onto the sample's retention times before subtraction. Where the sample
#' runs past either end of the blank, the blank's first or last value is
#' carried forward, with a warning if this spans more than 5 time points.
#' Blanks and samples must share the same wavelengths.
#'
#' @return A `chrom_list` with blank-subtracted samples. Blanks
#'   themselves are removed from the returned list, and so is any
#'   `"estimated_baseline"` attribute on the samples.
#'
#' @family preprocessing functions
#' @seealso [`preprocess`], [`correct_baseline`]
#' @examples
#' data(Sa_warp)
#' subtract_blanks(Sa_warp, blank_id = 1)
#' @export
subtract_blanks <- function(chrom_list, pattern = NULL, blank_id,
                           method = c("mean", "nearest", "preceding"),
                           zero_floor = FALSE) {
  method <- match.arg(method)
  if (!is.null(pattern)) {
    if (is.null(names(chrom_list)))
      stop("`chrom_list` must have names to use `pattern`.")
    blank_id <- grep(pattern, names(chrom_list))
    if (length(blank_id) == 0)
      stop(sprintf("No elements of `chrom_list` matched pattern %s.", sQuote(pattern)))
  } else if (missing(blank_id)) {
    stop("Either `pattern` or `blank_id` must be supplied.")
  }
  if (is.character(blank_id)) {
    blank_id <- match(blank_id, names(chrom_list))
    if (any(is.na(blank_id)))
      stop("Some `blank_id` values were not found in `names(chrom_list)`.")
  }
  
  n <- length(chrom_list)
  if (is.numeric(blank_id) && (any(blank_id < 1) | any(blank_id > n)))
    stop("Some `blank_id` values are out of range.")
  blank_id <- sort(blank_id)
  sample_id <- setdiff(seq_len(n), blank_id)
  
  pick <- switch(method,
    mean = function(i) blank_id,
    nearest = function(i) blank_id[which.min(abs(blank_id - i))],
    preceding = function(i) {
      preceding <- blank_id[blank_id <= i]
      if (length(preceding) == 0)
        stop(sprintf("No preceding blank found for sample %d.", i))
      max(preceding)
    })
  last <- NULL
  out <- lapply(sample_id, function(i) {
    x <- chrom_list[[i]]
    key <- list(pick(i), dim(x), dimnames(x))
    if (!identical(key, last$key)) {
      blanks <- lapply(chrom_list[key[[1]]], match_time_grid, target = x)
      last <<- list(key = key, blank = Reduce(`+`, blanks) / length(blanks))
    }
    x - last$blank
  })
  blank_range <- lapply(chrom_list, function(b) {
    bt <- rt_values(b)
    if (length(bt) > 0 && !anyNA(bt)) range(bt)
  })
  n_extrap <- vapply(sample_id, function(i) {
    tt <- rt_values(chrom_list[[i]])
    max(vapply(blank_range[pick(i)], function(r) {
      if (is.null(r) || length(tt) == 0 || anyNA(tt)) return(0)
      sum(tt < r[1] | tt > r[2])
    }, numeric(1)))
  }, numeric(1))
  if (any(n_extrap > 5)) {
    ids <- if (is.null(names(chrom_list))) sample_id else names(chrom_list)[sample_id]
    warning(sprintf(paste("Blanks were extrapolated over more than 5 time points",
                          "for: %s. Check that the blanks cover each sample's",
                          "retention time range."),
                    format_list(ids[n_extrap > 5], quote = FALSE)), call. = FALSE)
  }
  out <- Map(function(new, old){
    attr(new, "estimated_baseline") <- NULL
    transfer_metadata(new, old, transfer_class = FALSE)
  }, out, chrom_list[sample_id])
  if (zero_floor){
    floored <- vapply(out, function(x) mean(x < 0, na.rm = TRUE), numeric(1))
    out <- lapply(out, function(x) pmax(x, 0))
    if (any(floored > 0.05))
      warning(sprintf(paste("`zero_floor` set more than 5%% of values to zero",
                            "in %d of %d samples (up to %.0f%%). Set",
                            "`zero_floor = FALSE` to keep them for baseline",
                            "correction."),
                      sum(floored > 0.05), length(out), 100 * max(floored)),
              call. = FALSE)
  }
  names(out) <- names(chrom_list)[sample_id]
  class(out) <- class(chrom_list)
  out
}

#' @noRd
rt_values <- function(x) suppressWarnings(as.numeric(rownames(x)))

#' @noRd
match_time_grid <- function(blank, target) {
  if (!identical(colnames(blank), colnames(target)) ||
      ncol(blank) != ncol(target))
    stop("Blanks and samples must have the same wavelengths.")
  bt <- rt_values(blank)
  tt <- rt_values(target)
  if (nrow(blank) == nrow(target) &&
      (identical(bt, tt) || anyNA(bt) || anyNA(tt)))
    return(blank)
  if (length(bt) == 0 || length(tt) == 0 || anyNA(bt) || anyNA(tt))
    stop("Chromatograms of different lengths need numeric retention times as rownames.")
  if (length(bt) < 2)
    stop("Blanks need at least two time points to be matched to a sample.")
  j <- findInterval(tt, bt, all.inside = TRUE)
  w <- pmin(pmax((tt - bt[j]) / (bt[j + 1] - bt[j]), 0), 1)
  unname(blank[j, , drop = FALSE] * (1 - w) + blank[j + 1, , drop = FALSE] * w)
}
