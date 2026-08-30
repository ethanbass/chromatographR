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
#' @param zero_floor Logical. If `TRUE` (default), values below zero after
#' subtraction are set to zero.
#' @param method How to select the blank to subtract for each sample:
#'   `"mean"` subtracts the elementwise mean of all blanks from every
#'   sample; `"nearest"` subtracts the blank closest in run order
#'   (list position); `"preceding"` subtracts the closest blank
#'   occurring at or before the sample in run order.
#'
#' @return A `chrom_list` with blank-subtracted samples. Blanks
#'   themselves are removed from the returned list.
#'
#' @export
subtract_blanks <- function(chrom_list, pattern = NULL, blank_id,
                           method = c("mean", "nearest", "preceding"),
                           zero_floor = TRUE) {
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
  sample_id <- setdiff(seq_len(n), blank_id)
  
  if (method == "mean") {
    blank_mean <- Reduce(`+`, chrom_list[blank_id]) / length(blank_id)
    out <- lapply(chrom_list[sample_id], function(x) x - blank_mean)
  } else if (method == "nearest") {
    out <- lapply(sample_id, function(i) {
      nearest <- blank_id[which.min(abs(blank_id - i))]
      chrom_list[[i]] - chrom_list[[nearest]]
    })
  } else if (method == "preceding") {
    out <- lapply(sample_id, function(i) {
      preceding <- blank_id[blank_id <= i]
      if (length(preceding) == 0)
        stop(sprintf("No preceding blank found for sample %d.", i))
      closest <- max(preceding)
      chrom_list[[i]] - chrom_list[[closest]]
    })
  }
  
  out <- Map(function(new, old) transfer_metadata(new, old, transfer_class = FALSE),
             out, chrom_list[sample_id])
  if (zero_floor)
    out <- lapply(out, function(x) pmax(x, 0))
  names(out) <- names(chrom_list)[sample_id]
  class(out) <- class(chrom_list)
  out
}