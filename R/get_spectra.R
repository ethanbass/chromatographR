#' Extract a spectrum from a chromatogram
#'
#' Extracts the UV-Vis spectrum at a specified retention time from one or more
#' chromatograms in a `chrom_list`.
#'
#' @param loc Peak or retention time to extract the spectrum at. Interpretation
#' depends on `what`.
#' @param peak_table A `peak_table` object. Required if `what == "peak"`.
#' @param chrom_list A list of chromatographic matrices.
#' @param idx Indices of chromatograms to extract spectra from. Defaults to all
#' chromatograms in `chrom_list`.
#' @param lambda Wavelength at which to extract the spectrum. Use `"max"` to
#' select the wavelength with the highest absorbance.
#' @param scale_spectrum Logical. Whether to normalize absorbance values to
#' the range \[0, 1\]. Defaults to `FALSE`.
#' @param what How to interpret `loc`. One of `"peak"` (look up retention time
#' from `peak_table`), `"rt"` (treat `loc` as a retention time directly), or
#' `"idx"` (treat `loc` as a row index).
#' @param format One of `"wide"` or `"long"`. In wide format, wavelengths are
#' rows and each column is a sample. In long format, columns are `wavelength`,
#' `absorbance`, `sample`, and `rt`. Defaults to `"wide"`.
#'
#' @return A `data.frame` in the format specified by `format`.
#'
#' @export

get_spectra <- function(loc, peak_table, chrom_list = NULL, 
                        idx = seq_along(chrom_list),
                        lambda = "max",
                        scale_spectrum = FALSE,
                        what = c("peak", "rt", "idx"),
                        format = c("wide", "long")){
  if (is.null(chrom_list) & missing(peak_table))
    stop("Must provide either a peak_table or a chrom_list.")
  if (!missing(peak_table))
    check_peaktable(peak_table)
  chrom_list <- get_chrom_list(if (missing(peak_table)) NULL else peak_table, chrom_list)
  if (!(inherits(chrom_list, c("chrom_list", "list", "matrix"))))
    stop("The provided `chrom_list` does not appear to be valid. 
                          Please check `chrom_list` argument.")
  what <- match.arg(what)
  format <- match.arg(format)
  new.ts <- get_times(chrom_list, idx = idx[1])
  new.lambdas <- get_lambdas(chrom_list)
  sig <- max(nchar(gsub(".*\\.","",rownames(chrom_list[[1]]))))
  
  if (what == "peak"){
    if (is.null(peak_table))
      stop("A `peak_table` must be provided when `what == 'peak'`.")
    RT <- round(as.numeric(peak_table$pk_meta["rt", loc]), sig)
  } else if (what == "rt"){
    RT <- round(as.numeric(loc), sig)
  } else if (what == "idx"){
    check_idx(loc, chrom_list)
    RT <- new.ts[loc]
  }
  
  sp <- lapply(idx, function(chr){
    row.idx <- get_retention_idx(RT, times = get_times(chrom_list, idx = chr))
    y <- unlist(chrom_list[[chr]][row.idx, , drop = TRUE])
    if (all(is.na(y)))
      stop(sprintf("No data found at RT %g in chromatogram %d.", RT, chr))
    if (scale_spectrum)
      y <- scales::rescale(y)
    y
  })
  
  wide <- as.data.frame(do.call(cbind, sp))
  rownames(wide) <- new.lambdas
  colnames(wide) <- names(chrom_list)[idx]
  
  if (format == "wide"){
    return(wide)
  } else {
    long <- data.frame(
      wavelength = rep(new.lambdas, times = length(idx)),
      absorbance = unlist(sp),
      sample = rep(names(chrom_list)[idx], each = length(new.lambdas)),
      rt = RT
    )
    return(long)
  }
}
