#' Combine peaks in peak table
#' 
#' Utility function to combine duplicate peaks in peak table, i.e. peaks that
#' were integrated at more than one wavelength or component. Specify tolerance
#' (\code{tol}) for retention time matching and minimum spectral correlation
#' (\code{min.cor}) for a match.
#' 
#' @name combine_peaks
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param tol Tolerance for matching retention times (maximum retention time
#' difference). Defaults to \code{.01}.
#' @param min.cor Minimum spectral correlation to confirm a match. Defaults to 
#' \code{0.9}.
#' @param choose If \code{max} will retain peak with highest intensity. Otherwise,
#' the first column in the data.frame will be retained.
#' @param verbose Logical. Whether to print status to the console.
#' @return A peak table similar to the input peak table, but with duplicate
#' columns combined according to the specified criteria.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @examples
#' data(pk_tab)
#' data(Sa_warp)
#' pk_tab <- attach_ref_spectra(pk_tab)
#' combine_peaks(pk_tab, tol = .02, min.cor = .9)
#' @export combine_peaks

combine_peaks <- function(peak_table, tol = .01, min.cor = 0.9,
                          choose = 'max', verbose = getOption("verbose")){
  check_peaktable(peak_table)
  if (!(is.data.frame(peak_table$ref_spectra) | is.matrix(peak_table$ref_spectra))){
    stop("No reference spectra found. Use attach_ref_spectra function first.")
  }
  RTs <- as.numeric(peak_table$pk_meta['rt',])
  compare_rts <- function(rt1, rt2, tol){
    abs(rt1 - rt2) < tol
  }
  mat <- outer(RTs, RTs, compare_rts, tol = tol)
  # find columns with a retention time match
  d <- which(apply(mat, 2, sum) > 1)
  cors <- suppressWarnings(cor(peak_table$ref_spectra))
  if (verbose) warnings()
  # iterate over columns with retention time match
  # find columns to remove (sub)
  sub <- sapply(d, function(e){
    i <- which(mat[,e] == 1)
    # compare spectral correlation among rt matches
    j <- which(cors[i,e] > min.cor)
    k <- i[j]
    if (length(k) > 1){
      if (choose == 'max'){
        sub <- names(sort(colSums(peak_table$tab[,k]), decreasing = TRUE)[-1])
      }
      which(colnames(peak_table$tab) %in% sub)
    }
  })
  if (length(sub) > 0){
    peak_table$tab <- peak_table$tab[, -unlist(sub)]
    peak_table$pk_meta <- peak_table$pk_meta[, -unlist(sub)]
  }
  if (verbose){
    message(paste("Removed "), length(sub), " peaks from peak table.")
  }
  peak_table
}

#' Merge split peaks in peak table
#' 
#' Merges the specified peaks, by selecting the largest value from each column.
#' Utility function to combine split peaks into a single column of the peak table.
#'
#' @name merge_peaks
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param peaks A vector specifying the names or indices of peaks to be merged.
#' @param method Method to merge peaks. Either \code{max} to select the largest
#' peak from each sample or \code{sum} to sum the peaks together.
#' @return A peak table similar to the input peak table, but where the specified
#' columns are combined. 
#' @author Ethan Bass
#' @examples
#' data(pk_tab)
#' pk_tab <- merge_peaks(peak_table = pk_tab, peaks=c("V10","V11"))
#' @export

merge_peaks <- function(peak_table, peaks, method = c("max", "sum")){
  check_peaktable(peak_table)
  method <- match.arg(method, c("max", "sum"))
  if (is.character(peaks)){
    pks.idx <- which(colnames(peak_table$tab) %in% peaks)
  } else {
    pks.idx <- peaks
  }
  sel <- which.max(colMeans(peak_table$tab[, pks.idx], na.rm = TRUE))
  sel.idx <- which(colnames(peak_table$tab) == peaks[sel])
  if (method == "max"){
    peak_table$tab[[sel.idx]] <- do.call(pmax, peak_table$tab[, pks.idx])
  } else if (method == "sum"){
    peak_table$tab[[sel.idx]] <- apply(pk_tab$tab[,pks.idx], 1, sum)
  }
  peak_table$tab <- peak_table$tab[, -pks.idx[-sel]]
  peak_table$pk_meta <- peak_table$pk_meta[, -pks.idx[-sel]]
  if (inherits(peak_table$ref_spectra,"matrix")){
    peak_table$ref_spectra <- peak_table$ref_spectra[, -pks.idx[-sel]]
  }
  return(peak_table)
}
