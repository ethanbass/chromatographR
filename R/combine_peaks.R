#' Combine peaks in peak table
#' 
#' Utility function to combine duplicate peaks in peak table, i.e. peaks that
#' were integrated at more than one wavelength or component. Specify tolerance
#' (\code{tol}) for retention time matching and minimum spectral correlation
#' (\code{min.cor}) for a match.
#' 
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param tol Tolerance for matching retention times (maximum retention time
#' difference).
#' @param min.cor Minimum spectral correlation to confirm a match.
#' @param choose If "max" will retain peak with highest intensity. Otherwise,
#' the first column in the data.frame will be retained.
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
combine_peaks <- function(peak_table, tol=.01, min.cor=0.9, choose='max'){
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
  cors <- cor(peak_table$ref_spectra)
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
  if (length(sub > 0)){
  peak_table$tab <- peak_table$tab[, -unlist(sub)]
  peak_table$pk_meta <- peak_table$pk_meta[, -unlist(sub)]
  }
  message(paste("Removed "), length(sub), " peaks from peak table.")
  peak_table
}
