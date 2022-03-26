#' Combine peaks in peak table
#' 
#' Utility function to combine duplicate peaks in peak table, i.e. peaks that
#' were integrated at more than one wavelength or component. Specify number of
#' digits to match retention time and minimum spectral correlation for a match.
#' 
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param tol Tolerance for matching retention times.
#' @param min.cor Minimum spectral correlation to confirm a match.
#' @param choose If "max" will retain peak with highest intensity. Otherwise,
#' the first column in the data.frame will be retained.
#' @return A peak table similar to the input peak table, but with duplicate
#' columns combined according to the specified criteria.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @examples \dontrun{
#' combine_peaks(pk_tab, tol = .02, min.cor = .9)}
#' @export combine_peaks
combine_peaks <- function(peak_table, tol=.01, min.cor=0.9, choose='max'){
  RTs <- as.numeric(peak_table$pk_meta['RT',])
  mat <- outer(RTs, RTs, elementwise.all.equal, tolerance=0.01, scale=1)
  # diag(mat) <- 0
  d <- which(apply(mat,2, sum)>1)
  sub <- sapply(d, function(e){
    i <- which(mat[,e]==1)
    cors <- cor(peak_table$tab[,i])
    j <- which(cors[,1] > min.cor)
    if(length(j) > 1){
      if (choose == 'max'){
        sub <- names(sort(colSums(peak_table$tab[,i]),decreasing = TRUE)[-1])
      }
      which(colnames(peak_table$tab) %in% sub)
      # peak_table$tab <- peak_table$tab[,-sub]
    }
  })
  peak_table$tab <- peak_table$tab[, -unlist(sub)]
  peak_table$pk_meta <- peak_table$pk_meta[, -unlist(sub)]
  return(peak_table)
}

# sub <- rbind(
#   apply(peak_table$tab[,i], 2, max),
#   peak_table$pk_meta["Component",i])
# names(sort(colSums(peak_table$tab[,i]), decreasing = TRUE)[-1]),
# peak_table$pk_meta["Component",i])