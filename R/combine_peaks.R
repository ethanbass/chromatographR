#' Combine peaks in peak table
#' 
#' Utility function to combine duplicate peaks in peak table, i.e. peaks that
#' were integrated at more than one wavelength or component. Specify number of
#' digits to match retention time and minimum spectral correlation for a match.
#' 
#' 
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param digits Number of digits to use for retention time matching.
#' @param min.cor Minimum spectral correlation to confirm a match.
#' @param choose If "max" will retain peak with highest intensity. Otherwise,
#' the first column in the dataframe will be retained.
#' @return A peak table similar to the input peak table, but with duplicate
#' columns combined according to the specified criteria.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @export combine_peaks
combine_peaks <- function(peak_table, digits=2, min.cor=0.9, choose='max'){
  RTs <- round(as.numeric(peak_table['RT',]), digits = digits)
  d <- unique(RTs[duplicated(RTs)])
  peak_table.s <- peak_table
  for (rt in d){
    i <- which(RTs == rt)
    cors <- cor(peak_table[4:nrow(peak_table),i])
    j <- which(cors[,1] > min.cor)
    if(length(j)>1){
      if (choose == 'max'){
      sub <- names(sort(colSums(peak_table[,i]),decreasing = TRUE)[-1])
      }
      sub <- which(colnames(peak_table.s) %in% sub)
      peak_table.s <- peak_table.s[,-sub]
    }
  }
  peak_table.s
}
