## Function to combine peaks that have been integrated at multiple wavelengths,
## based on retention time and spectral similarity.

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
