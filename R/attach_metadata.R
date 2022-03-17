#' Attach experimental metadata
#' 
#' Attach experimental metadata to peak table.
#' 
#' @param peak_table Peak table
#' @param metadata Dataframe
#' @param column Name of the column with sample names.
#' @return Peak table with metadata added
#' @author Ethan Bass
#' @seealso \code{\link[ptw:ptw]{ptw}}, \code{\link{correctPeaks}}
#' @references Eilers, P.H.C. 2004.a
#' @export attach_metadata

attach_meta <- function(peak_table, metadata, column){
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[,column] <- as.character(metadata[,column])
  meta <- merge(meta, metadata, all.x=TRUE, all.y=FALSE, sort=FALSE)
  peak_table$sample_meta <- meta
  return(peak_table)
}
