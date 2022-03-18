#' Attach experimental metadata
#' 
#' Attach experimental meta-data to peak table.
#' 
#' @param peak_table Peak table
#' @param metadata Dataframe
#' @param column Name of the column with sample names.
#' @return Peak table with meta-data.
#' @author Ethan Bass
#' @seealso \code{\link[get_peaktable}}
#' @export attach_metadata

attach_metadata <- function(peak_table, metadata, column){
  if (!(column %in% colnames(metadata)))
    stop(paste0("Column, ", column, ", is not found."))
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[,column] <- as.character(metadata[,column])
  meta <- merge(meta, metadata, all.x=TRUE, all.y=FALSE, sort=FALSE)
  peak_table$sample_meta <- meta
  return(peak_table)
}
