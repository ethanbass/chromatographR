#' Import folders of chromatograms in CSV format.
#' 
#' Convenience function to import chromatograms from a list of folders.
#' Chromatograms must be in CSV format.
#' @importFrom utils read.csv
#' @param paths Paths of folders where chromatograms are stored.
#' @param dat Optional list of chromatograms. If list is provided, the function
#' will append chromatograms in path to the existing list of chromatograms.
#' @param ... Additional arguments to \code{\link{read.csv}}.
#' @return A list of chromatograms in matrix format.
#' @author Ethan Bass
#' @export load_chromes
load_chromes <- function(paths, dat=NULL, ...){
  if (is.null(dat)){
  dat<-list()
  }
  for (path in paths){
    files <- list.files(path=path, pattern = "CSV",full.names = TRUE)
    file_names <- gsub(pattern = ".CSV",x = basename(files), replacement = "")
    mydata <- lapply(X=files, FUN=read.csv, row.names = 1, header=TRUE,
                     fileEncoding="utf-16",check.names = FALSE, ...)
    mydata <- lapply(mydata, FUN=as.matrix)
    names(mydata) <- file_names
    dat <- append(dat,mydata)
  } 
  dat
}
