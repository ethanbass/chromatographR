#' Import folders of chromatograms in CSV format.
#' 
#' Convenience function to import chromatograms from a list of folders.
#' Chromatograms must be in CSV format.
#' 
#' 
#' @param paths Paths of folders where chromatograms are stored.
#' @param dat Optional list of chromatograms. If list is provided, the function
#' will append chromatograms in path to the existing list of chromatograms.
#' @return A list of chromatograms in matrix format.
#' @author Ethan Bass
#' @keywords manip
#' @export load_chromes
load_chromes <- function(paths, dat=NULL){
  if (is.null(dat)){
  dat<-list()
  }
  for (path in paths){
    files <- list.files(path=path, pattern = "CSV",full.names = T)
    file_names <- gsub(pattern = ".CSV",x = basename(files), replacement = "")
    mydata <- lapply(X=files, FUN=read.csv, row.names = 1, header=T,
                     fileEncoding="utf-16",check.names = F)
    mydata <- lapply(mydata, FUN=as.matrix)
    names(mydata) <- file_names
    dat <- append(dat,mydata)
  } 
  dat
}
