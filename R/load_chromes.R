#' Import folders of chromatograms in CSV format.
#' 
#' Convenience function to import chromatograms from a list of folders.
#' Chromatograms must be in CSV format.
#' @importFrom utils read.csv
#' @param paths Path(s) to folders where chromatograms are stored.
#' @param pattern Pattern to identify files (i.e. a file extension).
#' Defaults to 'CSV'.
#' @param dat Optional list of chromatograms. If list is provided, the function
#' will append newly imported chromatograms to the existing list.
#' @param ... Additional arguments to \code{\link{read.csv}}.
#' @return A list of chromatograms in matrix format.
#' @author Ethan Bass
#' @examples \dontrun{
#' ###  import from single folder
#' dat <- load_chromes(paths = path)
#' ### import from multiple folders
#' path = 'foo'
#' folders <- list.files(path = path, pattern = "EXPORT3D")
#' dat <- load_chromes(folders)
#' }
#' @export load_chromes
load_chromes <- function(paths, pattern="CSV", dat=NULL, ...){
  if (is.null(dat)){
  dat<-list()
  }
  exists <- dir.exists(paths)
  if (mean(exists) == 0){
    stop("None of the provided paths exist.")
  }
  for (path in paths){
    files <- list.files(path=path, pattern = pattern, full.names = TRUE)
    if (length(files)==0){
      if (!dir.exists(path)){
        warning(paste0("The directory '", basename(path), "' does not exist."))
    } else{
      warning(paste0("No files matching the pattern '", pattern, "' were found in '", basename(path), "'"))
    }
    }
    file_names <- gsub(pattern = ".CSV",x = basename(files), replacement = "")
    mydata <- lapply(X=files, FUN=read.csv, row.names = 1, header=TRUE,
                     fileEncoding="utf-16",check.names = FALSE, ...)
    mydata <- lapply(mydata, FUN=as.matrix)
    names(mydata) <- file_names
    dat <- append(dat,mydata)
  }
  dat
}
