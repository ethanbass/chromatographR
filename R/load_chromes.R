#' Import folders of chromatograms in CSV format.
#' 
#' Convenience function to import chromatograms from a list of folders.
#' Chromatograms must be in CSV format.
#' @importFrom utils read.csv
#' @param paths Path(s) to folders where chromatograms are stored.
#' @param format.in Format of files.
#' @param sep Argument provided to \code{read.csv}. Defaults to ",".
#' @param dat Optional list of chromatograms. If list is provided, the function
#' will append newly imported chromatograms to the existing list.
#' @param ... Additional arguments to \code{\link{read.csv}}.
#' @return A list of chromatograms in matrix format.
#' @author Ethan Bass
#' @note Relies on the file parsers from the
#' \url{https://github.com/bovee/aston}{Aston} package to import chemstation
#' \code{.uv} and masshunter \code{.sp} files.
#' @examples \dontrun{
#' ###  import from single folder
#' dat <- load_chromes(paths = path)
#' ### import from multiple folders
#' path = 'foo'
#' folders <- list.files(path = path, pattern = "EXPORT3D")
#' dat <- load_chromes(folders)
#' }
#' @export load_chromes

load_chromes <- function(paths,
                         format.in=c("csv", "chemstation", "masshunter"),
                         sep = ",", dat=NULL, ...){
  format.in = match.arg(format.in, c("csv", "chemstation", "masshunter"))
  if (is.null(dat)){
    dat<-list()
  }
  exists <- dir.exists(paths)
  if (mean(exists) == 0){
    stop("None of the provided paths exist.")
  }
  if (format.in=="csv"){
    pattern <- ".csv|.CSV"
    get_names <- function(files){
      gsub(pattern,"", basename(files))
      }
    converter <- function(file){
      read.csv(file, row.names = 1, header=TRUE,
                          fileEncoding="utf-16",check.names = FALSE, ...)}
  }
  if (format.in!="csv" & !requireNamespace("chromConverter", quietly = TRUE)){
      stop("To import binary formats, you must install chromConverter. You can use
      devtools::install_github(https://github.com/ethanbass/chromConverter).")
  }
  if (format.in %in% c("chemstation", "masshunter")){
    get_names <- function(files){
      file_names <- sapply(files, function(f){
      split_path <- strsplit(f,"/")[[1]]
      split_path[grep("\\.[Dd]", split_path)]
      })
      gsub("\\.[Dd]","", file_names)
      }
  }
  if (format.in == "chemstation"){
    pattern <- ".uv"
    converter <- chromConverter:::trace_converter
  }
  if (format.in == "masshunter"){
    pattern <- ".sp"
    converter <- chromConverter:::sp_converter
  }
  for (path in paths){
    files <- list.files(path=path, pattern = pattern, full.names = TRUE,
                        recursive=TRUE)
    file_names <- get_names(files)
    if (length(files)==0){
      if (!dir.exists(path)){
        warning(paste0("The directory '", basename(path), "' does not exist."))
    } else{
      warning(paste0("No files matching the pattern '", pattern, "' were found in '", basename(path), "'"))
    }
    }
    data <- lapply(X=files, FUN=converter)
    data <- lapply(data, FUN=as.matrix)
    names(data) <- file_names
    dat <- append(dat,data)
  }
  dat
}

