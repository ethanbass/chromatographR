#' Import chromatograms.
#' 
#' Convenience function to import chromatograms from a list of folders or paths.
#' 
#' Chromatograms may be CSVs, ChemStation \code{.uv} files, or MassHunter \code{
#' .sp} files. Parsers from the \href{https://github.com/bovee/Aston}{Aston}
#' package for python are used to load binary files.
#' 
#' @name load_chroms
#' @importFrom utils read.csv
#' @import chromConverter
#' @param paths Path(s) to chromatograms or the folders containing the files
#' @param find_files Logical. Set to \code{TRUE} (default) if you are providing
#' the function with a folder or vector of folders containing the files.
#' Otherwise, set to\code{FALSE}.
#' @param format.in Format of files.
#' @param sep Argument provided to \code{read.csv}. Defaults to ",".
#' @param dat Optional list of chromatograms. If provided, newly imported
#' chromatograms will be appended to the existing list.
#' @param ... Additional arguments to \code{\link{read.csv}}.
#' @return A list of chromatograms in matrix format.
#' @author Ethan Bass
#' @note Relies on the file parsers from the
#' \href{https://github.com/bovee/aston}{Aston} package to import ChemStation
#' \code{.uv} and MassHunter \code{.sp} files.
#' @examples \dontrun{
#' ###  import from single folder
#' dat <- load_chromes(paths = path)
#' ### import from multiple folders
#' path = 'foo'
#' folders <- list.files(path = path, pattern = "EXPORT3D")
#' dat <- load_chroms(folders)
#' }
#' @export load_chroms

load_chroms <- function(paths, find_files = TRUE,
                         format.in=c("csv", "chemstation", "masshunter"),
                         sep = ",", dat=NULL, ...){
  .Deprecated("read_chroms", package="chromConverter",
              msg="The `load_chroms` function is deprecated. It is suggested to use `read_chroms` in the chromConverter package instead.",
              old = as.character(sys.call(sys.parent()))[1L])
  format.in <- match.arg(format.in, c("csv", "chemstation", "masshunter"))
  exists <- dir.exists(paths) | file.exists(paths)
  if (mean(exists) == 0){
    stop("Cannot locate files. None of the provided paths exist.")
  }
  if (is.null(dat)){
    dat<-list()
  }
  # choose converter
  if (format.in == "csv"){
    pattern <- ".csv|.CSV"
    get_names <- function(files){
      gsub(pattern,"", basename(files))
      }
    converter <- function(file){
      read.csv(file, row.names = 1, header = TRUE,
                          fileEncoding="utf-16",check.names = FALSE, ...)}
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
    converter <- uv_converter
  }
  if (format.in == "masshunter"){
    pattern <- ".sp"
    converter <- sp_converter
  }
  if (find_files){
    files <- unlist(lapply(paths, function(path){
      files <- list.files(path = path, pattern = pattern,
                          full.names = TRUE, recursive = TRUE)
      if (length(files)==0){
        if (!dir.exists(path)){
          warning(paste0("The directory '", basename(path), "' does not exist."))
        } else{
          warning(paste0("No files matching the pattern '", pattern,
                         "' were found in '", basename(path), "'"))
        }
      }
      files
    }))
  } else{
    files <- paths
    match <- grep(pattern, files)
    if (length(match)==0){
      warning("The provided files do not match the expected file extension.
      Please confirm that the specified format ('format.in') is correct.",
              immediate. = TRUE)
    } else if (length(match) < length(files)){
      warning(paste("Some of the files do not have the expected file extension:",
                    files[match]), immediate. = TRUE)
    }
  }
  data <- lapply(X = files, function(f){
    as.matrix(converter(f))
  })
  names(data) <- lapply(X = files, function(f){
    get_names(f)
  })
  dat <- append(dat,data)
  dat
}

#' Read chromatograms.
#' @importFrom chromConverter read_chroms
#' @name read_chroms
#' @export read_chroms
chromConverter::read_chroms