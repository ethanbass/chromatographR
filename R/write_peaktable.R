#' Export peak table
#' @param peak_table Peak table object from \code{\link{get_peaktable}}.
#' @param path Path to write file.
#' @param filename File name. Defaults to "peak_table".
#' @param format File format to export. Either \code{csv} or \code{xlsx}.
#' @param what Which elements of the \code{peak_table} to export.
#' @importFrom utils write.csv
#' @return No return value. The function is called for its side effects.
#' @section Side effects:
#' Exports peak_table object as \code{.csv} or \code{.xlsx} file according to the value
#' of \code{format}.
#' @examples \donttest{
#' data(pk_tab)
#' path_out = tempdir()
#' write_peaktable(pk_tab, path = path_out, what = c("tab"))
#' }
#' @export

write_peaktable <- function(peak_table, path, filename = "peak_table",
                            format=c("csv", "xlsx"),
                             what = c("tab", "pk_meta",
                                      "sample_meta", "ref_spectra",
                                      "args")){
  check_peaktable(peak_table)
  if (missing(path)){
    stop("Please provide a path to write the files.")
  }
  if (!fs::dir_exists(path)){
    stop("The specified directory does not exist.")
  }
  what <- match.arg(what, c("tab", "pk_meta",
                            "sample_meta", "ref_spectra",
                            "args"), several.ok=TRUE)
  format <- match.arg(format, c("csv", "xlsx"))
  sheets <- list("tab" = peak_table$tab,
                 "pk_meta" = peak_table$pk_meta,
                 "sample_meta" = peak_table$sample_meta,
                 "ref_spectra" = peak_table$ref_spectra,
                 "args" = as.matrix(peak_table$args)
  )

  sheets <- sheets[what]
  
  if (format == "csv"){
    writer <- purrr::partial(write_csvs, filename = filename)
  } else if (format == "xlsx"){
      check_for_pkg("openxlsx")
      sheets$args[is.na(sheets$args)]<-''
      names(sheets) <- c(tab = "Peak Table",
                         pk_meta = "Peak Metadata",
                         sample_meta = "Sample Metadata",
                         ref_spectra = "Reference Spectra",
                         args = "Arguments")[what]
      path <- fs::path(path, filename, ext = "xlsx")
      writer <- openxlsx::write.xlsx
  } 
  # else if (format == "html" & !requireNamespace("ymlthis", quietly = TRUE)) {
  #   stop(
  #     "Package `ymlthis` must be installed to export html peak table:
  #     try `install.packages('ymlthis')`",
  #     call. = FALSE
  #   )
  #   suffix = "html"
  # }
  
  # write sheets
  writer(sheets, file = path, rowNames = TRUE)
}

#' Write CSV from peak_table
#' @noRd
write_csvs <- function(x, file, filename="peak_table", rowNames){
  # file <- gsub(suffix, "", file)
  x <- lapply(names(x), function(i){
    file <- fs::path(file, paste(filename, i, sep="-"), ext = "csv")
    write.csv(x = x[[i]],
              file = file,
              row.names = rowNames)
  })
}
