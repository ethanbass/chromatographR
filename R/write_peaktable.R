#' Export peak table
#' 
#' Exports peak table in `csv` or `xlsx` format.
#' 
#' @param peak_table Peak table object from [`get_peaktable`].
#' @param path Path to an existing directory in which to write the files.
#' @param filename File name. Defaults to `"peak_table"`.
#' @param format File format to export. Either `"csv"` (default) or `"xlsx"`.
#' @param what Which elements of the `peak_table` to export. One or more of
#' `"tab"`, `"pk_meta"`, `"sample_meta"`, `"instrument_meta"`,
#' `"ref_spectra"`, and `"args"`. Defaults to all six.
#' @importFrom utils write.csv
#' @return No return value. The function is called for its side effects.
#' @section Side effects:
#' With `format = "csv"`, writes one file per element of `what`, named
#' `<filename>-<element>.csv`. With `format = "xlsx"`, writes a single
#' `<filename>.xlsx` workbook with one sheet per element. Requires the
#' openxlsx package.
#' @examples \donttest{
#' data(pk_tab)
#' path_out = tempdir()
#' write_peaktable(pk_tab, path = path_out, what = c("tab"))
#' }
#' @export

write_peaktable <- function(peak_table, path, filename = "peak_table",
                            format=c("csv", "xlsx"),
                             what = c("tab", "pk_meta", "sample_meta",
                                      "instrument_meta", "ref_spectra", "args")){
  check_peaktable(peak_table)
  if (missing(path)){
    stop("Please provide a path to write the files.")
  }
  if (!fs::dir_exists(path)){
    stop("The specified directory does not exist.")
  }
  what <- match.arg(what, c("tab", "pk_meta", "sample_meta",
                            "instrument_meta", "ref_spectra",
                            "args"), several.ok = TRUE)
  format <- match.arg(format, c("csv", "xlsx"))
  sheets <- list("tab" = peak_table$tab,
                 "pk_meta" = peak_table$pk_meta,
                 "sample_meta" = peak_table$sample_meta,
                 "instrument_meta" = peak_table$instrument_meta,
                 "ref_spectra" = peak_table$ref_spectra,
                 "args" = as.matrix(vapply(peak_table$args, flatten_arg,
                                          character(1)))
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
                         instrument_meta = "Instrument Metadata",
                         ref_spectra = "Reference Spectra",
                         args = "Arguments")[what]
      path <- fs::path(path, filename, ext = "xlsx")
      writer <- openxlsx::write.xlsx
  } 
  # write sheets
  writer(sheets, file = path, rowNames = TRUE)
}

#' Flatten a peak table argument to one string
#'
#' @param a One element of a peak table's `args`.
#' @return Its values joined by `", "`, or `NA` if it has none.
#' @author Ethan Bass
#' @noRd
flatten_arg <- function(a){
  a <- unlist(a)
  if (!length(a) || all(is.na(a))) NA_character_ else paste(a, collapse = ", ")
}

#' Write CSV from peak_table
#' @noRd
write_csvs <- function(x, file, filename = "peak_table", rowNames){
  # file <- gsub(suffix, "", file)
  x <- lapply(names(x), function(i){
    file <- fs::path(file, paste(filename, i, sep = "-"), ext = "csv")
    write.csv(x = x[[i]],
              file = file,
              row.names = rowNames)
  })
}
