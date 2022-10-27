#' Export peak table
#' @param peak_table Peak table object from 
#' @param path path to write file
#' @param format format for export. Either \code{csv} or \code{xlsx}.
#' @param what which elements of the \code{peak_table} to export.
#' @importFrom utils write.csv
#' @export

write_peaktable <- function(peak_table, path, format=c("csv", "xlsx"),
                             what = c("tab", "pk_meta",
                                      "sample_meta", "ref_spectra",
                                      "args")){
  check_peaktable(peak_table)
  what <- match.arg(what, c("tab", "pk_meta",
                            "sample_meta", "ref_spectra",
                            "args"), several.ok=TRUE)
  format <- match.arg(format, c("csv", "xlsx"))
  
  sheets <- list("tab" = peak_table$tab,
                 "pk_meta" = peak_table$pk_meta,
                 "sample_meta" = peak_table$sample_meta,
                 "ref_spectra" = peak_table$ref_spectra,
                 "args" = data.frame(value=peak_table$args))
  sheets <- sheets[what]
  
  path <- check_path(path, format)
  
  if (format == "csv"){
    suffix <- "csv"
    writer <- write_csvs
  } else if (format == "xlsx"){
    check_for_pkg("openxlsx")
      names(sheets) <- c("Peak Table", "Peak Metadata", "Sample Metadata",
                         "Reference Spectra", "Arguments")[what]
    suffix = "xlsx"
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

#' @noRd
write_csvs <- function(x, file, rowNames){
  suffix <- ".csv"
  file <- gsub(suffix,"", file)
  lapply(names(x), function(i){
    file <- paste0(paste(file, i, sep="-"), ".csv")
    write.csv(x = x[[i]],
              file = file,
              row.names = rowNames)
  })
}

#' check path
#' check that path is properly formatted
#' @param path path as character string
#' @noRd
check_path <- function(path, format){
  if (.Platform$OS.type != "windows"){
    # check for leading slash
    if (!(substr(path, 1, 1) %in% c("/", "~"))){
      path <- paste0("/", path)
    }
  }
  no_suffix <- !grepl(format, path, ignore.case = TRUE)
  if (no_suffix){
    # check if path leads to directory
    if (dir.exists(path)){
      # check for trailing slash
      n <- nchar(path)
      if (substr(path, n, n) != "/"){
        path <- paste0(path, "/")
      }
      # add file name
      path <- paste0(path, "peak_table")
    }
    # add suffix
    path <- paste(path, format, sep = ".")
  }
  path
}
