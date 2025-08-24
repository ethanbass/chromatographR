#' Check that peak_table is of proper class
#' 
#' This function validates that the provided object is of the \code{peak_table}
#' class.
#' @author Ethan Bass
#' @noRd
check_peaktable <- function(peak_table){
  if (!inherits(peak_table, "peak_table"))
    stop("The provided peak_table must be of the `peak_table` class")
}

#' Retrieve chrom_list from peak_table object and check properties
#' @param peak_table Peak table object
#' @author Ethan Bass
#' @noRd
get_chrom_list <- function(x, chrom_list, verbose = FALSE){
  if (missing(chrom_list)){
    if (inherits(x, "peak_table")){
      string <- x$args[["chrom_list"]]
    } else if (inherits(x, "peak_list") | inherits(x, "ptw_list")){
      string <- attr(x, "chrom_list")
    }
    if (grepl("\\[*\\]", string)){
      subsetted <- TRUE
      idx <- extract_idx(string)
      string <- gsub("\\[\\[?(.*?)\\]?\\]", "", string)
    } else {
      subsetted <- FALSE
    }
    chrom_list <- try(get(string))
    if (inherits(chrom_list, "try-error")){
      stop("Chromatograms not found! Please make sure the appropriate chrom_list
             object has been loaded")
    }
    if (subsetted){
      chrom_list <- chrom_list[idx]
    }
  }
  check_chrom_list(x, chrom_list, verbose = verbose)
  chrom_list
}

#' Extracts idx from string
#' @noRd
extract_idx <- function(string, chrom_names){
  idx <- sub(".*?\\[(.*?)\\].*", "\\1", string)
  if (any(grepl(":", idx))){
    if (grepl("c(.*?)", idx)){
      matches <- regmatches(idx, gregexpr('\\(.*\\)', idx))[[1]]
      idx <- gsub('[()]', '', matches)
    }
    split <- as.numeric(strsplit(idx, ":")[[1]])
    idx <- split[1]:split[2]
  } else if (grepl("c(.*?)", idx)){
    if (grepl("[']|[\"]", idx)){
      matches <- regmatches(idx, gregexpr("'([^']*)'|\"([^\"]*)\"", idx))[[1]]
      idx <- gsub('[\'\"]', '', matches)
    } else{
      matches <- regmatches(idx, gregexpr('\\(.*\\)', idx))[[1]]
      idx <- gsub('[()]', '', matches)
      idx <- strsplit(idx,",")[[1]]
    }
  }
  if (any(!is.na(suppressWarnings(as.numeric(idx))))){
    idx <- as.numeric(idx)
  }
  idx
}

#' Check that chrom_list is matching a peak table or peak list.
#' @noRd
check_chrom_list <- function(x, chrom_list, verbose = FALSE){
  if (inherits(x, "peak_table")){
    if (length(chrom_list) != nrow(x$tab)){
      stop("Dimensions of chrom_list and peak_table do not match.")
    } else{
      if (verbose & any(names(chrom_list) != rownames(x$tab)))
        warning("Names of chromatograms do not match peak_table")
    }
  } else if (inherits(x, "peak_list")){
    if (length(chrom_list) != length(x)){
      stop("Dimensions of chrom_list and peak_list do not match.")
    } else{
      if (verbose & any(names(chrom_list) != names(x)))
        warning("Names of chromatograms do not match peak_list")
    }
  }
}

#' Get retention index
#' @author Ethan Bass
#' @noRd
get_retention_idx <- function(RT, times){
  if (!is.numeric(RT))
    stop("Retention time not found!")
  if (RT > tail(times, 1) | RT < head(times, 1))
    stop("The supplied retention time falls outside the bounds of the chromatogram.")
  idx <- which.min(abs(RT - times))
  idx
}

#' Check index
#' @author Ethan Bass
#' @noRd
check_idx <- function(idx, chrom_list){
  if (idx > nrow(chrom_list[[1]]))
    stop("The supplied index time exceeds the bounds of the chromatogram.")
}

#' Get wavelength index
#' @param lambda lambda to be found
#' @param lambdas lambdas to be matched
#' @param y Signal as numeric vector.
#' @param allow_max Logical. Whether to return lambda at maximum signal intensity
#' of vector \code{y}. Defaults to TRUE.
#' @author Ethan Bass
#' @noRd
get_lambda_idx <- function(lambda, lambdas, y, allow_max = TRUE){
  if (lambda == 'max'){
    if (allow_max){
      lambda.idx <- which.max(y)
    } else{
      stop("Wavelength (`lambda`) must be specified for interactive scanning.")
    }
  } else{
    lambda.idx <- ifelse((length(lambdas) == 1), 1,
                                which(lambdas == as.numeric(lambda))
    )
  }
  if (is.na(lambda.idx) | length(lambda.idx) == 0)
    stop("The specified wavelength (`lambda`) could not be found!")
  lambda.idx
}

#' Check chromatogram
#' @noRd
check_chr <- function(chr, loc = NULL, peak_table, chrom_list, allow_max = TRUE){
  if (chr[[1]] == 'max'){
    if (allow_max){
      chr <- which.max(peak_table$tab[, loc])
    } else{
      stop("Chromatogram must be specified for scan function.")
    }
  }
  if (is.character(chr)){
    if (!(chr %in% names(chrom_list))){
      stop("Chromatogram not found. Please check `idx` argument and try again!")
    } else{
        chr <- which(names(chrom_list) == chr)
      }
  }
  chr
}

#' Elementwise all equal function
#' @author Brian Diggs
#' @references https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
#' @noRd
elementwise.all.equal <- Vectorize(function(x, y, ...){
  isTRUE(all.equal(x, y, ...))
})

#' Get retention times
#' 
#' Get retention times from a list of chromatograms or a \code{peak_table} object.
#' 
#' @param x A list of chromatograms or \code{peak_table} object.
#' @param idx Index of chromatogram from which to extract times.
#' @return Numeric vector of retention times from the chromatogram specified by
#' \code{idx}.
#' @family utility functions
#' @export
get_times <- function(x, idx = 1){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    as.numeric(rownames(x[[idx]]))
  } else if (inherits(x, "matrix")){
    as.numeric(rownames(x))
  }
}

#' Get lambdas
#' 
#' Get wavelengths from a list of chromatograms or a \code{peak_table} object.
#' 
#' @param x A list of chromatograms or \code{peak_table} object.
#' @return A numeric vector of wavelengths.
#' @family utility functions
#' @export

get_lambdas <- function(x){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    as.numeric(colnames(x[[1]]))
  } else if (inherits(x, "matrix")){
    as.numeric(colnames(x))
  }
}

#' Get time resolution
#' @return Returns average gap between time points.
#' @noRd
get_time_resolution <- function(chrom_list, idx = 1){
  ts <- get_times(x = chrom_list, idx = idx)
  signif(median(diff(ts)))
}

#' Check for suggested package
#' 
#' This function checks for a suggested package and returns an error if the 
#' package is not installed (if \code{return_boolean} is FALSE. Otherwise, it 
#' returns a boolean value.
#' 
#' @noRd
check_for_pkg <- function(pkg, return_boolean = FALSE){
  pkg_exists <- requireNamespace(pkg, quietly = TRUE)
  if (return_boolean){
    return(pkg_exists)
  } else if (!pkg_exists) {
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
}

#' Extract variables from the left-hand-side of a formula.
#' @param formula A \code{\link{formula}} object.
#' @importFrom Formula Formula
#' @noRd
#' @note Adapted from https://github.com/adibender/pammtools/blob/master/R/formula-utils.R
#' (c) Copyright © 2017 Andreas Bender and Fabian Scheipl under MIT license:
#' Permission is hereby granted, free of charge, to any person obtaining a copy of
#' this software and associated documentation files (the “Software”), to deal in 
#' the Software without restriction, including without limitation the rights to use, 
#' copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
#' Software, and to permit persons to whom the Software is furnished to do so, 
#' subject to the following conditions:
#'   
#' The above copyright notice and this permission notice shall be included in all 
#' copies or substantial portions of the Software.
#' 
#' THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#' WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#' CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

get_lhs_vars <- function(formula) {
  if (is.character(formula) ) formula <- as.formula(formula)
  form <- formula(Formula::Formula(formula), lhs = TRUE, rhs = FALSE)
  all.vars(form)
}

#' Transfer metadata
#' Transfers metadata attributes between objects.
#' @noRd
transfer_metadata <- function(new_object, old_object, transfer_class = TRUE,
                              exclude = c('names','row.names','dim','dimnames')){
  if (!transfer_class){
    exclude <- c(exclude, "class")
  }
  a <- attributes(old_object)
  a[exclude] <- NULL
  attributes(new_object) <- c(attributes(new_object), a)
  new_object
}

#' Choose apply function
#' This function chooses an apply function based on arguments provided by the
#' user. The options are \code{lapply}, \code{pblapply} and \code{mclapply}.
#' @importFrom parallel mclapply
#' @return Returns \code{\link[pbapply]{pblapply}} if \code{progress_bar == TRUE},
#' otherwise returns \code{\link{lapply}}.
#' @noRd
choose_apply_fnc <- function(show_progress, parallel = NULL, cl = 2){
  pbapply_installed <- check_for_pkg("pbapply", return_boolean = TRUE)
  is_not_windows <- .Platform$OS.type != "windows"
  
  if (is.null(show_progress)){
    show_progress <- pbapply_installed
  }
  
  if (is.null(parallel)){
    if (!is.null(cl) && 
        (class(cl)[1] == "SOCKcluster" || (is_not_windows && cl > 1))){
      parallel <- TRUE
    } else {
      parallel <- FALSE
    }
  }
  
  if (!parallel){
    cl <- 1
  }

  if (.Platform$OS.type == "windows" && parallel && !pbapply_installed){
    stop("Please install the pbapply package to use parallel processing on Windows.")
  }
  
  if (show_progress){
    check_for_pkg("pbapply")
    pblapply <- pbapply::pblapply
    fn <- purrr::partial(pblapply, cl = cl)
  } else if (!show_progress && is_not_windows && is.numeric(cl) && cl > 1){
      fn <- purrr::partial(mclapply, mc.cores = cl)
  } else {
    fn <- lapply
  }
  fn
}

#' Get y bounds
#' @noRd
get_y_bounds <- function(x, idx, lambdas.idx, pad = 1.1){
  mx <- get_maximum(x, idx, lambdas.idx = lambdas.idx)*pad
  c(0, mx)
}

#' Get maximum
#' @noRd
get_maximum <- function(x, idx, lambdas.idx){
  max(sapply(x[idx], function(xx){
    max(xx[, lambdas.idx], na.rm = TRUE)
  }))
}

#' Get minimum
#' @noRd
get_minimum <- function(x, idx, lambdas.idx){
  min(sapply(x[idx], function(xx){
    min(xx[, lambdas.idx], na.rm = TRUE)
  }))
}

#' Simple cap
#' Copied from \code{toupper} examples.
#' @noRd
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
