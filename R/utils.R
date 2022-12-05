#' Check that peaktable is of proper class
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
  if (inherits(x, "peak_table")){
    if (missing(chrom_list)){
      chrom_list <- try(get(x$args[["chrom_list"]]))
      if (inherits(chrom_list, "try-error")){
        stop("Chromatograms not found! Please make sure the appropriate chrom_list
             object has been loaded")
      }
    }
    if (length(chrom_list) != nrow(x$tab)){
      stop("Dimensions of chrom_list and peak_table do not match.")
    } else{
      if (verbose & any(names(chrom_list) != rownames(x$tab)))
        warning("Names of chromatograms do not match peak_table")
    }
  } else if (inherits(x, "peak_list")){
    if (missing(chrom_list)){
      chrom_list <- try(get(attr(x, "chrom_list")))
      if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
    }
    if (length(chrom_list) != length(x)){
      stop("Dimensions of chrom_list and peak_list do not match.")
    } else{
      if (verbose & any(names(chrom_list) != names(x)))
        warning("Names of chromatograms do not match peak_list")
    }
  }
  chrom_list
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
    lambda.idx <- ifelse(length(lambdas == 1), 1,
                                which(lambdas == as.numeric(lambda))
    )
  }
  if (length(lambda.idx) == 0)
    stop("The specified wavelength (`lambda`) could not be found!")
  lambda.idx
}

check_chr <- function(chr, loc=NULL, peak_table, chrom_list, allow_max = TRUE){
  if (chr == 'max'){
    if (allow_max){
      chr <- which.max(peak_table$tab[,loc])
    } else{
      stop("Chromatogram must be specified for scan function.")
    }
  }
  if (is.character(chr) & !(chr %in% names(chrom_list))){
    stop("Chromatogram not found. Check `chr` argument!")
  }
  chr
}

#' Elementwise all equal function
#' @author Brian Diggs
#' @references https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
#' @noRd
elementwise.all.equal <- Vectorize(function(x, y, ...) {isTRUE(all.equal(x, y, ...))})

#' @noRd
get_times <- function(x, index=1){
  if (inherits(x, "chrom_list") | inherits(x, "list")){
    as.numeric(rownames(x[[index]]))
  } else if (inherits(x, "matrix")){
    as.numeric(rownames(x))
  }
}

#' @noRd
get_lambdas <- function(chrom_list){
  as.numeric(colnames(chrom_list[[1]]))
}

#' @noRd
get_time_resolution <- function(chrom_list){
  signif(median(diff(as.numeric(rownames(chrom_list[[1]])))))
}

#' @noRd
check_for_pkg <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
}

#' Extract variables from the left-hand-side of a formula
#'
#' @param formula A \code{\link{formula}} object.
#' @importFrom Formula Formula
#' @noRd
#' @note Adapted from https://github.com/adibender/pammtools/blob/master/R/formula-utils.R
get_lhs_vars <- function(formula) {
  if (is.character(formula) ) formula <- as.formula(formula)
  form <- formula(Formula::Formula(formula), lhs = TRUE, rhs = FALSE)
  all.vars(form)
}
