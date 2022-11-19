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
      if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
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

get_times <- function(chrom_list){
  as.numeric(rownames(chrom_list[[1]]))
}

get_lambdas <- function(chrom_list){
  as.numeric(colnames(chrom_list[[1]]))
}