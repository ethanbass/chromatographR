#' Reshape chromatograms
#' 
#' Converts a list of chromatographic matrices into a single long-format
#' data frame with one row per sample × retention time × wavelength combination.
#'
#' @details Each row corresponds to a single measurement of signal intensity at
#' a given retention time and wavelength for a specific sample.
#' @name reshape_chroms
#' @param x A list of chromatographic matrices in wide format.
#' @param idx Indices of chromatograms to convert.
#' @param time_resolution Time resolution used when reshaping. Can be used
#' to subsample the time axis. Defaults to full resolution.
#' @param sample_var Name of the new column containing sample identifiers
#' @param lambdas Wavelength(s) to include.
#' @param rts Retention times to include.
#' @param transfer_metadata Logical. If `TRUE`, metadata attributes are
#' transferred to the output. Defaults to `FALSE`.
#' @return A data frame in long format with columns `rt`, `lambda`, `absorbance`,
#' and a sample identifier column specified by `sample_var`.
#' @author Ethan Bass
#' @family utility functions
#' @examples
#' data(Sa_warp)
#' reshape_chroms(Sa_warp, lambdas = 210)
#' @export

reshape_chroms <- function(x, idx, time_resolution = NULL,  
                           sample_var = "sample", lambdas = NULL, rts = NULL, 
                           transfer_metadata = FALSE){
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (missing(lambdas)){
    lambdas <- colnames(x[[1]])
  } else{
    lambdas.idx <- sapply(lambdas, function(lambda){
      get_lambda_idx(lambda, lambdas = get_lambdas(x),
                     allow_max = FALSE)
    }) 
  }
  dat <- lapply(idx, function(i){
    xx <- reshape_chrom(x = x[[i]], time_resolution = time_resolution,
                        lambdas = lambdas.idx, rts = rts,
                        transfer_metadata)
    xx[, sample_var] <- names(x)[[i]]
    xx
  })
  df <- do.call(rbind, dat)
  df$sample <- as.factor(df$sample)
  df
}

#' Reshapes a single chromatogram from wide to long format
#' @name reshape_chrom
#' @param x A chromatographic matrix in wide format.
#' @param lambdas Vector specifying wavelength(s) to include.
#' @param rts Vector specifying retention times to include.
#' @return A chromatographic matrix in long format.
#' @author Ethan Bass
#' @noRd
reshape_chrom <- function(x, time_resolution = NULL, lambdas = NULL, rts = NULL, 
                          transfer_metadata = FALSE){
  xx <- as.data.frame(x)
  if (!is.null(lambdas)){
    xx <- xx[, lambdas, drop = FALSE]
  }
  if (!is.null(rts)){
    times <- get_times(x = x)
    if (is.character(rts)){
      rts <- as.numeric(rts)
    }
    rts.idx <- sapply(rts, function(rt){
      get_retention_idx(RT = rt, times = times)})
    xx <- xx[rts.idx, , drop = FALSE]
  }
  if (!is.null(time_resolution)){
    times <- get_times(xx)
    time_diff <- mean(diff(times[1:10]))
    thin_factor <- round(time_resolution / time_diff)
    keep_idx <- seq(1, length(times), by = thin_factor)
    xx <- xx[keep_idx, , drop = FALSE]
  }
  data <- data.frame(
    rt = as.numeric(rep(rownames(xx), ncol(xx))),
    lambda = as.numeric(rep(colnames(xx), each = nrow(xx))),
    absorbance = as.vector(as.matrix(xx)),
    row.names = NULL
  )
  if (transfer_metadata){
    data <- transfer_metadata(data, x, transfer_class = TRUE)
  }
 data
}
