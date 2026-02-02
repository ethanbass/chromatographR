#' Reshape chromatograms
#' 
#' Reshapes a list of chromatograms from wide to long format.
#' 
#' @name reshape_chroms
#' @param x A list of chromatographic matrices in wide format.
#' @param idx Indices of chromatograms to convert.
#' @param time_resolution Time resolution for plot. This argument can be used
#' to thin the time axis while reshaping. By default the time resoution is not
#' altered.
#' @param sample_var String with name of new column containing sample IDs.
#' @param lambdas Vector specifying wavelength(s) to include.
#' @param rts Vector specifying retention times to include.
#' @param transfer_metadata Logical. Whether to transfer metadata attributes or
#' not. Defaults to \code{FALSE}.
#' @return A list of chromatographic matrices in long format.
#' @author Ethan Bass
#' @family utility functions
#' @export

reshape_chroms <- function(x, idx, time_resolution = NULL,  
                           sample_var = "sample", lambdas = NULL, rts = NULL, 
                           transfer_metadata = FALSE){
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (missing(lambdas)){
    lambdas <- colnames(x[[1]])
  }
  dat <- lapply(idx, function(i){
    xx <- reshape_chrom(x = x[[i]], time_resolution = time_resolution,
                        lambdas = lambdas, rts = rts,
                        transfer_metadata)
    xx[, sample_var] <- names(x)[[i]]
    xx
  })
  do.call(rbind, dat)
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
