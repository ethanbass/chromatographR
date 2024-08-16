#' Reshape chromatograms
#' Reshapes a list of chromatograms from wide to long format.
#' @name reshape_chroms
#' @param x A list of chromatographic matrices in wide format.
#' @param idx Indices of chromatograms to convert
#' @param sample_var String with name of new column containing sample IDs.
#' @param lambdas Vector specifying wavelength(s) to include.
#' @param rts Vector specifying retention times to include.
#' @return A list of chromatographic matrices in long format.
#' @author Ethan Bass
#' @export

reshape_chroms <- function(x, idx, sample_var = "sample",
                            lambdas = NULL, rts = NULL){
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (missing(lambdas)){
    lambdas <- colnames(x[[1]])
  }
  dat <- lapply(idx, function(i){
    xx <- reshape_chrom(x = x[[i]], lambdas = lambdas, rts = rts)
    xx[, sample_var] <- names(x)[[i]]
    xx
  })
  do.call(rbind, dat)
}

#' Reshapes a single chromatogram from wide to long format
#' @name reshape_chrom
#' @importFrom stats reshape
#' @param x A chromatographic matrix in wide format.
#' @param lambdas Vector specifying wavelength(s) to include.
#' @param rts Vector specifying retention times to include.
#' @return A chromatographic matrix in long format.
#' @author Ethan Bass
#' @noRd
reshape_chrom <- function(x, lambdas = NULL, rts = NULL){
  if (ncol(x) == 1)
    stop("The provided data is already in long format!")
  times <- get_times(x = x)
  xx <- as.data.frame(x)
  if (!is.null(lambdas)){
    xx <- xx[, lambdas, drop = FALSE]
  }
  if (!is.null(rts)){
    if (is.character(rts)){
      rts <- as.numeric(rts)
    }
    rts.idx <- sapply(rts, function(rt){
      get_retention_idx(RT = rt, times = times)})
    xx <- xx[rts.idx, , drop = FALSE]
  }
  data <- reshape(as.data.frame(rt = rownames(xx), xx), direction = "long",
                  varying = list(seq_len(ncol(xx))), v.names = "absorbance",
                  times = colnames(xx), timevar = "lambda",
                  idvar = "rt", ids = rownames(xx))
  rownames(data) <- NULL
  data$rt <- as.numeric(data$rt)
  data$lambda <- as.numeric(data$lambda)
  data <- data[,c(3,2,1)]
  transfer_metadata(data, x, transfer_class = TRUE)
}

#' Reshapes peak table from wide to long format
#' @name reshape_peaktable
#' @param x A \code{peak_table} object.
#' @param peaks A character vector specifying the peaks to include. If the
#' character vector is named, the names of the vector elements will be used in
#' place of the original peak names.
#' @param metadata A character vector specifying the metadata fields to include.
#' @param fixed_levels Logical. Whether to fix factor levels of features in the
#' order provided. Defaults to \code{TRUE}.
#' @return A data.frame containing the information for the specified peaks in
#' long format.
#' @author Ethan Bass
#' @export

reshape_peaktable <- function(x, peaks, metadata, fixed_levels = TRUE){
  if (!missing(peaks)){
    if (is.numeric(peaks)){
      peaks <- colnames(x$tab)[peaks]
    }
    df <- x$tab[, match(peaks, colnames(x$tab)), drop = FALSE]
    if (!is.null(names(peaks))){
      colnames(df) <- names(peaks)
      peaks <- colnames(df)
    }
  } else {
    df <- x$tab
  }
  if (!missing(metadata)){
    meta_idx <- which(colnames(x$sample_meta) %in% metadata)
    x$sample_meta <- x$sample_meta[, meta_idx, drop = FALSE]
  }
  xx <- reshape(as.data.frame(chr = rownames(df), df), direction = "long",
    varying = list(seq_len(ncol(df))), v.names = x$args[["response"]],
    times = colnames(df), timevar = "peak",
    idvar = "sample", ids = rownames(df))
  rownames(xx) <- NULL
  xx <- xx[,c(3,1,2)]
  if (!is.null(dim(x$sample_meta))){
    xx <- merge(xx, data.frame(sample = row.names(df), x$sample_meta),
                by = "sample", all.x = TRUE)
  }
  if (fixed_levels){
    xx$peak <- factor(xx$peak, levels = peaks)
  }
  xx
}

# reshape <- function(x,...){
#   UseMethod("reshape")
# }

