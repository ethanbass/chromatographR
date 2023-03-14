#' Reshapes list of chromatograms from wide to long format
#' @name reshape_chroms
#' @param x A list of chromatographic matrices in wide format.
#' @param idx Indices of chromatograms to convert
#' @param sample_var String with name of new column containing sample IDs.
#' @param lambdas Wavelength(s) to include.
#' @return A list of chromatographic matrices in long format.
#' @author Ethan Bass
#' @export

reshape_chroms <- function(x, idx, sample_var = "sample", lambdas){
  if (missing(idx)){
    idx <- seq_along(x)
  }
  dat <- lapply(idx, function(i){
    xx <- reshape_chrom(x[[i]], lambdas)
    xx[,sample_var] <- names(x)[[i]]
    xx
  })
  do.call(rbind,dat)
}

#' Reshapes  a single chromatogram from wide to long format
#' @name reshape_chrom
#' @importFrom stats reshape
#' @param x A chromatographic matrix in wide format.
#' @param lambdas Wavelength(s) to include.
#' @return A chromatographic matrix in long format.
#' @author Ethan Bass
#' @noRd
reshape_chrom <- function(x, lambdas){
  if (ncol(x) == 1)
    stop("The provided data is already in long format!")
  x <- as.data.frame(x)
  if (!missing(lambdas)){
    x <- x[,lambdas,drop=FALSE]
  }
  data <- reshape(as.data.frame(rt=rownames(x),x), direction = "long",
                  varying = list(1:ncol(x)), v.names="absorbance",
                  times = colnames(x), timevar = "lambda",
                  idvar="rt", ids=rownames(x))
  rownames(data) <- NULL
  data$rt <- as.numeric(data$rt)
  data$lambda <- as.numeric(data$lambda)
  data[,c(3,2,1)]
}

#' Reshapes peak table from wide to long format
#' @name reshape_peaktable
#' @param x A \code{peak_table} object.
#' @param peaks A character vector specifying peaks to include.
#' @return A data.frame containing the information for the specified peaks in
#' long format.
#' @author Ethan Bass
#' @export

reshape_peaktable <- function(x, peaks){
  if (!missing(peaks)){
    x$tab <- x$tab[,which(colnames(x$tab) %in% peaks), drop=FALSE]
  }
  xx <- reshape(as.data.frame(chr=rownames(x$tab), x$tab), direction = "long",
    varying = list(1:ncol(x$tab)), v.names = x$args[["response"]],
    times = colnames(x$tab), timevar = "peak",
    idvar = "sample", ids = rownames(x$tab))
  rownames(xx) <- NULL
  xx <- xx[,c(3,1,2)]
  xx <- merge(xx, data.frame(sample=row.names(x$tab), x$sample_meta), by="sample", all.x = TRUE)
  xx
}

                  