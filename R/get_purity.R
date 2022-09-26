#' Define noise spectra based on specified threshold
#' @param x A chromatogram in matrix format 
#' @param noise_threshold Threshold to define noise. Highest proportion of maximum absorbance.
#' @param lambdas Wavelengths to include.
#' @author Ethan Bass
find_noise <- function(x, noise_threshold=0.01, lambdas){
  lambdas <- which(as.numeric(colnames(x)) %in% lambdas)
  max_abs <- apply(x[,lambdas], 1, max)
  which(max_abs < max(max_abs)*noise_threshold)
}

#alternatively, could define baseline as areas where no peak is detected,
# since we've already done peak detection.

#' Calculate variance of noise regions
#' @param x A chromatogram in matrix format
#' @param noise_threshold Threshold to define noise. Highest proportion of maximum absorbance.
#' @param lambdas Wavelengths to include
#' @importFrom stats var
#' @author Ethan Bass
get_noise_variance <- function(x, noise_threshold=.005, lambdas=c(210:400)){
  noise_idx <- find_noise(x = x, noise_threshold = noise_threshold, 
                      lambdas = lambdas)
  mean(apply(x[noise_idx,], 1, var))
  # if (plot_it){
  #   matplot(x,type='l')
  #   
  # }
}

# 
# matplot(chrom[,lambda],type='l')
# i<-find_noise(chrom, thresh=.005, lambdas=c(220:400))
# points(i, chrom[i,lambda],col="red",pch=20)

# noise_variance <- get_noise_variance(chrom)

#' calculate purity thresholds
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param weight scaling parameter affecting stringency of threshold
#' @param noise_threshold Threshold to define noise. Highest proportion of maximum absorbance.
#' @param lambdas Wavelengths to include
#' @author Ethan Bass
get_agilent_threshold <- function(x, pos, weight=1, noise_threshold=.005,
                                  lambdas=c(210:400)){
  var_noise <- get_noise_variance(x, noise_threshold = noise_threshold, lambdas=lambdas)
  idx <- seq(as.numeric(pos[2]), as.numeric(pos[3]))
  xx <- sapply(idx, function(i){
    (max(0, 1 - weight *
           (var_noise / var(x[i, ]) +
              var_noise / var(x[as.numeric(pos[1]), ]))))^2
  })
  xx
}

#' Calculate spectral similarity
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @author Ethan Bass
get_spectral_similarity <- function(x, pos){
  idx <- seq(as.numeric(pos[2]), as.numeric(pos[3]))
  suppressWarnings(cor(x[as.numeric(pos[1]),], t(x[idx,])))
}

# calc_purity_value <- function(x, pos){
#   get_spectral_similarity(x,pos) > get_agilent_threshold(x, pos)
# }

#' Calculate peak purity values
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param weight weight provided to \code{\link{get_agilent_threshold}}
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' /href{https://www.agilent.com/cs/library/applications/5988-8647EN.pdf}
#' @author Ethan Bass

get_purity_values <- function(x, pos, weight =1){
  ((1 - get_spectral_similarity(x, pos)))/
    (1 - get_agilent_threshold(x, pos, weight = weight))
}

#' Trim peak
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param cutoff Proportion of maximum absorbance to use as cutoff.
#' @author Ethan Bass
trim_peak <- function(x, pos, cutoff = 0.05){
  idx <- pos[2]:pos[3]
  which(x[idx] > cutoff*x[pos[1]])
}

#' Calculate mean peak purity
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param weight weight provided to \code{\link{get_agilent_threshold}}
#' @param cutoff Proportion of maximum absorbance to use as cutoff.
#' Argument to \code{\link{trim_peak}}.
#' @author Ethan Bass
get_mean_purity <- function(x, pos, weight=1, cutoff=0.05){
  p <- get_purity_values(x, pos, weight = weight)
  mean(p[trim_peak(x, pos, cutoff=cutoff)] < 1)
}
