#' Calculate mean peak purity
#' 
#' Estimates peak purity by assessing the dissimilarity of the spectra
#' comprising the peak, using the method described in Stahl 2003.
#' 
#' @param x A chromatogram in matrix format
#' @param pos A vector containing the center, lower and upper bounds of a peak
#' as numeric indices.
#' @param weight Weight provided to \code{\link{get_agilent_threshold}}.
#' @param cutoff Proportion of maximum absorbance to use as cutoff.
#' Argument to \code{\link{trim_peak}}. Defaults to \code{.05}.
#' @param noise_variance Variance of noise. Argument to 
#' \code{\link{get_agilent_threshold}}.
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance. Defaults to \code{.01}.
#' @param lambdas Wavelengths to include in calculations.
#' @param try Logical. Whether to estimate the purity or not. Defaults to TRUE.
#' @return Returns the mean purity of the peak specified by \code{pos}, defined
#' as the proportion of timepoints with purity values below 1.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' /href{https://www.agilent.com/cs/library/applications/5988-8647EN.pdf}
#' @author Ethan Bass
#' @keywords internal
#' @export

get_purity <- function(x, pos, weight = 1, cutoff = 0.05, 
                       noise_variance = NULL, 
                       noise_threshold = 0.01,
                       lambdas, try = TRUE){
  if (try){
    try({
      if (missing(lambdas)){
        lambdas <- seq_len(ncol(x))
      }
      if (is.character(lambdas)){
        lambdas <- which(as.numeric(colnames(x)) %in% lambdas) 
      }
      p <- get_purity_values(x, pos, weight = weight,
                             noise_variance = noise_variance,
                             lambdas = lambdas)
      mean(p[trim_peak(x, pos, cutoff = cutoff)] < 1, na.rm = TRUE)
      }, NA
    )
  } else NA
}

#' Calculate variance of noise regions
#' @param x A chromatogram in matrix format
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance. Defaults to \code{.01}.
#' @param lambdas Wavelengths to include.
#' @importFrom stats var
#' @return Returns the average variance of the signal over the retention times
#' defined as noise according to \code{\link{find_noise}}.
#' @author Ethan Bass
#' @noRd

get_noise_variance <- function(x, noise_threshold = .01, lambdas){
  if (missing(lambdas)){
    lambdas <- seq_len(ncol(x))
  }
  noise_idx <- find_noise(x = x, noise_threshold = noise_threshold, 
                      lambdas = lambdas)
  mean(apply(x[noise_idx,], 1, var), na.rm = TRUE)
  # if (plot_it){
  #   matplot(x,type='l')
  #   
  # }
}

#' Define noise spectra based on specified threshold
#' @param x A chromatogram in matrix format 
#' @param noise_threshold Threshold to define noise. Highest proportion of maximum absorbance.
#' @param lambdas Wavelengths to include.
#' @return Returns indices of retention times where the signal falls below the
#' specified noise threshold.
#' @author Ethan Bass
#' @noRd

find_noise <- function(x, noise_threshold = 0.01, lambdas){
  max_abs <- apply(x[,lambdas, drop = FALSE], 1, max)
  which(max_abs < max(max_abs, na.rm = TRUE) * noise_threshold)
}

#' Calculate purity thresholds
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param weight Scaling parameter affecting stringency of threshold. Defaults
#' to \code{1}.
#' @param noise_variance Variance of noise.
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance. Defaults to \code{.005}.
#' @param lambdas Wavelengths to include
#' @return Returns a vector of purity thresholds at each retention time index
#' within the peak specified by \code{pos}.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' /href{https://www.agilent.com/cs/library/applications/5988-8647EN.pdf}
#' @author Ethan Bass
#' @keywords internal
#' @export

get_agilent_threshold <- function(x, pos, weight = 1, noise_variance = NULL,
                                  noise_threshold = .005,
                                  lambdas){
  if (missing(lambdas)){
    lambdas <- seq_len(ncol(x))
  }
  if (is.null(noise_variance)){
    var_noise <- get_noise_variance(x, noise_threshold = noise_threshold,
                                    lambdas = lambdas)
  } else {
    var_noise <- noise_variance
  }
  idx <- seq(as.numeric(pos[2]), as.numeric(pos[3]))
  xx <- sapply(idx, function(i){
    (max(0, 1 - weight *
           (var_noise / var(x[i, ], na.rm = TRUE) +
              var_noise / var(x[as.numeric(pos[1]), ], na.rm = TRUE))))^2
  })
  xx
}

#' Calculate spectral similarity
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information.
#' @return Returns a vector of spectral similarities of the reference spectrum
#' to the spectrum at each timepoint within the peak specified by \code{pos}.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' /href{https://www.agilent.com/cs/library/applications/5988-8647EN.pdf}
#' @author Ethan Bass
#'@noRd

get_spectral_similarity <- function(x, pos){
  idx <- seq(as.numeric(pos[2]), as.numeric(pos[3]))
  suppressWarnings(cor(x[as.numeric(pos[1]),], t(x[idx,])))
}

#' Calculate peak purity values
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param weight weight provided to \code{\link{get_agilent_threshold}}.
#' Defaults to \code{1}.
#' @param noise_variance Variance of noise. Argument to 
#' \code{\link{get_agilent_threshold}}.
#' @param noise_threshold Threshold to define noise. Highest proportion of 
#' maximum absorbance. Defaults to \code{.005}.
#' @param lambdas Wavelengths to include in calculations.
#' @return Returns a vector of peak purity values at each timepoint within the
#' peak specified by \code{pos}.
#' @references
#' Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array Technology.”
#' Agilent Technologies, April 1, 2003, 16.
#' /href{https://www.agilent.com/cs/library/applications/5988-8647EN.pdf}
#' @author Ethan Bass
#' @noRd

get_purity_values <- function(x, pos, weight = 1, noise_variance = NULL, 
                              noise_threshold = 0.005,
                              lambdas){
  if (missing(lambdas)){
    lambdas <- seq_len(ncol(x))
  }
  ((1 - get_spectral_similarity(x, pos)))/
    (1 - get_agilent_threshold(x, pos, weight = weight,
                               noise_variance = noise_variance,
                               noise_threshold = noise_threshold,
                               lambdas = lambdas))
}

#' Trim peak
#' @param x A chromatogram in matrix format
#' @param pos A vector containing peak information
#' @param cutoff Proportion of maximum absorbance to use as cutoff. Defaults to
#' \code{.05}.
#' @return Returns indices within the peak specified by \code{pos} with a higher
#' signal intensity than the specified cutoff.
#' @author Ethan Bass
#' @keywords internal
#' @export

trim_peak <- function(x, pos, cutoff = 0.05){
  idx <- pos[2]:pos[3]
  which(x[idx] > cutoff*x[pos[1]])
}