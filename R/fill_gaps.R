#' Fill gaps in peak table
#' @param peak_table A \code{peak_table} object.
#' @param chrom_list A list of chromatograms in matrix format
#' @param ref How to select reference spectrum: either maximum correlation
#' (\code{max.cor}) or maximum intensity (\code{max.int}).
#' @param similarity_threshold Minimum spectral similarity.
#' @param rt_tolerance Maximum retention time difference.
#' @param spectral_weight Relative weight of spectral similarity and retention
#' time similarity.
#' @param only_zeros Only consider zeros in peak_table for filling.
#' @param peaks_only Only consider peaks found by \code{\link{find_peaks}}.
#' @param plot_it Whether to plot difference between filled and unfilled peak tables.
#' @author Ethan Bass
#' @export fill_gaps
fill_gaps <- function(peak_table, chrom_list, ref = c("max.cor","max.int"),
                      similarity_threshold=0.95, rt_tolerance=0.15,
                      spectral_weight=0.5, only_zeros = FALSE,
                      peaks_only=TRUE, plot_it = FALSE){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  }
  ref <- match.arg(ref, c("max.cor","max.int"))
  if (is.na(peak_table$ref_spectra)){
    peak_table <- attach_ref_spectra(peak_table, ref = ref)
  } else{
    if (peak_table$args["reference_spectra"] != ref)
      warning(paste("reference spectra in peak table are", peak_table$args["reference_spectra"],
                     "but selection is", ref))
  }
  for (j in colnames(peak_table$tab)){
    pktab <- fill_gap(j, peak_table, similarity_threshold = similarity_threshold,
             rt_tolerance = rt_tolerance, spectral_weight = spectral_weight,
             only_zeros = only_zeros, peaks_only = peaks_only, plot_it = plot_it)
  }
  pktab
}

#' Fill gap in one column of peaktable
#' @param peak_table A \code{peak_table} object.
#' @param chrom_list A list of chromatograms in matrix format
#' @param reference chromatogram to use as reference. Defaults to "max".
#' @param similarity_threshold Minimum spectral similarity.
#' @param rt_tolerance Maximum retention time difference.
#' @param only_zeros Only consider zeros in peak_table for filling.
#' @param peaks_only Only consider peaks found by \code{\link{find_peaks}}.
#' @param plot_it Whether to plot difference between filled and unfilled peak tables.
#' @author Ethan Bass
#' @noRd
fill_gap <- function(peak, peak_table, chrom_list, similarity_threshold=0.95,
                     rt_tolerance = 0.5, spectral_weight = 0.5, only_zeros = FALSE,
                     peaks_only = TRUE, plot_it = FALSE){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  }
  ref <- peak_table$ref_spectra[,peak]
  ts <- as.numeric(rownames(chrom_list[[1]]))
  lambdas <- as.numeric(colnames(chrom_list[[1]]))
  OG_areas <- peak_table$tab[,peak]
  # lambdas <- rownames(ref)
  ref.s <- rescale(ref)
  rt <- round(as.numeric(peak_table$pk_meta['rt', peak]), 2)
  pk.idx <- which(elementwise.all.equal(rt, ts))
  if (only_zeros){
    chrs <- which(peak_table$tab[,peak] == 0)
  } else {
    chrs <- (seq_along(chrom_list))
  }
  rad <- rt_tolerance/median(diff(ts))
    for (chr in chrs){
      start <- max(0, (pk.idx - rad))
      end <- min((pk.idx + rad), length(ts))
      spec <- t(chrom_list[[chr]][c(start:end),])
      spec.s <- rescale(spec)
      spec_cor <- as.numeric(suppressWarnings(cor(ref.s, spec.s, method='pearson')))
      lambda <- lambdas[which.max(ref.s)]
      lambda <- which(lambdas == lambda)
      if (peaks_only){
        idx <- find_peaks(spec[lambda,])$pos
      } else{
        idx <- seq_len(2*rad+1)
      }
      idx <- idx[spec_cor[idx] >= similarity_threshold]
      rt_diff <- abs(c(-rad:rad)[idx])
      score <- spectral_weight*spec_cor[idx] + (1-spectral_weight)*(1-(rt_diff/(2*rad)))
      idx_select <- idx[which.max(score)]
      if (length(idx_select) == 0){
        break
      } else{
        # replace values in peak_table
        peak_table$tab[chr,peak] <- spec[peak_table$pk_meta["lambda",peak],idx_select]
      }
    }
  if (plot_it){
    matplot(seq_dim(peak_table), cbind(OG_areas, peak_table$tab[,peak]), 
            pch=c(1,20), xlab='old',ylab='new')
    # legend()
    # plot(peak_table$tab[,peak] ~ OG_areas, pch=20)
  }
  peak_table
}


#' @noRd
seq_dim <- function(x, d = c("rows","cols")){
  d<-match.arg(d, c("rows","cols"))
  d <- switch(d, "rows" = 1, "cols" = 2)
  seq_len(dim(x)[d])
}

