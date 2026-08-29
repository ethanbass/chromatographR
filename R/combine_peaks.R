#' Combine duplicate peaks
#'
#' Identifies groups of duplicate peaks in a peak table and retains a single
#' representative peak from each group. Peaks are considered duplicates when
#' their retention times differ by less than `tol` and their spectral
#' correlation exceeds `min.cor`. This is useful for collapsing peaks that were
#' integrated at more than one wavelength or chromatographic component.
#' 
#' @name combine_peaks
#' @param peak_table Peak table from [`get_peaktable`].
#' @param tol Tolerance for matching retention times (maximum retention time
#' difference). Defaults to `0.01`.
#' @param min.cor Minimum spectral correlation to confirm a match. Defaults to 
#' `0.9`.
#' @param choose Method used to select a representative peak from each group of
#' duplicate peaks. `"max"` retains the peak with the greatest total intensity,
#' `"least_sparse"` retains the peak detected in the greatest number of samples,
#' and `"lambda"` retains peaks matching the preferred wavelength(s).
#' @param lambda Character vector of preferred wavelength(s) to retain when
#' `choose = "lambda"`. The function keeps the duplicate peak whose integration
#' wavelength matches one of the supplied values. An error is thrown if none of
#' the supplied values match any peak in the peak table.
#' @param verbose Logical. Whether to print status to the console.
#' @return A peak table derived from the original, but with columns corresponding
#' to duplicate peaks combined according to the specified criteria.
#' @author Ethan Bass
#' @seealso [`get_peaks`]
#' @examples
#' data(pk_tab)
#' data(Sa_warp)
#' pk_tab <- attach_ref_spectra(pk_tab)
#' combine_peaks(pk_tab, tol = 0.02, min.cor = 0.9)
#' @family utility functions
#' @export combine_peaks

combine_peaks <- function(peak_table, tol = .01, min.cor = 0.9,
                          choose = c("max", "least_sparse", "lambda"),
                          lambda = NULL,
                          verbose = getOption("verbose")){
  check_peaktable(peak_table)
  if (!(is.data.frame(peak_table$ref_spectra) | is.matrix(peak_table$ref_spectra))){
    stop("No reference spectra found. Use attach_ref_spectra function first.")
  }
  choose <- match.arg(choose, c("max", "least_sparse", "lambda"))
  if (choose == "lambda") {
    if (is.null(lambda))
      stop("Please select preferred lambda.")
    available <- unique(as.character(peak_table$pk_meta["lambda", ]))
    if (!any(lambda %in% available))
      stop(sprintf(
        "None of the specified lambda value(s) (%s) match any peak in the peak table.",
        paste(sQuote(lambda), collapse = ", ")
      ))
  }
  RTs <- as.numeric(peak_table$pk_meta['rt',])
  compare_rts <- function(rt1, rt2, tol){
    abs(rt1 - rt2) < tol
  }
  mat <- outer(RTs, RTs, compare_rts, tol = tol)
  # find columns with a retention time match
  rt_group_count <- colSums(mat)
  candidate_rt_duplicates <- which(rt_group_count > 1)
  spectral_cor_mat <- suppressWarnings(cor(peak_table$ref_spectra))
  if (verbose) warnings()
  # iterate over columns with retention time match
  # find columns to remove (sub)
  remove_idx <- lapply(candidate_rt_duplicates, function(idx){
    
    i <- which(mat[, idx] == 1)
    
    j <- which(spectral_cor_mat[i, idx] > min.cor)
    k <- i[j]
    
    if (length(k) <= 1) return(NULL)
    
    peak_names <- colnames(peak_table$tab)[k]
    tab_k <- peak_table$tab[, k, drop = FALSE]
    
    keep <- switch(
      choose,
      
      max = peak_names[
        which.max(colSums(tab_k, na.rm = TRUE))
      ],
      
      least_sparse = {
        coverage <- colSums(!is.na(tab_k))
        names(coverage)[which.max(coverage)]
      },
      
      lambda = {
        lambdas <- peak_table$pk_meta["lambda", k]
        lambda_idx <- which(lambdas %in% lambda)
        if (length(lambda_idx) > 0) {
          peak_names[lambda_idx[1]]
        } else {
          peak_names[which.max(colSums(tab_k, na.rm = TRUE))]
        }
      }
    )
    
    setdiff(peak_names, keep)
  })
  
  remove_idx <- unique(unlist(remove_idx))
  remove_idx <- match(remove_idx, colnames(peak_table$tab))
  
  if (length(remove_idx) > 0) {
    peak_table$tab <- peak_table$tab[, -remove_idx, drop = FALSE]
    peak_table$pk_meta <- peak_table$pk_meta[, -remove_idx, drop = FALSE]
  }
  if (verbose){
    message(paste("Removed "), length(remove_idx), " peaks from peak table.")
  }
  peak_table
}

#' Merge split peaks
#' 
#' Utility function to combine split peaks into a single column of the peak 
#' table.
#' 
#' Merges the specified peaks in peak table, by selecting the largest value from 
#' each column if `method` is `"max"`. If `method` is `"sum"`, peaks will be 
#' merged by summing their values.
#'
#' @name merge_peaks
#' @param peak_table A `peak_table` object from [`get_peaktable`].
#' @param peaks A vector specifying the names or indices of peaks to be merged.
#' @param method Method to merge peaks. Either `max` to select the largest
#' peak from each sample or `sum` to add the peaks together.
#' @return A peak table similar to the input peak table, but where the specified
#' columns are combined. 
#' @author Ethan Bass
#' @examples
#' data(pk_tab)
#' pk_tab <- merge_peaks(peak_table = pk_tab, peaks = c("V10","V11"))
#' @family utility functions
#' @export

merge_peaks <- function(peak_table, peaks, method = c("max", "sum")){
  check_peaktable(peak_table)
  method <- match.arg(method, c("max", "sum"))
  if (is.character(peaks)){
    pks.idx <- which(colnames(peak_table$tab) %in% peaks)
  } else {
    pks.idx <- peaks
  }
  sel <- which.max(colMeans(peak_table$tab[, pks.idx, drop=FALSE], na.rm = TRUE))
  sel.idx <- which(colnames(peak_table$tab) == peaks[sel])
  if (method == "max"){
    peak_table$tab[[sel.idx]] <- do.call(pmax, peak_table$tab[, pks.idx, 
                                                              drop = FALSE])
  } else if (method == "sum"){
    peak_table$tab[[sel.idx]] <- apply(peak_table$tab[,pks.idx], 1, sum)
  }
  peak_table$tab <- peak_table$tab[, -pks.idx[-sel], drop = FALSE]
  peak_table$pk_meta <- peak_table$pk_meta[, -pks.idx[-sel], drop = FALSE]
  if (inherits(peak_table$ref_spectra, "matrix")){
    peak_table$ref_spectra <- peak_table$ref_spectra[, -pks.idx[-sel]]
  }
  return(peak_table)
}
