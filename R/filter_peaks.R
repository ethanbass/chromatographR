filter_peaks <- function(peak_list, min_height, min_area, min_sd, max_sd) {
  if (missing(min_height) & missing(min_area) &
      missing(max_sd) & missing(min_sd)) {
    warning("Nothing to filter...")
    return(peak_list)
  }
  
  if (missing(max_sd)) ## find max value in peak_list and add 1 to be sure...
    max_sd <- 1 + max(unlist(sapply(peak_list,
                                    function(samp)
                                      sapply(samp,
                                             function(comp) comp[,"sd"]))))
  if (missing(min_sd)) min_sd <- 0
  if (missing(min_height)) min_height <- 0
  if (missing(min_area)) min_area <- 0
  
  lapply(peak_list,
         function(smpl)
           lapply(smpl,
                  function(comp)
                    comp[comp[,"sd"] < max_sd &
                           comp[,"sd"] > min_sd &
                           comp[,"height"] > min_height &
                           comp[,"area"] > min_area,,drop = FALSE]))
}
