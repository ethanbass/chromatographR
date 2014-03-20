correctPeaks <- function(peakList, modList) {
  mapply(function(samp, mod) {
    lapply(samp,
           function(profile) {
             if (nrow(profile) > 0) {
               cbind(profile,
                     rt.cor = c(predict(mod, profile[,1], what = "time")))
             } else {
               cbind(profile, rt.cor = rep(0, 0))
             }
           })},
         peakList, modList, SIMPLIFY = FALSE)
}
