correctRT <- function(CList, reference,
                      what = c("corrected.values", "models"),
                      init.coef = c(0, 1, 0), ...) {
  what <- match.arg(what)

  allmats <- sapply(CList, identity, simplify = "array")
 
  ptwmods <- lapply((1:dim(allmats)[3]),
                    function(ii)
                    ptw(t(allmats[,,reference]),
                        t(allmats[,,ii]),
                        init.coef = init.coef,
                        ...,
                        warp.type = "global"))
  
  if (what == "corrected.values") {
    result <- lapply(ptwmods, function(x) t(x$warped.sample))
    
    for (i in 1:length(result))
        dimnames(result[[i]]) <- dimnames(CList[[i]])
    names(result) <- names(CList)
    
    result
  } else {
    ptwmods
  }
}

