correctRT5 <- function(CList, warpings, reference=1, ...){
  msg<-"The correctRT5 function is deprecated. You can use the 'models' argument in correctRT to warp chromatograms according to an existing model"
  .Deprecated(correctRT, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  allmats <- sapply(CList, identity, simplify = "array")
  warp.coef <- lapply(warpings,FUN=function(x){
    x$warp.coef[1,]
  })
  
  ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
    ptw(t(allmats[,, reference]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
        try=TRUE, warp.type = "global")})
  
  result <- lapply(ptwmods, function(x) t(x$warped.sample))
  for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(CList[[i]])[[1]]
  names(result) <- names(CList)
  result
}