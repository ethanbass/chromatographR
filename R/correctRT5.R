correctRT5 <- function(CList, warpings, reference=1, ...){
  
  allmats <- sapply(CList, identity, simplify = "array")
  warpings2 <- lapply(warpings,FUN=function(x){
    x$warp.coef[1,]
  })
  
  ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
    ptw(t(allmats[,, reference]), t(allmats[, , ii]), init.coef=warpings2[[ii]],
        try=TRUE, warp.type = "global")})
  
  result <- lapply(ptwmods, function(x) t(x$warped.sample))
  for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(CList[[i]])[[1]]
  names(result) <- names(CList)
  result
}