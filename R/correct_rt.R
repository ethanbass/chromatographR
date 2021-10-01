correct_rt <- function(chrom_list, models=NULL, lambdas, reference='best', what = c("models", "corrected.values"), 
                       init.coef = c(0, 1, 0), n.traces=NULL, n.zeros=0, scale=T,
                       trwdth=200, plot_it=T, ...){
  what <- match.arg(what, c("models", "corrected.values"))
  if (is.null(models)){
    chrom_list <- lapply(chrom_list,function(x){
      apply(x, 2, padzeros, nzeros=n.zeros, side='both')
    })
    if (scale){
      chrom_list<-lapply(chrom_list, rescale)
    }
    allmats <- sapply(chrom_list, function(x)x[,lambdas], simplify = "array")
    allmats.t <- sapply(chrom_list, function(x) t(x[,lambdas]), simplify = "array")
    if (is.null(n.traces)){
      traces=lambdas
    } else{
      traces <- select.traces(X=allmats.t, criterion='coda')
      traces <- traces$trace.nrs[1:n.traces]
    }
    if (reference=='best'){
      best<-bestref(allmats.t)
      reference <- as.numeric(names(sort(table(best$best.ref),decreasing=TRUE))[1])
    } else{
      reference <- reference
    }
    ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
      ptw(allmats.t[,, reference],
          allmats.t[, , ii], selected.traces = traces, init.coef=init.coef, ..., warp.type = "global")})
  } else {
    allmats <- sapply(chrom_list, identity, simplify = "array")
    warp.coef <- lapply(models,FUN=function(x){
      x$warp.coef[1,]
    })
    ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
      ptw(t(allmats[,,1]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
          try=TRUE, warp.type = "global")})
  }
  if (what == "corrected.values" || !is.null(models)) {
    result <- lapply(ptwmods, function(x) t(x$warped.sample))
    for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(chrom_list[[i]])[[1]]
    names(result) <- names(chrom_list)
    result
  } else {
    ptwmods
  }
}

# if (plot_it==T){
#   sapply(ptwmods$
# }

correctRT <- function(chrom_list, models=NULL, reference='best', what = c("corrected.values", "models"), 
                      init.coef = c(0, 1, 0), n.traces=NULL, n.zeros=0, lambdas, scale=T,
                      trwdth=200, plot_it=T, ...){
  msg<-"The function `correctRT` is deprecated. Please use `correct_rt` instead"
  .Deprecated(correct_rt, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  what <- match.arg(what, c("models", "corrected.values"))
  if (is.null(models)){
  chrom_list<-lapply(chrom_list,function(x){
    apply(x,2,padzeros, nzeros=n.zeros, side='both')
  })
  if (scale){
    chrom_list<-lapply(chrom_list,rescale)
  }
  allmats <- sapply(chrom_list, function(x)x[,lambdas], simplify = "array")
  allmats.t <- sapply(chrom_list, function(x) t(x[,lambdas]), simplify = "array")
  if (is.null(n.traces)){
    traces=lambdas
  } else{
    traces <- select.traces(X=allmats.t,criterion='coda')
    traces <- traces$trace.nrs[1:n.traces]
  }
  if (reference=='best'){
    best<-bestref(allmats.t)
    reference <- as.numeric(names(sort(table(best$best.ref),decreasing=TRUE))[1])
  } else{
    reference <- reference
  }
  ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
    ptw(allmats.t[,, reference],
        allmats.t[, , ii], selected.traces = traces, init.coef=init.coef, ..., warp.type = "global")})
  } else {
    allmats <- sapply(chrom_list, identity, simplify = "array")
    warp.coef <- lapply(models,FUN=function(x){
      x$warp.coef[1,]
    })
    ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
      ptw(t(allmats[,,1]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
          try=TRUE, warp.type = "global")})
  }
  if (what == "corrected.values" || !is.null(models)) {
    result <- lapply(ptwmods, function(x) t(x$warped.sample))
    for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(chrom_list[[i]])[[1]]
    names(result) <- names(chrom_list)
    result
  } else {
    ptwmods
  }
}

correctRT5 <- function(chrom_list, warpings, reference=1, ...){
  msg<-"The `correctRT5` function is deprecated. You can use the 'models' argument in `correct_rt` to warp chromatograms according to an existing model"
  .Deprecated(correctRT, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  allmats <- sapply(chrom_list, identity, simplify = "array")
  warp.coef <- lapply(warpings,FUN=function(x){
    x$warp.coef[1,]
  })
  
  ptwmods <- lapply((1:dim(allmats)[3]), function(ii){
    ptw(t(allmats[,, reference]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
        try=TRUE, warp.type = "global")})
  
  result <- lapply(ptwmods, function(x) t(x$warped.sample))
  for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(chrom_list[[i]])[[1]]
  names(result) <- names(chrom_list)
  result
}