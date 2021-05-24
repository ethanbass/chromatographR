correctRT <- function(CList, reference='best', what = c("corrected.values", "models"), 
                      init.coef = c(0, 1, 0), n.traces=NULL, n.zeros, selected.traces, scale=T,
                      trwdth=200, ...) 
{
  what <- match.arg(what)
  CList<-lapply(CList,function(x){
    apply(x,2,padzeros,nzeros=n.zeros,side='both')
  })
  if (scale){
    CList<-lapply(CList,scale)
  }
  allmats <- sapply(CList, function(x)x[,selected.traces], simplify = "array")
  allmats.t <- sapply(CList, function(x) t(x[,selected.traces]), simplify = "array")
  if (is.null(n.traces)){
    traces=selected.traces
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
  if (what == "corrected.values") {
    result <- lapply(ptwmods, function(x) t(x$warped.sample))
    for (i in 1:length(result)) dimnames(result[[i]])[[1]] <- dimnames(CList[[i]])[[1]]
    names(result) <- names(CList)
    result
  } else {
    ptwmods
  }
}

