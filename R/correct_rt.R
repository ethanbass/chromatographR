#' Correct retention time
#' 
#' Corrects retention time differences using parametric time warping as 
#' implemented in \code{\link[ptw:ptw]{ptw}}.
#' 
#' @aliases correct_rt
#' @import ptw
#' @importFrom scales rescale
#' @param chrom_list List of matrices containing concentration profiles.
#' @param models List of models to warp by.
#' @param lambdas Select wavelengths to use by name.
#' @param reference Index of the sample that is to be considered the reference
#' sample.
#' @param what What to return: either the 'corrected.values' (useful for visual
#' inspection) or the warping 'models' (for further programmatic use).
#' @param init.coef Starting values for the optimization.
#' @param n.traces Number of traces to use.
#' @param n.zeros Number of zeros to add.
#' @param scale Logical. If true, scale chromatograms before warping.
#' @param trwdth width of the triangle in the WCC criterion.
#' @param \dots Optional arguments for the \code{\link[ptw:ptw]{ptw}} function.
#' The only argument that cannot be changed is \code{warp.type}: this is always
#' equal to \code{"global"}.
#' @return A list of `ptw` objects or a list of warped absorbance profiles.
#' @author Ethan Bass
#' @note Adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/correctRT.R}{correctRT}
#' function in the alsace package by Ron Wehrens.
#' @seealso \code{\link[ptw:ptw]{ptw}}, \code{\link{correct_peaks}}
#' @references Eilers, P.H.C. 2004. Parametric Time Warping.
#' \emph{Anal. Chem.} \bold{76}:404-411. \doi{10.1021/ac034800e}.
#' 
#' Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. \emph{Bioinformatics}
#' \bold{31}:3063-3065. \doi{10.1093/bioinformatics/btv299}.
#' 
#' Wehrens, R., Carvalho, E., Fraser, P.D. 2015.
#' Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11:1}:143-154. \doi{10.1007/s11306-014-0683-5}
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' warping.models <- correct_rt(Sa_pr, what = "models", lambdas=c("210"))
#' warp <- correct_rt(chrom_list = Sa_pr, models = warping.models)
#' @export correct_rt

correct_rt <- function(chrom_list, lambdas, models=NULL, reference='best',
                       alg = c("ptw", "sptw"), what = c("models", "corrected.values"), 
                       init.coef = c(0, 1, 0), n.traces=NULL, n.zeros=0, 
                       scale=TRUE, trwdth=200, ndx = 40, ...){
  what <- match.arg(what, c("models", "corrected.values"))
  alg <- match.arg(alg, c("ptw", "sptw"))
  if (is.null(models)){
    if (is.null(lambdas) & is.null(n.traces)){
      stop("Must specify wavelengths ('lambdas') or number of traces ('n.traces')
           to use for alignment.")
    }
    if (is.null(lambdas)){
      lambdas <- colnames(chrom_list[[1]])
    }
    lambdas <- as.character(lambdas)
    if (n.zeros > 0){
    chrom_list <- lapply(chrom_list,function(x){
      apply(x, 2, padzeros, nzeros=n.zeros, side='both')
    })
    }
    if (scale){
      chrom_list <- lapply(chrom_list, rescale)
    }
    allmats <- sapply(chrom_list, function(x) x[,lambdas,drop=FALSE], simplify = "array")
    allmats.t <- sapply(chrom_list, function(x) t(x[,lambdas,drop=F]), simplify = "array")
    if (is.null(n.traces)){
      traces <- ifelse(length(lambdas) == 1, 1, lambdas)
    } else{
      traces <- select.traces(X=allmats.t, criterion='coda')
      traces <- traces$trace.nrs[1:n.traces]
    }
    if (reference == 'best'){
      best <- bestref(allmats.t)
      reference <- as.numeric(names(sort(table(best$best.ref),decreasing=TRUE))[1])
    } else{
      reference <- reference
    }
    ptwmods <- lapply(seq_len(dim(allmats)[3]), function(ii){
      ptw(allmats.t[,, reference],
          allmats.t[,, ii], selected.traces = traces, init.coef=init.coef,
          alg = alg, warp.type = "global", ndx=ndx, ...)})
  } else {
    allmats <- sapply(chrom_list, identity, simplify = "array")
    warp.coef <- lapply(models,FUN=function(x){
      x$warp.coef[1,]
    })
    ptwmods <- lapply(seq_len(dim(allmats)[3]), function(ii){
      ptw(t(allmats[,,1]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
          try=TRUE, alg = models[[1]]$alg, warp.type = "global", ...)})
  }
  if (what == "corrected.values" || !is.null(models)) {
    result <- lapply(ptwmods, function(x) t(x$warped.sample))
    for (i in seq_along(result)) rownames(result[[i]]) <- rownames(chrom_list[[i]])
    names(result) <- names(chrom_list)
    result
  } else {
    ptwmods
  }
}

