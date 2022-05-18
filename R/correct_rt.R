#' Correct retention time
#' 
#' Corrects retention time differences using parametric time warping, as 
#' implemented in \code{\link[ptw]{ptw}}, or variable penalty dynamic 
#' time warping, as implemented in \code{\link[VPdtw]{VPdtw}}.
#' 
#' To use variable penalty dynamic time warping, the \code{VPdtw} package must
#' be manually installed since it's no longer available from CRAN: \code{
#' install.packages('VPdtw', repos='https://ethanbass.github.io/drat/')}.
#'
#' @aliases correct_rt
#' @import ptw
#' @importFrom scales rescale
#' @importFrom stats approx
#' @param chrom_list List of matrices containing concentration profiles.
#' @param lambdas Select wavelengths to use by name.
#' @param models List of models to warp by.
#' @param reference Index of the sample that is to be considered the reference
#' sample.
#' @param alg algorithm to use: parametric time warping(\code{ptw}) or variable
#' penalty dynamic time warping \code{vpdtw}.
#' @param what What to return: either the 'corrected.values' (useful for visual
#' inspection) or the warping 'models' (for further programmatic use).
#' @param init.coef Starting values for the optimization.
#' @param n.traces Number of traces to use.
#' @param n.zeros Number of zeros to add.
#' @param scale Logical. If true, scale chromatograms before warping.
#' @param trwdth width of the triangle in the WCC criterion.
#' @param plot Logical. Whether to plot alignment.
#' @param penalty Divisor for dilation calculated by \code{\link[VPdtw]{dilation}}.
#' Adjusts penalty for variable penalty dynamic time warping.
#' @param maxshift Integer. Maximum allowable shift for \code{\link[VPdtw]{VPdtw}}.
#' @param verbose Whether to be verbose.
#' @param \dots Optional arguments for the \code{\link[ptw:ptw]{ptw}} function.
#' The only argument that cannot be changed is \code{warp.type}: this is always
#' equal to \code{"global"}.
#' @return A list of `ptw` objects or a list of warped absorbance profiles,
#' depending on the value of the \code{what} argument.
#' @author Ethan Bass
#' @note Adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/correctRT.R}{correctRT}
#' function in the alsace package by Ron Wehrens.
#' @seealso \code{\link[ptw:ptw]{ptw}}, \code{\link{correct_peaks}},
#' \code{\link[VPdtw]{VPdtw}}
#' @references 
#' * Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy, P.,
#' ... & Kochhar, S. 2009. Alignment using variable penalty dynamic time warping.
#' \emph{Analytical chemistry}, \bold{81(3)}:1000-1007. \doi{10.1021/ac802041e}.
#'
#' * Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping Code
#' for Aligning Mass Spectrometry Chromatograms in R. \emph{Journal of
#' Statistical Software}, \bold{47(8)}:1-17. \doi{10.18637/jss.v047.i08}.
#' 
#' * Eilers, P.H.C. 2004. Parametric Time Warping.
#' \emph{Anal. Chem.}, \bold{76}:404-411. \doi{10.1021/ac034800e}.
#' 
#' * Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. \emph{Bioinformatics},
#' \bold{31}:3063-3065. \doi{10.1093/bioinformatics/btv299}.
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics}, \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' warping.models <- correct_rt(Sa_pr, what = "models", lambdas=c("210"))
#' warp <- correct_rt(chrom_list = Sa_pr, models = warping.models)
#' @md
#' @export correct_rt
correct_rt <- function(chrom_list, lambdas, models=NULL, reference='best',
                       alg = c("ptw", "vpdtw"), what = c("models", "corrected.values"), 
                       init.coef = c(0, 1, 0), n.traces=NULL, n.zeros=0, 
                       scale=FALSE, trwdth=200, plot=FALSE,
                       penalty=5, maxshift=50,
                       verbose = FALSE, ...){
  what <- match.arg(what, c("models", "corrected.values"))
  alg <- match.arg(alg, c("ptw", "vpdtw"))
  if (alg == "vpdtw" & !requireNamespace("VPdtw", quietly = TRUE)) {
    stop(
      "Package VPdtw must be installed to use VPdtw alignment:
      install.packages('VPdtw', repos='https://ethanbass.github.io/drat/')",
      call. = FALSE
    )
  }
  if (missing(lambdas)){
    if (is.null(models) & is.null(n.traces)){
        stop("Must specify wavelengths ('lambdas') or number of traces ('n.traces')
             to use for alignment.")
      } else lambdas <- colnames(chrom_list[[1]])
    }
    lambdas <- as.character(lambdas)
    if (scale){
      chrom_list <- lapply(chrom_list, rescale)
    }
    chrom_list_og <- chrom_list
    if (n.zeros > 0){
    chrom_list <- lapply(chrom_list,function(x){
      apply(x, 2, padzeros, nzeros=n.zeros, side='both')
    })
    }
    allmats <- sapply(chrom_list, function(x) x[,lambdas,drop=FALSE], simplify = "array")
    allmats.t <- sapply(chrom_list, function(x) t(x[,lambdas,drop=FALSE]), simplify = "array")
    if (is.null(n.traces)){
      traces <- ifelse(length(lambdas) == 1, 1, lambdas)
    } else{
      traces <- select.traces(X = allmats.t, criterion='coda')
      traces <- traces$trace.nrs[1:n.traces]
    }
    if (reference == 'best'){
      best <- bestref(allmats.t)
      reference <- as.numeric(names(sort(table(best$best.ref),decreasing=TRUE))[1])
      if (verbose) message(paste("Selected chromatogram", reference, "as best reference."))
    } else{
      reference <- reference
    }
    if (alg == "ptw"){
      if (!is.null(models)){
        ptwmods <- models
      } else{
        ptwmods <- lapply(seq_len(dim(allmats)[3]), function(ii){
          ptw(allmats.t[,, reference],
              allmats.t[,, ii], selected.traces = traces, init.coef = init.coef,
              warp.type = "global", ...)})
      }
      if (what == "corrected.values"){
        allmats <- sapply(chrom_list_og, identity, simplify = "array")
        warp.coef <- lapply(ptwmods,FUN=function(x){
          x$warp.coef[1,]
        })
        ptwmods <- lapply(seq_len(dim(allmats)[3]), function(ii){
          ptw(t(allmats[,,1]), t(allmats[, , ii]), init.coef=warp.coef[[ii]],
              try=TRUE, alg = ptwmods[[1]]$alg, warp.type = "global", ...)})
        result <- lapply(ptwmods, function(x) t(x$warped.sample))
        for (i in seq_along(result)) rownames(result[[i]]) <- rownames(chrom_list[[i]])
        names(result) <- names(chrom_list)
        result
      } else {
        ptwmods
      }
  } else{
    allmats <- sapply(chrom_list_og, function(x) x[,lambdas,drop=FALSE])
    if (length(lambdas) > 1)
      stop("VPdtw only supports warping by a single wavelength")
    if (is.null(models)){
      penalty <- VPdtw::dilation(allmats[,reference], 350) / penalty
      models <- VPdtw::VPdtw(query=allmats, reference=allmats[,reference],
                             penalty = penalty, maxshift = maxshift)
    }
    if (plot)
      VPdtw::plot.VPdtw(models)
    if (what == "corrected.values"){
      jset <- models$xVals + models$shift
      iset <- models$query
      jmax <- nrow(jset)
      short <- jmax - nrow(iset)
      res <- median(diff(as.numeric(rownames(chrom_list_og[[1]]))))
      result <- lapply(1:ncol(allmats), function(samp){
        x<-as.data.frame(apply(chrom_list_og[[samp]], 2, function(j){
          iset <- c(rep(NA, short), j)
          suppressWarnings(stats::approx(x = jset[,samp], y = iset, 1:jmax)$y)
        }))
        # fix times
        old_ts <- c(rep(NA,short), rownames(chrom_list_og[[samp]]))
        times <- suppressWarnings(stats::approx(x = jset[,samp],
                                                y = old_ts, 1:jmax)$y)
        mi<-min(which(!is.na(times)))
        if (mi>1){
          beg<-sort(seq(from = times[mi]-res, by=-res, length.out = mi-1),
                    decreasing=FALSE)
        } else beg<-NULL
        ma<-max(which(!is.na(times)))
        if (ma < length(times)){
          end <- seq(from = times[ma]+res, length.out = length(times)-ma, by=res)
        } else end <- NULL
        new.times <- c(beg,times[!is.na(times)], end)
        rownames(x) <- new.times
        x
      })
      names(result) <- names(chrom_list)
      result
    } else{
      models
    }
  }
}

#' Correct peak positions according to a ptw warping model
#' 
#' Corrects retention time differences using parametric time warping as 
#' implemented in \code{\link[ptw]{ptw}}.
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in a separate column
#' in the list of peak tables.
#' 
#' @importFrom stats predict
#' @param peak_list A nested list of peak tables: the first level is the sample,
#' and the second level is the component. Every component is described by a
#' matrix where every row is one peak, and the columns contain information on
#' retention time, full width at half maximum (FWHM), peak width, height, and
#' area.
#' @param mod_list A list of ptw models.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention time.
#' @author Ron Wehrens
#' @seealso \code{\link{correct_rt}}
#' @export correct_peaks
correct_peaks <- function(peak_list, mod_list){
  mapply(function(samp, mod){
    lapply(samp,
           function(profile){
             if (nrow(profile) > 0) {
               cbind(profile,
                     rt.cor = c(predict(mod, profile[,1], what = "time")))
             } else {
               cbind(profile, rt.cor = rep(0, 0))
             }
           })},
    peak_list, mod_list, SIMPLIFY = FALSE)
}
