#' Correct retention time
#' 
#' Aligns chromatograms using one of two algorithms, according to the value of 
#' \code{alg}: either parametric time warping, as implemented in 
#' \code{\link[ptw]{ptw}}, or variable penalty dynamic time warping, as 
#' implemented in \code{\link[VPdtw]{VPdtw}}. The \code{init.coef} and 
#' \code{n.traces} arguments apply only to \code{ptw} warping, while 
#' \code{penalty} and \code{maxshift} apply only to \code{vpdtw}
#' warping.
#' 
#' @aliases correct_rt
#' @import ptw
#' @importFrom scales rescale
#' @importFrom stats approx
#' @param chrom_list List of chromatograms in matrix format.
#' @param lambdas Select wavelengths to use by name.
#' @param models List of models to warp by. The models provided here (if any)
#' must match the algorithm selected in \code{alg}.
#' @param reference Index of the sample that is to be considered the reference
#' sample.
#' @param alg algorithm to use: parametric time warping (\code{ptw}) or variable
#' penalty dynamic time warping (\code{vpdtw}).
#' @param what What to return: either the 'corrected.values' (useful for visual
#' inspection) or the warping 'models' (for further programmatic use).
#' @param init.coef Starting values for the optimization.
#' @param n.traces Number of traces to use.
#' @param n.zeros Number of zeros to add.
#' @param scale Logical. If true, scale chromatograms before warping.
#' @param trwdth width of the triangle in the WCC criterion.
#' @param plot_it Logical. Whether to plot alignment.
#' @param penalty The divisor used to calculate the penalty for
#' \code{\link[VPdtw]{VPdtw}}. The warping penalty is calculated by dividing the
#' \code{\link[VPdtw]{dilation}} by this number. Thus, a higher number will
#' produce a lower penalty and be more permissive, while a lower number will 
#' produce a higher penalty and allow less warping. Defaults to 5.
#' @param maxshift Integer. Maximum allowable shift for \code{\link[VPdtw]{VPdtw}}.
#' Defaults to 50.
#' @param verbose Whether to print verbose output.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' \code{TRUE} if \code{\link[pbapply]{pbapply}} is installed. Currently works 
#' only for \code{ptw} alignments.
#' @param cl Argument to \code{\link[pbapply]{pblapply}} or \code{\link[parallel]{mclapply}}.
#' Either an integer specifying the number of clusters to use for parallel
#' processing or a cluster object created by \code{\link[parallel]{makeCluster}}.
#' Defaults to 2. On Windows integer values will be ignored.
#' @param \dots Optional arguments for the \code{\link[ptw:ptw]{ptw}} function.
#' The only argument that cannot be changed is \code{warp.type}: this is always
#' equal to \code{"global"}.
#' @return A list of warping models or a list of warped absorbance profiles,
#' depending on the value of the \code{what} argument.
#' @author Ethan Bass
#' @note Adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/correctRT.R}{correctRT}
#' function in the alsace package by Ron Wehrens.
#' @seealso \code{\link[ptw:ptw]{ptw}}, \code{\link{correct_peaks}},
#' \code{\link[VPdtw]{VPdtw}}
#' @references 
#' * Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy, P.,
#' Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty dynamic time
#' warping. \emph{Analytical chemistry}, \bold{81(3)}:1000-1007. \doi{10.1021/ac802041e}.
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
#' Metabolomics}, \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, lambdas=210)
#' @md
#' @export correct_rt
correct_rt <- function(chrom_list, lambdas, models = NULL, reference = 'best',
                       alg = c("ptw", "vpdtw"), what = c("corrected.values", "models"), 
                       init.coef = c(0, 1, 0), n.traces = NULL, n.zeros = 0, 
                       scale = FALSE, trwdth = 200, plot_it = FALSE,
                       penalty = 5, maxshift = 50,
                       verbose = getOption("verbose"), show_progress = NULL, 
                       cl = 2, ...){
  what <- match.arg(what, c("corrected.values", "models"))
  alg <- match.arg(tolower(alg), c("ptw", "vpdtw"))
  
  if (!is.null(models)){
    model_class <- switch(alg, ptw = "ptw_list", "vpdtw" = "VPdtw")
    if (!inherits(models, model_class)){
      stop("The supplied models do not match the selected algorithm. Please check arguments
           and try again.")
    }
  }
  if (missing(lambdas)){
    if (ncol(chrom_list[[1]]) != 1 & is.null(models) & is.null(n.traces)){
      stop("Must specify wavelengths ('lambdas') or number of traces ('n.traces')
           to use for alignment.")
    } else lambdas <- colnames(chrom_list[[1]])
  }
  lambdas <- as.character(lambdas)
  if (!all(lambdas %in% colnames(chrom_list[[1]]))){
    stop("Lambdas not found!")
  }
  if (scale){
    chrom_list <- lapply(chrom_list, rescale)
  }
  chrom_list_og <- chrom_list
  if (n.zeros > 0){
    chrom_list <- lapply(chrom_list, function(x){
      apply(x, 2, padzeros, nzeros = n.zeros, side = 'both')
    })
  }
  allmats <- sapply(chrom_list, function(x){
    x[, lambdas, drop = FALSE]}, simplify = "array")
  allmats.t <- sapply(chrom_list, function(x){
    t(x[, lambdas, drop = FALSE])}, simplify = "array")
  if (is.null(n.traces)){
    traces <- ifelse(length(lambdas) == 1, 1, list(lambdas))[[1]]
  } else {
    traces <- select.traces(X = allmats.t, criterion = 'coda')
    traces <- traces$trace.nrs[seq_len(n.traces)]
  }
  # choose reference chromatogram
  if (reference == 'best'){
    best <- ptw::bestref(allmats.t)
    reference <- as.numeric(names(sort(table(best$best.ref), 
                                       decreasing = TRUE))[1])
    if (verbose) message(paste("Selected chromatogram", reference, 
                               "as best reference."))
  } else {
    reference <- reference
  }
  args <- substitute(list(lambdas = lambdas, models = models, 
                          reference = reference, alg = alg, 
                          init.coef = init.coef, n.traces = n.traces,
                          n.zeros = n.zeros, scale = scale, trwdth = trwdth,
                          penalty = penalty, maxshift = maxshift))
  if (alg == "ptw"){
    if (is.null(models)){
      laplee <- choose_apply_fnc(show_progress, cl = cl)
      models <- laplee(seq_len(dim(allmats)[3]), function(ii){
        pw <- ptw::ptw(allmats.t[,, reference],
            allmats.t[,, ii], selected.traces = traces, init.coef = init.coef,
            warp.type = "global", ...)
        if (nrow(pw$reference) == 1){
          colnames(pw$reference) <- get_times(chrom_list)
          colnames(pw$sample) <- get_times(chrom_list)
          colnames(pw$warped.sample) <- get_times(chrom_list)
          rownames(pw$reference) <- lambdas
          rownames(pw$sample) <- lambdas
          rownames(pw$warped.sample) <- lambdas
        }
        pw
        })
      names(models) <- names(chrom_list)
      models <- structure(models, chrom_list = deparse(substitute(chrom_list)), 
                          reference = reference, init.coef = init.coef,
                          n.traces = n.traces, n.zeros=n.zeros, scale = scale,
                          trwdth = trwdth, penalty = penalty, 
                          maxshift = maxshift, class = "ptw_list")
      if (plot_it){
        plot(models)
      }
    }
    if (what == "corrected.values"){
      allmats <- sapply(chrom_list_og, identity, simplify = "array")
      warp.coef <- lapply(models, FUN = function(x){
        x$warp.coef[1,]
      })
      models <- lapply(seq_len(dim(allmats)[3]), function(ii){
        ptw::ptw(t(allmats[,, 1]), t(allmats[,, ii]), init.coef = warp.coef[[ii]],
            try = TRUE, alg = models[[1]]$alg, warp.type = "global", ...)})
      result <- lapply(seq_along(models), function(i){
        x <- t(models[[i]]$warped.sample)
        rownames(x) <- rownames(chrom_list[[i]])
        colnames(x) <- colnames(chrom_list[[i]])
        x <- transfer_metadata(x, chrom_list_og[[i]])
      })
      result <- structure(result, warped = TRUE, warping_args = args)
      names(result) <- names(chrom_list)
      result
    } else {models}
  } else if (alg == "vpdtw"){
    if (length(lambdas) > 1)
      stop("VPdtw only supports warping by a single wavelength")
    allmats <- sapply(chrom_list_og, function(x) x[, lambdas, drop = FALSE])
    if (is.null(models)){
      penalty <- VPdtw::dilation(allmats[,reference], 350) / penalty
      models <- VPdtw::VPdtw(query = allmats, reference = allmats[, reference],
                             penalty = penalty, maxshift = maxshift)
    }
    if (plot_it){
      VPdtw::plot.VPdtw(models)
    }
    if (what == "corrected.values"){
      jset <- models$xVals + models$shift
      iset <- models$query
      jmax <- nrow(jset)
      short <- jmax - nrow(iset)
      res <- get_time_resolution(chrom_list_og, idx = reference)
      result <- lapply(seq_len(ncol(allmats)), function(samp){
        # warp retention times
        x <- apply(chrom_list_og[[samp]], 2, function(j){
          iset <- c(rep(NA, short), j)
          suppressWarnings(stats::approx(x = jset[,samp], y = iset, 
                                         seq_len(jmax))$y)
        })
      })
        # fix times
        old_ts <- c(rep(NA, short), get_times(chrom_list_og, idx = reference))
        times <- suppressWarnings(stats::approx(x = jset[, reference],
                                                y = old_ts, seq_len(jmax))$y)
        idx_start <- which.min(times)
        if (idx_start > 1){
          beg <- sort(seq(from = times[idx_start] - res, by = -res,
                          length.out = idx_start - 1), decreasing = FALSE)
        } else beg <- NULL
        idx_end <- which.max(times)
        if (idx_end < length(times)){
          end <- seq(from = times[idx_end] + res,
                     length.out = length(times) - idx_end, by = res)
        } else end <- NULL
        new.times <- c(beg, times[!is.na(times)], end)
        result <- mapply(function(x,idx){
          rownames(x) <- new.times
          x <- transfer_metadata(x, chrom_list_og[[idx]])
          x
        }, result, seq_along(result), SIMPLIFY = FALSE)
      names(result) <- names(chrom_list)
      # replace NAs with 0s and add additional metadata
      result <- lapply(result, function(xx){
          if(any(is.na(xx))){
            xx[which(is.na(xx))] <- 0
          }
          xx <- structure(xx, warped = TRUE, warping_args = args)
          xx
        })
      result
    } else {
      models
    }
  }
}
