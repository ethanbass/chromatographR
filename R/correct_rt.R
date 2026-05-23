#' Correct retention time
#' 
#' Aligns chromatograms using one of two algorithms, according to the value of
#' `alg`: parametric time warping (`"ptw"`), as implemented in [`ptw`][ptw::ptw]; 
#' variable penalty dynamic time warping (`"vpdtw"`), as implemented in 
#' [`VPdtw`][VPdtw::VPdtw].
#' 
#' Some arguments are specific to particular warping functions. For example 
#' the `init.coef` and `n.traces` arguments apply only to `"ptw"` warping, while 
#' `penalty` and `maxshift` apply only to  `"vpdtw"` warping.
#' 
#' @aliases correct_rt
#' @import ptw
#' @importFrom scales rescale
#' @importFrom stats approx
#' @param chrom_list List of chromatograms in matrix format.
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' use for alignment. Only one wavelength should be specified for VPdtw warping.
#' For one-dimensional chromatograms, this argument can be ignored.
#' @param models List of models to warp by. The models provided here (if any)
#' must match the algorithm selected in `alg`.
#' @param reference Index of the sample to be used as the reference. If no 
#' reference is specified, the reference will be chosen algorithmically from a
#' similarity matrix of the supplied chromatograms using the 
#' [`bestref`][ptw::bestref] function from `ptw`.
#' @param alg Alignment algorithm to use: parametric time warping (`"ptw"`), 
#' or variable penalty dynamic time warping (`"vpdtw"`).
#' @param what What to return: either the `"corrected.values"` (useful for 
#' visual inspection and downstream analysis) or the warping `"models"`
#' (for further programmatic use).
#' @param init.coef Starting values for the optimization.
#' @param n.traces Number of traces to use.
#' @param fill_zeros Logical. If `TRUE`, out-of-bounds regions produced by 
#' warping are filled with zeros. If `FALSE` (default), these regions are 
#' returned as `NA`.
#' @param n.zeros Number of zeros to add for padding chromatograms at the edges.
#' @param scale Logical. If `TRUE`, scale chromatograms before warping.
#' @param trwdth Argument to [`ptw`][ptw::ptw]. Width of the triangle in 
#' the WCC criterion. Defaults to `200`.
#' @param plot_it Logical. Whether to plot alignment. Defaults to `FALSE`.
#' @param penalty The divisor used to calculate the penalty for
#' [`VPdtw`][VPdtw::VPdtw]. The warping penalty is calculated by dividing the
#' [`dilation`][VPdtw::dilation] by this number. Thus, a higher number will
#' produce a lower penalty and be more permissive, while a lower number will 
#' produce a higher penalty and allow less warping. Defaults to `5`.
#' @param maxshift Integer. Maximum allowable shift for `VPdtw` warping.
#' Defaults to `50`.
#' @param verbose Logical. Whether to print verbose output.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed. Currently works 
#' only for `ptw` alignments.
#' @param cl Argument to `pbapply` or [`mclapply`][parallel::mclapply]. Either 
#' an integer specifying the number of clusters to use for parallel processing 
#' or a cluster object created by [`makeCluster`][parallel::makeCluster]. 
#' Defaults to `2`. On Windows systems, integer values will be ignored.
#' @param ... Optional additional arguments to `ptw`. The only argument that 
#' cannot be changed is `warp.type` which is hard-coded to `"global"` to permit
#' warping on multiple wavelengths.
#' @return A list of warping models or a list of warped absorbance profiles,
#' according to the value of the `what` argument.
#' @author Ethan Bass
#' @note Adapted from the `correctRT` function in the alsace package by Ron 
#' Wehrens (<https://github.com/rwehrens/alsace/blob/master/R/correctRT.R>).
#' @seealso [`correct_rt_group`], [`correct_peaks`], [`ptw`][ptw::ptw], 
#' [`VPdtw`][VPdtw::VPdtw]
#' @references 
#' * Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy, P.,
#' Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty dynamic time
#' warping. *Analytical chemistry*, 81(3):1000-1007. \doi{10.1021/ac802041e}.
#'
#' * Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping Code
#' for Aligning Mass Spectrometry Chromatograms in R. *Journal of
#' Statistical Software*, 47(8):1-17. \doi{10.18637/jss.v047.i08}.
#' 
#' * Eilers, P.H.C. 2004. Parametric Time Warping.
#' *Anal. Chem.*, 76:404-411. \doi{10.1021/ac034800e}.
#' 
#' * Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. *Bioinformatics*,
#' 31:3063-3065. \doi{10.1093/bioinformatics/btv299}.
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R.
#' *Metabolomics*, 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' 
#' @examplesIf interactive()
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, lambdas=210)
#' @md
#' @export correct_rt
correct_rt <- function(chrom_list, lambdas, models = NULL, reference = 'best',
                       alg = c("ptw", "vpdtw"),
                       what = c("corrected.values", "models"), 
                       init.coef = c(0, 1, 0), n.traces = NULL, 
                       fill_zeros = FALSE, n.zeros = 0, 
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
    chrom_list <- lapply(chrom_list, scales::rescale)
  }
  chrom_list_og <- chrom_list
  if (n.zeros > 0){
    chrom_list <- lapply(chrom_list, function(x){
      x_p <- pad_zeros(x, n_zeros = n.zeros, side = "both")
      transfer_metadata(x_p, x)
    })
  }
  allmats.t <- sapply(chrom_list, function(x){
    t(x[, lambdas, drop = FALSE])}, simplify = "array")
  if (is.null(n.traces)){
    traces <- ifelse(length(lambdas) == 1, 1, list(lambdas))[[1]]
  } else {
    traces <- select.traces(X = allmats.t, criterion = 'coda')
    traces <- traces$trace.nrs[seq_len(n.traces)]
  }
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
    laplee <- choose_apply_fnc(show_progress, cl = cl)
    if (is.null(models)){
      if (verbose) message("Fitting PTW warping models.")
      models <- laplee(seq_len(dim(allmats.t)[3]), function(ii){
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
      if (verbose) message("Applying PTW warping models to chromatograms.")
      result <- laplee(seq_along(models), function(ii){
        coef <- models[[ii]]$warp.coef[1, ]
        x <- warp_sample_ptw(chrom_list_og[[ii]], coef, 
                             fill_zeros = fill_zeros,
                             mode = models[[ii]]$mode)
        rownames(x) <- rownames(chrom_list_og[[ii]])
        colnames(x) <- colnames(chrom_list_og[[ii]])
        transfer_metadata(x, chrom_list_og[[ii]])
      })
      result <- structure(result, warped = TRUE, warping_args = args)
      names(result) <- names(chrom_list)
      result
    } else {
      return(models)
    }
  } else if (alg == "vpdtw"){
    if (length(lambdas) > 1)
      stop("VPdtw only supports warping by a single wavelength")
    if (is.null(models)){
      if (verbose) message("Fitting VPdtw warping models.")
      allmats <- sapply(chrom_list_og, function(x) x[, lambdas, drop = FALSE])
      penalty <- VPdtw::dilation(allmats[,reference], 350) / penalty
      models <- VPdtw::VPdtw(query = allmats, reference = allmats[, reference],
                             penalty = penalty, maxshift = maxshift)
      attr(models, "parameters") <- list(reference = reference,
                                         maxshift = maxshift)
    }
    if (plot_it){
      VPdtw::plot.VPdtw(models)
    }
    if (what == "corrected.values"){
      if (verbose) message("Applying VPdtw warping models to chromatograms.")
      jset <- models$xVals + models$shift
      iset <- models$query
      jmax <- nrow(jset)
      short <- jmax - nrow(iset)
      res <- get_time_resolution(chrom_list_og, idx = reference)
      result <- lapply(seq_along(chrom_list), function(samp){
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
      result <- lapply(result, function(xx){
        if (fill_zeros && any(is.na(xx))){
          xx[which(is.na(xx))] <- 0
        }
        xx <- structure(xx, warped = TRUE, warping_args = args)
        xx
      })
      result
    } else {
      return(models)
    }
  }
}

#' Apply VPdtw
#' @noRd
apply_vpdtw <- function(chrom_list, models, reference, args = NULL, 
                        fill_zeros = FALSE) {
  jset <- models$xVals + models$shift
  iset <- models$query
  jmax <- nrow(jset)
  short <- jmax - nrow(iset)
  res <- get_time_resolution(chrom_list, idx = reference)
  
  result <- lapply(seq_along(chrom_list), function(ii) {
    warp_sample_vpdtw(chrom_list[[ii]], jset[, ii], short, 
                      fill_zeros = fill_zeros)
  })
  
  old_ts <- c(rep(NA, short), get_times(chrom_list, idx = reference))
  times <- suppressWarnings(stats::approx(x = jset[, reference],
                                          y = old_ts, seq_len(jmax))$y)
  idx_start <- which.min(times)
  beg <- if (idx_start > 1) {
    sort(seq(from = times[idx_start] - res, by = -res,
             length.out = idx_start - 1), decreasing = FALSE)
  } else NULL
  idx_end <- which.max(times)
  end <- if (idx_end < length(times)) {
    seq(from = times[idx_end] + res,
        length.out = length(times) - idx_end, by = res)
  } else NULL
  new.times <- c(beg, times[!is.na(times)], end)
  
  # assign times and metadata
  result <- mapply(function(x, idx) {
    rownames(x) <- new.times
    transfer_metadata(x, chrom_list[[idx]])
  }, result, seq_along(result), SIMPLIFY = FALSE)
  names(result) <- names(chrom_list)
  
  # replace NAs with 0s and attach metadata
  lapply(result, function(xx) {
    if (fill_zeros && any(is.na(xx))) xx[which(is.na(xx))] <- 0
    structure(xx, warped = TRUE, warping_args = args)
  })
}

#' Warp sample VPdtw
#' @noRd
warp_sample_vpdtw <- function(mat, jset_col, short, fill_zeros = FALSE) {
  padded <- rbind(matrix(NA, nrow = short, ncol = ncol(mat)), mat)
  tp <- seq_len(nrow(padded))
  valid_range <- range(jset_col, na.rm = TRUE)
  jset_col[is.na(jset_col) & seq_along(jset_col) > valid_range[2]] <- Inf
  jset_col[is.na(jset_col)] <- -Inf
  lo <- findInterval(tp, jset_col)
  hi <- lo + 1
  lo_c <- pmax(1, lo)
  hi_c <- pmin(length(jset_col), hi)
  frac <- (tp - jset_col[lo_c]) / (jset_col[hi_c] - jset_col[lo_c])
  frac[is.na(frac)] <- 0
  frac_mat <- matrix(frac, nrow = length(tp), ncol = ncol(padded))
  x <- padded[lo_c, ] + frac_mat * (padded[hi_c, ] - padded[lo_c, ])
  x[which(tp < valid_range[1] | 
            tp > valid_range[2]), ] <- ifelse(fill_zeros, 0, NA)
  x
}

#' Warp sample PTW
#' @noRd
warp_sample_ptw <- function(mat, coef, fill_zeros = FALSE,
                            mode = c("forward", "backward")){
  mode <- match.arg(mode)
  tp <- seq_len(nrow(mat))
  w <- ptw::warp.time(tp, coef)
  if (mode == "forward"){
    lo <- findInterval(tp, w)
    hi <- lo + 1
    oob <- lo < 1 | hi > length(w)
    lo_c <- pmax(1, lo)
    hi_c <- pmin(length(w), hi)
    frac <- (tp - w[lo_c]) / (w[hi_c] - w[lo_c])
  } else {
    stop("Backward mode is not yet supported.")
  }
  frac[oob] <- 0
  frac_mat <- matrix(frac, nrow = length(tp), ncol = ncol(mat))
  x <- mat[lo_c, ] + frac_mat * (mat[hi_c, ] - mat[lo_c, ])
  x[which(tp < min(w) | tp > max(w)), ] <- ifelse(fill_zeros, 0, NA)
  x
}


#' Correct retention time with group-based warping
#'
#' Aligns chromatograms using parametric time warping (`"ptw"`), as
#' implemented in [`ptw`][ptw::ptw], with warping functions estimated
#' from within-group averages. This is useful when samples fall into batches
#' with shared retention time shifts, or when individual samples contain peaks
#' absent in other groups that would otherwise confound alignment. In addition
#' to potentially being more accurate, this should also be faster than computing
#' warping functions individually on each sample.
#'
#' @aliases correct_rt_group
#' @import ptw
#' @param chrom_list A list of chromatograms in matrix format.
#' @param lambdas A character or numeric vector specifying the wavelengths to
#' use for alignment.
#' @param groups A vector of group assignments for each chromatogram, or a
#' single string naming a metadata attribute from which to extract the group
#' assignments (e.g. `"batch"`). If a vector is provided, it must be the same 
#' length as `chrom_list` and may be named (matching `names(chrom_list)`) or
#' positional. All samples must have a group assignment; `NA` values will
#' trigger an error.
#' @param reference Index or name of the group average to use as the alignment
#' reference. Defaults to `"best"`, which selects the reference
#' automatically using [`bestref`][ptw::bestref] applied to the group averages.
#' @param reference_group Name of the group to use as the reference group. If
#' supplied, overrides `reference`. Defaults to `NULL`.
#' @inheritParams correct_rt
#' @return A list of warped chromatogram matrices in the same order as
#' `chrom_list`, with each sample warped using the warp coefficients
#' estimated from its group average.
#' @author Ethan Bass
#' @seealso [`correct_rt`], [`ptw`][ptw::ptw]
#' @examplesIf interactive()
#' \dontrun{
#' data(Sa_pr)
#' groups <- c("a", "a", "b", "b")
#' warp <- correct_rt_group(chrom_list = Sa_pr, lambdas = 210, groups = groups)
#' }
#' @md
#' @export

correct_rt_group <- function(chrom_list, lambdas, groups, reference = 'best',
                             reference_group = NULL,
                             init.coef = c(0, 1, 0), n.traces = NULL,
                             fill_zeros = FALSE, n.zeros = 0, scale = FALSE, 
                             trwdth = 200, plot_it = FALSE, penalty = 5, 
                             maxshift = 50, verbose = getOption("verbose"),
                             show_progress = NULL, cl = 2,
                             ...) {
  
  if (length(groups) == 1 && is.character(groups)) {
    groups_vec <- sapply(chrom_list, function(x) attr(x, groups))
    if (any(is.na(groups_vec))) {
      stop(sprintf("Metadata attribute '%s' not found on all chromatograms.", 
                   groups))
    }
  } else {
    groups_vec <- groups
    groups_vec <- as.factor(groups_vec)
    if (is.null(names(groups_vec))) {
      if (length(groups_vec) != length(chrom_list)){
        stop("Length of `groups` must match the number of samples in `chrom_list`.")
      }
      names(groups_vec) <- names(chrom_list)
    } else {
      groups_vec <- groups_vec[names(chrom_list)] 
    }
  }
  if (any(is.na(groups_vec))) {
    na_idx <- which(is.na(groups_vec))
    na_samples <- names(groups_vec)[na_idx]
    stop(sprintf("All samples must have a group assignment. NA%s found for:\n%s",
                 ifelse(length(na_idx)>1, "s", ""),
                 paste(
                   sprintf("         %d (%s)", na_idx, sQuote(na_samples)),
                   collapse = ", \n")))
  }
  
  group_levels <- levels(groups_vec)
  
  if (verbose) message("Computing within-group averages.")
  group_avg_list <- lapply(setNames(group_levels, group_levels), function(g) {
    members <- names(groups_vec)[which(groups_vec == g)]
    mats <- lapply(chrom_list[members], as.matrix)
    avg <- Reduce("+", mats) / length(mats)
    dimnames(avg) <- dimnames(mats[[1]])
    avg
  })
  
  if (is.null(reference_group)) {
    if (reference == 'best') {
      reference_group <- NULL  # passed through, correct_rt will resolve
    } else {
      ref_name <- if (is.numeric(reference)) names(chrom_list)[reference] else reference
      reference_group <- groups_vec[ref_name]
      reference <- which(names(group_avg_list) == reference_group)
    }
  } else {
    reference <- which(names(group_avg_list) == reference_group)
  }
  
  if (verbose) message("Fitting PTW warp models on group averages.")
  group_models <- correct_rt(
    chrom_list  = group_avg_list,
    lambdas     = lambdas,
    reference   = reference,
    alg         = "ptw",
    what        = "models",
    init.coef   = init.coef,
    n.traces    = n.traces,
    n.zeros     = n.zeros,
    scale       = scale,
    trwdth      = trwdth,
    plot_it     = plot_it,
    penalty     = penalty,
    maxshift    = maxshift,
    verbose     = verbose,
    show_progress = show_progress,
    cl          = cl,
    ...
  )
  laplee <- choose_apply_fnc(show_progress, cl = cl)
  if (verbose) message("Applying group warp models to individual chromatograms.")
  result <- laplee(group_levels, function(g) {
    members <- names(groups_vec)[which(groups_vec == g)]
    member_list <- chrom_list[members]
    correct_rt(
      chrom_list  = member_list,
      lambdas     = lambdas,
      models      = structure(rep(group_models[g], 
                                  length(member_list)), class="ptw_list"), 
      alg         = "ptw",
      what        = "corrected.values",
      init.coef   = init.coef,
      n.zeros     = n.zeros,
      fill_zeros = fill_zeros,
      scale       = scale,
      trwdth      = trwdth,
      verbose     = verbose,
      show_progress = FALSE,
      cl          = 1,
      ...
    )
  })
  
  result_flat <- do.call(c, result)
  result_flat[names(chrom_list)]
}

#' Pad zeros
#' Helper for [`correct_rt`].
#' @noRd
pad_zeros <- function(x, n_zeros, side = "both") {
  if (is.vector(x)) x <- as.matrix(x)
  zeros <- matrix(0, nrow = n_zeros, ncol = ncol(x))
  switch(side,
         both = rbind(zeros, x, zeros),
         left = rbind(zeros, x),
         right = rbind(x, zeros)
  )
}

#' Strip zeros
#' @noRd
strip_zeros <- function(x, n_zeros, side = "both") {
  n <- nrow(x)
  switch(side,
         both = x[(n_zeros + 1):(n - n_zeros), ],
         left = x[(n_zeros + 1):n, ],
         right = x[1:(n - n_zeros), ]
  )
}
