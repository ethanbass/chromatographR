#' Correct peak positions according to a PTW warping model
#' 
#' Corrects retention time differences using parametric time warping as 
#' implemented in \code{\link[ptw]{ptw}}.
#' 
#' Once an appropriate warping model has been established, corrected retention
#' times can be predicted for each peak. These are stored in a separate column
#' in the list of peak tables.
#' 
#' @importFrom stats predict
#' @param peak_list A `peak_list` object created by \code{\link{get_peaks}},
#' containing a nested list of peak tables where the first level is the sample,
#' and the second level is the spectral wavelength. Every component is described
#' by a matrix where each row corresponds to a feature, and the columns contain 
#' information on that feature (e.g., retention time, peak width (FWHM), height,
#' area, etc.)
#' @param mod_list A list of \code{ptw} models.
#' @param chrom_list List of chromatograms from which the ptw models are derived.
#' @param match_names Logical. Whether to actively match the names of the 
#' \code{peak_list} to the list of models (\code{mod_list}). Defaults to 
#' \code{TRUE}.
#' @return The input list of peak tables is returned with extra columns
#' containing the corrected retention time.
#' @author Ron Wehrens, Ethan Bass
#' @note This function is adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/correctPeaks.R}{getPeakTable}
#' function in the alsace package by Ron Wehrens.
#' @seealso \code{\link{correct_rt}}
#' @export correct_peaks

correct_peaks <- function(peak_list, mod_list, chrom_list, match_names = TRUE){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(mod_list)
  }
  ref_times <- get_times(chrom_list, idx = attr(mod_list, "reference"))
  if (match_names){
    mod_list <- mod_list[match(names(peak_list), names(mod_list))]
  } else {
    if (length(peak_list) != length(mod_list)){
      stop("The length of the provided peaklist does not match the length of the model list.")
    }
  }
  if (length(ref_times) != length(mod_list[[1]]$warp.fun)){
    stop("Dimensions of the warping models and chromatograms do not match.")
  }
  corrected_pks <- mapply(function(samp, mod){
    lapply(samp, function(profile){
      if (nrow(profile) > 0){
        profile <- cbind(profile, rt.cor = c(predict.ptw(mod, profile[, "rt"],
                                                         what = "time",
                                                         RTref = ref_times)))
        
        if (all(c("start", "end") %in% colnames(profile))){
          profile <- cbind(profile,
                           start.cor = c(predict.ptw(mod, profile[, "start"],
                                                     what = "time",
                                                     RTref = ref_times)),
                           end.cor = c(predict.ptw(mod, profile[, "end"],
                                                   what = "time",
                                                   RTref = ref_times)))
        }
      } else {
        profile <- cbind(profile, rt.cor = rep(0, 0))
      }
      profile
    }
    )}, peak_list, mod_list, SIMPLIFY = FALSE)
  corrected_pks <- transfer_metadata(corrected_pks, peak_list, 
                                     transfer_class = TRUE)
  corrected_pks
}

#' Plot PTW alignments
#' 
#' Plots \code{\link{ptw}} alignments.
#' 
#' @importFrom graphics matplot
#' @param x A \code{ptw_list} object created by \code{\link{correct_rt}}.
#' @param what What type of plot to return. Either \code{traces} or 
#' \code{heatmap}.
#' @param engine What plotting engine to use. Either \code{base}, \code{ggplot}
#' or \code{plotly}
#' @param lambdas Which lambdas to plot.
#' @param show_legend Logical. Whether to include sample legend.
#' @param ... Additional arguments (placeholder).
#' @author Ethan Bass
#' @examplesIf interactive()
#' data(Sa_pr)
#' warp <- correct_rt(chrom_list = Sa_pr, what = "models", lambdas = 210)
#' plot(warp)
#' @export

plot.ptw_list <- function(x, what = c("traces", "heatmap"),
                          engine = c("base", "ggplot", "plotly"),
                          lambdas, show_legend = TRUE, ...){
  what <- match.arg(what, c("traces", "heatmap"))
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  all.lambdas <- as.numeric(rownames(x[[1]]$warped.sample))
  ts <- as.numeric(colnames(x[[1]]$sample))
  
  if (missing(lambdas)){
    lambdas <- all.lambdas
  }
  if (any(!(lambdas %in% all.lambdas))){
    stop("Lambdas not found. Please check argument and try again")
  }
  raw_chroms <- lapply(x, function(samp) t(samp$sample))
  warped_chroms <- lapply(x, function(samp){
    xx <- t(samp$warped.sample)
    rownames(xx) <- colnames(samp$sample)
    xx
  })
  lambda.idx <- which(lambdas %in% all.lambdas)
  if (what == "traces"){
    if (engine == "base"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mfrow=c(2, 1), mar=c(2, 3, 2, 3))
      plot_chroms(raw_chroms, lambdas = lambdas, title = "Raw data",
                  show_legend = show_legend)
      plot_chroms(warped_chroms, lambdas = lambdas, title = "PTW alignment", 
                  show_legend = show_legend)
    } else {
      p1 <- plot_chroms(raw_chroms, lambdas = lambdas, title = "Raw data",
                        engine = engine, show_legend = show_legend)
      p2 <- plot_chroms(warped_chroms, lambdas = lambdas, title = "PTW alignment",
                        engine = engine, show_legend = show_legend)
      if (engine == "plotly"){
        p2 <- p2 |> plotly::style(showlegend = FALSE) |> plotly::layout(title="")
      }
      combine_plots <- switch(engine,
                              plotly = purrr::partial(plotly::subplot, nrows = 2),
                              ggplot = purrr::partial(cowplot::plot_grid, 
                                                       nrow = 2))
      p_c <- combine_plots(p1, p2)
      if (engine == "plotly"){
        p_c <- p_c |> plotly::layout(
          annotations = list(
            list(x = 0.5, y = 1, text = "Raw data", showarrow = F, 
                 xref = 'paper', yref = 'paper', font = list(size = 16)),
            list(x = 0.5, y = 0.45, text = "PTW alignments", showarrow = F, 
                 xref = 'paper', yref = 'paper', font = list(size = 16))
          )
        )
      }
      p_c
    }
  } else if (what == "heatmap"){
    original <- sapply(x, function(x) x$sample[1,])
    warped <- sapply(x, function(x) x$warped.sample[1,])
    times <- colnames(x[[1]]$reference)
    rownames(original) <- colnames(x[[1]]$reference)
    rownames(warped) <- colnames(x[[1]]$reference)
    if (engine == "base"){
      plot_chroms_heatmap_base(original, show_legend = FALSE,
                               title = "Raw data")
      plot_chroms_heatmap_base(warped, show_legend = FALSE,
                               title = "PTW alignment")
    } else{
        if (engine == "ggplot"){
        p_o <- plot_chroms_heatmap_ggplot(original, title = "Raw data")
        p_w <- plot_chroms_heatmap_ggplot(warped, title = "PTW alignment")
        p <- cowplot::plot_grid(p_o + ggplot2::theme(legend.position = "none"),
                                p_w + ggplot2::theme(legend.position = "none"), 
                                nrow = 2)
        if (show_legend){
            legend <- suppressWarnings(cowplot::get_legend(p_o + 
                        ggplot2::theme(legend.box.margin = ggplot2::margin(0, 6, 0, 0))))
            p <- cowplot::plot_grid(p, legend, ncol = 2, rel_widths  = c(1, 0.2))
        }
      } else if (engine == "plotly"){
        p_o <- plot_chroms_heatmap_plotly(original)
        p_w <- plot_chroms_heatmap_plotly(warped, show_legend = FALSE)
        p <- plotly::subplot(p_o, p_w, nrows = 2, margin = 0.1)
        p <- p |> plotly::layout(
          margin = list(t = 50),
          annotations = list(
            list(
              x = 0.5, y = 1.01, 
              text = "Raw data", 
              showarrow = FALSE,
              xref = "paper", yref = "paper",
              xanchor = "center", yanchor = "bottom",
              font = list(size = 14)
            ),
            list(
              x = 0.5, y = 0.41, 
              text = "PTW alignment", 
              showarrow = FALSE,
              xref = "paper", yref = "paper",
              xanchor = "center", yanchor = "bottom",
              font = list(size = 14)
            )
          )
        )
      }
      p
    }
  }
}

#' Predict PTW
#' @note This is the function from the ptw package, reproduced here because it
#' isn't exported from ptw.
#' @noRd
predict.ptw <- function (object, newdata, what = c("response", "time"), 
                         RTref = NULL, 
                         ...){
  what <- match.arg(what)
  switch(what, response = {
    if (missing(newdata)) return(object$warped.sample)
    if (!is.matrix(newdata)) newdata <- matrix(newdata, nrow = 1)
    if (object$warp.type == "individual" & nrow(newdata) > 
        1 & nrow(newdata) != nrow(object$warp.fun))
      stop("Incorrect number of rows in newdata")
    if (object$warp.type == "individual") {
      WF <- object$warp.fun
    } else{
      WF <- matrix(object$warp.fun, nrow(object$sample), 
                   ncol(object$warp.fun), byrow = TRUE)
    }
    if (object$mode == "backward") {
      t(sapply(seq_len(nrow(newdata)), function(i){
        approx(x = seq_len(ncol(newdata)), y = newdata[i, ], xout = WF[i, ])$y
        }))
    } else{
      t(sapply(seq_len(nrow(newdata)), function(i){
        approx(x = WF[i,], y = newdata[i, ], xout = seq_len(ncol(newdata)))$y
      }))
    }
  }, time = {
    correctedTime <- switch(object$mode, 
                            backward = -sweep(object$warp.fun, 2, 2 * 
                                                (seq_len(ncol(object$ref))), 
                                              FUN = "-"),
                            object$warp.fun)
    if (is.null(RTref)) {
      if (is.null(colnames(object$ref))) {
        RTref <- seq_len(ncol(object$ref))
      } else {
        RTref <- as.numeric(colnames(object$ref))
      }
    }
    if (missing(newdata)){
      newdata <- RTref
      newdataIndices <- seq_len(length(RTref))
    } else{
      newdataIndices <- round((newdata - min(RTref)) * 
                                (length(RTref) - 1)/diff(range(RTref)) + 1)
    }
    t(sapply(seq_len(nrow(correctedTime)), function(i){
      approx(RTref, NULL, correctedTime[i, newdataIndices])$y
      }))
  })
}
