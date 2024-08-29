#' Plot traces from list of chromatograms.
#' 
#' Plots the specified traces from a list of chromatograms.
#' 
#' @importFrom graphics matplot axis box
#' @param x A list of chromatograms in matrix format (timepoints x wavelengths).
#' @param lambdas The wavelength(s) you wish to plot the trace at.
#' @param idx A vector representing the names or numerical indices of the 
#' chromatograms to plot.
#' @param xlim Range of x axis.
#' @param ylim Range of y axis.
#' @param ylab Y label. Defaults to "Absorbance".
#' @param xlab X label.
#' @param engine Plotting engine. Either \code{base} (\code{\link[graphics]{matplot}}), 
#' \code{\link[plotly]{plotly}}, or \code{\link[ggplot2]{ggplot2-package}}.
#' @param linewidth Line width.
#' @param show_legend Logical. Whether to display legend or not. Defaults to TRUE.
#' @param legend_position Position of legend.
#' @param ... Additional arguments to plotting function specified by \code{engine}.
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots the traces of the specified chromatograms \code{idx} at the specified
#' wavelengths \code{lambdas}. Plots can be produced using base graphics, ggplot2,
#' or plotly, according to the value of \code{engine}.
#' @author Ethan Bass
#' @examples 
#' data(Sa_pr)
#' plot_chroms(Sa_pr, idx = c(1:2), lambdas = c(210))
#' @concept Visualization
#' @export

plot_chroms <- function(x, lambdas, idx, xlim, ylim, xlab = "", ylab = "Absorbance",
                        engine = c("base", "ggplot", "plotly"), linewidth = 1, 
                        show_legend = TRUE, legend_position = "topright", ...){
  engine <- match.arg(engine, c("base","ggplot","plotly"))
  if (!class(x) %in% c("list", "chrom_list")){
    stop("Please supply list of chromatograms.")
  }
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (ncol(x[[1]]) == 1){
    lambdas.idx <- 1
  } else {
    lambdas.idx <- sapply(lambdas, function(lambda){
      get_lambda_idx(lambda, lambdas = get_lambdas(x),
                     allow_max = FALSE)
    }) 
  }
  zoom_x <- TRUE
  zoom_y <- TRUE
  if (missing(xlim)){
    zoom_x <- FALSE
    xlim <- c(head(get_times(x), 1), tail(get_times(x), 1))
  }
  if (missing(ylim)){
    zoom_y <- FALSE
    ylim <- get_y_bounds(x, idx, lambdas.idx)
  }
  if (engine == "base"){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    title(xlab = xlab, ylab = ylab)
    axis(1)
    axis(2)
    box()
    for (i in seq_along(idx)){
        matplot(get_times(x, idx = idx[i]), x[[idx[i]]][,lambdas.idx], type = 'l',
                add = TRUE, col = i, lwd = linewidth, ...)
    }
    if (show_legend)
      legend(x = legend_position, legend = names(x)[idx], fill = seq_along(x[idx]))
  } else {
    xx <- reshape_chroms(x, idx = idx, lambdas = lambdas.idx)
    if (engine == "ggplot"){
      check_for_pkg("ggplot2")
      .data <- ggplot2::.data
      p <- ggplot2::ggplot(xx[which(xx$lambda %in% lambdas),],
                           ggplot2::aes(x = .data$rt, y = .data$absorbance,
                                        color = .data$sample)) +
        ggplot2::geom_line(linewidth = linewidth*0.5, ...) +
        ggplot2::ylab(ylab) + ggplot2::xlab(xlab)
      if (!show_legend){
        p <- p + ggplot2::theme(legend.position = "none")
      }
      if (legend_position != "topright"){
        p <- p + ggplot2::theme(legend.position = legend_position)
      }
      if (zoom_x){
        p <- p + ggplot2::xlim(xlim)
      }
      if (zoom_y){
        p <- p + ggplot2::ylim(ylim)
      }
    } else if (engine == "plotly"){
      check_for_pkg("plotly")
      p <- plotly::plot_ly(data = xx[which(xx$lambda %in% lambdas),],
                           x = ~rt, y = ~absorbance, color = ~sample,
                           type='scatter', mode = 'lines',
                           line = list(width = linewidth, ...))
      p <-  plotly::layout(p, xaxis = list(title = xlab),
                           yaxis = list(title = ylab)
      )
      if (!show_legend){
        p <- plotly::hide_legend(p)
      }
      if (legend_position != "topright"){
        p <- plotly::layout(p, legend = position_plotly_legend(legend_position))
      }
      if (zoom_x || zoom_y){
        p <- plotly::layout(p, xaxis = list(range = xlim), yaxis = list(range = ylim))
      }
    }
    suppressWarnings(p)
  }
}

#' Position plotly legend
#' @author Ethan Bass
#' @noRd
position_plotly_legend <- function(pos){
  switch(pos, 
         bottomright = list("x" = 100, "y" = 0),
         topright = list("x" = 100, "y"= 1),
         right = list("x" = 100, "y" = .5),
         bottom = list("y" = -100, orientation = "h"),
         top = list("y" = 100, orientation = "h"),
         bottomleft = list("x" = 0.1, "y" = 0),
         topleft = list("x" = 0.1, "y"= 1),
         left = list("x" = 0.1, "y" = .5),
  )
}

