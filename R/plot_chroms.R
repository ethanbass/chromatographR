#' Plot traces from list of chromatograms.
#' 
#' Plots the specified traces from a list of chromatograms. 
#' @importFrom graphics matplot
#' @param x A list of chromatograms in matrix format (timepoints x wavelengths).
#' @param lambdas The wavelength(s) you wish to plot the trace at.
#' @param idx A vector representing the names or numerical indices of the chromatograms to plot.
#' @param ylab Y label. Defaults to "abs".
#' @param xlab X label.
#' @param engine Plotting engine. Either \code{\link[graphics]{matplot}}, 
#' \code{\link[plotly]{plotly}}, or \code{\link[ggplot2]{ggplot2-package}}.
#' @param linewidth Line width.
#' @param show_legend Logical. Whether to display legend or not. Defaults to TRUE.
#' @param ... Additional arguments to plotting function specified by \code{engine}.
#' @examples 
#' data(Sa_pr)
#' plot_chroms(Sa_pr, idx=1:2, lambdas=c("210"))
#' @export

plot_chroms <- function(x, lambdas, idx, ylab = "Absorbance", xlab = "",
                        engine = c("base", "ggplot", "plotly"), linewidth = 1, 
                        show_legend = TRUE, ...){
  engine <- match.arg(engine, c("base","ggplot","plotly"))
  if (!class(x) %in% c("list", "chrom_list")){
    stop("Please supply list of chromatograms.")
  }
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (engine == "base"){
    first <- which.max(sapply(x[idx], function(xx) max(xx[,lambdas])))
    matplot(get_times(x[[idx[first]]]), x[[idx[first]]][,lambdas], type = 'l', lwd = linewidth,
            ylab = ylab, xlab = xlab, ...)
    if (length(idx) > 1){
      count=1
      for (i in idx[-first]){
        count <- count+1
        matplot(get_times(x, index=i), x[[i]][,lambdas], type = 'l', add = TRUE,
                col=count, lwd = linewidth, ...)
      }
    }
    if (show_legend)
      legend("topright", legend = names(x[idx]), fill=seq_along(x[idx]))
  } else {
    xx <- reshape_chroms(x, idx = idx, lambdas = lambdas)
    if (engine == "ggplot"){
      check_for_pkg("ggplot2")
      .data <- ggplot2::.data
      p <- ggplot2::ggplot(xx[which(xx$lambda %in% lambdas),],
                  ggplot2::aes(x = .data$rt, y = .data$absorbance,
                               color=.data$sample)) +
        ggplot2::geom_line(linewidth = linewidth*0.5, ...) +
        ggplot2::ylab(ylab) + ggplot2::xlab(xlab)
      if (!show_legend)
        p <- p + ggplot2::theme(legend.position = "none")
    } else if (engine == "plotly"){
      check_for_pkg("plotly")
      p <- plotly::plot_ly(data = xx[which(xx$lambda %in% lambdas),],
                           x = ~rt,y = ~absorbance, color = ~sample,
                           type='scatter', mode = 'lines', line = list(width = linewidth, ...))
      p <-  plotly::layout(p, xaxis = list(title = xlab),
                           yaxis = list(title = ylab)
      )
      if (!show_legend)
        p <- plotly::hide_legend(p)
    }
    p
  }
}

