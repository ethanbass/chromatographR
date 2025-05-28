#' Plot traces from list of chromatograms.
#' 
#' Plots the specified traces from a list of chromatograms.
#' 
#' @importFrom graphics matplot axis box
#' @param x A list of chromatograms in matrix format (timepoints x wavelengths).
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' plot.
#' @param idx A vector representing the names or numerical indices of the 
#' chromatograms to plot.
#' @param xlim Range of x axis.
#' @param ylim Range of y axis.
#' @param ylab Y label. Defaults to "Absorbance".
#' @param xlab X label.
#' @param engine Plotting engine. Either \code{base} (\code{\link[graphics]{matplot}}), 
#' \code{\link[plotly]{plotly}}, or \link[ggplot2:ggplot2-package]{ggplot}.
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
#' data(Sa_warp)
#' plot_chroms(Sa_warp, lambdas = 210)
#' @family visualization functions
#' @export

plot_chroms <- function(x, lambdas, idx, xlim = NULL, ylim = NULL, 
                        xlab = "", ylab = "Absorbance",
                        engine = c("base", "ggplot", "plotly"), linewidth = 1, 
                        show_legend = TRUE, legend_position = "topright", 
                        title = "", ...){
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  if (!class(x) %in% c("list", "chrom_list")){
    stop("Please supply list of chromatograms.")
  }
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (is.character(idx)){
    idx <- match(idx, names(x))
  }
  if (ncol(x[[1]]) == 1){
    lambdas.idx <- 1
  } else {
    lambdas.idx <- sapply(lambdas, function(lambda){
      get_lambda_idx(lambda, lambdas = get_lambdas(x),
                     allow_max = FALSE)
    }) 
  }
  zoom_x <- zoom_y <- TRUE
  if (is.null(xlim)){
    zoom_x <- FALSE
    xlim <- c(head(get_times(x), 1), tail(get_times(x), 1))
  }
  if (is.null(ylim)){
    zoom_y <- FALSE
    ylim <- get_y_bounds(x, idx, lambdas.idx)
  }
  if (engine == "base"){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    title(main = title, xlab = xlab, ylab = ylab)
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
      p <- ggplot2::ggplot(xx, ggplot2::aes(x = .data$rt, y = .data$absorbance,
                                        color = .data$sample)) +
        ggplot2::geom_line(linewidth = linewidth*0.5, ...) +
        ggplot2::ylab(ylab) + 
        ggplot2::xlab(xlab) + 
        ggplot2::ggtitle(title) +
        ggplot2::theme_light()
      if (!show_legend){
        p <- p + ggplot2::theme(legend.position = "none")
      } else if (legend_position != "topright"){
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
      p <- plotly::plot_ly(data = xx, x = ~rt, y = ~absorbance, color = ~sample,
                           type = 'scatter', mode = 'lines',
                           line = list(width = linewidth, ...))
      p <-  plotly::layout(p, xaxis = list(title = xlab),
                           yaxis = list(title = ylab),
                           title = title
      )
      if (!show_legend){
        p <- plotly::hide_legend(p)
      } else if (legend_position != "topright"){
        p <- plotly::layout(p, legend = position_plotly_legend(legend_position))
      }
      if (zoom_x){
        p <- plotly::layout(p, xaxis = list(range = xlim))
      }
      if (zoom_y){
        p <- plotly::layout(p, yaxis = list(range = ylim))
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

#' Plot chromatograms as heatmap
#' 
#' Plots the specified traces from a list of chromatograms as a heatmap.
#' 
#' Adapted from \code{\link[VPdtw]{plot.VPdtw}}.
#' 
#' @param chrom_list List of chromatograms to plot
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' plot.
#' @param idx A vector representing the names or numerical indices of the 
#' chromatograms to plot.
#' @param xlim Range of x axis.
#' @param engine Plotting engine. Either \code{base} (\code{\link[graphics]{matplot}}), 
#' \code{\link[plotly]{plotly}}, or \link[ggplot2:ggplot2-package]{ggplot}.
#' @param show_legend Logical. Whether to display legend or not. Defaults to
#' \code{TRUE}.
#' @param legend_position Position of legend.
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots the traces of the specified chromatograms \code{idx} at the specified
#' wavelengths \code{lambdas} as a heatmap. Plots can be produced using base 
#' graphics engine, \code{ggplot2}, or \code{plotly}, according to the value of 
#' \code{engine}.
#' @author Ethan Bass
#' @examples 
#' data(Sa_warp)
#' plot_chroms_heatmap(Sa_warp, lambdas = 210)
#' @family visualization functions
#' @export
plot_chroms_heatmap <- function(chrom_list, idx = NULL, lambdas, 
                                engine = c("base", "ggplot", "plotly"),
                                show_legend = TRUE, xlim = NULL,
                                legend_position = "topright", title = "") {
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  if (ncol(chrom_list[[1]]) == 1){
    lambdas.idx <- 1
  } else {
    lambdas.idx <- get_lambda_idx(lambdas, lambdas = get_lambdas(chrom_list))
  }
  fn <- switch(engine, base = plot_chroms_heatmap_base,
               ggplot = plot_chroms_heatmap_ggplot,
               plotly = plot_chroms_heatmap_plotly)
  fn(chrom_list = chrom_list, lambdas.idx = lambdas.idx, idx = idx,
     show_legend = show_legend, xlim = xlim,
     legend_position = legend_position, title = title)
}

#' Plot chromatograms as heatmap using base graphics
#' Adapted from \code{VPdtw::plot.VPdtw}
#' @param ... Additional arguments (currently unused).
#' @noRd
plot_chroms_heatmap_base <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                     title = "", xlim = NULL, show_legend = TRUE,
                                     ...){
  viridis_base <- hcl.colors(100, "viridis")
  if (is.null(idx)) idx <- seq_along(chrom_list)
  if (inherits(chrom_list, "list")){
    x <- sapply(chrom_list[idx], function(x)x[, lambdas.idx])
  } else{
    x <- chrom_list
  }
  
  bgcol <- grey(0.9)
  if (is.null(xlim)){
    # xlim <- c(1, nrow(x))
    xlim <- c(head(get_times(chrom_list), 1), tail(get_times(chrom_list), 1))
  }
  if (show_legend){
    old_par <- par(no.readonly = TRUE)
    layout(matrix(c(1, 2), ncol = 2), widths = c(8, 1))
    par(mar = c(5, 4, 4, 1))  # c(bottom, left, top, right)
  }
  plot(xlim, c(0.8, ncol(x) + 0.2), type = "n", 
       xlab = "", ylab = "Sample", main = title, 
       xlim = xlim, yaxt = "n")
  axis(2, at = seq_len(ncol(x)), labels = colnames(x))
  rect(xlim[1] - 1000, -1000,
       xlim[2] + 1000, ncol(x) + 1000, 
       col = bgcol)
  times <- seq(xlim[1], xlim[2], by = get_time_resolution(chrom_list))
  image(times, seq_len(ncol(x)), 
        x[seq(which.min(abs(get_times(chrom_list) - head(times, 1))),
              which.min(abs(get_times(chrom_list) - tail(times, 1)))),], 
        col = viridis_base, add=TRUE)
  axis(2, at = seq_len(ncol(x)), labels = colnames(x))
  box()
  if (show_legend){
    par(mar = c(5, 0, 4, 3))  # reset margins
    
    plot(0, 0, type = "n", xlim = c(0, 1), 
         ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    
    legend_height <- 0.6
    y_start <- (1 - legend_height) / 2
    y_end <- y_start + legend_height
    
    # Generate scaled values for the legend
    n_colors <- length(viridis_base)
    y_vals <- seq(y_start, y_end, length.out = n_colors + 1)
    
    # Draw the rectangles
    for (i in 1:n_colors) {
      rect(0, y_vals[i], 1, y_vals[i + 1], col = viridis_base[i], border = NA)
    }
    
    # Add axis labels that correspond to data values
    data_range <- range(x[xlim[1]:xlim[2],], na.rm = TRUE)
    tick_vals <- seq(y_start, y_end, length.out = 5)
    tick_labels <- round(seq(data_range[1], data_range[2], length.out = 5) / 100) * 100
    axis(4, at = tick_vals, labels = tick_labels, las = 1, lwd = 0, hadj = 0.2)
    mtext("Intensity", side = 3, line = 1, padj = 5.2, adj = .1)
    par(old_par)
  }
}

#' Plot chromatograms as heatmap using ggplot2
#' @author Ethan Bass
#' @noRd
plot_chroms_heatmap_ggplot <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                       show_legend = TRUE, xlim = NULL,
                                       legend_position = "topright", title = ""){
  check_for_pkg("ggplot2")
  if (is.null(idx)) idx <- seq_along(chrom_list)
  if (inherits(chrom_list, "list")){
    x <- sapply(chrom_list[idx], function(x)x[, lambdas.idx])
  } else{
    x <- chrom_list
  }
  # Format data for ggplot
  df <- expand.grid(Index = get_times(x), 
                    Sample = colnames(x))
  df$Intensity <- as.vector(x)
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$Sample, 
                                        fill = .data$Intensity)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = "", y = "Sample", 
                  title = title) +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (!show_legend){
    p <- p + ggplot2::theme(legend.position = "none")
  } else if (legend_position != "topright"){
    p <- p + ggplot2::theme(legend.position = legend_position)
  }
  if (!is.null(xlim)){
    p <- p + ggplot2::xlim(xlim)
  }
  p
}

#' Plot chromatograms as heatmap using plotly
#' @author Ethan Bass
#' @noRd
plot_chroms_heatmap_plotly <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                       show_legend = TRUE, xlim = NULL,
                                       legend_position = "topright", title = "") {
  check_for_pkg("plotly")
  if (is.null(idx)) idx <- seq_along(chrom_list)
  if (inherits(chrom_list, "list")){
    x <- sapply(chrom_list[idx], function(x)x[, lambdas.idx])
  } else{
    x <- chrom_list
  }
  # Format data for plotly
  x_vals <- get_times(x)
  y_vals <- colnames(x)
  
  p <- plotly::plot_ly(
    x = x_vals,
    y = y_vals,
    z = t(x), 
    type = "heatmap") |>
    plotly::layout(
      title = title,
      xaxis = list(
        title = "",
        showgrid = FALSE
      ),
      yaxis = list(
        title = "Sample",
        showgrid = FALSE,
        autorange = "reversed",  # Reverses the y-axis to match the ggplot ordering
        type = "category"
      )
    )
  if (!show_legend){
    p <-  plotly::hide_colorbar(p)
  } else if (legend_position != "topright"){
    p <- plotly::layout(p, legend = position_plotly_legend(legend_position))
  }
  if (!is.null(xlim)){
    p <- plotly::layout(p, xaxis = list(range = xlim))
  }
  return(p)
}
