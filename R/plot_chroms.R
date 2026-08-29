#' Plot traces from list of chromatograms.
#' 
#' Visualizes absorbance traces from a list of one- or multi-wavelength
#' chromatograms using base R, ggplot2, or plotly graphics. For
#' multi-wavelength chromatograms, one or more wavelengths can be selected
#' using `lambdas`. When plotting many chromatograms, traces are automatically
#' thinned according to `time_resolution` to improve rendering performance. 
#' Each chromatogram specified by `idx` is plotted as a separate trace and
#' assigned a distinct color.
#' 
#' @importFrom graphics matplot axis box
#' @param x A list of chromatograms in matrix format (timepoints × wavelengths)
#' or a `peak_table` containing a pointer to a list of chromatograms accessible
#' in the current environment.
#' @param lambdas A character or numeric vector specifying the wavelengths to 
#' plot. Ignored for one-dimensional chromatograms.
#' @param idx A vector representing the names or numerical indices of the 
#' chromatograms to plot.
#' @param time_resolution Desired temporal resolution for plotting (in minutes).
#' Chromatograms are thinned to approximately this interval prior to plotting,
#' which can substantially improve performance for large datasets.  Defaults to 
#' `0.01`.
#' @param time_unit Time units of the provided chromatograms: either `min`, `s`,
#' or `ms`. If possible, units will be detected automatically from chromatogram
#' metadata. If `time_unit` attribute is not present and no argument is 
#' provided, the time units will default to minutes.
#' @param xlim Range of x axis values.
#' @param ylim Range of y axis values.
#' @param ylab Y label. Defaults to `"Absorbance"`.
#' @param xlab X label.
#' @param engine Plotting backend to use. One of [`"base"`][graphics::matplot], 
#' [`"ggplot"`][ggplot2::ggplot2-package], or [`"plotly"`][plotly::plotly].
#' @param linewidth Line width. Defaults to `1`.
#' @param show_legend Logical. Whether to display legend or not. Defaults to 
#' `FALSE`.
#' @param legend_position Position of legend. Defaults to `"topright"`
#' @param title Title for plot.
#' @param ... Additional arguments to plotting function specified by `engine`.
#' @return A `plotly` or `ggplot` object when `engine = "plotly"` or
#' `engine = "ggplot"`, respectively. No value is returned when
#' `engine = "base"`.
#' @section Side effects:
#' If `engine == "base"`, the plot is rendered to the active graphics device.
#' @details
#' For multi-wavelength chromatograms, wavelength indices are resolved
#' internally using `get_lambda_idx()`. Time values are converted according
#' to `time_unit`, either supplied explicitly or inferred from chromatogram
#' metadata.
#'
#' When using the `"ggplot"` or `"plotly"` engines, chromatograms are first
#' reshaped into long format using `reshape_chroms()`.
#'
#' @author Ethan Bass
#' @examples 
#' data(Sa_warp)
#' plot_chroms(Sa_warp, lambdas = 210)
#' @family visualization functions
#' @export

plot_chroms <- function(x, lambdas, idx, time_resolution = 0.01,
                        time_unit = NULL,
                        xlim = NULL, ylim = NULL, 
                        xlab = "", ylab = "Absorbance",
                        engine = c("base", "ggplot", "plotly"), linewidth = 1, 
                        show_legend = FALSE, legend_position = "topright", 
                        title = "", ...){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (!inherits(x, c("list", "chrom_list"))){
    stop("Please supply a list of chromatograms.")
  }
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  if (is.null(time_unit)){
    time_unit <- get_time_unit(x, na_value = "min")
  }
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
  time_resolution <- tfac*time_resolution
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
      times <- get_times(x, idx = idx[i])
      time_diff <- mean(diff(times[1:10]))
      thin_factor <- round(time_resolution/time_diff)
      keep_idx <- seq(1, length(times), by = thin_factor)
      matplot(times[keep_idx], x[[idx[i]]][keep_idx, lambdas.idx], type = 'l',
              add = TRUE, col = i, lwd = linewidth, ...)
    }
    if (show_legend)
      legend(x = legend_position, legend = names(x)[idx], 
             fill = seq_along(x[idx]))
  } else {
    xx <- reshape_chroms(x, idx = idx, lambdas = lambdas, 
                         time_resolution = time_resolution)
    if (engine == "ggplot"){
      check_for_pkg("ggplot2")
      .data <- ggplot2::.data
      p <- ggplot2::ggplot(xx, ggplot2::aes(x = .data$rt, y = .data$absorbance,
                                        color = .data$sample)) +
        ggplot2::geom_line(linewidth = linewidth*0.5, na.rm = TRUE, ...) +
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
      if (nlevels(xx$sample) == 1){
        p <- p + ggplot2::scale_color_manual(values = "black")
      }
    } else if (engine == "plotly"){
      check_for_pkg("plotly")
      plotly_colors <- if (nlevels(xx$sample) == 1) "black" else NULL
      p <- plotly::plot_ly(data = xx, x = ~rt, y = ~absorbance, color = ~sample,
                           type = 'scatter', mode = 'lines', 
                           colors = plotly_colors,
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
    attr(p, "lambdas") <- lambdas
    attr(p, "idx") <- idx
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
#' Plots the traces of the chromatograms specified by `idx` at the specified
#' wavelengths (`lambdas`) as a heatmap. Plots can be produced using base 
#' graphics engine, `ggplot2`, or `plotly`, according to the value of the 
#' `engine` argument. Adapted from [`VPdtw::plot.VPdtw`].
#' 
#' @importFrom grDevices grey hcl.colors
#' @importFrom graphics image layout mtext rect
#' @inheritParams plot_chroms
#' @param show_legend Logical. Whether to display legend or not. Defaults to
#' `TRUE`.
#' @param show_ylabs Logical. Whether to show y labels. Defaults to `FALSE`.
#' @return A `plotly` or `ggplot` object when `engine = "plotly"` or
#' `engine = "ggplot"`, respectively. No value is returned when
#' `engine = "base"`.
#' @section Side effects:
#' If `engine == "base"`, the plot is rendered to the active graphics device.
#' @author Ethan Bass
#' @examples 
#' data(Sa_warp)
#' plot_chroms_heatmap(Sa_warp, lambdas = 210)
#' @family visualization functions
#' @export
plot_chroms_heatmap <- function(x, idx = NULL, lambdas, 
                                engine = c("base", "ggplot", "plotly"),
                                show_legend = TRUE, xlim = NULL,
                                legend_position = "topright", title = "", 
                                show_ylabs = FALSE) {
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  if (!inherits(x, c("list", "chrom_list"))){
    stop("Please supply a list of chromatograms.")
  }
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  if (ncol(x[[1]]) == 1){
    lambdas.idx <- 1
  } else {
    lambdas.idx <- get_lambda_idx(lambdas, lambdas = get_lambdas(x))
  }
  fn <- switch(engine, base = plot_chroms_heatmap_base,
               ggplot = plot_chroms_heatmap_ggplot,
               plotly = plot_chroms_heatmap_plotly)
  fn(chrom_list = x, lambdas.idx = lambdas.idx, idx = idx,
     show_legend = show_legend, xlim = xlim,
     legend_position = legend_position, title = title, 
     show_ylabs = show_ylabs)
}

#' Plot chromatograms as heatmap using base graphics
#' Adapted from `VPdtw::plot.VPdtw`
#' @param ... Additional arguments (currently unused).
#' @noRd
plot_chroms_heatmap_base <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                     title = "", xlim = NULL, show_legend = TRUE,
                                     show_ylabs = FALSE, ...){
  viridis_base <- hcl.colors(100, "viridis")
  if (is.null(idx)) idx <- seq_along(chrom_list)
  if (inherits(chrom_list, "list")){
    x <- sapply(chrom_list[idx], function(x)x[, lambdas.idx])
  } else{
    x <- chrom_list
  }
  
  bgcol <- grey(0.9)
  if (is.null(xlim)){
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
  # if (show_ylabs){
  #   axis(2, at = seq_len(ncol(x)), labels = colnames(x))
  # }
  rect(xlim[1] - 1000, -1000,
       xlim[2] + 1000, ncol(x) + 1000, 
       col = bgcol)
  times <- seq(xlim[1], xlim[2], by = get_time_resolution(chrom_list))
  image(times, seq_len(ncol(x)), 
        x[seq(which.min(abs(get_times(chrom_list) - head(times, 1))),
              which.min(abs(get_times(chrom_list) - tail(times, 1)))),], 
        col = viridis_base, add = TRUE)
  if (show_ylabs){
    axis(2, at = seq_len(ncol(x)), labels = colnames(x))
  }
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
                                       legend_position = "topright", title = "",
                                       show_ylabs = FALSE){
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
  if (!show_ylabs){
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
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
                                       legend_position = "topright", title = "",
                                       show_ylabs = FALSE) {
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
  if (!show_ylabs)
    p <- p |> plotly::layout(yaxis = list(showticklabels = FALSE, ticks=""))
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

#' Annotate peaks on a chromatogram
#'
#' Adds text labels to peaks on a chromatogram plot produced by
#' [plot_chroms()]. Peak retention times and intensities are looked up
#' automatically from the peak table.
#'
#' @param p A `ggplot` object produced by [plot_chroms()].
#' @param loc A character vector of peak names (e.g. `c("V8", "V15")`).
#'   If named, the names are used as labels (e.g.
#'   `c(Sinigrin = "V8", `1ME` = "V15")`).
#' @param peak_table A `peak_table` object created by [`get_peaktable`].
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find
#' the `chrom_list` object using the pointer in the `peak_table`.
#' @param label Labels to display at each peak. Can be a character vector of
#'   labels (one per peak), `"rt"` to use retention times from the peak table.
#'   If labels are not specified here, the names of `loc` will be applied, or
#'   if the `loc` vector is not named, the peak names will be used directly.
#' @param lambda Wavelength(s) to use for locating the peak apex. Inherited
#'   from `p` if `NULL` (default). If multiple wavelengths are provided, the
#'   one with the highest absorbance is used.
#' @param idx Index of chromatogram(s) to use for locating the peak apex.
#'   Inherited from `p` if `NULL` (default).
#' @param vjust Vertical justification of the label relative to the peak apex.
#'   Defaults to `-0.5`.
#' @param ... Additional arguments passed to [ggplot2::annotate()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_chroms()], [plot_spectrum_inset()]
#'
#' @examples
#' \dontrun{
#' plot_chroms(dat, idx = 1:5, lambda = 228, engine = "ggplot") |>
#'   annotate_peaks(c(Sinigrin = "V8", `1ME` = "V15"), peak_table = pktab)
#' }
#' @export
annotate_peaks <- function(p, loc, peak_table, chrom_list=NULL, label = NULL, 
                           lambda = NULL, idx = NULL, vjust = -0.5, ...) {
  if (!inherits(p, "ggplot")){
    stop("Please supply a ggplot object to `p`.")
  }
  lambda <- if (is.null(lambda)) attr(p, "lambdas") else lambda
  idx    <- if (is.null(idx))    attr(p, "idx")    else idx

  if (is.null(label) && !is.null(names(loc))) {
    label <- names(loc)
  } else if (identical(label, "rt")){
    label <- peak_table$pk_meta["rt", loc]
  } else if (is.null(label)) {
    label <- loc
  }
  
  for (i in seq_along(loc)) {
    peak_data <- get_peak_coordinates(peak_table, chrom_list = chrom_list, 
                                      loc[i], lambda = lambda, idx = idx)
    coords <- peak_data$coordinates
    p <- p + ggplot2::annotate("text", x = coords[1], y = coords[2], 
                               label = label[i], vjust = vjust, ...)
  }
  p
}
  