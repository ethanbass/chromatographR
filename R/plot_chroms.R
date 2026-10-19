#' Plot traces from list of chromatograms
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
#' plot. Required for multi-wavelength chromatograms and optional for
#' one-dimensional ones, where a supplied wavelength must match the column name.
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
#' @param xlim Range of x axis values. Unless `ylim` is given, the y axis is
#' then fitted to the traces within this range.
#' @param ylim Range of y axis values.
#' @param ylab Y label. Defaults to `"Absorbance"`.
#' @param xlab X label.
#' @param engine Plotting backend to use. One of [`"base"`][graphics::matplot], 
#' [`"ggplot"`][ggplot2::ggplot2-package], or [`"plotly"`][plotly::plotly].
#' @param linewidth Line width. Defaults to `1`.
#' @param show_legend Logical. Whether to display legend or not. Defaults to 
#' `FALSE`.
#' @param legend_position Position of legend. Defaults to `"topright"`. For
#' `engine = "base"`, a top or bottom keyword spreads the entries over columns
#' and, unless `ylim` is given, widens the y-range so the legend clears the
#' traces. For `engine = "ggplot"`, any value accepted by the `legend.position`
#' argument of [ggplot2::theme()]; the default leaves ggplot2's own placement.
#' @param title Title for plot.
#' @param ... Additional arguments to plotting function specified by `engine`.
#' @param what Which element to plot when each chromatogram in `x` is itself a
#' list, as returned by [chromConverter::read_chroms()] with several `what`
#' values. A single name (e.g. `"dad"`), or a vector of names for a nested
#' element (e.g. `c("instrument", "PMP1, Pressure")`). A nested element can
#' also be written as one string split at the first colon, as in
#' `"instrument:PMP1, Pressure"`, unless that string is itself an element name.
#' Defaults to `"dad"` for such lists; ignored when `x` is a list of matrices.
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

plot_chroms <- function(x, lambdas = NULL, idx, time_resolution = 0.01,
                        time_unit = NULL,
                        xlim = NULL, ylim = NULL, 
                        xlab = "", ylab = "Absorbance",
                        engine = c("base", "ggplot", "plotly"), linewidth = 1, 
                        show_legend = FALSE, legend_position = "topright", 
                        title = "", what = NULL, ...){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  x <- unnest_chrom_list(x, what)
  if (!inherits(x, c("list", "chrom_list"))){
    stop("Please supply a list of chromatograms.")
  }
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  if (is.null(time_unit)){
    time_unit <- get_time_unit(x, na_value = "min")
  }
  time_unit <- match.arg(time_unit, c("min", "s", "ms"))
  tfac <- time_factor(time_unit)
  if (missing(idx)){
    idx <- seq_along(x)
  }
  if (is.character(idx)){
    orig_idx <- idx
    idx <- match(idx, names(x))
    if (all(is.na(idx))) stop("Chromatogram could not be found. Please check names and try again.")
    if (any(is.na(idx))){
      warning(sprintf("The following chromatograms could not be found: %s",
                      format_list(orig_idx[is.na(idx)])))
      idx <- idx[!is.na(idx)]
    }
  }
  if (is.null(lambdas) && ncol(x[[1]]) == 1){
    lambdas.idx <- 1
  } else if (is.null(lambdas)){
    stop("`lambdas` must be supplied for multi-wavelength chromatograms.")
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
    ylim <- get_y_bounds(x, idx, lambdas.idx, xlim = if (zoom_x) xlim)
  }
  if (engine == "base"){
    if (show_legend){
      band <- legend_band(names(x)[idx], legend_position, ylim)
      if (!zoom_y) ylim <- band$ylim
    }
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    title(main = title, xlab = xlab, ylab = ylab)
    axis(1)
    axis(2)
    box()
    for (i in seq_along(idx)){
      times <- get_times(x, idx = idx[i])
      keep_idx <- thin_idx(times, time_resolution * tfac)
      matplot(times[keep_idx], x[[idx[i]]][keep_idx, lambdas.idx], type = 'l',
              add = TRUE, col = i, lwd = linewidth, ...)
    }
    if (show_legend)
      legend(x = legend_position, legend = names(x)[idx], 
             fill = seq_along(x[idx]), ncol = band$ncol, cex = band$cex)
  } else {
    xx <- reshape_chroms(x, idx = idx, lambdas = lambdas,
                         time_resolution = time_resolution,
                         time_unit = time_unit)
    if (engine == "ggplot"){
      check_for_pkg("ggplot2")
      .data <- ggplot2::.data
      p <- ggplot2::ggplot(xx, ggplot2::aes(x = .data$rt, y = .data$absorbance,
                                        color = .data$sample,
                                        group = interaction(.data$sample,
                                                            .data$lambda))) +
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
      lams <- unique(xx$lambda)
      p <- plotly::plot_ly(colors = plotly_colors)
      for (l in lams){
        p <- plotly::add_lines(p, data = xx[xx$lambda %in% l, ], x = ~rt,
                               y = ~absorbance, color = ~sample,
                               legendgroup = ~sample,
                               showlegend = identical(l, lams[1]),
                               line = list(width = linewidth, ...))
      }
      p <-  plotly::layout(p, xaxis = list(title = xlab),
                           yaxis = list(title = ylab),
                           title = title
      )
      if (!show_legend){
        p <- plotly::hide_legend(p)
      } else {
        pos <- if (legend_position != "topright")
          position_plotly_legend(legend_position)
        p <- plotly::layout(p, legend = c(list(tracegroupgap = 0), pos))
      }
      if (zoom_x){
        p <- plotly::layout(p, xaxis = list(range = xlim))
      }
      if (zoom_x || zoom_y){
        p <- plotly::layout(p, yaxis = list(range = ylim))
      }
    }
    attr(p, "lambdas") <- lambdas
    attr(p, "idx") <- idx
    attr(p, "time_unit") <- time_unit
    suppressWarnings(p)
  }
}

#' Lay out a top or bottom legend so it clears the traces
#'
#' Spreads the entries over as many columns as fit across the plot region and
#' widens `ylim` by a band tall enough for the resulting rows, shrinking the
#' text when the band would exceed half the plot and widening the band if it
#' still does not fit, so a keyword
#' legend at the top or bottom sits beside the data rather than on it. Other
#' positions come back as a single column with `ylim` unchanged.
#' @importFrom graphics strwidth
#' @noRd
legend_band <- function(labels, position, ylim, cex = 1){
  if (!is.character(position) || !grepl("^(top|bottom)", position))
    return(list(ncol = 1L, ylim = ylim, cex = cex))
  pin <- par("pin")
  for (size in cex * c(1, 0.9, 0.8, 0.7)){
    entry <- max(strwidth(labels, "inches", cex = size)) +
      4 * strwidth("M", "inches", cex = size)
    ncol <- max(1L, min(length(labels), floor(pin[1] / entry)))
    rows <- ceiling(length(labels) / ncol)
    frac <- (rows + 1.3) * par("csi") * size / pin[2]
    if (frac <= 0.5) break
  }
  frac <- min(0.8, frac)
  pad <- diff(ylim) * frac / (1 - frac)
  ylim <- if (startsWith(position, "top")) ylim + c(0, pad) else ylim - c(pad, 0)
  list(ncol = ncol, ylim = ylim, cex = size)
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
#' @param xlim Range of retention times to plot.
#' @param legend_position Position of legend. Defaults to `"topright"`. Cannot
#' be repositioned when using the `base` plotting engine.
#' @param zlim Range of intensities spanned by the color scale. Values outside
#' it are drawn in the end colors. Defaults to the range within `xlim`. Pass
#' the same `zlim` to several plots to put them on one scale.
#' @param lambdas A single wavelength to plot. Required for multi-wavelength
#' chromatograms and optional for one-dimensional ones, where a supplied
#' wavelength must match the column name.
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
#'
#' # zoomed in, on the color scale of the whole run
#' zlim <- range(sapply(Sa_warp, function(m) m[, "210"]), na.rm = TRUE)
#' plot_chroms_heatmap(Sa_warp, lambdas = 210, xlim = c(13, 15), zlim = zlim)
#' @family visualization functions
#' @export
plot_chroms_heatmap <- function(x, idx = NULL, lambdas = NULL, 
                                engine = c("base", "ggplot", "plotly"),
                                show_legend = TRUE, xlim = NULL, zlim = NULL,
                                legend_position = "topright", title = "", 
                                show_ylabs = FALSE) {
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  x <- unnest_chrom_list(x)
  if (!inherits(x, c("list", "chrom_list"))){
    stop("Please supply a list of chromatograms.")
  }
  engine <- match.arg(engine, c("base", "ggplot", "plotly"))
  check_zlim(zlim)
  if (is.null(lambdas) && ncol(x[[1]]) == 1){
    lambdas.idx <- 1
  } else if (is.null(lambdas)){
    stop("`lambdas` must be supplied for multi-wavelength chromatograms.")
  } else {
    if (length(lambdas) > 1){
      stop(sprintf(
        "`plot_chroms_heatmap` displays one wavelength at a time, but `lambdas` has %d.",
        length(lambdas)))
    }
    lambdas.idx <- get_lambda_idx(lambdas, lambdas = get_lambdas(x),
                                  allow_max = FALSE)
  }
  fn <- switch(engine, base = plot_chroms_heatmap_base,
               ggplot = plot_chroms_heatmap_ggplot,
               plotly = plot_chroms_heatmap_plotly)
  fn(chrom_list = x, lambdas.idx = lambdas.idx, idx = idx,
     show_legend = show_legend, xlim = xlim, zlim = zlim,
     legend_position = legend_position, title = title, 
     show_ylabs = show_ylabs)
}

#' Check a color-scale range
#' @noRd
check_zlim <- function(zlim){
  if (!is.null(zlim) && (!is.numeric(zlim) || length(zlim) != 2 ||
                         !all(is.finite(zlim)) || zlim[1] >= zlim[2])){
    stop("`zlim` must be two finite numbers in increasing order.")
  }
  invisible(zlim)
}

#' Color-scale range of a heatmap: `zlim` if given, else the plotted window
#' @noRd
heatmap_zlim <- function(x, times, xlim, zlim){
  if (!is.null(zlim)) return(zlim)
  if (!is.null(xlim)) x <- x[times >= xlim[1] & times <= xlim[2], , drop = FALSE]
  range(x, na.rm = TRUE)
}

#' Plot chromatograms as heatmap using base graphics
#' Adapted from `VPdtw::plot.VPdtw`
#' @param ... Additional arguments (currently unused).
#' @noRd
plot_chroms_heatmap_base <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                     title = "", xlim = NULL, show_legend = TRUE,
                                     show_ylabs = FALSE, zlim = NULL, ...){
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
    old_par <- colorbar_layout()
    par(mar = c(5, 4, 4, 2))  # c(bottom, left, top, right)
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
  zlim <- heatmap_zlim(x, get_times(chrom_list), xlim, zlim)
  z <- x[seq(which.min(abs(get_times(chrom_list) - head(times, 1))),
             which.min(abs(get_times(chrom_list) - tail(times, 1)))), ]
  image(times, seq_len(ncol(x)), pmin(pmax(z, zlim[1]), zlim[2]),
        col = viridis_base, add = TRUE, zlim = zlim)
  if (show_ylabs){
    axis(2, at = seq_len(ncol(x)), labels = colnames(x))
  }
  box()
  if (show_legend){
    draw_colorbar(viridis_base, zlim = zlim, label = "Intensity")
    par(old_par)
  }
}

#' Plot chromatograms as heatmap using ggplot2
#' @author Ethan Bass
#' @noRd
plot_chroms_heatmap_ggplot <- function(chrom_list, lambdas.idx = 1, idx = NULL,
                                       show_legend = TRUE, xlim = NULL, zlim = NULL,
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
    ggplot2::scale_fill_viridis_c(limits = heatmap_zlim(x, get_times(x), xlim, zlim),
                                  oob = scales::squish) +
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
                                       show_legend = TRUE, xlim = NULL, zlim = NULL,
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
  zl <- heatmap_zlim(x, x_vals, xlim, zlim)
  p <- plotly::plot_ly(
    x = x_vals,
    y = y_vals,
    z = t(x), 
    zmin = zl[1], zmax = zl[2],
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

#' Plot a chromatogram as a contour plot
#'
#' Plots a single multi-wavelength chromatogram as an image, with retention
#' time on the x axis, wavelength on the y axis, and color encoding absorbance.
#'
#' With `scale = "asinh"` (the default), color follows `asinh(x / s)`, where
#' `s` is 1/100 of the largest absolute value in the plotted region, or in
#' `zlim` when it is given. The scale is linear near zero and logarithmic
#' above `s`, so peaks up to about 100 times smaller than the tallest one stay
#' visible. The color bar is labeled in absorbance units.
#'
#' By default each plot is scaled to its own data, so colors cannot be
#' compared between plots. Pass the same `zlim` to every plot to put them on
#' one scale (see examples).
#'
#' @importFrom grDevices hcl.colors
#' @importFrom graphics image
#' @inheritParams plot_chroms
#' @param idx Name or numerical index of the chromatogram to plot. Defaults
#' to the first.
#' @param xlim Range of retention times to plot.
#' @param ylim Range of wavelengths to plot.
#' @param scale Color scale: `"asinh"` (the default) or `"linear"`.
#' @param zlim Range of absorbance values spanned by the color scale. Values
#' outside it are drawn in the end colors. Defaults to the range of the
#' plotted region.
#' @param show_legend Logical. Whether to display a color bar. Defaults to
#' `TRUE`.
#' @param title Title for plot. Defaults to the name of the chromatogram;
#' `""` omits it.
#' @param ylab Y label. Defaults to `"Wavelength"`.
#' @param ... Additional arguments to [graphics::image()], e.g.
#' `useRaster = TRUE` to speed up rendering of large chromatograms on a
#' regular grid.
#' @return No return value.
#' @section Side effects:
#' The plot is rendered to the active graphics device.
#' @author Ethan Bass
#' @examples
#' data(Sa_pr)
#' plot_chrom_contour(Sa_pr, idx = 1)
#'
#' # one color scale for every sample
#' zlim <- range(unlist(Sa_pr))
#' plot_chrom_contour(Sa_pr, idx = 1, zlim = zlim)
#' plot_chrom_contour(Sa_pr, idx = 2, zlim = zlim)
#'
#' # the same, zoomed in
#' zoomed <- crop_chroms(Sa_pr, rt_range = c(12, 13))
#' plot_chrom_contour(zoomed, idx = 1, zlim = range(unlist(zoomed)))
#' @family visualization functions
#' @export
plot_chrom_contour <- function(x, idx = 1, xlim = NULL, ylim = NULL,
                               scale = c("asinh", "linear"), zlim = NULL,
                               show_legend = TRUE, title = NULL, xlab = "",
                               ylab = "Wavelength", ...){
  if (inherits(x, "peak_table")){
    x <- get_chrom_list(x)
  }
  x <- unnest_chrom_list(x)
  if (!inherits(x, c("list", "chrom_list"))){
    stop("Please supply a list of chromatograms.")
  }
  scale <- match.arg(scale, c("asinh", "linear"))
  check_zlim(zlim)
  if (length(idx) != 1){
    stop("`plot_chrom_contour` displays one chromatogram at a time, but `idx` has ",
         length(idx), " elements.")
  }
  if (is.character(idx) && !(idx %in% names(x))){
    stop("Chromatogram could not be found. Please check names and try again.")
  }
  m <- x[[idx]]
  title <- title %||% if (is.character(idx)) idx else names(x)[idx] %||% ""
  if (ncol(m) < 2){
    stop("A contour plot needs a multi-wavelength chromatogram.")
  }
  times <- get_times(x, idx = idx)
  lambdas <- as.numeric(colnames(m))
  xlim <- xlim %||% range(times)
  ylim <- ylim %||% range(lambdas)
  rows <- times >= xlim[1] & times <= xlim[2]
  cols <- lambdas >= ylim[1] & lambdas <= ylim[2]
  if (sum(rows) < 2 || sum(cols) < 2){
    stop("`xlim` and `ylim` must each span at least two points.")
  }
  z <- m[rows, cols, drop = FALSE]
  zlim <- zlim %||% range(z, na.rm = TRUE)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  tick_labels <- function(ticks) format(ticks, trim = TRUE)
  if (scale == "asinh"){
    s <- max(abs(zlim)) / 100
    if (s > 0){
      z <- asinh(z / s)
      zlim <- asinh(zlim / s)
      tick_labels <- function(ticks) as.character(signif(s * sinh(ticks), 2))
    }
  }
  pal <- hcl.colors(100, "viridis")
  if (show_legend){
    old_par <- colorbar_layout()
    on.exit(par(old_par), add = TRUE)
    par(mar = c(5, 4, 4, 2))
  }
  image(times[rows], lambdas[cols], z, col = pal, zlim = zlim,
        xlab = xlab, ylab = ylab, main = title, ...)
  box()
  if (show_legend){
    draw_colorbar(pal, zlim = zlim, label = "Absorbance",
                  tick_labels = tick_labels)
  }
  invisible(NULL)
}

#' Annotate peaks on a chromatogram
#'
#' Adds text labels to peaks on a chromatogram plot produced by
#' [plot_chroms()]. Peaks named in a peak table are labeled at their
#' retention time, above the highest point of the chromatograms in `idx`
#' between the peak's start and end. Peaks given as retention times are
#' labeled at that time, above the highest plotted point within `tol` of it,
#' and need no peak table; the label then marks the tallest trace in that
#' window, whichever peak it belongs to.
#'
#' @param p A `ggplot` object produced by [plot_chroms()].
#' @param loc Either a character vector of peak names from `peak_table`
#'   (e.g. `c("V8", "V15")`) or a numeric vector of retention times, in the
#'   units of the plot's time axis. If named, the names are used as labels.
#' @param peak_table A `peak_table` object created by [`get_peaktable`].
#'   Required when `loc` gives peak names.
#' @param chrom_list A list of chromatograms in matrix format (timepoints ×
#' wavelengths). If no argument is provided here, the function will try to find
#' the `chrom_list` object using the pointer in the `peak_table`.
#' @param label Labels to display at each peak. Can be a character vector of
#'   labels (one per peak), or `"rt"` to use the retention times.
#'   If labels are not specified here, the names of `loc` will be applied, or
#'   if the `loc` vector is not named, the peak names will be used directly.
#' @param lambda Wavelength(s) to use for locating the peak apex. Inherited
#'   from `p` if `NULL` (default). If multiple wavelengths are provided, the
#'   one with the highest absorbance at the peak's retention time is used.
#' @param idx Index of chromatogram(s) to use for locating the peak apex.
#'   Inherited from `p` if `NULL` (default).
#'
#'   `lambda` and `idx` apply only to peak names; retention times use the
#'   plotted traces.
#' @param tol Half-width of the window around each retention time in which the
#'   highest plotted point is found, in minutes. Defaults to `0.05`. Ignored for
#'   peak names.
#' @param vjust Vertical justification of the label relative to the peak apex.
#'   Defaults to `-0.5`.
#' @param ... Additional arguments passed to [ggplot2::annotate()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_chroms()], [plot_spectrum_inset()]
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(Sa_warp)
#'   data(pk_tab)
#'   plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot") |>
#'     annotate_peaks(c(C1="V7", C2="V9", C3="V15"), peak_table = pk_tab)
#'   plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot") |>
#'     annotate_peaks(c(C1 = 12.18, C2 = 12.70))
#' }
#' @export
annotate_peaks <- function(p, loc, peak_table = NULL, chrom_list = NULL,
                           label = NULL, lambda = NULL, idx = NULL, tol = 0.05,
                           vjust = -0.5, ...) {
  if (!inherits(p, "ggplot")){
    stop("Please supply a ggplot object to `p`.")
  }
  by_rt <- is.numeric(loc)
  if (!by_rt && is.null(peak_table)){
    stop("Please supply a `peak_table`, or give `loc` as retention times.")
  }
  lambda <- if (is.null(lambda)) attr(p, "lambdas") else lambda
  idx    <- if (is.null(idx))    attr(p, "idx")    else idx

  if (is.null(label) && !is.null(names(loc))) {
    label <- names(loc)
  } else if (identical(label, "rt")){
    label <- if (by_rt) as.character(loc) else peak_table$pk_meta["rt", loc]
  } else if (is.null(label)) {
    label <- if (by_rt) as.character(loc) else loc
  }

  for (i in seq_along(loc)) {
    if (by_rt){
      near <- abs(p$data$rt - loc[[i]]) <= tol * rt_factor(attr(p, "time_unit")) &
        !is.na(p$data$absorbance)
      if (!any(near)){
        stop("No plotted data within `tol` of retention time ", loc[[i]], ".",
             call. = FALSE)
      }
      coords <- c(loc[[i]], max(p$data$absorbance[near]))
    } else {
      coords <- get_peak_coordinates(peak_table, chrom_list = chrom_list,
                                     loc[i],
                                     lambda = lambda %||%
                                       peak_table$pk_meta["lambda", loc[i]],
                                     idx = idx)$coordinates
    }
    p <- p + ggplot2::annotate("text", x = coords[1], y = coords[2], 
                               label = label[i], vjust = vjust, ...)
  }
  p
}
  