#' Plot PTW alignments
#' 
#' Plots [`ptw`][ptw::ptw] alignments.
#' 
#' Plots raw and warped chromatograms, one above the other, either as
#' individual traces or as a heatmap, according to the value of `what`. The
#' plot can be produced using either base R graphics, `ggplot2`, or `plotly`,
#' according to the value of the `engine` argument.
#' 
#' @importFrom graphics matplot
#' @param x A `ptw_list` object created by [`correct_rt`].
#' @param what What type of plot to return. Either `"traces"` (default) or
#' `"heatmap"`.
#' @param engine What plotting engine to use. One of `"base"` (default),
#' `"ggplot"`, or `"plotly"`.
#' @param lambdas Wavelengths to plot, from among those used for warping.
#' Defaults to all of them. The heatmap shows only the first.
#' @param show_legend Logical. Whether to include sample legend.
#' @param ... Ignored.
#' @return A `plotly` or `ggplot` object when `engine = "plotly"` or
#' `engine = "ggplot"`, respectively. No value is returned when
#' `engine = "base"`.
#' @section Side effects:
#' If `engine == "base"`, plots are rendered to the active graphics device.
#' @author Ethan Bass
#' @examples
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
  
  lambdas_missing <- missing(lambdas)
  if (lambdas_missing){
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
        p2 <- p2 |> plotly::style(showlegend = FALSE) |> 
          plotly::layout(title = "")
      }
      combine_plots <- switch(engine,
                              plotly = purrr::partial(plotly::subplot, nrows = 2),
                              ggplot = purrr::partial(cowplot::plot_grid, 
                                                      nrow = 2))
      p_c <- combine_plots(p1, p2)
      if (engine == "plotly"){
        p_c <- p_c |> plotly::layout(
          annotations = list(
            list(x = 0.5, y = 1, text = "Raw data", showarrow = FALSE, 
                 xref = 'paper', yref = 'paper', font = list(size = 16)),
            list(x = 0.5, y = 0.45, text = "PTW alignments", showarrow = FALSE, 
                 xref = 'paper', yref = 'paper', font = list(size = 16))
          )
        )
      }
      p_c
    }
  } else if (what == "heatmap"){
    if (length(lambdas) > 1 && !lambdas_missing)
      warning("The heatmap shows one wavelength; using ", lambdas[1], ".")
    li <- match(lambdas[1], all.lambdas)
    original <- sapply(x, function(s) s$sample[li, ])
    warped <- sapply(x, function(s) s$warped.sample[li, ])
    times <- colnames(x[[1]]$reference)
    rownames(original) <- colnames(x[[1]]$reference)
    rownames(warped) <- colnames(x[[1]]$reference)
    zlim <- range(original, warped, na.rm = TRUE)
    if (engine == "base"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      if (show_legend){
        layout(matrix(c(1, 3, 2, 3), 2, byrow = TRUE),
               widths = c(1, lcm(2.5)))
      } else {
        par(mfrow = c(2, 1))
      }
      par(mar = c(3, 4, 3, 2))
      plot_chroms_heatmap_base(original, show_legend = FALSE,
                               title = "Raw data", zlim = zlim)
      plot_chroms_heatmap_base(warped, show_legend = FALSE,
                               title = "PTW alignment", zlim = zlim)
      if (show_legend)
        draw_colorbar(hcl.colors(100, "viridis"), zlim = zlim,
                      label = "Intensity")
    } else{
      if (engine == "ggplot"){
        p_o <- plot_chroms_heatmap_ggplot(original, title = "Raw data", zlim = zlim)
        p_w <- plot_chroms_heatmap_ggplot(warped, title = "PTW alignment",
                                          zlim = zlim)
        p <- cowplot::plot_grid(p_o + ggplot2::theme(legend.position = "none"),
                                p_w + ggplot2::theme(legend.position = "none"), 
                                nrow = 2)
        if (show_legend){
          legend <- suppressWarnings(cowplot::get_legend(
            p_o + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 6, 0, 0))
          ))
          p <- cowplot::plot_grid(p, legend, ncol = 2, rel_widths  = c(1, 0.2))
        }
      } else if (engine == "plotly"){
        p_o <- plot_chroms_heatmap_plotly(original, zlim = zlim)
        p_w <- plot_chroms_heatmap_plotly(warped, show_legend = FALSE, zlim = zlim)
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

#' @method [ ptw_list
#' @export
#' @keywords internal
`[.ptw_list` <- function(x, i) {
  structure(NextMethod(), class = "ptw_list")
}  
