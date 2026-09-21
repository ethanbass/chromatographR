#' Plot an estimated baseline
#'
#' Plots the baseline retained by [`preprocess`] or [`correct_baseline`] over
#' the trace it was estimated from, annotated with the smoothing parameter that
#' produced it.
#'
#' Only simple diagnostic plots using base graphics are supported currently. 
#' For more complex plots, the baseline can be extracted using [`get_baseline`].
#'
#' @importFrom graphics lines mtext
#' @param x A matrix with an attached `"baseline"`, or a list of such matrices,
#' as returned by [`preprocess`] or [`correct_baseline`] with
#' `attach_baseline = TRUE`.
#' @param lambda Wavelength to plot, as a string or number. Optional for a
#' single-wavelength chromatogram.
#' @param idx Index or name of the chromatogram to plot, if `x` is a list.
#' Defaults to `1`.
#' @param zoom Logical. Whether to scale the y axis to the baseline, letting
#' the peaks run off the top of the panel. A baseline is usually a small
#' fraction of the signal, so at full scale it is pressed against the axis and
#' the failure mode worth looking for -- the estimate bending up into a peak --
#' is not visible. Defaults to `TRUE`; set `FALSE` to see the whole trace.
#' @param col,lty,lwd Graphical parameters for the baseline itself. The trace
#' beneath it is always drawn in gray.
#' @param add Logical. Whether to draw the baseline over an existing plot
#' rather than starting a new one, for comparing two estimates on the same
#' trace. The title and annotation are drawn only by the first call, and `...`
#' is ignored. Defaults to `FALSE`.
#' @param annotate Logical. Whether to annotate the plot with the estimator and
#' its smoothing parameter. Set `FALSE` when overlaying two estimates with
#' `add`, since a single annotation cannot describe both; label them with
#' [`legend`][graphics::legend] instead. Defaults to `TRUE`.
#' @param ... Additional arguments to [`plot`][graphics::plot.default], which
#' override the defaults chosen here -- `main = ""` suppresses the title, and
#' an explicit `ylim` overrides `zoom`.
#' @return The annotation label, invisibly, so it can be reused in a
#' [`legend`][graphics::legend]. The plot is rendered to the active graphics
#' device.
#' @author Ethan Bass
#' @seealso [`get_baseline`] [`preprocess`] [`correct_baseline`]
#' @examples \donttest{
#' data(Sa)
#' x <- preprocess(Sa[1], dim1 = seq(10, 18, by = .01),
#'                 dim2 = seq(200, 318, by = 2), attach_baseline = TRUE,
#'                 show_progress = FALSE, cl = 1)
#' plot_baseline(x, lambda = 210)
#' }
#' @md
#' @export
plot_baseline <- function(x, lambda = NULL, idx = 1, zoom = TRUE,
                          col = "red", lty = 1, lwd = 2, add = FALSE,
                          annotate = TRUE, ...){
  samp <- if (inherits(x, "list")) names(x)[[if (is.character(idx))
    match(idx, names(x)) else idx]] else NULL
  if (inherits(x, "list")) x <- x[[idx]]
  bl <- get_baseline(x)
  lambdas <- colnames(x)
  if (ncol(x) == 1){
    j <- 1L
  } else if (is.null(lambda)){
    stop("`lambda` must be supplied for multi-wavelength chromatograms.")
  } else {
    j <- get_lambda_idx(lambda, lambdas, y = x, allow_max = FALSE)
  }
  tt <- as.numeric(rownames(x))
  y <- x[, j] + bl[, j]
  wl <- if (!is.null(lambdas)) paste(format_axis(lambdas[j]), "nm") else NULL
  main <- if (!is.null(samp) && !is.null(wl)) sprintf("%s (%s)", samp, wl)
    else samp %||% wl %||% ""
  if (add){
    lines(tt, bl[, j], col = col, lty = lty, lwd = lwd)
    return(invisible(baseline_param_label(
      attr(x, "baseline_param", exact = TRUE))))
  }
  dots <- list(...)
  defaults <- list(type = "l", col = "darkgray", main = main,
                   xlab = add_unit("Time", list(x), idx = 1),
                   ylab = "Absorbance")
  if (zoom)
    defaults$ylim <- range(bl[, j]) + c(-0.1, 1) * diff(range(bl[, j]))
  defaults <- defaults[setdiff(names(defaults), names(dots))]
  do.call(plot, c(list(tt, y), defaults, dots))
  lines(tt, bl[, j], col = col, lty = lty, lwd = lwd)
  label <- baseline_param_label(attr(x, "baseline_param", exact = TRUE))
  if (annotate && !is.null(label))
    mtext(label, side = 3, line = 0.25, cex = 0.85)
  invisible(label)
}

#' Format a baseline smoothing parameter for annotation
#' @param param A named numeric, as attached by `subtract_baseline`.
#' @return A single string, or `NULL` if `param` is `NULL`.
#' @author Ethan Bass
#' @noRd
baseline_param_label <- function(param){
  if (is.null(param)) return(NULL)
  dt <- attr(param, "dt")
  paste(c(attr(param, "method"),
          paste0(names(param), " = ", format(param, digits = 4)),
          if (!is.null(dt) && isTRUE(is.finite(dt)))
            paste0("dt = ", format(dt, digits = 3))),
        collapse = ", ")
}
