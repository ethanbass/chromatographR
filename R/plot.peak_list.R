#' Overlay fitted peak shapes on chromatograms
#' 
#' Visually assess peak integration accuracy by overlaying fitted peak shapes
#' over chromatographic traces.
#'
#' @importFrom stats median
#' @importFrom graphics polygon arrows
#' @importFrom scales alpha
#' @param x A `peak_list` object. Output from the `get_peaks` function.
#' @param ... Additional arguments to main plot function.
#' @param chrom_list List of chromatograms (retention time x wavelength
#' matrices). If missing, extracted from environment using the pointer in `x`.
#' @param idx Index or name of chromatogram to plot from `chrom_list`.
#' @param lambda Wavelength (column) to use for plotting.
#' @param points Logical. If `TRUE`, display peak apex locations as points. 
#' Defaults to `FALSE`.
#' @param ticks Logical. If `TRUE`, mark peak boundaries with tick marks. 
#' Defaults to `FALSE`.
#' @param alpha Transparency of fitted peak shapes. Defaults to `0.5`.
#' @param color Color used to fill fitted peak shapes. If `NULL`, a default 
#' color is chosen based on the fitted model type.
#' @param cex.points Size of points. Defaults to `0.5`.
#' @param numbers If `TRUE`, label peaks with numeric identifiers. Defaults to 
#' `FALSE`.
#' @param cex.font Font size if peaks are numbered. Defaults to `0.5`.
#' @param y.offset Y offset for peak numbers. Defaults to `25`.
#' @param plot_purity Logical. If `TRUE`, overlays peak purity traces based on
#' peak boundaries. Defaults to `FALSE`.
#' @param res Time resolution for peak fitting. If missing, inferred from 
#' `chrom_list`.
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots a chromatographic trace from the specified chromatogram (`chr`)
#' at the specified wavelength (`lambda`) with fitted peak shapes from the
#' provided `peak_list` drawn underneath the curve. 
#' @details
#' The appearance of fitted peaks depends on the `"fit"` attribute of `x`, which
#' may be `"gaussian"`, `"egh"`, `"bemg"`, or `"raw"`.
#'
#' Peak lists are expected to contain columns such as `rt`, `height`, `start`, 
#' `end`, and `sd`, with additional parameters depending on the fit type.
#'
#' Time units in `x` are used to rescale width parameters for plotting.
#'
#' Peak rendering errors are silently ignored.
#' 
#' @author Ethan Bass
#' @examples 
#' data(Sa_warp)
#' pks <- get_peaks(chrom_list = Sa_warp[1], lambdas = 210)
#' plot(pks, points = TRUE, ticks = TRUE)
#' @seealso [`get_peaks`]
#' @rdname plot.peak_list
#' @family visualization functions
#' @export

plot.peak_list <- function(x, ..., chrom_list, idx = 1, lambda = NULL,
                           points = FALSE, ticks = FALSE, alpha = 0.5, color = NULL,
                           cex.points = 0.5, numbers = FALSE, cex.font = 0.5, 
                           y.offset = 25, plot_purity = FALSE, res){
  time.units <- attributes(x)$time.units
  time.units <- ifelse(is.null(time.units), "", time.units)
  tfac <- switch(time.units, "min" = 1, "s" = 1/60, "ms" = 1/60000, 1)
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(x)
  }
  if (is.null(lambda)){
    lambda <- names(x[[1]])[1]
  }
  if (!(lambda %in% names(x[[1]]))){
    stop('Wavelength (`lambda`) must match one of the wavelengths in your peak list.')
  }
  if (is.numeric(lambda)){
    lambda <- as.character(lambda)
  }
  new.ts <- get_times(x = chrom_list, idx = idx)
  y <- chrom_list[[idx]][,lambda]
  pks <- data.frame(x[[idx]][[lambda]])
  fit <- attr(x,"fit")
  plot(new.ts, y, type = 'l', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ...)
  if (points){
    points(pks$rt, pks$height, pch = 20, cex = cex.points, col = 'red')
  }
  if (ticks){
    arrows(pks$start, y[which(new.ts %in% pks$start)] - 5,
           pks$start, y[which(new.ts %in% pks$start)] + 5,
           col = "blue", length = 0)
    arrows(pks$end, y[which(new.ts %in% pks$end)] - 5,
           pks$end,y[which(new.ts %in% pks$end)] + 5,
           col = "blue", length = 0)
  }
  if (numbers){
    text(pks$rt, y[pks$rt] + y.offset, labels = seq_len(nrow(pks)), 
         cex = cex.font)
  }
  if (missing(res))
    res <- get_time_resolution(chrom_list)
  for (i in seq_len(nrow(pks))){
    try({
      peak.loc <- seq.int((pks$start[i]),(pks$end[i]), by = res)
      if (fit == "gaussian"){
        yvals <- gaussian(peak.loc, center=pks$rt[i],
                          width=pks$sd[i]*tfac, height = pks$height[i])
        if (is.null(color)) color <- "red"
      }
      else if (fit == "egh"){
        yvals <- egh(x = peak.loc, center = pks$rt[i],
                     width=pks$sd[i]*tfac, height = pks$height[i],
                     tau = pks$tau[i]*tfac)
        if (is.null(color)) color <- "purple"
      } else if (fit == "bemg"){
        yvals <- bemg(x = peak.loc, center = pks$rt[i], width = pks$sd[i]*tfac,
                      height = pks$h[i], a = pks$tau_right[i]*tfac, 
                      b = pks$tau_left[i]*tfac)
        if (is.null(color)) color <- "darkgreen"
      }
      else if (fit == "raw"){
        yvals <- chrom_list[[idx]][as.character(peak.loc), lambda]
        if (is.null(color)) color <- "hotpink"
      }
      draw_trapezoid(peak.loc, yvals, color, alpha = alpha)
    }, silent = TRUE)
  }
  if (plot_purity){
    try({
      peaks <- x[[idx]][[lambda]][,3:5]
      # color <- "#FFB000"
      color="black"
      p <- apply(peaks, 1, function(pos){
        pos[1] <- which(new.ts %in% pos[[1]])
        pos[2] <- which(new.ts %in% pos[[2]])
        pos[3] <- which(new.ts %in% pos[[3]])
        pk_indices <- seq(pos[[2]], pos[[3]])
        yvals <- chrom_list[[idx]][,lambda][pk_indices]
        p <- get_purity_values(chrom_list[[idx]], pos)
        lim = -20
        draw_trapezoid(new.ts[pk_indices], scales::rescale(p, c(0, lim)), 
                       color = "black", alpha = 0.6)
        abline(h = lim, lty = 3, col = "darkgray")
      })
    }, silent = TRUE)
  }
}

#' @noRd
draw_trapezoid <- function(peak.loc, yvals, color, alpha){
  sapply(seq_len(length(peak.loc) - 1), function(i){
    polygon(peak.loc[c(i, i, (i + 1), (i + 1))], c(0, yvals[i:(i + 1)], 0),
            col = scales::alpha(color, alpha), lty = 3, border = NA)
  })
}
