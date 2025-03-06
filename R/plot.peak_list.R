#' Plot fitted peak shapes.
#' 
#' Visually assess integration accuracy by plotting fitted peaks over trace.
#'
#' @importFrom stats median
#' @importFrom graphics polygon arrows
#' @importFrom scales alpha
#' @param x A \code{peak_list} object. Output from the \code{get_peaks} function.
#' @param ... Additional arguments to main plot function.
#' @param chrom_list List of chromatograms (retention time x wavelength
#' matrices)
#' @param idx Index or name of chromatogram to be plotted.
#' @param index This argument is deprecated. Please use \code{idx} instead.
#' @param lambda Wavelength for plotting.
#' @param points Logical. If TRUE, plot peak maxima. Defaults to FALSE.
#' @param ticks Logical. If TRUE, mark beginning and end of each peak. Defaults
#' to FALSE.
#' @param a Alpha parameter controlling the transparency of fitted shapes.
#' @param color The color of the fitted shapes.
#' @param cex.points Size of points. Defaults to 0.5
#' @param numbers Whether to number peaks. Defaults to FALSE.
#' @param cex.font Font size if peaks are numbered. Defaults to 0.5.
#' @param y.offset Y offset for peak numbers. Defaults to 25.
#' @param plot_purity Whether to add visualization of peak purity.
#' @param res time resolution for peak fitting
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots a chromatographic trace from the specified chromatogram (\code{chr})
#' at the specified wavelength (\code{lambda}) with fitted peak shapes from the
#' provided \code{peak_list} drawn underneath the curve. 
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @rdname plot.peak_list
#' @family visualization functions
#' @export

plot.peak_list <- function(x, ..., chrom_list, idx = 1, lambda = NULL,
                           points = FALSE, ticks = FALSE, a = 0.5, color = NULL,
                           cex.points = 0.5, numbers = FALSE, cex.font = 0.5, 
                           y.offset = 25, plot_purity = FALSE, res, index = NULL){
  if (!is.null(index)){
    idx <- index
    message("The `index` argument is deprecated. Please use `idx` instead")
  }
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
  if ("r.squared" %in% colnames(pks)){
    fit <- ifelse("tau" %in% colnames(pks), "egh", "gaussian")
  } else{
    fit <- "raw"
  }
  plot(new.ts, y, type = 'l', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ...)
  if (points){
    points(pks$rt, pks$height, pch = 20, cex=cex.points, col = 'red')
  }
  if (ticks){
    arrows(pks$start, y[which(new.ts %in% pks$start)] - 5,
           pks$start, y[which(new.ts %in% pks$start)] + 5,
           col="blue", length = 0)
    arrows(pks$end, y[which(new.ts %in% pks$end)] - 5,
           pks$end,y[which(new.ts %in% pks$end)] + 5,
           col="blue", length = 0)
  }
  if (numbers){
    text(pks$rt, y[pks$rt] + y.offset, labels = seq_len(nrow(pks)), 
         cex = cex.font)
  }
  if (missing(res))
    res <- get_time_resolution(chrom_list)
  for (i in seq_len(nrow(pks))){
    try({
      peak.loc<-seq.int((pks$start[i]),(pks$end[i]), by = res)
      if (fit == "gaussian"){
        yvals <- gaussian(peak.loc, center=pks$rt[i],
                          width=pks$sd[i]*tfac, height = pks$height[i])
        if (is.null(color))
          color <- "red"
      }
      else if (fit == "egh"){
        yvals <- egh(x = peak.loc, center = pks$rt[i],
                     width=pks$sd[i]*tfac, height = pks$height[i],
                     tau = pks$tau[i]*tfac)
        if (is.null(color))
          color <- "purple"
      }
      else if (fit == "raw"){
        yvals <- chrom_list[[idx]][as.character(peak.loc), lambda]
        if (is.null(color))
          color <- "hotpink"
      }
      draw_trapezoid(peak.loc, yvals, color, a)
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
                       color = "black", a = 0.6)
        abline(h = lim, lty = 3, col = "darkgray")
      })
    }, silent = TRUE)
  }
}

#' @noRd
draw_trapezoid <- function(peak.loc, yvals, color, a){
  sapply(seq_len(length(peak.loc) - 1), function(i){
    polygon(peak.loc[c(i, i, (i + 1), (i + 1))], c(0, yvals[i:(i + 1)], 0),
            col = scales::alpha(color, a), lty = 3, border = NA)
  })
}
