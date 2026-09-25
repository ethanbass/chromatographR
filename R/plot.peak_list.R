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
#' @param plot_purity Whether and how to draw peak purity traces below the
#' chromatogram. `"none"` (default) draws nothing. `"ratio"` plots the purity
#' ratio as lollipops on a log scale from 0.01 to 100, clamping values beyond
#' it: blue for pure (ratio < 1), red for impure, with a dotted reference line
#' at ratio = 1. `"curves"` draws spectral dissimilarity (1 - r^2) against the
#' noise threshold on a fixed \[0, 1\] scale. The threshold is a gray dotted
#' line, and the dissimilarity line is blue where it falls below the threshold
#' (pure) and red where it exceeds it (impure). `"curves_local"` is the same,
#' but rescales each peak to its own maximum. `TRUE` maps to `"curves"` and
#' `FALSE` to `"none"`.
#' @param purity_args A named list of additional arguments controlling the
#' purity calculation, e.g. `list(background = "none")`. See [`get_purity`].
#' Arguments that do not apply to the per-peak purity trace drawn here
#' (`metric`, `min_run` and `try`) are ignored, so the same list can be passed
#' to [`get_peaks`] and to this function.
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
                           y.offset = 25,
                           plot_purity = c("none", "ratio", "curves", "curves_local"),
                           purity_args = list(), res){
  if (isTRUE(plot_purity))  plot_purity <- "curves"
  if (isFALSE(plot_purity)) plot_purity <- "none"
  plot_purity <- match.arg(plot_purity)
  baseline <- attr(x, "baseline")
  if (is.null(baseline)) baseline <- "none"
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
  purity_lim <- NULL
  plot_args <- list(new.ts, y, type = 'l', xlab = '', ylab = '',
                    xaxt = 'n', yaxt = 'n', ...)
  yr <- range(y, na.rm = TRUE)
  if (plot_purity != "none"){
    purity_lim <- -0.1 * diff(yr)
    bands <- 2
    if (is.null(plot_args$ylim)){
      plot_args$ylim <- c(min(yr[1], purity_lim * bands), yr[2])
    }
  }
  do.call(plot, plot_args)
  if (points){
    apex_y <- if (baseline == "flat" && !is.null(pks$floor)) {
      pks$height + pks$floor
    } else if (baseline == "sloped" && !is.null(pks$floor_start) && !is.null(pks$floor_end)) {
      # apex is at rt; interpolate baseline linearly between start and end
      frac <- (pks$rt - pks$start) / (pks$end - pks$start)
      pks$height + pks$floor_start + frac * (pks$floor_end - pks$floor_start)
    } else {
      pks$height
    }
    points(pks$rt, apex_y, pch = 20, cex = cex.points, col = 'red')
  }
  if (ticks){
    draw_peak_ticks(pks, new.ts, y)
  }
  if (numbers){
    text(pks$rt, y[pks$rt] + y.offset, labels = seq_len(nrow(pks)), 
         cex = cex.font)
  }
  if (missing(res))
    res <- get_time_resolution(chrom_list)
  if (is.null(color)) {
    color <- switch(fit,
                    gaussian = "red",
                    egh = "purple",
                    bemg = "darkgreen",
                    raw = "hotpink")
  }
  fit_fun <- switch(fit,
                    gaussian = gaussian_from_peak,
                    egh      = egh_from_peak,
                    bemg     = bemg_from_peak,
                    raw      = raw_from_peak
  )
  for (i in seq_len(nrow(pks))){
    tryCatch({
      peak.loc <- seq.int((pks$start[i]),(pks$end[i]), by = res)
      yvals <- fit_fun(peak.loc, pks[i,], tfac, chrom_list, idx, lambda)
      floor_vec <- if (baseline == "flat" && !is.null(pks$floor)) {
        pks$floor[i]
      } else if (baseline == "sloped" && !is.null(pks$floor_start)) {
        seq(pks$floor_start[i], pks$floor_end[i], length.out = length(peak.loc))
      } else {
        0
      }
      draw_trapezoid(peak.loc = peak.loc, yvals = yvals,
                     color = color, alpha = alpha, floor = floor_vec)
    }, error = function(e){
      warning(
        sprintf(
          "Peak %d failed (rt=%.3f): %s",
          i, pks$rt[i], conditionMessage(e)
        ),
        call. = FALSE
      )
    })
  }
  if (plot_purity != "none"){
    plot_peak_purity(x, chrom_list, idx, lambda, new.ts, y,
                     purity_args = purity_args, lim = purity_lim,
                     purity_display = plot_purity)
  }
}

#' Draw trapezoid
#' @noRd
draw_trapezoid <- function(peak.loc, yvals, color, alpha, floor = 0){
  floor_vec <- if (length(floor) == 1) rep(floor, length(peak.loc)) else floor
  sapply(seq_len(length(peak.loc) - 1), function(i){
    polygon(peak.loc[c(i, i, (i + 1), (i + 1))],
            c(floor_vec[i], yvals[i:(i + 1)], floor_vec[i + 1]),
            col = scales::alpha(color, alpha), lty = 3, border = NA)
  })
}

#' Gaussian from peak
#' Helper for `plot.peak_list`
#' @noRd
gaussian_from_peak <- function(x, peak, tfac = 1, ...) {
  if (anyNA(c(peak$rt, peak$sd, peak$height))) {
    warning("gaussian: NA parameters for peak at rt=", peak$rt)
    return(rep(NA_real_, length(x)))
  }
  floor_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start)) {
    peak$floor_start
  } else if (!is.null(peak$floor) && !is.na(peak$floor)) peak$floor else 0
  slope_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start) &&
                   length(x) > 1 && (x[length(x)] != x[1])) {
    (peak$floor_end - peak$floor_start) / (x[length(x)] - x[1])
  } else 0
  gaussian(
    x,
    center = peak$rt,
    width  = peak$sd * tfac,
    height = peak$height,
    floor  = floor_val,
    slope  = slope_val
  )
}

#' EGH from peak
#' Helper for `plot.peak_list`
#' @noRd
egh_from_peak <- function(x, peak, tfac = 1, ...) {
  if (anyNA(c(peak$rt, peak$sd, peak$height, peak$tau))) {
    warning("egh: NA parameters for peak at rt=", peak$rt)
    return(rep(NA_real_, length(x)))
  }
  floor_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start)) {
    peak$floor_start
  } else if (!is.null(peak$floor) && !is.na(peak$floor)) peak$floor else 0
  slope_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start) &&
                   length(x) > 1 && (x[length(x)] != x[1])) {
    (peak$floor_end - peak$floor_start) / (x[length(x)] - x[1])
  } else 0
  egh(
    x,
    center = peak$rt,
    width  = peak$sd * tfac,
    height = peak$height,
    tau    = peak$tau * tfac,
    floor  = floor_val,
    slope  = slope_val
  )
}

#' BEMG from peak
#' Helper for `plot.peak_list`
#' @noRd
bemg_from_peak <- function(x, peak, tfac = 1, ...) {
  if (anyNA(c(peak$rt, peak$sd, peak$h, peak$tau_left, peak$tau_right))) {
    warning("bemg: NA parameters for peak at rt=", peak$rt)
    return(rep(NA_real_, length(x)))
  }
  center_val <- if (!is.null(peak$center) && !is.na(peak$center)) peak$center else peak$rt
  floor_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start)) {
    peak$floor_start
  } else if (!is.null(peak$floor) && !is.na(peak$floor)) peak$floor else 0
  slope_val <- if (!is.null(peak$floor_start) && !is.na(peak$floor_start) &&
                   length(x) > 1 && (x[length(x)] != x[1])) {
    (peak$floor_end - peak$floor_start) / (x[length(x)] - x[1])
  } else 0
  bemg(
    x,
    center = center_val,
    width  = peak$sd * tfac,
    height = peak$h,
    a      = peak$tau_right * tfac,
    b      = peak$tau_left * tfac,
    floor  = floor_val,
    slope  = slope_val
  )
}

#' Raw from peak
#' Helper for `plot.peak_list`
#' @noRd
raw_from_peak <- function(x, peak, tfac, chrom_list, idx, lambda) {
  chrom_list[[idx]][as.character(x), lambda]
}

#' Draw peak boundaries
#' Helper for `plot.peak_list`
#' @noRd
draw_peak_ticks <- function(pks, new.ts, y) {
  
  arrows(pks$start, y[match(pks$start, new.ts)] - 5,
         pks$start, y[match(pks$start, new.ts)] + 5,
         col = "blue", length = 0)
  
  arrows(pks$end, y[match(pks$end, new.ts)] - 5,
         pks$end, y[match(pks$end, new.ts)] + 5,
         col = "blue", length = 0)
}

## Half-span of the purity band's log scale: ratios from 0.01 to 100 are
## drawn, anything beyond is clamped to the edge.
PURITY_LOG_SPAN <- 2

#' Map purity ratios onto the plot
#' 
#' The band runs from `0` (the trace's baseline) down to `2 * lim`, with the
#' reference line at `lim`. A ratio of one lands on the reference line, `100`
#' or more on `2 * lim`, `0.01` or less on `0`. `lim` is negative, so larger
#' ratios — worse purity — draw further below the trace. A purity that could
#' not be computed stays `NA` so `polygon` skips it, matching
#' `summarize_purity`, which drops non-finite values from the reported number.
#' @noRd

scale_purity_ratio <- function(p, lim){
  lim * pmin(pmax(1 + log10(p) / PURITY_LOG_SPAN, 0), 2)
}

#' Plot peak purity
#'
#' Helper for `plot.peak_list`. See `plot.peak_list` for display modes.
#' @noRd
plot_peak_purity <- function(x, chrom_list, idx, lambda, new.ts, y,
                             purity_args = list(), lim = NULL,
                             purity_display = "ratio") {
  peaks <- x[[idx]][[lambda]][, 3:5]
  noise_variance <- hoist_noise_variance(chrom_list[[idx]], purity_args)
  purity_args <- purity_args[setdiff(
    intersect(names(purity_args), names(formals(get_purity_values))),
    "noise_variance")]
  usr <- par("usr")
  yr  <- usr[3:4]
  visible <- sum(new.ts >= usr[1] & new.ts <= usr[2])
  zoom_scale <- (length(new.ts) / max(visible, 1L))^(1/3)
  do_ratio  <- purity_display == "ratio"
  do_curves <- purity_display %in% c("curves", "curves_local")
  global_scale <- purity_display == "curves"
  drew_ratio_ref  <- FALSE
  drew_curves_ref <- FALSE
  invisible(apply(peaks, 1, function(pos) {
    pos <- match(as.numeric(pos), new.ts)
    if (anyNA(pos)) return(NULL)
    p <- do.call(get_purity_values,
                 c(list(x = chrom_list[[idx]], pos = pos,
                        noise_variance = noise_variance), purity_args))
    if (is.null(p)) return(NULL)
    keep <- attr(p, "keep")
    if (length(keep) < 2) return(NULL)
    ts_win <- new.ts[seq(pos[2], pos[3])][keep]
    if (do_curves) {
      sim  <- attr(p, "similarity")[keep]
      thr  <- attr(p, "threshold")[keep]
      num  <- 1 - sim
      den  <- 1 - thr
      scale_max <- if (global_scale) 1 else max(c(num, den), na.rm = TRUE)
      if (scale_max < .Machine$double.eps) return(NULL)
      y_num <- 2 * lim * num / scale_max
      y_den <- 2 * lim * den / scale_max
      n <- length(ts_win)
      ok <- which(!is.na(y_num) & !is.na(y_den))
      for (j in ok[ok < n]) {
        if ((j + 1L) %in% ok) {
          pure_seg <- num[j] <= den[j] && num[j + 1L] <= den[j + 1L]
          polygon(ts_win[c(j, j + 1L, j + 1L, j)],
                  c(y_num[j], y_num[j + 1L], y_den[j + 1L], y_den[j]),
                  col = scales::alpha(if (pure_seg) "steelblue" else "firebrick", 0.3),
                  border = NA)
        }
      }
      line_colors <- ifelse(num <= den,
                            scales::alpha("steelblue", 0.8),
                            scales::alpha("firebrick", 0.7))
      for (j in seq_len(n - 1L)) {
        if (!is.na(y_num[j]) && !is.na(y_num[j + 1L]))
          lines(ts_win[j:(j + 1L)], y_num[j:(j + 1L)],
                col = line_colors[j], lwd = 1)
      }
      lines(ts_win, y_den, col = "gray40", lty = 3, lwd = 1.5)
      if (!drew_curves_ref) {
        legend("topright",
               legend = c("Dissimilarity (1\u2212r\u00b2)", "Threshold"),
               col    = c("steelblue", "gray40"),
               lty    = c(1, 1), lwd = c(1, 1),
               bty = "n", cex = 0.7)
        drew_curves_ref <<- TRUE
      }
    }
    if (do_ratio) {
      if (!drew_ratio_ref) {
        abline(h = lim, lty = 3, col = "darkgray")
        legend("topright",
               legend = c("Pure (ratio < 1)", "Impure (ratio \u2265 1)"),
               col    = c(scales::alpha("steelblue", 0.8),
                          scales::alpha("firebrick",  0.8)),
               pch = 16, bty = "n", cex = 0.7)
        drew_ratio_ref <<- TRUE
      }
      yvals  <- scale_purity_ratio(p[keep], lim = lim)
      colors <- ifelse(p[keep] < 1, "steelblue", "firebrick")
      segments(ts_win, 0, ts_win, yvals, col = "gray70", lwd = 0.8 * zoom_scale)
      points(ts_win, yvals, pch = 16, cex = 0.6 * zoom_scale,
             col = scales::alpha(colors, 0.8))
    }
  }))
}