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
#' @param targets Optional target retention times, as a numeric vector or a
#' table of targets (see [`targeted_analysis`]). If supplied, the plot is
#' zoomed to the span of every peak whose apex falls within `tol` of a target,
#' or, for a target with no such peak, whose boundaries (widened by `tol`)
#' enclose it. The vertical axis is rescaled to the contents of
#' that window, so that small peaks are not flattened against the axis.
#' Defaults to `NULL`.
#' @param tol Numeric. Half-width of the window around each target, in
#' minutes, whatever the time unit of `chrom_list`. Either a single value applied to
#' every target, or a vector giving a separate tolerance for each target. A
#' `tol` column in `targets` takes precedence. Defaults to `0.1`. Ignored if
#' `targets` is `NULL`.
#' @param pad Fraction of the selected peak span to add on either side when
#' zooming. Defaults to `0.5`, which keeps neighboring peaks in view. An
#' explicit `xlim` or `ylim` passed through `...` takes precedence over the
#' zoom; if `xlim` is given, the vertical axis is rescaled to that window
#' instead.
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
#' purity calculation, e.g. `list(background_correction = "none")`. See
#' [`get_purity`].
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
#' plot(pks, targets = 11.12, ticks = TRUE)
#' op <- par(mfrow = c(1, 2))
#' for (t in c(12.19, 14.41)) plot(pks, targets = t)
#' par(op)
#' @seealso [`get_peaks`]
#' @rdname plot.peak_list
#' @family visualization functions
#' @export

plot.peak_list <- function(x, ..., chrom_list, idx = 1, lambda = NULL,
                           targets = NULL, tol = 0.1, pad = 0.5,
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
  axis_unit <- get_time_unit(chrom_list, idx = idx, na_value = "min")
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
  zoom <- peak_zoom_xlim(pks, targets, tol, pad, new.ts, rt_factor(axis_unit))
  plot_args$xlim <- plot_args$xlim %||% zoom
  win <- if (!is.null(zoom)) zoom_ylim(plot_args$xlim, pks, new.ts, y, baseline)
  yr <- win %||% range(y, na.rm = TRUE)
  if (plot_purity != "none"){
    purity_lim <- -0.1 * diff(yr)
    bands <- 2
    if (is.null(plot_args$ylim)){
      plot_args$ylim <- c(min(yr[1], purity_lim * bands), yr[2])
    }
  } else if (!is.null(win) && is.null(plot_args$ylim)){
    plot_args$ylim <- win
  }
  do.call(plot, plot_args)
  if (points){
    points(pks$rt, peak_apex_y(pks, baseline), pch = 20, cex = cex.points,
           col = 'red')
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

#' Plot matched peak pairs
#'
#' Shows how peaks in each sample were matched to the reference by
#' [`match_peaks`].
#'
#' With `type = "ladder"` (the default), each sample gets one row, with
#' retention time on the x-axis. Dotted vertical lines mark the reference
#' peaks, and a short tick marks each one in every row. A matched sample peak
#' is drawn as a point filled by its spectral correlation with the reference
#' peak, joined to the reference tick by a horizontal bar, so the length and
#' direction of the bar show the retention time shift. A tick without a bar is
#' a reference peak left unmatched in that sample, and an open circle is a
#' sample peak left unmatched. Where a bar would overlap another bar or run
#' through an unmatched peak, it is lifted slightly above its row.
#'
#' With `type = "scatter"`, each sample gets its own panel, plotting every
#' candidate pair (a pair within the retention time window and above the
#' spectral cutoff) at its reference and sample retention times, so retention
#' time drift shows as distance from the \eqn{y = x} diagonal. Candidate pairs
#' are colored by spectral correlation, matched pairs are outlined in black,
#' and unmatched peaks are marked on the axes.
#'
#' @importFrom grDevices hcl.colors
#' @importFrom graphics layout image axis box rect mtext lines points rug legend abline segments lcm grconvertX grconvertY
#' @param x A `matched_peak_list` object from [`match_peaks`], or a
#'   `matched_peak_list_set`, whose wavelengths are plotted in turn.
#' @param idx Integer or character vector of sample indices or names to plot.
#'   Defaults to all non-reference samples.
#' @param type Plot type, either `"ladder"` (default) or `"scatter"`. See
#'   *Details*.
#' @param ... Additional arguments passed to [`plot`][graphics::plot] for
#'   `type = "ladder"`, such as `xlim` to zoom in on part of the run.
#' @return Invisibly returns `x`.
#' @author Ethan Bass
#' @seealso [`match_peaks`], [`get_peaktable`]
#' @rdname plot.matched_peak_list
#' @family visualization functions
#' @export
plot.matched_peak_list <- function(x, idx = NULL, type = c("ladder", "scatter"),
                                   ...) {
  type <- match.arg(type)
  ref_idx <- reference_position(x)
  if (is.null(idx)) idx <- seq_along(x)[-ref_idx]
  if (is.character(idx)) idx <- match(idx, names(x))
  zlim <- c(attr(x, "spectral_cutoff"), 1)
  pal <- hcl.colors(100, "viridis")
  cor_col <- function(r){
    if (diff(zlim) == 0) return(rep(pal[100], length(r)))
    pal[pmax(1L, pmin(100L, as.integer((r - zlim[1]) / diff(zlim) * 99) + 1L))]
  }
  view <- if (type == "ladder") plot_matched_ladder else plot_matched_scatter
  view(x, idx, attr(x, "rt_ref"), names(x)[ref_idx], cor_col, pal, zlim, ...)
  invisible(x)
}

#' @rdname plot.matched_peak_list
#' @export
plot.matched_peak_list_set <- function(x, ...){
  for (lam in names(x)){
    args <- list(x[[lam]], ...)
    if (is.null(args$main) && !identical(args$type, "scatter"))
      args$main <- paste(lam, "nm")
    do.call(plot, args)
  }
  invisible(x)
}

#' Ladder view for `plot.matched_peak_list`
#' @noRd
plot_matched_ladder <- function(x, idx, rt_ref, ref_name, cor_col, pal, zlim,
                                xlim = NULL, ylim = NULL, xlab = NULL,
                                ylab = "", ...){
  old_par <- colorbar_layout()
  on.exit(par(old_par))
  par(mar = c(5, 6, 4, 2))
  y <- rev(seq_along(idx))
  xlim <- xlim %||%
    range(c(rt_ref, unlist(lapply(x[idx], `[[`, "rt_samp"))), na.rm = TRUE)
  ylim <- ylim %||% c(0.5, length(idx) + 0.5)
  xlab <- xlab %||% sprintf("Retention time (reference: %s)", ref_name)
  plot(NA, xlim = xlim, ylim = ylim, yaxt = "n", ylab = ylab, xlab = xlab, ...)
  axis(2, at = y, labels = names(x)[idx], las = 1, tick = FALSE)
  shade <- y[c(FALSE, TRUE)]
  if (length(shade) > 0){
    rect(par("usr")[1], shade - 0.5, par("usr")[2], shade + 0.5,
         col = "gray96", border = NA)
  }
  abline(v = rt_ref, col = "gray85", lty = 3)
  in_y <- function(inches) diff(grconvertY(c(0, inches), "inches", "user"))
  dot_radius <- 0.375 * 1.3 * par("cex") * par("ps") / 72
  tick_half <- in_y(0.5625 * dot_radius)
  lane_step <- in_y(2 * dot_radius)
  pad <- diff(grconvertX(c(0, dot_radius), "inches", "user"))
  empty <- list(matched_rows = integer(0), matched_cols = integer(0),
                rt_samp = numeric(0), cor_scores = numeric(0))
  for (k in seq_along(idx)){
    m <- x[[idx[k]]] %||% empty
    mr <- m$matched_rows
    mc <- m$matched_cols
    un_ref <- rt_ref[!seq_along(rt_ref) %in% mr]
    un_samp <- m$rt_samp[!seq_along(m$rt_samp) %in% mc]
    lane <- pack_lanes(pmin(rt_ref[mr], m$rt_samp[mc]),
                       pmax(rt_ref[mr], m$rt_samp[mc]), pad,
                       obstacles = c(un_ref, un_samp))
    yl <- y[k] + lane * lane_step
    tick_y <- rep(y[k], length(rt_ref))
    tick_y[mr] <- yl
    segments(rt_ref, tick_y - tick_half, rt_ref, tick_y + tick_half,
             col = "gray45")
    segments(rt_ref[mr], yl, m$rt_samp[mc], yl, col = "gray30", lwd = 1.5)
    points(un_samp, rep(y[k], length(un_samp)), pch = 1, col = "gray55",
           cex = 0.8)
    points(m$rt_samp[mc], yl, pch = 21, bg = cor_col(m$cor_scores), cex = 1.3)
  }
  draw_colorbar(pal, zlim)
}

#' Scatter view for `plot.matched_peak_list`
#' @noRd
plot_matched_scatter <- function(x, idx, rt_ref, ref_name, cor_col, pal, zlim,
                                 ...){
  rt_window <- attr(x, "rt_window")
  for (i in idx){
    m <- x[[i]]
    if (is.null(m)){
      message(sprintf("Skipping '%s': no matches.", names(x)[i]))
      next
    }
    old_par <- colorbar_layout()
    par(mar = c(5, 4, 4, 2))
    rt_samp <- m$rt_samp
    cand <- m$candidates
    rt_range <- range(c(rt_ref, rt_samp))
    plot(rt_range, rt_range, type = "n", xlim = rt_range, ylim = rt_range,
         xlab = paste("Reference RT:", ref_name),
         ylab = paste("Sample RT:", names(x)[i]))
    abline(0, 1, lty = 2, col = "gray50")
    if (!is.null(rt_window)){
      abline(rt_window, 1, lty = 3, col = "gray70")
      abline(-rt_window, 1, lty = 3, col = "gray70")
    }
    points(rt_ref[cand[, "ref"]], rt_samp[cand[, "samp"]],
           col = cor_col(cand[, "cor"]), pch = 19, cex = 0.9)
    points(rt_ref[m$matched_rows], rt_samp[m$matched_cols],
           pch = 21, cex = 1.4, col = "black", bg = NA, lwd = 1.5)
    rug(rt_ref[setdiff(seq_along(rt_ref), m$matched_rows)], side = 1,
        col = "gray40", lwd = 1.5)
    rug(rt_samp[setdiff(seq_along(rt_samp), m$matched_cols)], side = 2,
        col = "gray40", lwd = 1.5)
    draw_colorbar(pal, zlim)
    par(old_par)
  }
}

#' Assign overlapping intervals to lanes
#'
#' Each interval goes in the first lane (0, then 1, -1, 2, -2, ...) where it
#' clears the previous interval in that lane by `pad`. Lane 0 is also refused
#' when an obstacle falls strictly inside the interval.
#' @param lo,hi Interval ends.
#' @param pad Minimum gap between intervals sharing a lane.
#' @param obstacles Positions that intervals in lane 0 must not contain.
#' @return An integer lane for each interval.
#' @noRd
pack_lanes <- function(lo, hi, pad, obstacles = numeric(0)){
  lanes <- c(0L, 1L, -1L, 2L, -2L, 3L, -3L)
  last_end <- rep(-Inf, length(lanes))
  lane <- integer(length(lo))
  for (i in order(lo)){
    inside <- any(obstacles > lo[i] & obstacles < hi[i])
    free <- which(last_end < lo[i] - pad & !(lanes == 0L & inside))[1]
    if (is.na(free)) free <- 1L
    lane[i] <- lanes[free]
    last_end[free] <- hi[i]
  }
  lane
}

#' @noRd
colorbar_layout <- function() {
  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(1, 2), ncol = 2), widths = c(1, lcm(2.5)))
  old_par
}

#' @noRd
draw_colorbar <- function(colors, zlim, label = "Correlation") {
  par(mar = c(2, 0, 2, 3.5))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")
  legend_height <- 0.6
  y_start <- (1 - legend_height) / 2
  y_end   <- y_start + legend_height
  y_vals  <- seq(y_start, y_end, length.out = length(colors) + 1)
  rect(0, y_vals[-length(y_vals)], 1, y_vals[-1], col = colors, border = NA)
  ticks <- colorbar_ticks(zlim)
  tick_pos <- if (diff(zlim) == 0) (y_start + y_end) / 2 else
    y_start + (ticks - zlim[1]) / diff(zlim) * legend_height
  axis(4, at = tick_pos, labels = format(ticks, trim = TRUE), las = 1,
       lwd = 0, lwd.ticks = 0.75, tcl = -0.15, mgp = c(3, 0.3, 0),
       cex.axis = 0.8)
  text(0.5, y_end + 0.04, label, adj = c(0.5, 0), xpd = NA)
}

#' Round tick values inside `zlim`, in steps of 1, 2, 2.5 or 5 times a power of 10
#' @param zlim Range of the color scale.
#' @param n Target number of ticks.
#' @return Numeric vector of tick values.
#' @noRd
colorbar_ticks <- function(zlim, n = 5) {
  if (diff(zlim) == 0) return(zlim[1])
  mag <- 10^floor(log10(diff(zlim) / n) + 1e-9)
  steps <- c(1, 2, 2.5, 5, 10) * mag
  eps <- steps * 1e-9
  counts <- floor((zlim[2] + eps) / steps) - ceiling((zlim[1] - eps) / steps) + 1
  i <- which.min(abs(counts - n))
  step <- steps[i]
  eps <- eps[i]
  seq(ceiling((zlim[1] - eps) / step), floor((zlim[2] + eps) / step)) * step
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
    # `a`/`b` are rates, the reciprocals of the reported time constants
    a      = 1 / (peak$tau_left * tfac),
    b      = 1 / (peak$tau_right * tfac),
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

#' Apex height of each peak, including any fitted baseline
#' Helper for `plot.peak_list`
#' @noRd
peak_apex_y <- function(pks, baseline) {
  if (baseline == "flat" && !is.null(pks$floor)) {
    pks$height + pks$floor
  } else if (baseline == "sloped" && !is.null(pks$floor_start) &&
             !is.null(pks$floor_end)) {
    # apex is at rt; interpolate baseline linearly between start and end
    frac <- (pks$rt - pks$start) / (pks$end - pks$start)
    pks$height + pks$floor_start + frac * (pks$floor_end - pks$floor_start)
  } else {
    pks$height
  }
}

#' Horizontal limits spanning a subset of peaks
#' Helper for `plot.peak_list`
#' @noRd
peak_zoom_xlim <- function(pks, rts, tol, pad, new.ts, rtf = 1) {
  if (is.null(rts)) return(NULL)
  if (nrow(pks) == 0){
    stop("Cannot zoom on a peak: the peak list is empty.")
  }
  targets <- check_targets(rts, tol)
  rts <- targets$rt
  tol <- targets$tol
  win <- tol * rtf
  near <- function(k){
    hits <- which(abs(pks$rt - rts[k]) <= win[k])
    if (length(hits)) hits
    else which(rts[k] >= pks$start - win[k] & rts[k] <= pks$end + win[k])
  }
  hits <- lapply(seq_along(rts), near)
  sel <- sort(unique(unlist(hits)))
  if (length(sel) == 0){
    stop(sprintf(
      "No peak within %s of %s. The peaks in this chromatogram run from %s to %s.",
      paste(unique(tol), collapse = "/"),
      paste(rts, collapse = ", "),
      format(min(pks$rt)), format(max(pks$rt))))
  }
  missed <- which(lengths(hits) == 0)
  if (length(missed) > 0){
    warning(sprintf("No peak within %s of %s.",
                    paste(tol[missed], collapse = ", "),
                    paste(rts[missed], collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(pad) || length(pad) != 1 || is.na(pad) || pad < 0){
    stop("`pad` must be a single non-negative number.")
  }
  xlim <- range(c(pks$start[sel], pks$end[sel]))
  xlim + c(-1, 1) * pad * max(diff(xlim), diff(range(new.ts)) * 1e-3)
}

#' Vertical limits covering the trace, peaks and baseline within a time window
#' Helper for `plot.peak_list`
#' @noRd
zoom_ylim <- function(xlim, pks, new.ts, y, baseline) {
  w <- which(new.ts >= xlim[1] & new.ts <= xlim[2])
  sel <- which(pks$rt >= xlim[1] & pks$rt <= xlim[2])
  floors <- if (baseline == "flat" && !is.null(pks$floor)){
    pks$floor[sel]
  } else if (baseline == "sloped" && !is.null(pks$floor_start)){
    c(pks$floor_start[sel], pks$floor_end[sel])
  } else 0
  ylim <- suppressWarnings(
    range(c(y[w], floors, peak_apex_y(pks, baseline)[sel]), na.rm = TRUE))
  if (!all(is.finite(ylim)) || diff(ylim) == 0){
    ylim <- range(y, na.rm = TRUE)
  }
  ylim
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
    pos <- get_retention_idx(as.numeric(pos), new.ts, clamp = TRUE)
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