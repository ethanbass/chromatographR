#' Plot spectrum from peak table
#' 
#' Plots the trace and/or spectrum for a given peak or retention time in a 
#' `peak_table` object or a list of chromatograms.
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Retention times can
#' also be selected by clicking on the plotted trace if `what == 'click'`.
#' Plots can be produced using either base R graphics, 
#' [`ggplot2`][ggplot2::ggplot2], or `plotly`, according to the value
#' of the `engine` argument.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text
#' @importFrom utils head tail
#' @param loc The peak, retention time, or index for which you wish to
#' extract spectral data.
#' @param peak_table A `peak_table` object created by [`get_peaktable`].
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find
#' the `chrom_list` object using the pointer in the `peak_table`.
#' @param idx Numerical index of chromatogram you wish to plot, or "max" to
#' automatically plot the chromatogram with the largest signal at the given peak
#' or retention time.
#' @param lambda The wavelength you wish to plot the trace at if `plot_trace` is
#' `TRUE` and/or the wavelength to be used for the determination of signal
#' abundance.
#' @param plot_spectrum Logical. If `TRUE`, plots the spectrum of the chosen
#' peak. Defaults to `TRUE`.
#' @param plot_trace Logical. If `TRUE`, plots the trace of the chosen peak at
#' lambda. Defaults to `TRUE`.
#' @param spectrum_labels Logical. If `TRUE`, plots labels on maxima in spectral
#' plot. Defaults to `TRUE`.
#' @param scale_spectrum Logical. If `TRUE`, scales spectrum to unit height.
#' Defaults to `FALSE`.
#' @param export_spectrum Logical. If `TRUE`, exports spectrum to console.
#' Defaults to `FALSE`.
#' @param verbose Logical. If `TRUE` (default), prints verbose output to console.
#' @param what What to look for. Either `peak` to extract spectral 
#' information for a certain peak, `rt` to scan by retention time, 
#' `idx` to scan by numeric index, or `click` to manually select 
#' retention time by clicking on the chromatogram. Defaults to `"peak"` mode.
#' @param engine Which plotting engine to use: `base`, `ggplot2`, or
#' `plotly`.
#' @param ... Additional arguments.
#' @return If `export_spectrum` is `TRUE`, returns the spectrum as a 
#' `data.frame` with wavelengths as rows and a single column encoding the
#' absorbance (or normalized absorbance, if `scale_spectrum` is `TRUE`)
#' at each wavelength. If `export_spectrum` is `FALSE`, the output depends on
#' the plotting `engine`. If `engine == "plotly"`, returns a `plotly` object 
#' containing the specified plots. Otherwise, if `engine == "base"`, there is no
#' return value.
#' @section Side effects:
#' * If `plot_trace` is `TRUE`, plots the chromatographic trace of the
#' specified chromatogram (`idx`), at the specified wavelength (`lambda`) with a
#' dotted red line to indicate the retention time given by `loc`. The trace is a
#' single column from the chromatographic matrix.
#' * If `plot_spectrum` is `TRUE`, plots the spectrum for the specified
#' chromatogram at the specified retention time. The spectrum is a single row
#' from the chromatographic matrix.
#' @author Ethan Bass
#' @examplesIf interactive()
#' data(Sa)
#' pks <- get_peaks(Sa, lambda = "220.00000")
#' pk_tab <- get_peaktable(pks)
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 1))
#' plot_spectrum(loc = "V10", peak_table = pk_tab, what = "peak")
#' par(oldpar)
#' @family visualization functions
#' @export plot_spectrum
#' @md

plot_spectrum <- function(loc = NULL, peak_table, chrom_list,
                          idx = 'max', lambda = 'max',
                          plot_spectrum = TRUE, plot_trace = TRUE,
                          spectrum_labels = TRUE, scale_spectrum = FALSE, 
                          export_spectrum = FALSE, verbose = TRUE, 
                          what = c("peak", "rt", "idx", "click"),
                          engine = c('base', "plotly", "ggplot2"), ...){
  if (missing(chrom_list) & missing(peak_table))
    stop("Must provide either a peak_table or a chrom_list.")
  if (!missing(peak_table))
    check_peaktable(peak_table)
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else{
    if (!missing(peak_table)) get_chrom_list(peak_table, chrom_list)
  }
  if (!(inherits(chrom_list, c("chrom_list", "list", "matrix"))))
    stop("The provided `chrom_list` does not appear to be valid. 
                            ......Please check `chrom_list` argument")
  if (is.matrix(chrom_list)){
    chrom_list <- list(chrom_list)
    idx <- 1
  }
  what <- match.arg(what, c("peak", "rt", "idx", "click"))
  engine <- match.arg(engine, c("base", "plotly", "ggplot2"))
  if (what %in% c("peak", "rt", "idx") && is.null(loc)){
    stop("Please supply argument to `loc`")
  }
  if ((plot_spectrum | export_spectrum) & ncol(chrom_list[[1]]) == 1)
    stop("Spectral data is unidimensional.")
  if (what %in% c("rt", "idx", "click")){
    if (idx == "max")
      stop("Chromatogram must be specified for scan function.")
    if (is.null(chrom_list))
      stop("List of chromatograms must be provided for scan function.")
  } else if (what == "peak"){
    if (missing(peak_table)){
      stop("Peak table must be provided to locate peak.")}
    if (!(loc %in% colnames(peak_table$tab))){
      stop(paste0("No match found for peak \'", loc, "\' in peak table."))}
  }
  chr_idx <- check_chr(idx, loc, peak_table, chrom_list)
  plt <- switch(engine,
                "base" = plot_spectrum_base,
                "plotly" = plot_spectrum_ggpl,
                "ggplot2" = plot_spectrum_ggpl)
  
  plt(loc = loc, peak_table = peak_table, chrom_list = chrom_list,
                       chr = chr_idx, lambda = lambda,
                       plot_spectrum = plot_spectrum, plot_trace = plot_trace,
                       spectrum_labels = spectrum_labels, scale_spectrum = scale_spectrum,
                       export_spectrum = export_spectrum, verbose = verbose, 
                       what = what, engine = engine, ...)
}

#' Plot trace and/or spectrum with ggplot2 or plotly
#' @author Ethan Bass
#' @noRd
plot_spectrum_ggpl <- function(loc, peak_table, chrom_list,
                               chr = 'max', lambda = 'max',
                               plot_spectrum = TRUE, plot_trace = TRUE,
                               spectrum_labels = TRUE, scale_spectrum = FALSE,
                               export_spectrum = FALSE, verbose = TRUE, 
                               what = c("peak", "rt", "idx", "click"), 
                               zoom = FALSE, engine = c("plotly", "ggplot2"),
                               ...){
  check_for_pkg(engine)
  if (what == "click"){
    stop(paste0("The ", engine, " engine does not currently support clicking."))
  }
  new.ts <- get_times(chrom_list, idx = chr)
  new.lambdas <- get_lambdas(chrom_list)
  sig <- max(nchar(gsub(".*\\.","",rownames(chrom_list[[1]]))))
  if (what == "peak"){
    RT <- round(as.numeric(peak_table$pk_meta['rt', loc]), sig)
  } else if (what == "rt"){
    RT <- round(as.numeric(loc), sig)
  } else if (what == "idx"){
    idx <- loc
    check_idx(idx, chrom_list)
    RT <- new.ts[idx]
  }
  idx <- get_retention_idx(RT, times = new.ts)
  y <- unlist(chrom_list[[chr]][idx, , drop = TRUE])
  lambda.idx <- get_lambda_idx(lambda, lambdas = new.lambdas, y = y)
  if (plot_trace){
    plot_fn <- switch(engine, plotly = plotly_trace,
           ggplot2 = ggplot_trace)
    trace <- plot_fn(chrom_list = chrom_list, chr = chr, 
                            lambda.idx = lambda.idx, line.idx = idx)
  }
  if (verbose){
    message(sprintf("chrome no. %d (`%s`); \n RT = %g; \n lambda = %g nm; \n abs = %g",
                    chr, names(chrom_list)[chr],
                    round(RT, 2), new.lambdas[lambda.idx], 
                    round(chrom_list[[chr]][,lambda.idx][idx], 2)))

    ### report closest match ###
    if (!missing(peak_table) & what != "peak"){
      pk <- names(which.min(abs(peak_table$pk_meta["rt",] - RT)))
      message(paste("nearest peak:", pk))
    }
  }
  y <- unlist(chrom_list[[chr]][idx, , drop = TRUE])
  if (all(is.na(y))){
    stop("No data was found at the specified retention time.")
  }
  if (scale_spectrum){
    y <- scales::rescale(y)
  }
  if (plot_spectrum){
    plot_fn <- switch(engine, plotly = plotly_spec,
                        ggplot2 = ggplot_spec)
    spectrum <- plot_fn(y = y, chr = chr, RT = RT, 
                        spectrum_labels = spectrum_labels, ...)
  }
  if (plot_spectrum & plot_trace){
    combine_plots <- switch(engine,
                            plotly = purrr::partial(plotly::subplot, nrows = 2),
                            ggplot2 = purrr::partial(cowplot::plot_grid, 
                                                     nrow = 2))
    sub <- combine_plots(trace, spectrum)
  } else if (plot_spectrum & !plot_trace){
    sub <- spectrum
  } else if (!plot_spectrum & plot_trace){
    sub <- trace
  }
  if (export_spectrum){
    y <- data.frame(y)
    colnames(y) <- names(chrom_list)[chr]
    a <- attributes(chrom_list[[chr]])
    a$sample_name <- ifelse(is.null(a$sample_name), names(chrom_list)[chr], 
                            a$sample_name)
    a <- c(a, rt = RT, loc=loc)
    attr(y, "meta") <- a[-which(names(a) %in% c(
      "time_range", "time_interval", "dimnames", "row.names",
      "class", "dim", "format", "names"))]
    print(sub)
    y
  } else sub
}

#' Plot trace and/or spectrum with base R plotting engine
#' @author Ethan Bass
#' @noRd
plot_spectrum_base <- function(loc, peak_table, chrom_list,
                               chr = 'max', lambda = 'max',
                               plot_spectrum = TRUE, plot_trace = TRUE,
                               spectrum_labels = TRUE, scale_spectrum = FALSE,
                               export_spectrum = FALSE, verbose = TRUE, 
                               what=c("peak", "rt", "idx", "click"),
                               zoom = FALSE, engine="base", ...){
  new.ts <- get_times(chrom_list, idx = chr)
  new.lambdas <- get_lambdas(chrom_list)
  sig <- max(nchar(gsub(".*\\.","",rownames(chrom_list[[1]]))))
  if (what == "peak"){
    RT <- round(as.numeric(peak_table$pk_meta['rt', loc]), sig)
  } else if (what == "rt"){
    RT <- round(as.numeric(loc), sig)
  } else if (what == "idx"){
    idx <- loc
    check_idx(idx, chrom_list)
    RT <- new.ts[idx]
  } else{
    idx <- scan_chrom(chrom_list = chrom_list, peak_table = peak_table,
                      idx = chr, lambda = lambda,
                       plot_spectrum = FALSE)
    RT <- new.ts[idx]
    plot_trace <- FALSE
  }
  idx <- get_retention_idx(RT, times = new.ts)
  chr <- check_chr(chr, loc, peak_table, chrom_list)
  y <- unlist(chrom_list[[chr]][idx, , drop=TRUE])
  if (all(is.na(y))){
    stop("The peak does not exist in the specified chromatogram")
  }
  lambda.idx <- get_lambda_idx(lambda, lambdas = new.lambdas, y = y)
  if (plot_trace){
    idx <- plot_trace(chrom_list, chr, lambda.idx = lambda.idx, 
                      line.idx = idx, what = what)
  }
  if (verbose){
    message(sprintf("chrome no. %d (`%s`); \n RT = %g; \n lambda = %g nm; \n abs = %g",
                    chr, names(chrom_list)[chr],
                    round(RT, 2), new.lambdas[lambda.idx], 
                    round(chrom_list[[chr]][,lambda.idx][idx], 2)))
    
    ### report closest match ###
    if (!missing(peak_table) & what != "peak"){
      pk <- names(which.min(abs(peak_table$pk_meta["rt",] - RT)))
      message(paste("nearest peak:", pk))
    }
  }
  y <- unlist(chrom_list[[chr]][idx, , drop=TRUE])
  if (all(is.na(y))){
    stop("No data was found at the specified retention time.")
  }
  if (scale_spectrum){
    y <- scales::rescale(y)
  }
  if (plot_spectrum){
    plot_spec(y = y, spectrum_labels = spectrum_labels, main = loc, ...)
  }
  if (export_spectrum){
    y <- data.frame(y)
    colnames(y) <- names(chrom_list)[chr]
    a <- attributes(chrom_list[[chr]])
    a$sample_name <- ifelse(is.null(a$sample_name), names(chrom_list)[chr], 
                            a$sample_name)
    a <- c(a, rt = RT, loc = loc)
    attr(y, "meta") <- a[-which(names(a) %in% c(
      "time_range", "time_interval", "dimnames", "row.names",
      "class", "dim", "format", "names"))]
    y
  }
}


#' Plot spectra by clicking on the chromatogram.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text abline
#' @inheritParams plot_spectrum
#' @param ... Additional arguments to `plot_spectrum`.
#' @inherit plot_spectrum return
#' @section Side effects:
#' Plots a chromatographic trace from the specified chromatogram (`idx`),
#' at the specified wavelength (`lambda`) with a dotted red line to indicate
#' the user-selected retention time. The trace is a single column from the
#' chromatographic matrix.
#' 
#' If `plot_spectrum` is `TRUE`, plots the spectrum for the specified
#' chromatogram at the user-specified retention time. The spectrum is a single
#' row from the chromatographic matrix.
#' 
#' @author Ethan Bass
#' @examplesIf interactive()
#' data(Sa_pr)
#' scan_chrom(Sa_pr, lambda = "210", idx = 2, export_spectrum = TRUE)
#' @export scan_chrom
#' @family visualization functions
#' @md

scan_chrom <- function(chrom_list, idx, lambda,
                        plot_spectrum = TRUE, peak_table=NULL,
                        scale_spectrum = FALSE, spectrum_labels = TRUE,
                        export_spectrum = FALSE, ...){
  # check chrom_list
  if (missing(chrom_list))
    stop("List of chromatograms must be provided for scan function.")
  if (!inherits(chrom_list, c("list", "chrom_list")) | inherits(chrom_list, "peak_table"))
    stop("`chrom_list` argument should be a list of chromatograms in matrix format")
  if (missing(idx)){
    idx <- as.numeric(readline(
      prompt = "Which chromatogram do you wish to plot? \n"))
  }
  idx <- check_chr(idx, loc = NULL, peak_table, chrom_list, allow_max = FALSE)
  
  new.ts <- get_times(chrom_list, idx = idx)
  new.lambdas <- get_lambdas(chrom_list)
  sig <- max(nchar(gsub(".*\\.","",rownames(chrom_list[[1]]))))

  #check lambdas
  if (missing(lambda))
    stop("Please specify wavelength (`lambda`) to proceed with scanning.")
  # y <- unlist(chrom_list[[chr]][idx, , drop = TRUE])
  lambda.idx <- get_lambda_idx(lambda, lambdas = new.lambdas, allow_max=FALSE)
  
  line.idx <- plot_trace(chrom_list = chrom_list, chr = idx, lambda.idx = lambda.idx, 
                         what = "click")
  
  if (plot_spectrum){
    plot_spectrum(loc = line.idx, chrom_list = chrom_list,
                  idx = idx, lambda = lambda, what = "idx",
                  scale_spectrum = scale_spectrum,
                  spectrum_labels = spectrum_labels, 
                  export_spectrum = export_spectrum,
                  engine = "base", ...)
  }
  line.idx
}

#' Plot all spectra for chosen peak.
#' 
#' Plot multiple for a given peak in peak table. Wrapper for [`plot_spectrum`].
#' 
#' @inheritParams plot_spectrum
#' @param peak The name of a peak to plot (in character format).
#' @param idx Vector of chromatograms to plot.
#' @param overlapping Logical. If `TRUE`, plot spectra in single plot.
#' @param verbose Logical. If `TRUE`, prints verbose output to console. Defaults
#' to `FALSE`.
#' @param what What to look for. Either `"peak"` to extract spectral 
#' information for a certain peak, `"rt"` to scan by retention time, or 
#' `"idx"` to scan by numeric index. Defaults to "peak" mode.
#' @param ... Additional arguments to `plot_spectrum.`
#' @return If `export_spectrum` is `TRUE`, invisibly returns the 
#' spectra as a `data.frame` with wavelengths as rows and one column for 
#' each sample in the `chrom_list` encoding the absorbance (or normalized 
#' absorbance if `scale_spectrum` is `TRUE`) at each wavelength. 
#' Otherwise, there is no return value.
#' @section Side effects:
#' If `plot_spectrum` is `TRUE`, plots the spectra for the specified 
#' chromatogram (`idx`) of the given peak. Each spectrum represents a single 
#' row from the chromatographic matrix.
#' @author Ethan Bass
#' @examplesIf interactive()
#' data(Sa_warp)
#' pks <- get_peaks(Sa_warp, lambda = "220")
#' pk_tab <- get_peaktable(pks)
#' plot_all_spectra(peak = "V13", peak_table = pk_tab, overlapping = TRUE)
#' @family visualization functions
#' @export

plot_all_spectra <- function(loc, peak_table, chrom_list, idx = "all",
                             engine = c("base", "ggplot2", "plotly"),
                             plot_spectrum = TRUE, export_spectrum = TRUE,
                             scale_spectrum = TRUE, overlapping = TRUE,
                             verbose = FALSE, what = c("peak", "rt", "idx"), 
                             peak = NULL, ...){
  engine <- match.arg(engine, c("base","ggplot2","plotly"))
  what <- match.arg(what, c("peak", "rt", "idx"))
  check_peaktable(peak_table)
  if (!is.null(peak)){
    message("The `peak` argument is deprecated. Please us `loc` instead.")
    loc <- peak
  }
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  if (!(inherits(chrom_list, "list") | inherits(chrom_list, "chrom_list")))
    stop("chrom_list is not a list")
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  if ("all" %in% idx)
    idx <- seq_along(chrom_list)
  sp <- sapply(idx, function(chr){
    tryCatch({
      plot_spectrum(loc = loc, peak_table = peak_table, chrom_list = chrom_list,
                  idx = chr, plot_spectrum = FALSE, plot_trace = FALSE, 
                  export_spectrum = TRUE, scale_spectrum = scale_spectrum,
                  verbose = verbose, what = what, engine = "base")
      }, error = function(e) NA
    )
  })
  if (engine == "base"){
    sp <- as.data.frame(do.call(cbind, sp))
    colnames(sp) <- names(chrom_list)[idx]
    rownames(sp) <- colnames(chrom_list[[1]])
    if (plot_spectrum){
      if (overlapping){
        matplot(new.lambdas, sp, type = 'l', 
                xlab = 'wavelength', ylab = 'intensity', las = 2,
                main = peak)
      } else {
        apply(sp, 2, function(spp){
          plot(new.lambdas, spp, type = 'l', xlab = '', ylab = '', las = 2)
        })
      }
    }
    if (export_spectrum)
      return(invisible(sp))
  } else {
    sp <- lapply(seq_along(sp), function(i){
      data.frame(lambda = get_lambdas(chrom_list), absorbance = sp[[i]], 
                 sample = names(sp)[i])
    })
    sp <- do.call(rbind, sp)
    plot_fn <- switch(engine, plotly = plotly_spec,
                      ggplot2 = ggplot_spec)
    p <- plot_fn(sp)
    if (export_spectrum){
      if (plot_spectrum){
      print(p)
      }
      return(sp)
    } else {
      return(p)}
  }
}

#' Plot spectrum
#' @importFrom scales rescale
#' @param y Numeric vector containing spectral data.
#' @param spectrum_labels Logical. Whether to label peaks in spectrum.
#' @author Ethan Bass
#' @noRd
plot_spec <- function(y, spectrum_labels = TRUE, ...){
  matplot(x = as.numeric(names(y)), y = as.numeric(y), type = 'l',
          ylab = 'Intensity', xlab = 'Wavelength (nm)',
          ylim = c(0, max(y, na.rm = TRUE)*1.2), ...)
  if (spectrum_labels){
    suppressWarnings(pks <- find_peaks(y, slope_thresh = .00001, bounds = FALSE, 
                                       smooth_type = "none", amp_thresh = .05*max(y)))
    if (length(pks) > 0){
      pks <- data.frame(round(as.numeric(names(y)[pks]), 0), y[pks],
                        stringsAsFactors = FALSE)
      text(pks[,1], pks[,2], pks[,1], pos = 3, offset = .3, cex = .8)
    }
  }
}

#' Plot spectrum with ggplot2
#' @param y Numeric vector containing spectral data.
#' @param spectrum_labels Logical. Whether to label peaks in spectrum.
#' @author Ethan Bass
#' @noRd
ggplot_spec <- function(y, chr, RT, spectrum_labels = TRUE, color="black", 
                        width = 1.2, hide_legend = TRUE, group = TRUE, ...){
  check_for_pkg("ggplot2")
  .data <- ggplot2::.data
  if (inherits(y, "numeric")){
    df <- data.frame(lambda = as.numeric(names(y)), absorbance = y)
    group <- FALSE
  } else if (inherits(y, "list")){
    df <- reshape_chroms(x = y, idx = chr, rts = RT)
  } else if (inherits(y, "data.frame")){
    df <- y
  }
  p <- ggplot2::ggplot(data = df, 
                       ggplot2::aes(x = .data$lambda, y = .data$absorbance)) +
    ggplot2::geom_line()
  if ("sample" %in% colnames(df)){
    p <- p + ggplot2::aes(group = .data$sample, color = .data$sample)
  }
  if (spectrum_labels && inherits(y, "numeric")) {
    suppressWarnings(pks <- find_peaks(df$absorbance, slope_thresh = .00001, bounds = FALSE, 
                                       smooth_type = "none", 
                                       amp_thresh = .05*max(df$absorbance)))
    if (length(pks) > 0) {
      p <- p + ggplot2::geom_text(
        data = df[pks,],
        ggplot2::aes(x = .data$lambda, y = .data$absorbance, label = .data$lambda),
        nudge_y = .06*max(y), nudge_x=.7,
        size = 3
      )
    }
  }
  p <- p + ggplot2::theme_light() + 
    ggplot2::xlab("Wavelength") + 
    ggplot2::ylab("Absorbance")
  if (hide_legend)
    p <- p + ggplot2::theme(legend.position = "none")
  p
}


#' Plot spectrum with plotly
#' @param y Numeric vector containing spectral data.
#' @param spectrum_labels Logical. Whether to label peaks in spectrum.
#' @author Ethan Bass
#' @noRd

plotly_spec <- function(y, chr, RT, reshape = TRUE, color = "black",
                        spectrum_labels = TRUE, width = 1.2, 
                         hide_legend = TRUE, ...){
  check_for_pkg("plotly")
  if (inherits(y, "numeric")){
    df <- data.frame(lambda = as.numeric(names(y)), absorbance = y)
    group <- FALSE
  } else if (inherits(y, "list")){
    df <- reshape_chroms(x = y, idx = chr, rts = RT)
  } else if (inherits(y, "data.frame")){
    df <- y
  }
  
  if ("sample" %in% colnames(df)){
    p <- plotly::plot_ly(data = df, x = ~lambda, y = ~absorbance, 
                         color = ~sample, type = 'scatter',
                         mode = 'lines', line = list(width = width))
  } else{
    p <- plotly::plot_ly(data = df, x = ~lambda, y = ~absorbance, type = 'scatter', 
                          mode = 'lines', line = list(width = width, 
                                                      color = color))
  }
  p <-  plotly::layout(p,
                       # title = list(text=sprintf("Chr %d;   RT: %g",
                       #                           as.integer(chr), RT)
                       # ),
                       xaxis = list(title = "Wavelength"),
                       yaxis = list(title= "Absorbance (mAU)")
  )
  if (spectrum_labels) {
    if ("sample" %in% colnames(df)) {
      warning("Peak labeling for multiple samples is not implemented yet.")
    } else {
      y <- df$absorbance
      names(y) <- as.character(df$lambda)
      
      # Find peaks using the same function as in plot_spec
      suppressWarnings(pks <- find_peaks(y, slope_thresh = .00001, bounds = FALSE, 
                                         smooth_type = "none", amp_thresh = .05*max(y)))
      
      if (length(pks) > 0) {
        # Create annotations for each peak
        annotations <- list()
        
        for (i in seq_along(pks)) {
          peak_x <- as.numeric(names(y)[pks[i]])
          peak_y <- y[pks[i]]
          peak_label <- round(peak_x, 0)
          
          # Create an annotation for this peak
          ann <- list(
            x = peak_x,
            y = peak_y,
            text = as.character(peak_label),
            showarrow = FALSE,
            arrowhead = 0,
            ax = 0,
            yshift = 20
          )
          
          annotations[[i]] <- ann
        }
        
        # Add all annotations to the plot
        p <- plotly::layout(p, annotations = annotations)
      }
    }
  }
  if (hide_legend)
    p <- plotly::hide_legend(p)
  p
}

#' Plot trace
#' @param chrom_list Numeric vector containing spectral data.
#' @param line.idx What index to plot line at.
#' @author Ethan Bass
#' @noRd

plot_trace <- function(chrom_list, chr, lambda.idx, line.idx = NULL, what = ""){
  new.ts <- get_times(chrom_list, idx = chr)
  lambda <- colnames(chrom_list[[1]])[lambda.idx]
  y_trace <- chrom_list[[chr]][, lambda.idx]
  matplot(x = new.ts, y = y_trace, type = 'l', ylab = '', xlab = '')
  if (what == "click"){
    message("Click trace to select timepoint")
    line.idx <- identify(new.ts, y_trace, n = 1, plot = FALSE)
  }
  if (!is.null(line.idx))
  RT <- new.ts[line.idx]
  abline(v = RT, col = 'red', lty=3)
  title(bquote(paste("Chr ", .(chr),  " ;   RT: ", .(RT), " ;  ",
                     lambda, ": ", .(lambda), " nm")
  ))
  line.idx
}

#' Plot trace with plotly
#' @param y chrom_list A list of chromatograms in matrix format
#' @param chr Index of chromatogram to plot.
#' @param lambda.idx Index of wavelength to plot
#' @param line.idx Index to plot vertical line.
#' @author Ethan Bass
#' @noRd
plotly_trace <- function(chrom_list, chr, lambda.idx, line.idx = NULL,
                         color="black", width=1.2, hide_legend = TRUE, ...){
  check_for_pkg("plotly")
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  lambda <- colnames(chrom_list[[1]])[lambda.idx]
  RT <- new.ts[line.idx]
  y_trace <- chrom_list[[chr]][, lambda.idx]
  df <- data.frame(rt = new.ts, abs = y_trace)
  p <- plotly::plot_ly(data = df, x = ~rt, y = ~abs, 
                       type='scatter', mode = 'lines',
                       line = list(color = color, width = width, ...))
  if (!is.null(line.idx)){
    p <- plotly::add_trace(p, x = ~RT, mode = "lines",
                         line = list(dash = 3, color = "red"))
  }
  p <- plotly::layout(p,
                      title = list(text=sprintf("Chr %d;   RT: %g;  lambda: %s nm",
                                                as.integer(chr), RT, lambda)
                      ),
                      xaxis = list(title = "Wavelength"),
                      yaxis = list(title = "Absorbance (mAU)")
  )
  if (hide_legend){
    p <- plotly::hide_legend(p)
  }
  p
}

#' Plot trace with ggplot2
#' @author Ethan Bass
#' @noRd
ggplot_trace <- function(chrom_list, chr, lambda.idx, line.idx = NULL){
  check_for_pkg("ggplot2")
  .data <- ggplot2::.data
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  RT <- new.ts[line.idx]
  lambda <- colnames(chrom_list[[1]])[lambda.idx]
  p <- ggplot2::ggplot(reshape_chroms(x = chrom_list, idx = chr, 
                                      lambdas = lambda),
                       ggplot2::aes(x = .data$rt, y = .data$absorbance)) +
    ggplot2::geom_line(na.rm = TRUE)
  if (!is.null(line.idx)){
    p <- p + ggplot2::geom_vline(xintercept = RT, color = "red", 
                                 linetype = "dashed") +
      ggplot2::ggtitle(sprintf("Chr %d;   RT: %g;  \u03bb: %s nm",
                                                 as.integer(chr), RT, lambda)) +
      ggplot2::theme_light() + 
      ggplot2::xlab(add_unit("Retention time", chrom_list, unit = "time_unit", 
                             idx = chr)) +
      ggplot2::ylab(add_unit("Absorbance", chrom_list, unit = "detector_y_unit", 
                             idx = chr))
  }
  p 
}

#' Plot a chromatogram or chromatograms with an inset UV spectrum
#'
#' Plots a chromatogram (either a single trace or multiple traces) with an
#' inset UV spectrum for a specified peak. Dotted lines connect the peak apex
#' to the corners of the inset box.
#'
#' @inheritParams plot_spectrum
#' @param loc Peak name or vector of names in the peak table (e.g. `"V10"`).
#' @param idx Index of chromatogram(s) to plot. If `"max"` (default) or a
#'   scalar, plots a single trace using [plot_spectrum()]. If a vector,
#'   plots multiple traces using [plot_chroms()].
#' @param lambda Wavelength to use for the chromatographic trace and for
#'   locating the peak apex. Either a numeric value or `"max"` (default)
#'   to use the wavelength of maximum absorbance.
#' @param position Position of the inset. One of `"topright"` (default),
#'   `"topleft"`, `"top"`, `"left"`, or `"right"`. Alternatively, a numeric
#'   vector of length 2 giving the proportional coordinates of the bottom-left
#'   corner, e.g., `c(0.6, 0.6)`, used together with `inset_width` and
#'   `inset_height`.
#' @param inset_width Width of the inset as a proportion of the plot width.
#'   Defaults to `0.35`.
#' @param inset_height Height of the inset as a proportion of the plot height.
#'   Defaults to `0.35`.
#' @param ylim Y axis limits for the main plot. Passed to [plot_chroms()] or
#'   [plot_spectrum()].
#' @param xlim X axis limits for the main plot. Passed to [plot_chroms()] or
#'   [plot_spectrum()].
#' @param engine Plotting engine. Currently only `"ggplot"` is supported.
#' @param inset_text_size Base size for text elements in inset. Defaults to `8`. 
#' @param ... Additional arguments passed to [plot_chroms()] or [plot_spectrum()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_spectrum()], [plot_chroms()]
#'
#' @examples
#' \dontrun{
#' # single trace with inset spectrum
#' plot_spectrum_inset("V10", peak_table = pktab)
#'
#' # multiple traces, inset on the left
#' plot_spectrum_inset("V15", peak_table = pk_tab)
#'
#' # manual inset placement using bottom-left corner
#' plot_spectrum_inset("V15", peak_table = pk_tab, 
#'                    position=c(0.6,0.3),
#'                    inset_width = 0.35, inset_height = 0.25)
#' }
#' @export
plot_spectrum_inset <- function(loc, peak_table, chrom_list=NULL,
                                idx = "max", lambda = "max",
                                position = "topright",
                                inset_width = 0.35, inset_height = 0.35,
                                ylim = NULL, xlim = NULL,
                                engine = "ggplot", scale_spectrum = TRUE,
                                verbose = getOption("verbose"),
                                inset_text_size = 8,
                                ...) {
  check_for_pkg("ggplot2")
  chrom_list <- get_chrom_list(peak_table, chrom_list)
  peak_data <- resolve_idx_lambda(idx=idx, lambda = lambda, loc = loc, 
                                  peak_table=peak_table, chrom_list = chrom_list)
  idx <- peak_data$idx
  lambda <- peak_data$lambda
  named <- !is.null(names(loc))
  if (length(loc) > 1 && length(position) != length(loc))
    stop("When annotating multiple peaks, `position` must have the same length as `loc`.")
  
  p1 <- plot_chroms(chrom_list, idx = idx, lambdas = lambda,
                    ylim = ylim, xlim = xlim, engine = "ggplot", ...) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  # build inset
  for (i in seq_along(loc)) {
    if (length(idx) == 1){
      p2 <- plot_spectrum(loc = loc[i], peak_table = peak_table,
                          chrom_list = chrom_list, idx = idx,
                          lambda = lambda, plot_spectrum = TRUE,
                          plot_trace = FALSE, engine = engine, verbose = verbose)
    } else {
      p2 <- plot_all_spectra(loc = loc[i], peak_table = peak_table, 
                             idx = idx, engine = engine, 
                             scale_spectrum = scale_spectrum,
                             export_spectrum = FALSE, verbose = verbose)
    }
    if (named){
      p2 <- p2 + ggplot2::ggtitle(names(loc)[i]) + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    p2 <- p2 +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = "transparent", 
                                                 colour = "black"),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        axis.text.y      = ggplot2::element_blank(),
        axis.ticks.y     = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(size = inset_text_size),
        axis.text = ggplot2::element_text(size = inset_text_size - 2),
        axis.title = ggplot2::element_text(size = inset_text_size),
        plot.title = ggplot2::element_text(size = inset_text_size + 2, 
                                           hjust = 0.5)
      )
    
    # get peak coordinates
    peak_data <- get_peak_coordinates(peak_table, chrom_list = chrom_list, 
                                      loc = loc[i], idx = idx, lambda = lambda)
    coords <- peak_data$coordinates
    
    # get axis ranges
    built   <- ggplot2::ggplot_build(p1)
    x_range <- built$layout$panel_params[[1]]$x.range
    y_range <- built$layout$panel_params[[1]]$y.range
    
    peak_x <- coords[1]
    peak_y <- min(coords[2], y_range[2])
    
    # get inset position
    pos <- get_inset_position(position[[i]], inset_width, inset_height)
    
    # convert proportional to data coordinates
    xmin <- x_range[1] + pos$inset_left   * diff(x_range)
    xmax <- x_range[1] + pos$inset_right  * diff(x_range)
    ymin <- y_range[1] + pos$inset_bottom * diff(y_range)
    ymax <- y_range[1] + pos$inset_top    * diff(y_range)
    
    closest <- get_connector_corners(peak_x, peak_y, xmin, xmax, ymin, ymax, 
                                     x_range, y_range)
    
    p1 <- p1 +
      ggplot2::annotation_custom(
        grob = ggplot2::ggplotGrob(p2),
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax
      ) +
      ggplot2::annotate("segment",
                        x = peak_x, y = peak_y, xend = closest[[1]][1], 
                        yend = closest[[1]][2],
                        linetype = "dotted") +
      ggplot2::annotate("segment",
                        x = peak_x, y = peak_y, xend = closest[[2]][1], 
                        yend = closest[[2]][2],
                        linetype = "dotted")
  }
  suppressWarnings(p1)
}

#' Get inset position
#' @noRd
get_inset_position <- function(position, inset_width, inset_height) {
  if (is.numeric(position)) {
    if (length(position) == 2) {
      return(list(inset_left   = position[1],
                  inset_right  = position[1] + inset_width,
                  inset_bottom = position[2],
                  inset_top    = position[2] + inset_height))
    } else {
      stop("`position` must be a string or a numeric vector of length 2 c(left, bottom).")
    }
  }
  position <- match.arg(position, c("topright", "topleft", "top", "left", "right"))
  left <- switch(position,
                 topright = ,
                 right    = 1 - inset_width - 0.05,
                 topleft  = ,
                 left     = 0.05,
                 top      = 0.5 - inset_width / 2)
  bottom <- switch(position,
                   topright = ,
                   topleft  = ,
                   top      = 1 - inset_height - 0.05,
                   left     = ,
                   right    = 0.5 - inset_height / 2)
  list(inset_left   = left,
       inset_right  = left   + inset_width,
       inset_bottom = bottom,
       inset_top    = bottom + inset_height)
}

#' Get connector corners
#' Helper function for `plot_spectrum_inset`.
#' @noRd
get_connector_corners <- function(peak_x, peak_y, xmin, xmax, ymin, ymax, 
                                  x_range, y_range) {
  xmid <- mean(c(xmin, xmax))
  ymid <- mean(c(ymin, ymax))
  
  ndist <- function(x1, y1, x2, y2) {
    sqrt(((x1 - x2) / diff(x_range))^2 + ((y1 - y2) / diff(y_range))^2)
  }
  
  edge_dists <- c(
    left   = ndist(peak_x, peak_y, xmin, ymid),
    right  = ndist(peak_x, peak_y, xmax, ymid),
    bottom = ndist(peak_x, peak_y, xmid, ymin),
    top    = ndist(peak_x, peak_y, xmid, ymax)
  )
  
  switch(names(which.min(edge_dists)),
         left   = list(c(xmin, ymin), c(xmin, ymax)),
         right  = list(c(xmax, ymin), c(xmax, ymax)),
         bottom = list(c(xmin, ymin), c(xmax, ymin)),
         top    = list(c(xmin, ymax), c(xmax, ymax))
  )
}
