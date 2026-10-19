#' Plot spectrum from peak table
#' 
#' Plots the trace and/or spectrum for a given peak in peak table. 
#' 
#' Use it to confirm the identity of a peak, or to check that a particular
#' column in the peak table represents a single compound. It also draws simple
#' box plots of the distribution of a peak against variables defined in sample
#' metadata.
#' 
#' When `plot_trace` is `TRUE`, plots the trace of chromatogram `idx` at
#' wavelength `lambda`, with a dotted red line at the retention time given by
#' `loc`. When `plot_spectrum` is `TRUE`, plots the spectrum of the same
#' chromatogram at that retention time. When `box_plot` is `TRUE`, draws a
#' [`boxplot`] of the peak with groups defined by `vars`.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text boxplot
#' @importFrom stats as.formula
#' @inheritParams plot_spectrum
#' @param x The peak table (output from [`get_peaktable`]).
#' @param loc A vector specifying the peak(s) or retention time(s) to plot.
#' If missing when `what = "peak"`, you are prompted for a peak name. When a
#' value is returned (see below), only the first element of `loc` is used, with
#' a warning.
#' @param idx Numerical index of chromatogram you wish to plot; `"max"`
#' (default) to plot the chromatogram with the largest value for the peak in
#' the peak table; or `"all"` to plot spectra for all chromatograms.
#' @param lambda The wavelength at which to plot the trace. `"max"` (default)
#' selects the wavelength at which the spectrum at `loc` is highest.
#' @param box_plot Logical. If `TRUE`, plots box plot using categories
#' defined by `vars`.
#' @param vars Character vector of column names in the sample metadata,
#' used as the right-hand side of the boxplot formula.
#' @param verbose Logical. If `TRUE` (default), prints the index and name of
#' each chromatogram plotted, with its retention time, wavelength and
#' absorbance.
#' @param ... Additional arguments to [boxplot].
#' @return
#' * If `export_spectrum = FALSE` (default), a `plotly` or `ggplot` object, or
#' nothing if `engine == "base"`.
#' * If `export_spectrum = TRUE`, the spectrum as a `data.frame` with
#' wavelengths as rows and one column per sample giving absorbance at each
#' wavelength (normalized if `scale_spectrum = TRUE`).
#' @section Side effects:
#' If `engine == "base"`, plots are rendered to the active graphics device.
#' @author Ethan Bass
#' @rdname plot.peak_table
#' @family visualization functions
#' @export

plot.peak_table <- function(x, loc, chrom_list = NULL, what = "peak",
                            idx = 'max', lambda = 'max',
                            plot_spectrum = TRUE, plot_trace = TRUE,
                            box_plot = FALSE, vars = NULL,
                            spectrum_labels = TRUE, scale_spectrum = FALSE,
                            export_spectrum = FALSE, verbose = TRUE,
                            engine = c("base", "plotly", "ggplot"), ...){
  engine <- match.arg(engine, c("base", "plotly", "ggplot"))
  if (what == "peak" & missing(loc)){
    loc <- readline(prompt="Which peak would you like to plot? \n")
    loc <- gsub('\\"', '', loc)
    loc <- gsub("\\'", "", loc)
    if (!all(loc %in% colnames(x$tab)))
      stop("Peak(s) not found.")
  }
  if (length(loc) > 1 && (export_spectrum || engine != "base")){
    warning("Only the first element of `loc` is used when a value is returned.")
    loc <- loc[1]
  }
  out <- NULL
  if (plot_spectrum | plot_trace)
    chrom_list <- get_chrom_list(x, chrom_list)
  for (loc in loc){
    if (plot_spectrum | plot_trace){
      if (identical(idx, "all")){
        out <- plot_all_spectra(loc, x, chrom_list,
                                plot_spectrum = plot_spectrum,
                                export_spectrum = export_spectrum,
                                verbose = verbose, what = what)
      } else{
        out <- plot_spectrum(loc, x, chrom_list, idx = idx,
                             lambda = lambda, plot_spectrum = plot_spectrum,
                             plot_trace = plot_trace, spectrum_labels = spectrum_labels,
                             scale_spectrum = scale_spectrum,
                             export_spectrum = export_spectrum,
                             verbose = verbose, what = what, engine = engine)
      }
    }
    if (box_plot){
      if (!is.data.frame(x$sample_meta)){
        stop("Attach metadata to `peak_table` to make a boxplot.")
      }
      if (is.null(vars)){
        stop("Must provide independent variable(s) (`vars`) to make a boxplot.")
      }
      if (what != "peak"){
        stop("A peak name must be provided to `loc` to make a boxplot.")
      }
      boxplot(x, formula = reformulate(vars, response = loc), ...)
    }
    if (export_spectrum | engine %in% c("plotly", "ggplot")){
      return(out)
    } 
  }
}

#' Make boxplot from peak table
#' 
#' The function can take multiple response variables on the left hand side of the
#' formula (separated by `+`). In this case, a separate boxplot will be
#' produced for each response variable.
#' 
#' @param x A `peak_table` object.
#' @param formula A [formula] object.
#' @param ... Additional arguments to [`boxplot`].
#' @importFrom stats reformulate terms
#' @importFrom graphics boxplot
#' @return No return value, called for side effects.
#' @section Side effects:
#' Creates a boxplot according to the provided formula, using data from the
#' supplied `peak_table` object.
#' @author Ethan Bass
#' @examples
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' boxplot(pk_tab, formula=V9 ~ trt)
#' @family visualization functions
#' @export

boxplot.peak_table <- function(x, formula, ...){
  if (missing(formula)){
    stop("Please provide a `formula` to make a boxplot.")
  }
  lhs <- get_lhs_vars(formula)
  for (li in lhs){
    response <- paste0("x[['tab']][,'",li, "']")
    form <- reformulate(labels(terms(formula)), response = response)
    boxplot(form,
            data = x$sample_meta,
            main = paste(li, '\n',
                         'RT = ', round(as.numeric(x$pk_meta['rt', li]), 2)),
            ylab = simpleCap(x$args[["response"]]), xlab = "", ...)
  }
}

#' Make mirror plot from peak table
#' 
#' Plots chromatograms as a mirror plot.
#'
#' Compares the chromatograms in two groups, defined by the levels of `var`.
#' When `mirror` is `TRUE` both are drawn in one panel, with the second
#' inverted below the axis; when it is `FALSE` they are drawn in two stacked
#' panels.
#'
#' @importFrom graphics matplot legend plot.new plot.window par
#' @importFrom utils head tail
#' @param x The peak table (output from [`get_peaktable`]).
#' @param chrom_list A list of chromatograms in matrix format (timepoints ×
#' wavelengths). If no argument is provided here, the function will try to find the
#' `chrom_list` object used to create the `peak_table`.
#' @param lambdas The wavelength(s) at which to plot the traces. Required for
#' multi-wavelength chromatograms; ignored for one-dimensional chromatograms.
#' @param var Name of the column in the sample metadata that defines the
#' two groups.
#' @param subset Character vector of the two levels of `var` to compare.
#' Required if `var` has more than two levels.
#' @param print_legend Logical. Whether to print legend. Defaults to `TRUE`.
#' @param legend_txt Character vector containing labels for legend.
#' @param legend_pos Legend position, used only when `mirror = FALSE`. The
#' mirror layout always places the legends at the top and bottom right.
#' @param legend_size Legend size (`cex` argument). Default is `1`.
#' @param mirror Logical. Whether to plot as mirror or stacked plots.
#' Defaults to `TRUE`.
#' @param xlim Numerical vector specifying limits for x axis.
#' @param ylim Numerical vector specifying limits for y axis.
#' @param ... Additional arguments to [`matplot`].
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots are rendered to the active graphics device.
#' @author Ethan Bass
#' @examples
#' data(Sa_warp)
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' mirror_plot(pk_tab,lambdas = c("210","260"), var = "trt", mirror = TRUE, 
#'   col = c("green","blue"))
#' @family visualization functions
#' @export

mirror_plot <- function(x, chrom_list, lambdas = NULL, var, subset = NULL,
                        print_legend = TRUE, legend_txt = NULL,
                        legend_pos = "topright", legend_size = 1, mirror = TRUE, 
                        xlim = NULL, ylim = NULL, ...){
  check_peaktable(x)
  meta <- x$sample_meta
  if (missing(var))
    stop("Must provide independent variable or variables for mirror plot.")
  if (!is.data.frame(meta)){
    stop("Sample metadata must be attached to peak table to make mirror plot.")
  }
  if (!(var %in% colnames(meta))){
    stop(paste0("`", var, "`", " could not be found in sample metadata."))
  }
  chrom_list <- get_chrom_list(x, if (!missing(chrom_list)) chrom_list)
  new.ts <- round(as.numeric(rownames(chrom_list[[1]])), 2)
  if (ncol(chrom_list[[1]]) == 1){
    lambdas.idx <- 1
  } else if (is.null(lambdas)){
    stop("`lambdas` must be supplied for multi-wavelength chromatograms.")
  } else {
    lambdas.idx <- sapply(lambdas, function(lambda){
      get_lambda_idx(lambda, lambdas = get_lambdas(chrom_list),
                     allow_max = FALSE)
    }) 
  }
  fac <- factor(meta[,var])
  if (is.null(subset) & length(levels(fac)) > 2)
    stop("The mirror plot can only compare two factor levels. Please use the
    `subset` argument to specify which levels to use.")
  if (!is.null(subset)){
    if (mean(subset %in% levels(fac)) < 1){
      stop(paste0("The specified factor levels are not present in `", var, "`"))
    }
  }
  if (is.null(subset)){
    trt1 <- levels(fac)[1]
    trt2 <- levels(fac)[2]
  } else {
    trt1 <- subset[1]
    trt2 <- subset[2]
  }
  set1 <- which(meta[,var] == trt1)
  set2 <- which(meta[,var] == trt2)
  if (print_legend & is.null(legend_txt))
    legend_txt <- c(trt1,trt2)
  if (is.null(xlim))
    xlim <- c(head(new.ts, 1), tail(new.ts, 1))
  y_max <- max(sapply(chrom_list, function(xx){
    max(xx[, lambdas.idx], na.rm = TRUE)
  }))
  if (mirror){
    if (is.null(ylim))
      ylim <- c(-y_max, y_max)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    for (i in set1){
      matplot(new.ts, chrom_list[[i]][,lambdas.idx], type = 'l', 
              add = TRUE, ...)
    }
    if (print_legend)
      legend("topright", legend=legend_txt[[1]], cex = legend_size, bty = "n")
    for (i in set2){
      matplot(new.ts, -chrom_list[[i]][,lambdas.idx], type = 'l', 
              add = TRUE, ...)
    }
    if (print_legend)
      legend("bottomright", legend = legend_txt[[2]], 
             cex = legend_size, bty = "n")
  } else {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(2,1))
    if (is.null(ylim))
      ylim <- c(0, y_max)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    for (i in set1){
      matplot(new.ts, chrom_list[[i]][,lambdas.idx], type = 'l', 
              add = TRUE, ...)
    }
    if (print_legend)
      legend(legend_pos, legend = legend_txt[[1]], cex = legend_size, bty = "n")
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    for (i in set2){
      matplot(new.ts, chrom_list[[i]][,lambdas.idx], type = 'l', 
              add = TRUE, ...)
    }
    if (print_legend)
      legend(legend_pos, legend = legend_txt[[2]], cex = legend_size, bty = "n")
    par(oldpar)
  }
}
