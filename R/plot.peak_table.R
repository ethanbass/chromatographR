#' Plot spectrum from peak table
#' 
#' Plots the trace and/or spectrum for a given peak in peak table. 
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Can also be used
#' to create simple box-plots to examine the distribution of a peak with respect
#' to variables defined in sample metadata.
#' 
#' When `plot_trace` is `TRUE`, plots the chromatographic trace of the 
#' specified chromatogram (`idx`), at the specified wavelength 
#' (`lambda`) with a dotted red line to indicate the retention time given 
#' by `loc`. The trace is a single column from the chromatographic matrix. When
#' `plot_spectrum` is `TRUE`, plots the spectrum for the specified 
#' chromatogram at the specified retention time. The spectrum represents a single 
#' row from the chromatographic matrix. When `box_plot` is `TRUE`, produces a 
#' [`boxplot`] from the specified peak with groups defined by the `vars` argument.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text boxplot
#' @importFrom stats as.formula
#' @inheritParams plot_spectrum
#' @param x The peak table (output from [`get_peaktable`]).
#' @param loc A vector specifying the peak(s) or retention time(s) to plot.
#' @param idx Numerical index of chromatogram you wish to plot; "max" to
#' plot the chromatogram with the largest signal; or "all" to plot spectra
#' for all chromatograms.
#' @param lambda The wavelength you wish to plot the trace at (if `plot_chrom` 
#' is `TRUE`. Otherwise, the wavelength to be used for the determination
#' of signal abundance.
#' @param box_plot Logical. If `TRUE`, plots box plot using categories
#' defined by `vars`.
#' @param vars Independent variables for boxplot. Righthand side of formula.
#' @param verbose Logical. If `TRUE` (default), prints verbose output to the
#' console.
#' @param ... Additional arguments to [boxplot].
#' @return 
#' * If `export_spectrum = FALSE` (default), a `plotly` or `ggplot` object, or
#' nothing if `engine == "base"`. 
#' * If `export_spectrum = TRUE`, invisibly returns the spectrum as a
#' `data.frame` with wavelengths as rows and a single column per sample encoding absorbance at each 
#' wavelength (normalized if `scale_spectrum = TRUE`).
#' Otherwise, if `engine == "base"`, there is no 
#' return value. If `export_spectrum` is `TRUE`, returns the spectrum as a 
#' `data.frame` with wavelengths as rows and columns encoding the
#' absorbance (or normalized absorbance, if `scale_spectrum` is `TRUE`) for 
#' the specified sample(s).
#' @section Side effects:
#' If `engine == "base"`, plots are rendered to the active graphics device.
#' @author Ethan Bass
#' @rdname plot.peak_table
#' @family visualization functions
#' @export

plot.peak_table <- function(x, loc, chrom_list, what = "peak",
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
        stop("Must provide independent variable(s) (`var`) to make a boxplot.")
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

#' Make boxplot from peak table.
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
#' boxplot(pk_tab, formula=V11 ~ trt)
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

#' Make mirror plot from peak table.
#' 
#' Plots chromatograms as a mirror plot.
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Can also be used
#' to create simple box-plots to examine the distribution of a peak with respect
#' to variables defined in sample metadata.
#' 
#' @importFrom graphics matplot legend plot.new plot.window par
#' @importFrom utils head tail
#' @param x The peak table (output from [`get_peaktable`].
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find the
#' `chrom_list` object used to create the `peak_table`.
#' @param lambdas The wavelength you wish to plot the traces at.
#' @param var Variable to index chromatograms.
#' @param subset Character vector specifying levels to use (if more than two
#' levels are present in `var`).
#' @param print_legend Logical. Whether to print legend. Defaults to `TRUE`.
#' @param legend_txt Character vector containing labels for legend.
#' @param legend_pos Legend position.
#' @param legend_size Legend size (`cex` argument). Default is `1`.
#' @param mirror Logical. Whether to plot as mirror or stacked plots.
#' Defaults to `TRUE`.
#' @param xlim Numerical vector specifying limits for x axis.
#' @param ylim Numerical vector specifying limits for y axis.
#' @param ... Additional arguments to `matplot` function.
#' @return No return value, called for side effects.
#' @section Side effects:
#' Plots are rendered to the active graphics device.
#' If `mirror_plot` is `TRUE`, plots a mirror plot comparing two treatments
#' defined by `var` and `subset` (if more than two factors are present
#' in `var`).
#' Otherwise, if `mirror_plot` is `FALSE`, the treatments are plotted in two
#' separate panes.
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
  if (!exists("var"))
    stop("Must provide independent variable or variables for mirror plot.")
  if (!is.data.frame(meta)){
    stop("Sample metadata must be attached to peak table to make mirror plot.")
  }
  if (!(var %in% colnames(meta))){
    stop(paste0("`", var, "`", " could not be found in sample metadata."))
  }
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(x)
  } else {
    get_chrom_list(x, chrom_list)
  }
  new.ts <- round(as.numeric(rownames(chrom_list[[1]])), 2)
  if (ncol(chrom_list[[1]]) == 1){
    lambdas.idx <- 1
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
    if (nlevels(fac) > 2){
      warning(paste0("Selecting first two levels of ", sQuote(var)),
      ". To choose different levels, use the `subset` argument 
        to specify the desired levels.")
    }
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
