
test_that("plot_all_spectra works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  plot_spectra_scaled <- function(){
    plot_all_spectra("V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export = TRUE, overlapping = TRUE, scale = TRUE)
  }
  vdiffr::expect_doppelganger("plot_spectra_scaled", plot_spectra_scaled)
  
  plot_spectra_raw <- function(){
    plot_all_spectra("V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export = TRUE, overlapping = TRUE, scale = FALSE)
  }
  vdiffr::expect_doppelganger("plot_spectra_raw", plot_spectra_raw)
  
  x <- plot_all_spectra("V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                        export = TRUE, overlapping = TRUE, scale = TRUE)
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(colnames(x), rownames(pk_tab$tab))
})

test_that("plot_all_spectra works with ggplot2", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  plot_spectra <- function(){
    plot_all_spectra("V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export = FALSE, overlapping = TRUE, engine = "ggplot")
  }
  x <- plot_all_spectra("V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                        export=FALSE, overlapping=TRUE, engine = "ggplot")
  vdiffr::expect_doppelganger("plot_all_spectra_ggplot", plot_spectra)
})

test_that("plot_spectrum works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  plot_spec <- function(){
    par(mfrow=c(2,1))
    plot_spectrum("13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                  export=TRUE, what="rt", idx = 1, verbose = FALSE)
  }
  vdiffr::expect_doppelganger(title = "plot_spectrum", plot_spec)
  
  x <- plot_spectrum(loc = "V15", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
})

test_that("plot_spectrum works", {
  data(pk_tab); data(Sa_warp)
  x <- plot_spectrum(loc = "V14", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what = "click"), 
               regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="click", idx = 1), 
               regexp = "must be specified")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt", lambda="210"),
               regexp = "Please supply argument")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt", idx = 1), regexp = "Please supply argument")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt"), regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt"), regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(loc = 12, chrom_list = pk_tab, what = "rt", idx = 1))
  expect_error(plot_spectrum(loc = 12, what="rt", idx = 1), 
               regexp = "Must provide either a peak_table or a chrom_list")
  expect_error(plot_spectrum(loc = 12, chrom_list = Sa_warp, what="peak", idx = 1),
               regexp = "Peak table must be provided")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="peak", idx = 1), regexp = "No match found for peak")
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", engine = "plotly"),
               regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", idx = 1, engine = "plotly"),
               regexp = "does not currently support clicking")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", idx = 1, lambda = "210",
                             engine = "plotly"),
               regexp = "does not currently support clicking")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click",lambda = "210", engine = "plotly"),
               regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", lambda="210", engine="plotly"),
               regexp = "Please supply argument")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", idx = 1, engine = "plotly"),
               regexp = "Please supply argument")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", engine = "plotly"), 
               regexp = "Chromatogram must be specified")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", engine = "plotly"), 
               regexp = "Chromatogram must be specified")
})


test_that("plot_spectrum works with ggplot2", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  data(pk_tab); data(Sa_warp)
  
  p1 <- plot_spectrum(loc = "13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="ggplot")
  vdiffr::expect_doppelganger(title = "plot_both_ggplot", fig = p1)
  
  p2 <- plot_spectrum(loc = "13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="ggplot", plot_trace = FALSE)
  vdiffr::expect_doppelganger(title = "plot_spectrum_ggplot", fig = p2)
  
  p3 <- plot_spectrum(loc = "13.75", peak_table = pk_tab, chrom_list = Sa_warp, 
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine = "ggplot", plot_spectrum = FALSE)
  vdiffr::expect_doppelganger(title = "plot_trace_ggplot", fig = p3)
})

test_that("plot_spectrum_inset works correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  
  data(Sa_warp)
  data(pk_tab)
  plot_inset <- function(){
    plot_spectrum_inset(c(Unk1 = "V9", Unk2 = "V15", Unk3 = "V26"), 
                        peak_table = pk_tab, chrom_list = Sa_warp,
                        position=list(c(0.35, 0.7), c(0.45, 0.4),
                                      c(0.72,0.4)), 
                        inset_width = 0.25, inset_height = 0.25)
  }
  vdiffr::expect_doppelganger("plot_inset",
                              plot_inset)
})

test_that("annotate_peaks works correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  
  data(Sa_warp)
  data(pk_tab)
  annotated_chrom <- function(){
    p <- plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot")
    print(annotate_peaks(p, loc = c(A = "V9", B = "V15", C = "V26"), 
                   peak_table = pk_tab, vjust=-1))
  }
  vdiffr::expect_doppelganger("annotated_chrom",
                              annotated_chrom)
  p <- plot_chroms(Sa_warp, lambdas = 210)
  expect_error(annotate_peaks(p, "V9", peak_table = pk_tab), 
               regexp = c("Please supply a ggplot"))
})

test_that("plot_chroms_heatmap works to plot alignments with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  data(Sa_warp)
  data(pk_tab)
  alignment_ggplot_heatmap <- function(){
    plot_chroms_heatmap(Sa_warp, lambdas = 210, engine = "ggplot")
  }
  vdiffr::expect_doppelganger("alignment_ggplot_heatmap",
                              alignment_ggplot_heatmap)
  
  alignment_ggplot_heatmap_zoom <- function(){
    plot_chroms_heatmap(Sa_warp, lambdas = 210, engine = "ggplot", 
                        show_legend = TRUE, xlim = c(10, 15))
  }
  suppressWarnings(vdiffr::expect_doppelganger("alignment_ggplot_heatmap_zoom",
                                               alignment_ggplot_heatmap_zoom))
  expect_error(plot_chroms_heatmap(Sa_warp), 
               regexp = 'argument "lambdas" is missing')
  expect_error(plot_chroms_heatmap(pk_tab$tab, lambdas = 210), 
               regexp = "Please supply a list of chromatograms.")
})

test_that("plot_chroms can subset chromatograms correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  numeric_indices <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c(1, 3), lambdas = 210)
  }
  vdiffr::expect_doppelganger("heatmap_plot_trace", numeric_indices)
  
  numeric_indices_rev <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c(3, 1), lambdas = 210)
  }
  vdiffr::expect_doppelganger("heatmap_plot_trace_rev", numeric_indices_rev)
})

# run the following line to activate plotly tests:
# Sys.setenv("VISUAL_TESTS" = "true")

test_that("plot_chroms_heatmap works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  data(Sa_warp)
  
  p <- plot_chroms_heatmap(Sa_warp, lambdas = "210", engine = "plotly")
  expect_doppelganger_plotly(name = "heatmap_plotly", 
                             p = p)
  # plot retention times instead of indices
  p2 <- plot_chroms_heatmap(Sa_warp, lambdas = "210", engine = "plotly", 
                            show_legend = FALSE, xlim = c(10, 15))
  expect_doppelganger_plotly(name = "heatmap_plotly_zoom", p = p2)
})

test_that("plot_chroms works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  data(Sa_warp)
  
  p <- plot_chroms(Sa_warp, lambdas = "210", engine = "plotly", show_legend = TRUE)
  expect_doppelganger_plotly(name = "alignment_plotly", p = p)
  
  p2 <- plot_chroms(Sa_warp, lambdas="210", engine = "plotly", show_legend = FALSE)
  expect_doppelganger_plotly(name = "alignment_plotly_no_legend", p = p2)
  
  p3 <- plot_chroms(Sa_warp, lambdas="210", engine = "plotly", show_legend = TRUE,
                    legend_position = "topleft", 
                    xlim=c(15, 18), ylim=c(0,400))
  expect_doppelganger_plotly(name = "alignment_plotly_zoom", p = p3)
})


test_that("plot_spectrum works with plotly engine", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  data(pk_tab); data(Sa_warp)
  
  p1 <- plot_spectrum("13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                      export_spectrum = FALSE, what = "rt", idx = 1,
                      verbose = FALSE, engine = "plotly")
  expect_doppelganger_plotly("plot_both_plotly", p = p1)
  
  p2 <- plot_spectrum("13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                      export_spectrum = FALSE, what = "rt", idx = 1,
                      verbose = FALSE, engine = "plotly", plot_trace = FALSE)
  expect_doppelganger_plotly("plot_trace_plotly", p = p2)
  
  p3 <- plot_spectrum("13.75", peak_table = pk_tab, chrom_list = Sa_warp,
                      export_spectrum = FALSE, what = "rt", idx = 1,
                      verbose = FALSE, engine = "plotly", plot_spectrum = FALSE)
  expect_doppelganger_plotly("plot_spectrum_plotly", p = p3)
})

test_that("boxplot works as expected", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  
  data("pk_tab"); data("Sa_warp")
  path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
  meta <- read.csv(path)
  pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")
  
  boxplot1 <- function(){
    boxplot(pk_tab, V11 ~ trt)
  }
  vdiffr::expect_doppelganger("boxplot1", boxplot1)
  
  boxplot2 <- function(){
    boxplot(pk_tab, V11~trt, las=2)
  }
  vdiffr::expect_doppelganger("boxplot2", boxplot2)
})

test_that("mirror_plot works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  
  data("pk_tab"); data("Sa_warp")
  path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
  meta <- read.csv(path)
  pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")
  
  mirror1 <- function(){
    mirror_plot(pk_tab, chrom_list = Sa_warp, lambdas = c("210","260"),
                var = "trt", legend_size = 2)
  }
  vdiffr::expect_doppelganger("mirror1", mirror1)
  
  mirror2 <- function(){
    mirror_plot(pk_tab, chrom_list = Sa_warp, lambdas = c("210","260"),
                var = "trt", legend_size = 2, mirror = FALSE)
  }
  vdiffr::expect_doppelganger("mirror2", mirror2)
  expect_error(mirror_plot(pk_tab, chrom_list = Sa_warp), regexp = 'argument "var" is missing')
  expect_error(mirror_plot(pk_tab, chrom_list = Sa_warp, var = "invalid_variable"), 
               regexp = "could not be found")
})

test_that("plot.peak_table works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
  meta <- read.csv(path)
  pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")
  
  plot_peak_table <- function(){
    par(mfrow=c(3, 1))
    plot(pk_tab, loc = "V9", chrom_list = Sa_warp, box_plot = TRUE,
         vars = "trt", verbose = FALSE, spectrum_labels = TRUE)
  }
  vdiffr::expect_doppelganger("plot.peak_table", plot_peak_table)
  expect_error(plot(pk_tab, chrom_list = Sa_warp, loc = "V13", idx = 1,
                    lambda = "210",  box_plot = TRUE, verbose = FALSE),
               regexp="Must provide independent variable")
  expect_error(plot(pk_tab, chrom_list = Sa_warp, loc = "15", what = "rt",
                    idx = 1, lambda = "210",  
                    box_plot = TRUE, var = "trt", verbose = FALSE),
               regexp = "A peak name must be provided")
})