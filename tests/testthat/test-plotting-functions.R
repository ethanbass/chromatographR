
test_that("plot_all_spectra works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  plot_spectra_scaled <- function(){
    plot_all_spectra("V11", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export = TRUE, overlapping = TRUE, scale = TRUE)
  }
  vdiffr::expect_doppelganger("plot_spectra_scaled", plot_spectra_scaled)
  
  plot_spectra_raw <- function(){
    plot_all_spectra("V11", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export = TRUE, overlapping = TRUE, scale = FALSE)
  }
  vdiffr::expect_doppelganger("plot_spectra_raw", plot_spectra_raw)
  
  x <- plot_all_spectra("V11", peak_table = pk_tab, chrom_list = Sa_warp, 
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
    print(plot_all_spectra("V11", peak_table = pk_tab, chrom_list = Sa_warp,
                           export_spectrum = FALSE, overlapping = TRUE,
                           engine = "ggplot2"))
  }
  vdiffr::expect_doppelganger("plot_all_spectra_ggplot", plot_spectra)
})

test_that("plot_all_spectra ggplot2 has correct wavelengths on x-axis", {
  skip_if_not_installed("ggplot2")
  data(pk_tab); data(Sa_warp)
  p <- plot_all_spectra("V11", peak_table = pk_tab, chrom_list = Sa_warp,
                        export_spectrum = FALSE, engine = "ggplot2")
  expect_s3_class(p, "gg")
  expect_equal(sort(unique(p$data$wavelength)), sort(get_lambdas(Sa_warp)))
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
  
  x <- plot_spectrum(loc = "V13", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
})

test_that("plot_spectrum works", {
  data(pk_tab); data(Sa_warp)
  x <- plot_spectrum(loc = "V12", peak_table = pk_tab, chrom_list = Sa_warp, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what = "click"), 
               regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="click", idx = 1), 
               regexp = "specify a wavelength")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt", lambda="210"),
               regexp = "Please supply argument to `loc`")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt", idx = 1), regexp = "Please supply argument to `loc`")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt"), regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="rt"), regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(loc = 12, chrom_list = pk_tab, what = "rt", idx = 1),
               regexp = "must be a list of chromatograms or a single chromatogram matrix, not a `peak_table`")
  expect_error(attach_ref_spectra(pk_tab, chrom_list = pk_tab),
               regexp = "not a `peak_table`")
  expect_error(plot_spectrum(loc = 12, what="rt", idx = 1), 
               regexp = "Must provide either a peak_table or a chrom_list")
  expect_error(plot_spectrum(loc = 12, chrom_list = Sa_warp, what="peak", idx = 1),
               regexp = "Peak table must be provided")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp, 
                             what="peak", idx = 1), regexp = "No match found for peak")
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", engine = "plotly"),
               regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", idx = 1, engine = "plotly"),
               regexp = "does not currently support clicking")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click", idx = 1, lambda = "210",
                             engine = "plotly"),
               regexp = "does not currently support clicking")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "click",lambda = "210", engine = "plotly"),
               regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", lambda="210", engine="plotly"),
               regexp = "Please supply argument to `loc`")
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", idx = 1, engine = "plotly"),
               regexp = "Please supply argument to `loc`")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", engine = "plotly"), 
               regexp = "specify a chromatogram with `idx`")
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = Sa_warp,
                             what = "rt", engine = "plotly"), 
               regexp = "specify a chromatogram with `idx`")
})


test_that("get_spectra returns correct structure (wide format)", {
  data(pk_tab); data(Sa_warp)
  x <- get_spectra("V12", peak_table = pk_tab, chrom_list = Sa_warp, what = "peak")
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
  expect_equal(ncol(x), length(Sa_warp))
  expect_equal(colnames(x), names(Sa_warp))
})

test_that("get_spectra works with what = 'rt'", {
  data(Sa_warp)
  x <- get_spectra("13.75", chrom_list = Sa_warp, what = "rt")
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
})

test_that("get_spectra works with what = 'idx'", {
  data(Sa_warp)
  x <- get_spectra(100, chrom_list = Sa_warp, what = "idx")
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(get_lambdas(Sa_warp)))
})

test_that("get_spectra format = 'long' returns correct columns", {
  data(pk_tab); data(Sa_warp)
  x <- get_spectra("V12", peak_table = pk_tab, chrom_list = Sa_warp,
                   what = "peak", format = "long")
  expect_equal(class(x), "data.frame")
  expect_equal(colnames(x), c("wavelength", "absorbance", "sample", "rt"))
  expect_equal(nrow(x), length(get_lambdas(Sa_warp)) * length(Sa_warp))
})

test_that("get_spectra subset of idx returns correct number of columns", {
  data(Sa_warp)
  x <- get_spectra("13.75", chrom_list = Sa_warp, what = "rt", idx = c(1, 2))
  expect_equal(ncol(x), 2)
  expect_equal(colnames(x), names(Sa_warp)[c(1, 2)])
})

test_that("get_spectra scale_spectrum normalizes to [0, 1]", {
  data(Sa_warp)
  x <- get_spectra("13.75", chrom_list = Sa_warp, what = "rt",
                   scale_spectrum = TRUE)
  expect_true(all(x >= 0 & x <= 1))
})

test_that("get_spectra errors when neither peak_table nor chrom_list provided", {
  expect_error(get_spectra("13.75", what = "rt"),
               regexp = "Must provide either")
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

test_that("plot_spectrum multi-sample works with base engine", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(pk_tab); data(Sa_warp)
  plot_multi_base <- function(){
    par(mfrow = c(2, 1))
    plot_spectrum(loc = "13.75", chrom_list = Sa_warp, what = "rt",
                 idx = c(1, 2), lambda = "210", verbose = FALSE)
  }
  vdiffr::expect_doppelganger("plot_spectrum_multi_base", plot_multi_base)
})

test_that("plot_spectrum multi-sample works with ggplot2 engine", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  data(pk_tab); data(Sa_warp)
  p <- plot_spectrum(loc = "13.75", chrom_list = Sa_warp, what = "rt",
                    idx = c(1, 2), lambda = "210", verbose = FALSE,
                    engine = "ggplot2")
  vdiffr::expect_doppelganger("plot_spectrum_multi_ggplot", p)
})

test_that("plot_spectrum_inset works correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  
  data(Sa_warp)
  data(pk_tab)
  plot_inset <- function(){
    print(plot_spectrum_inset(c(Unk1 = "V7", Unk2 = "V13", Unk3 = "V25"), idx=c(1:3),
                              peak_table = pk_tab, chrom_list = Sa_warp,
                              scale_spectrum = FALSE,
                              position = list(c(0.35, 0.7), c(0.45, 0.4),
                                              c(0.72, 0.4)),
                              inset_width = 0.25, inset_height = 0.25,
                              verbose = FALSE))
  }
  vdiffr::expect_doppelganger("plot_inset", plot_inset)
})

test_that("annotate_peaks works correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  
  data(Sa_warp)
  data(pk_tab)
  annotated_chrom <- function(){
    p <- plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot")
    print(annotate_peaks(p, loc = c(A = "V7", B = "V13", C = "V25"), 
                   peak_table = pk_tab, vjust=-1))
  }
  vdiffr::expect_doppelganger("annotated_chrom",
                              annotated_chrom)
  p <- plot_chroms(Sa_warp, lambdas = 210)
  expect_error(annotate_peaks(p, "V7", peak_table = pk_tab), 
               regexp = c("Please supply a ggplot"))
})

test_that("plot_chroms plots an element of nested chromatograms", {
  skip_if_not_installed("ggplot2")
  data(Sa_warp)
  nested <- lapply(Sa_warp, function(m){
    list(dad = m, instrument = list(pressure = m[, "210", drop = FALSE]))
  })
  flat <- plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot")
  expect_equal(plot_chroms(nested, lambdas = 210, engine = "ggplot")$data,
               flat$data)
  expect_equal(plot_chroms(nested, engine = "ggplot",
                           what = c("instrument", "pressure"))$data$absorbance,
               flat$data$absorbance)
  expect_equal(plot_chroms(nested, engine = "ggplot",
                           what = "instrument: pressure")$data$absorbance,
               flat$data$absorbance)
  expect_error(plot_chroms(nested, what = "instrument:flow"),
               "`what = \"instrument:flow\"` was not found")
  colon_named <- lapply(nested, function(el) c(el, list("dad:raw" = el$dad)))
  expect_equal(plot_chroms(colon_named, lambdas = 210, engine = "ggplot",
                           what = "dad:raw")$data, flat$data)
  expect_error(plot_chroms(nested, lambdas = 210, what = "uv"),
               "not found. Elements available: 'dad', 'instrument'")
  expect_error(plot_chroms(nested, what = "instrument"),
               "is a list, not a matrix.*c\\(\"instrument\", \"pressure\"\\)")
  expect_error(plot_chroms(nested, what = c("instrument", "flow")),
               "Elements available in .*: 'pressure'")
})

test_that("get_element extracts nested elements as a list of matrices", {
  data(Sa_pr)
  nested <- lapply(Sa_pr, function(m){
    list(dad = m, instrument = list(pressure = m[, "210", drop = FALSE]))
  })
  class(nested) <- class(Sa_pr)
  expect_equal(get_element(nested, "dad"), Sa_pr)
  p <- get_element(nested, "instrument:pressure")
  expect_equal(names(p), names(Sa_pr))
  expect_equal(p[[2]], Sa_pr[[2]][, "210", drop = FALSE])
})

test_that("annotate_peaks labels retention times without a peak table", {
  skip_if_not_installed("ggplot2")
  data(Sa_warp)
  p <- plot_chroms(Sa_warp, lambdas = 210, engine = "ggplot")
  top <- function(rt, tol = 0.05){
    k <- abs(p$data$rt - rt) <= tol
    max(p$data$absorbance[k], na.rm = TRUE)
  }
  a <- annotate_peaks(p, c(C1 = 12.18, C2 = 12.70))
  n <- length(a$layers)
  expect_equal(a$layers[[n - 1]]$data$x, 12.18)
  expect_equal(a$layers[[n - 1]]$data$y, top(12.18))
  expect_equal(a$layers[[n]]$data$y, top(12.70))
  expect_equal(a$layers[[n]]$aes_params$label, "C2")
  b <- annotate_peaks(p, 12.18, tol = 0.2)
  expect_equal(b$layers[[length(b$layers)]]$data$y, top(12.18, 0.2))
  expect_equal(b$layers[[length(b$layers)]]$aes_params$label, "12.18")
  expect_error(annotate_peaks(p, 100), "No plotted data within `tol`")
  expect_error(annotate_peaks(p, "V7"), "Please supply a `peak_table`")
})

test_that("plot_chroms_heatmap works to plot alignments with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  data(Sa_warp)
  data(pk_tab)
  alignment_ggplot_heatmap <- function(){
    print(plot_chroms_heatmap(Sa_warp, lambdas = 210, engine = "ggplot"))
  }
  vdiffr::expect_doppelganger("alignment_ggplot_heatmap",
                              alignment_ggplot_heatmap)

  alignment_ggplot_heatmap_zoom <- function(){
    print(plot_chroms_heatmap(Sa_warp, lambdas = 210, engine = "ggplot",
                              show_legend = TRUE, xlim = c(10, 15)))
  }
  suppressWarnings(vdiffr::expect_doppelganger("alignment_ggplot_heatmap_zoom",
                                               alignment_ggplot_heatmap_zoom))
  expect_error(plot_chroms_heatmap(Sa_warp), 
               regexp = "`lambdas` must be supplied")
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

test_that("plot_chroms handles unrecognized character idx", {
  data(Sa_warp)
  expect_error(plot_chroms(Sa_warp, lambdas = "210", idx = "bad_name"),
               regexp = "could not be found")
  expect_warning(plot_chroms(Sa_warp, lambdas = "210",
                             idx = c(names(Sa_warp)[1], "bad_name")),
                 regexp = "bad_name")
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

test_that("a top or bottom base legend gets its own band above or below the traces", {
  pdf(NULL); on.exit(dev.off())
  labels <- c("119", "121", "122", "458")

  top <- legend_band(labels, "topright", c(0, 100))
  expect_equal(top$ylim[1], 0)
  expect_gt(top$ylim[2], 100)
  expect_equal(top$ncol, 4)

  bottom <- legend_band(labels, "bottom", c(0, 100))
  expect_lt(bottom$ylim[1], 0)
  expect_equal(bottom$ylim[2], 100)

  expect_equal(legend_band(labels, "right", c(0, 100)),
               list(ncol = 1L, ylim = c(0, 100), cex = 1))

  data(Sa_pr)
  plot_chroms(Sa_pr, lambdas = "210", show_legend = TRUE)
  r <- legend("topright", legend = names(Sa_pr), fill = 1:4, ncol = 4,
              plot = FALSE)$rect
  expect_gt(r$top - r$h, max(sapply(Sa_pr, function(m) max(m[, "210"]))))
})

test_that("a legend too tall for a small panel shrinks to fit its band", {
  pdf(NULL, width = 5, height = 5); on.exit(dev.off())
  par(mfrow = c(2, 1)); plot.new()
  labels <- c("460", "464", "467", "228", "224", "230", "406", "450", "385", "355")
  band <- legend_band(labels, "topright", c(0, 100))
  expect_lt(band$cex, 1)
  plot.window(xlim = c(0, 1), ylim = band$ylim)
  r <- legend("topright", legend = labels, fill = seq_along(labels), ncol = band$ncol,
              cex = band$cex, plot = FALSE)$rect
  expect_gt(r$top - r$h, 100)
})

test_that("plot_chroms fits the y axis to the traces within xlim", {
  data(Sa_pr)
  lam <- which(colnames(Sa_pr[[1]]) == "210")
  win <- c(12.8, 13.4)
  in_win <- max(sapply(Sa_pr, function(m){
    t <- as.numeric(rownames(m)); max(m[t >= win[1] & t <= win[2], lam])
  }))
  whole <- max(sapply(Sa_pr, function(m) max(m[, lam])))
  expect_lt(in_win, whole / 2)
  expect_equal(get_y_bounds(Sa_pr, seq_along(Sa_pr), lam, xlim = win)[2], in_win * 1.1,
               tolerance = 0.05)
  expect_equal(get_y_bounds(Sa_pr, seq_along(Sa_pr), lam, xlim = c(30, 40)),
               get_y_bounds(Sa_pr, seq_along(Sa_pr), lam))
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
    boxplot(pk_tab, V9 ~ trt)
  }
  vdiffr::expect_doppelganger("boxplot1", boxplot1)
  
  boxplot2 <- function(){
    boxplot(pk_tab, V9~trt, las=2)
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
    plot(pk_tab, loc = "V7", chrom_list = Sa_warp, box_plot = TRUE,
         vars = "trt", verbose = FALSE, spectrum_labels = TRUE)
  }
  vdiffr::expect_doppelganger("plot.peak_table", plot_peak_table)
  expect_error(plot(pk_tab, chrom_list = Sa_warp, loc = "V11", idx = 1,
                    lambda = "210",  box_plot = TRUE, verbose = FALSE),
               regexp="Must provide independent variable")
  expect_error(plot(pk_tab, chrom_list = Sa_warp, loc = "15", what = "rt",
                    idx = 1, lambda = "210",  
                    box_plot = TRUE, var = "trt", verbose = FALSE),
               regexp = "A peak name must be provided")
})

test_that("purity plot works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  pks <- get_peaks(Sa_pr[1], lambdas = 210, fit = "egh", estimate_purity = TRUE,
                   cl = 1, show_progress = FALSE, smooth_window = 2)
  plot_purity <- function(display){
    function() plot(pks, chrom_list = Sa_pr, points = TRUE, ticks = TRUE,
                    lambda = 210, plot_purity = display)
  }
  vdiffr::expect_doppelganger("plot_purity", plot_purity(TRUE))
  vdiffr::expect_doppelganger("plot_purity_ratio", plot_purity("ratio"))
  vdiffr::expect_doppelganger("plot_purity_curves_local", plot_purity("curves_local"))
})

test_that("plot.peak_list rejects a peak table passed as `chrom_list`", {
  data(Sa_pr); data(pk_tab)
  pks <- get_peaks(Sa_pr[1], lambdas = 210, show_progress = FALSE)
  expect_error(plot(pks, chrom_list = pk_tab), regexp = "not a `peak_table`")
})

test_that("plot_chroms checks a supplied wavelength on single-column data", {
  data(Sa_pr)
  x <- lapply(Sa_pr[1:2], function(m) m[, "210", drop = FALSE])
  pdf(NULL); on.exit(dev.off())
  expect_error(plot_chroms(x, lambdas = "280"), "could not be found")
  expect_error(plot_chroms_heatmap(x, lambdas = "280"), "could not be found")
  expect_no_error(plot_chroms(x, lambdas = "210"))
  expect_no_error(plot_chroms(x))
  bl <- correct_baseline(x[[1]], attach_baseline = TRUE)
  expect_error(plot_baseline(bl, lambda = "280"), "could not be found")
  expect_no_error(plot_baseline(bl, lambda = "210"))
  colnames(x[[1]]) <- colnames(x[[2]]) <- "intensity"
  expect_no_error(plot_chroms(x, lambdas = "280"))
})

test_that("plot_chroms does not thin below the data's own resolution", {
  skip_on_cran()
  data(Sa_pr)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_chroms(Sa_pr, lambdas = 210))
  expect_silent(plot_chroms(Sa_pr, lambdas = 210, time_resolution = 0.005))
  expect_silent(plot_chroms(Sa_pr, lambdas = 210, time_resolution = 0.001))
  coarse <- lapply(Sa_pr, function(x){
    y <- x[seq(1, nrow(x), by = 3), , drop = FALSE]
    rownames(y) <- as.character(seq(10, by = 0.06, length.out = nrow(y)))
    y
  })
  expect_silent(plot_chroms(coarse, lambdas = 210))
  tiny <- lapply(Sa_pr, function(x) x[1:5, , drop = FALSE])
  expect_silent(plot_chroms(tiny, lambdas = 210))
  expect_silent(plot_chroms(Sa_pr, lambdas = 210, time_resolution = 0.5))
})

test_that("plot_chroms works across dilation_span values for vpdtw output", {
  skip_on_cran()
  data(Sa_pr)
  pdf(NULL)
  on.exit(dev.off())
  for (sp in c(5, 2, 0.6, 0.3, 0.12)){
    warp <- suppressWarnings(
      correct_rt(Sa_pr, lambdas = 210, alg = "vpdtw", dilation_span = sp,
                 show_progress = FALSE))
    expect_silent(plot_chroms(warp, lambdas = 210))
  }
})

test_that("reshape_chroms defaults to every wavelength", {
  skip_on_cran()
  data(Sa_pr)
  all_lambdas <- reshape_chroms(Sa_pr)
  expect_equal(nrow(all_lambdas),
               nrow(Sa_pr[[1]]) * ncol(Sa_pr[[1]]) * length(Sa_pr))
  expect_equal(length(unique(all_lambdas$lambda)), ncol(Sa_pr[[1]]))
  expect_equal(reshape_chroms(Sa_pr, lambdas = NULL), all_lambdas)
  one <- reshape_chroms(Sa_pr, lambdas = 210)
  expect_equal(nrow(one), nrow(Sa_pr[[1]]) * length(Sa_pr))
  coarse <- lapply(Sa_pr, function(x){
    y <- x[seq(1, nrow(x), by = 3), , drop = FALSE]
    rownames(y) <- as.character(seq(10, by = 0.06, length.out = nrow(y)))
    y
  })
  expect_equal(nrow(reshape_chroms(coarse, lambdas = 210,
                                   time_resolution = 0.01)),
               nrow(coarse[[1]]) * length(coarse))
})

test_that("plot_chroms works on one-dimensional chromatograms in every engine", {
  skip_on_cran()
  data(Sa_warp)
  x1 <- structure(lapply(Sa_warp, function(m) m[, "210", drop = FALSE]),
                  names = names(Sa_warp))
  expect_no_error(plot_chroms(x1))
  expect_no_error(plot_chroms(x1, engine = "ggplot"))
  check_for_pkg("plotly")
  expect_no_error(plot_chroms(x1, engine = "plotly"))
})

test_that("plot_chroms asks for lambdas on multi-wavelength data", {
  data(Sa_warp)
  expect_error(plot_chroms(Sa_warp, engine = "ggplot"),
               regexp = "`lambdas` must be supplied")
})

test_that("purity_args accepted by get_peaks also work at plot time", {
  # `get_purity_values`' formals are a strict subset of `get_purity`'s
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr[1], lambdas = 210, cl = 1, show_progress = FALSE)
  pdf(NULL)
  on.exit(dev.off())
  for (pa in list(list(metric = "max"), list(min_run = 5), list(try = TRUE))){
    expect_no_error(suppressWarnings(
      plot(pks, chrom_list = Sa_pr, idx = 1, lambda = "210",
           plot_purity = TRUE, purity_args = pa)))
  }
})

test_that("plot_baseline draws the attached baseline", {
  skip_on_cran()
  data(Sa)
  args <- list(dim1 = seq(10, 12, by = .01), dim2 = seq(200, 318, by = 12),
               show_progress = FALSE, cl = 1)
  x <- do.call(preprocess, c(list(Sa[1]), args, baseline = "snip", m = 10,
                             attach_baseline = TRUE))
  expect_silent(plot_baseline(x, lambda = 212))
  expect_silent(plot_baseline(x[[1]], lambda = 212))
  expect_error(plot_baseline(x), "must be supplied")
  bare <- do.call(preprocess, c(list(Sa[1]), args))
  expect_error(plot_baseline(bare, lambda = 212), "No baseline found")
})

test_that("plot_baseline annotates the smoothing parameter it was given", {
  expect_equal(baseline_param_label(c(m = 9)), "m = 9")
  expect_equal(baseline_param_label(structure(c(m = 9), method = "snip")),
               "snip, m = 9")
  expect_equal(baseline_param_label(
    structure(c(stiffness = 18738.87), method = "arpls")),
    "arpls, stiffness = 18739")
  expect_null(baseline_param_label(NULL))
})

test_that("plot_baseline renders", {
  skip_if_not_installed("vdiffr")
  tt <- seq(0, 30, by = 0.05)
  bl <- 50 + 30 * exp(tt / 12)
  y <- bl + 1000 * dnorm(tt, 10, 0.4) + 60 * dnorm(tt, 20, 0.4)
  x <- matrix(y - bl, ncol = 1, dimnames = list(tt, "210"))
  attr(x, "estimated_baseline") <- matrix(bl, ncol = 1, dimnames = dimnames(x))
  attr(x, "baseline_param") <- structure(c(stiffness = 18738.87),
                                         method = "arpls")
  vdiffr::expect_doppelganger("plot-baseline", function() plot_baseline(x))
})

test_that("plot_baseline yields its defaults to `...` and returns its label", {
  tt <- seq(0, 30, by = 0.01)
  bl <- 50 + 30 * exp(tt / 12)
  x <- matrix(1000 * dnorm(tt, 10, 0.4), ncol = 1, dimnames = list(tt, "210"))
  attr(x, "estimated_baseline") <- matrix(bl, ncol = 1, dimnames = dimnames(x))
  attr(x, "baseline_param") <- structure(c(stiffness = 18738.87),
                                         method = "arpls", dt = 0.01)
  pdf(NULL)
  on.exit(dev.off())
  for (arg in list(list(main = "mine"), list(xlab = "t"), list(ylab = "A"),
                   list(ylim = c(0, 60)), list(col = "blue")))
    expect_no_error(do.call(plot_baseline, c(list(x), arg)))
  expect_equal(plot_baseline(x), "arpls, stiffness = 18739, dt = 0.01")
  expect_equal(plot_baseline(x, add = TRUE), plot_baseline(x))
})

test_that("plot_chroms draws one line per sample and wavelength", {
  skip_if_not_installed("ggplot2")
  data(Sa_pr)
  p <- plot_chroms(Sa_pr, idx = 1:2, lambdas = c(210, 260), engine = "ggplot")
  expect_equal(length(unique(ggplot2::layer_data(p)$group)), 4)
})

test_that("plot_chroms plotly splits traces by wavelength", {
  skip_if_not_installed("plotly")
  data(Sa_pr)
  p <- plotly::plotly_build(plot_chroms(Sa_pr, idx = 1:2, lambdas = c(210, 260),
                                        engine = "plotly"))
  expect_length(p$x$data, 4)
})

test_that("plot_baseline accepts a vector with an attached baseline", {
  data(Sa_pr)
  x <- correct_baseline(Sa_pr[[1]][, "210"], attach_baseline = TRUE)
  pdf(NULL); on.exit(dev.off())
  expect_match(plot_baseline(x), "als")
})
test_that("annotate_peaks labels peaks on a single-wavelength plot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data(Sa_pr)
  one <- lapply(Sa_pr, function(m) m[, "210", drop = FALSE])
  pt <- get_peaktable(get_peaks(one, lambdas = 210, cl = 1, show_progress = FALSE))
  expect_s3_class(annotate_peaks(plot_chroms(one, engine = "ggplot"),
                                 colnames(pt$tab)[3], peak_table = pt,
                                 chrom_list = one), "ggplot")
})

test_that("annotate_peaks reads `tol` in minutes on a seconds axis", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data(Sa_pr)
  s <- lapply(Sa_pr, function(m){
    rownames(m) <- as.numeric(rownames(m)) * 60
    attr(m, "time_unit") <- "s"
    m
  })
  expect_s3_class(annotate_peaks(plot_chroms(s, lambdas = 210, engine = "ggplot"),
                                 c(A = 724.5)), "ggplot")
})

test_that("plot.peak_list zooms on a targets data.frame", {
  skip_on_cran()
  data(Sa_pr)
  pm <- get_peaks(Sa_pr, lambdas = 210, cl = 1, show_progress = FALSE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(pm, chrom_list = Sa_pr, targets = data.frame(rt = 12.07)))
})
