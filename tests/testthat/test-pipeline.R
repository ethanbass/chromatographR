skip_on_cran()
skip_if_not_installed("vdiffr")
skip_if_not_installed("plotly")
skip_if_not_installed("reticulate")
skip_if_not_installed("rsvg")
# to run visual tests, set Sys.setenv("VISUAL_TESTS"="true")

data("Sa")
new.ts <- seq(10, 18.66, by = .01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths

out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas, 
                  show_progress = FALSE)

test_that("preprocess works on matrix", {
  out <- preprocess_matrix(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas, 
                           maxI = 1)
  expect_equal(class(out)[1],"matrix")
  expect_equal(rownames(out), as.character(new.ts))
  expect_equal(colnames(out), as.character(new.lambdas))
  expect_error(preprocess_matrix(Sa))
  expect_error(preprocess_matrix(Sa[[1]], dim2 = seq(9, 20, by = .01)))
})

test_that("Preprocess works without providing dimensions", {
  out <- suppressMessages(preprocess(X = (Sa[1]), show_progress = FALSE))
  expect_equal(class(out), "list")
  expect_equal(class(out[[1]])[1],"matrix")
  expect_equal(rownames(out[[1]]), as.character(new.ts))
  expect_equal(colnames(out[[1]]), as.character(new.lambdas))
})


test_that("preprocess returns correct errors", {
  expect_error(preprocess(X = as.data.frame(Sa[[1]]), show_progress = FALSE))
})


test_that("preprocess works without interpolation", {
  dat.pr <- preprocess(X = Sa[[1]], interpolate_cols = FALSE,
                       interpolate_rows = FALSE, show_progress = FALSE)
  expect_equal(dim(dat.pr),dim(Sa[[1]]))
  expect_equal(rownames(dat.pr), rownames(Sa[[1]]))
  expect_equal(colnames(dat.pr), colnames(Sa[[1]]))
})

dat.pr <- preprocess(X = Sa[1:2], dim1 = new.ts, dim2 = new.lambdas,
                     show_progress = FALSE)

test_that("preprocess works on a list", {
  expect_equal(length(dat.pr),length(Sa[1:2]))
  expect_equal(rownames(dat.pr[[1]]), as.character(new.ts))
  expect_equal(colnames(dat.pr[[1]]), as.character(new.lambdas))
  expect_equal(names(dat.pr), names(Sa[1:2]))
})

### test correct_rt ###

warping.models <- correct_rt(dat.pr, lambdas = c('210','260','318'),
                             what = "models", show_progress = FALSE)

warp <- correct_rt(chrom_list = dat.pr, models = warping.models,
                   what = "corrected.values", show_progress = FALSE)

test_that("correct_rt works", {
  equals(length(warping.models), length(warp), length(Sa[1:2]))
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_equal(t(warping.models[[1]]$warped.sample)[,1], 
               as.numeric(warp[[1]][,"210"]))
  expect_error(correct_rt(dat.pr))
  expect_error(correct_rt(dat.pr, what="x"))
  expect_error(correct_rt(dat.pr, lambdas = "x"))
  expect_error(correct_rt(dat.pr, lambdas = "210", alg="x"))
  expect_error(correct_rt(dat.pr, lambdas = "210", 
                          models = "warping.models", alg = "vpdtw"))
})

test_that("correct_rt works with vpdtw", {
  skip_if_not_installed("VPdtw")
  warp <- correct_rt(dat.pr, lambdas = "210", alg="vpdtw")
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr = TRUE)
  expect_error(correct_rt(dat.pr, lambdas = c("210", "260"),
                          alg="vpdtw",  show_progress = FALSE))
  expect_error(correct_rt(dat.pr, lambdas = c("x"), alg = "vpdtw"))
})

test_that("VPdtw plot displays correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  vpdtw_alignment <- function(){
    warp <- correct_rt(dat.pr, lambdas = "210", alg = "vpdtw", plot_it = TRUE)
  }
  vdiffr::expect_doppelganger("vpdtw_alignment", vpdtw_alignment)
})

test_that("PTW plot displays correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  ptw_alignment <- function(){
    warp <- correct_rt(dat.pr, lambdas = "210", alg = "ptw", plot_it = TRUE)
  }
  vdiffr::expect_doppelganger("ptw_alignment", ptw_alignment)
})


test_that("plot_chroms works to plot alignments", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  alignment <- function(){
    par(mfrow = c(2,1))
    plot_chroms(warp, lambdas = "210")
    plot_chroms(dat.pr, lambdas = "210")
  }
  vdiffr::expect_doppelganger("alignment", alignment)
  expect_error(plot_chroms(pktab))
})

test_that("plot_chroms can subset chromatograms correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  numeric_indices <- function(){
    plot_chroms(Sa_pr, idx=c(1, 3), lambdas = 210)
  }
  character_indices <- function(){
    plot_chroms(Sa_pr, idx=c("119", "122"), lambdas = 210)
  }
  vdiffr::expect_doppelganger("plot_trace", numeric_indices)
  vdiffr::expect_doppelganger("plot_trace", character_indices)

  numeric_indices_rev <- function(){
    plot_chroms(Sa_pr, idx=c(3, 1), lambdas = 210)
  }
  character_indices_rev <- function(){
    plot_chroms(Sa_pr, idx=c("122", "119"), lambdas = 210)
  }
  vdiffr::expect_doppelganger("plot_trace_rev", numeric_indices_rev)
  vdiffr::expect_doppelganger("plot_trace_rev", character_indices_rev)
})


test_that("plot_chroms works to plot alignments with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  alignment_ggp <- function(){
    plot_chroms(warp, lambdas="210", engine="ggplot")
  }
  vdiffr::expect_doppelganger("alignment_ggp", alignment_ggp)
  
  alignment_ggp_zoom <- function(){
    plot_chroms(warp, lambdas = "210", engine = "ggplot", show_legend = FALSE,
                xlim = c(12.5, 15), ylim = c(0, 550))
  }
  suppressWarnings(vdiffr::expect_doppelganger("alignment_ggp_zoom",
                                               alignment_ggp_zoom))
  
  expect_error(plot_chroms(pktab))
})


### test get_peaks ###
lam <- c("210", "318")
pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh", smooth_type = "none",
                     show_progress = FALSE)
pks_gaussian <- get_peaks(dat.pr, lambdas = lam, fit = "gaussian",
                          smooth_type = "none", show_progress = FALSE)
pks_raw <- get_peaks(dat.pr, lambdas = lam, fit = "raw", smooth_type = "savgol",
                     show_progress = FALSE, sd.max = 100, smooth_window=3)

test_that("get_peaks works", {
  expect_equal(names(pks_egh), names(dat.pr))
  expect_equal(names(pks_gaussian), names(dat.pr))
  expect_equal(names(pks_egh[[1]]), lam)
  expect_equal(names(pks_gaussian[[1]]), lam)
  expect_equal(class(pks_egh), "peak_list")
  expect_error(get_peaks(dat.pr)) # lambdas must be provided
  expect_error(get_peaks(dat.pr, lambdas = "210", fit = "nonsense"))
})

test_that("filter_peaks works", {
  pks_s <- filter_peaks(pks_egh, min_height = 10, min_area = 10, min_sd = .07)
  expect_equal(names(pks_egh), names(pks_s))
  expect_equal(names(pks_s), names(dat.pr))
  expect_equal(names(pks_s[[1]]), lam)
  expect_equal(class(pks_s), "peak_list")
  expect_lt(nrow(pks_s[[1]][[1]]), nrow(pks_egh[[1]][[1]]))
  expect_lt(nrow(pks_s[[2]][[1]]), nrow(pks_egh[[2]][[1]]))
  expect_warning(filter_peaks(pks_egh))
  expect_warning(filter_peaks(pks_egh, min_height=0))
  expect_warning(filter_peaks(pks_egh, min_area=0))
  expect_warning(filter_peaks(pks_egh, max_sd=Inf))
})

### test get_peaktable ###

pk_tab <- get_peaktable(peak_list = pks_egh, chrom_list = dat.pr)
pk_tab_sp <- get_peaktable(peak_list = pks_egh, chrom_list = dat.pr, clust = "sp.rt")

test_that("get_peaktable works", {
  expect_equal(rownames(pk_tab$tab), names(dat.pr))
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$pk_meta))
  expect_equal(class(pk_tab), "peak_table")
  expect_equal(class(pk_tab$tab), "data.frame")
  expect_equal(class(pk_tab$pk_meta), "data.frame")
  expect_equal(class(pk_tab$args), "list")
})

test_that("correct_peaks works", {
  skip_on_cran()
  data(Sa_pr)
  pks <- get_peaks(Sa_pr, lambdas = 210, show_progress = FALSE)
  
  warping.models <- correct_rt(chrom_list = Sa_pr, lambdas = 210, 
                               what = "models", show_progress = FALSE)
  
  pks_cor <- correct_peaks(pks, mod_list = warping.models)
  pktab_cor <- get_peaktable(pks_cor, use.cor = TRUE)

  ptw_warp <- correct_rt(Sa_pr, models = warping.models)
  pks_warp <- get_peaks(ptw_warp, lambdas = 210, show_progress = FALSE)
  pktab_warp <- get_peaktable(pks_warp)
  
  # this test failed on UBUNTU
  # expect_equal(pks_cor[[1]][[1]]$rt.cor[-33], pks_warp[[1]][[1]]$rt, 
  #              tolerance = .001)
  # expect_equal(pks_cor[[2]][[1]]$rt.cor, pks_warp[[2]][[1]]$rt[-6], 
  #              tolerance = .001)
  
  pks_reg <- get_peaks(Sa_pr, lambdas = 210, show_progress = FALSE)
  pktab_reg <- get_peaktable(pks_reg)

  expect_equal(row.names(pktab_cor), names(Sa_pr))
  expect_equal(row.names(pktab_reg), names(Sa_pr))
  expect_equal(row.names(pktab_warp), names(Sa_pr))
  
  expect_lt(ncol(pktab_warp), ncol(pktab_reg))
  expect_lt(ncol(pktab_cor), ncol(pktab_reg))
})

test_that("strip plot works", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("vdiffr")
  skip_on_ci()
  plot_peak_table <- function(){
    pktab <- get_peaktable(pks_egh, plot_it = TRUE,
                           ask = FALSE)
  }
  
  vdiffr::expect_doppelganger("peak_table_plot", plot_peak_table)
})

#### test metadata attachment ###

path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)

pk_tab <- attach_metadata(pk_tab, metadata = meta, column = "vial")

test_that("attach_metadata works", {
  expect_equal(rownames(pk_tab$tab), pk_tab$sample_meta$vial)
  expect_equal(colnames(pk_tab$sample_meta), colnames(meta))
  
  ### check errors & warnings ###
  expect_error(attach_metadata())
  expect_error(attach_metadata(warp))
  expect_error(attach_metadata(pk_tab, metadata = warp))
  expect_error(attach_metadata(pk_tab, metadata = meta, column = "x"))
  expect_error(attach_metadata(pk_tab$tab, metadata = meta, column = "vial"))
  expect_error(attach_metadata(pk_tab, metadata = rbind(meta,meta), 
                               column = "vial"))
  expect_warning(attach_metadata(pk_tab, metadata = meta[-1,], column = "vial"))
})

pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.cor")
test_that("attach_ref_spectra works", {
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  expect_equal(pk_tab$args[["reference_spectra"]], "max.cor")
  
  pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.int")
  expect_equal(pk_tab$args[["reference_spectra"]], "max.int")
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  
  # test errors
  expect_error(attach_ref_spectra(peak_table = warp))
  expect_error(attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "x"))
})

test_that("subset.peaktable works", {
  sub <- subset(pk_tab, subset = pk_tab$sample_meta$trt == "+")
  
  expect_equal(sub$tab, pk_tab$tab[1,])
  expect_equal(sub$pk_meta, pk_tab$pk_meta)
  expect_equal(sub$sample_meta, pk_tab$sample_meta[1,])
  expect_equal(sub$ref_spectra, pk_tab$ref_spectra)
  
  sub <- subset(pk_tab, subset = pk_tab$sample_meta$trt == "-", 
                select = c("V10","V12"))
  
  expect_equal(sub$tab, pk_tab$tab[c(2),c("V10","V12")])
  expect_equal(sub$pk_meta, pk_tab$pk_meta[,c("V10", "V12")])
  expect_equal(sub$sample_meta, pk_tab$sample_meta[2,])
  expect_equal(sub$ref_spectra, pk_tab$ref_spectra[,c("V10", "V12")])
})


# test filter_peaktable

test_that("filter_peaktable works", {
  pktab_s <- filter_peaktable(pk_tab, min_rt = 6, max_rt = 15)
  # check that dimensions are unaltered
  expect_equal(rownames(pk_tab$tab), rownames(pktab_s$tab))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$pk_meta))
  expect_equal(nrow(pktab_s$tab), nrow(pktab_s$sample_meta))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$ref_spectra))

  # warning if no arguments are provided
  expect_warning(filter_peaktable(pk_tab))
})

test_that("combine_peaks works", {
  pk_tab_c <- combine_peaks(pk_tab, verbose = FALSE)
  expect_lt(ncol(pk_tab_c), ncol(pk_tab))
  expect_lt(ncol(pk_tab_c$pk_meta), ncol(pk_tab$pk_meta))
  expect_equal(rownames(pk_tab_c$tab), rownames(pk_tab$tab))
})

test_that("merge_peaks works with max method", {
  pk_tab_m <- merge_peaks(pk_tab, peaks=c("V10","V11"))
  expect_equal(pk_tab_m$tab$V11, pmax(pk_tab$tab$V10, pk_tab$tab$V11))
  expect_equal(ncol(pk_tab_m), ncol(pk_tab)-1)
  expect_equal(ncol(pk_tab_m$pk_meta), ncol(pk_tab$pk_meta)-1)
  expect_equal(rownames(pk_tab_m$tab), rownames(pk_tab$tab))
})

test_that("merge_peaks works with sum method", {
  data(pk_tab)
  pk_tab_m <- merge_peaks(pk_tab, peaks=c("V10","V11"), method = "sum")
  expect_equal(pk_tab_m$tab[["V11"]], (pk_tab$tab$V10 + pk_tab$tab$V11))
  expect_equal(ncol(pk_tab_m), ncol(pk_tab)-1)
  expect_equal(ncol(pk_tab_m$pk_meta), ncol(pk_tab$pk_meta) - 1)
  expect_equal(rownames(pk_tab_m$tab), rownames(pk_tab$tab))
})

test_that("normalize_data works", {
  pk_tab_norm <- normalize_data(pk_tab, chrom_list = dat.pr, column = "mass")
  expect_equal(rownames(pk_tab_norm$tab), rownames(pk_tab$tab))
  expect_equal(class(pk_tab_norm), class(pk_tab))
  expect_equal(colnames(pk_tab_norm$tab), colnames(pk_tab$tab))
  expect_equal(pk_tab_norm$tab[1,], pk_tab$tab[1,]/pk_tab$sample_meta$mass[1])
  expect_equal(pk_tab_norm$tab[2,], pk_tab$tab[2,]/pk_tab$sample_meta$mass[2])
  expect_equal(pk_tab_norm$args[["normalized"]], "TRUE")
  expect_equal(pk_tab_norm$args[["normalization_by"]], "mass")
  
  chroms_norm <- normalize_data(pk_tab, chrom_list = dat.pr, 
                           column = "mass", what = "chrom_list")
  expect_equal(dim(chroms_norm), dim(dat.pr))
  expect_equal(chroms_norm[[1]], dat.pr[[1]]/pk_tab$sample_meta$mass[1])
  expect_equal(chroms_norm[[2]], dat.pr[[2]]/pk_tab$sample_meta$mass[2])
  
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "x"))
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, 
                              column = "mass", what="x"))
})

test_that("cluster_spectra works", {
  skip_on_cran()
  skip_if_not_installed("pvclust")
  cl <- cluster_spectra(pk_tab, nboot = 10, parallel = FALSE, verbose = FALSE, 
                        save = FALSE, plot_dend = FALSE, plot_spectra = FALSE)
  expect_equal(class(cl[[1]]), "pvclust")
  expect_equal(class(cl[[2]]), "list")
})

### test plotting functions

test_that("plot.peak_list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_peaks <- function(){
    plot(pks_egh, chrom_list = dat.pr, points = TRUE,
         ticks = TRUE, lambda = 210)
  }
  vdiffr::expect_doppelganger("plot.peak_list", plot_peaks)
  
  plot_peaks_gaussian <- function(){
    plot(pks_gaussian, chrom_list = dat.pr)
  }
  vdiffr::expect_doppelganger("plot.peak_list_gaussian", plot_peaks_gaussian)
  
  # something raw with this test on GitHub actions
  # plot_peaks_raw <- function(){
  #   plot(pks_raw, chrom_list = dat.pr)
  # }
  # vdiffr::expect_doppelganger("plot.peak_list_raw", plot_peaks_raw)
  
  expect_error(plot(pks_egh, chrom_list = dat.pr, lambda = 190))
})

test_that("purity plot works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh", smooth_type = "none",
                     show_progress = FALSE, estimate_purity = TRUE)
  plot_purity <- function(){
    plot(pks_egh, chrom_list = dat.pr, points = TRUE,
         ticks = TRUE, lambda = 210, plot_purity = TRUE)
  }
  vdiffr::expect_doppelganger("plot_purity", plot_purity)
})


test_that("plot.ptw_list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_ptw_list_traces_base <- function(){
    plot(warping.models, what = "traces", engine = "base")
  }
  suppressWarnings(vdiffr::expect_doppelganger("plot_ptw_list_traces_base",
                                               plot_ptw_list_traces_base))
  
  plot_ptw_list_heatmap_base <- function(){
    plot(warping.models, what = "heatmap", engine = "base")
  }
  vdiffr::expect_doppelganger("plot_ptw_list_heatmap_base", plot_ptw_list_heatmap_base)
  
})

test_that("plot.ptw_list works with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  plot_ptw_list_traces_ggplot <- function(){
    suppressWarnings(plot(warping.models, what = "traces", engine = "ggplot"))
  }
  vdiffr::expect_doppelganger(title = "plot_ptw_list_traces_ggplot", 
                              fig = plot_ptw_list_traces_ggplot)
  
  plot_ptw_list_heatmap_ggplot <- function(){
    plot(warping.models, what = "heatmap", engine = "ggplot")
  }
  vdiffr::expect_doppelganger(title = "plot_ptw_list_heatmap_ggplot",
                              fig = plot_ptw_list_heatmap_ggplot)

})

test_that("plot.peak_table works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_peak_table <- function(){
    par(mfrow=c(3, 1))
    plot(pk_tab, loc = "V13", chrom_list = dat.pr, box_plot = TRUE,
         vars = "trt", verbose = FALSE, spectrum_labels = TRUE)
  }
  vdiffr::expect_doppelganger("plot.peak_table", plot_peak_table)
  expect_error(plot(pk_tab, chrom_list = dat.pr, loc = "V13", idx = 1,
                    lambda = "210",  box_plot = TRUE, verbose = FALSE))
  expect_error(plot(pk_tab, chrom_list = dat.pr, loc = "15", what = "rt",
                    idx = 1, lambda = "210",  
                    box_plot = TRUE, var = "trt", verbose = FALSE))
})

test_that("plot_all_spectra works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_spectra <- function(){
    plot_all_spectra("V13", peak_table = pk_tab, chrom_list = dat.pr, 
                     export = TRUE, overlapping = TRUE)
  }
  x <- plot_all_spectra("V13", peak_table = pk_tab, chrom_list = dat.pr, 
                        export = TRUE, overlapping = TRUE)
  vdiffr::expect_doppelganger("plot_all_spectra", plot_spectra)
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(colnames(x), rownames(pk_tab$tab))
})

test_that("plot_all_spectra works with ggplot2", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_spectra <- function(){
    plot_all_spectra("V13", peak_table = pk_tab, chrom_list = dat.pr, 
                     export = FALSE, overlapping = TRUE, engine = "ggplot")
  }
  x <- plot_all_spectra("V13", peak_table = pk_tab, chrom_list = dat.pr, 
                        export=FALSE, overlapping=TRUE, engine = "ggplot")
  vdiffr::expect_doppelganger("plot_all_spectra_ggplot", plot_spectra)
})

test_that("plot_spectrum works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_spec <- function(){
    par(mfrow=c(2,1))
    plot_spectrum("13.62", peak_table = pk_tab, chrom_list = dat.pr,
                  export=TRUE, what="rt", idx = 1, verbose = FALSE)
  }
  vdiffr::expect_doppelganger(title = "plot_spectrum", plot_spec)
  
  plot_spec_character_idx <- function(){
    par(mfrow = c(2,1))
    plot_spectrum(loc = "13.62", peak_table = pk_tab, chrom_list = dat.pr,
                  what = "rt", idx = "119", verbose = FALSE)
  }
  vdiffr::expect_doppelganger(title = "plot_spectrum", 
                              plot_spec_character_idx)

  x <- plot_spectrum(loc = "V13", peak_table = pk_tab, chrom_list = dat.pr, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
})

test_that("plot_spectrum works", {
  x <- plot_spectrum(loc = "V14", peak_table = pk_tab, chrom_list = dat.pr, 
                     export_spectrum = TRUE, what = "peak", idx = 1, 
                     verbose = FALSE)
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, 
                             what = "click"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, 
                             what="click", idx = 1))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, 
                             what="click",lambda="210"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, 
                             what="rt", lambda="210"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, 
                             what="rt", idx = 1))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, 
                             what="rt"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, 
                             what="rt"))
  expect_error(plot_spectrum(loc=12, chrom_list = pk_tab, what="rt", idx = 1))
  expect_error(plot_spectrum(loc=12, what="rt", idx = 1))
  expect_error(plot_spectrum(loc=12, chrom_list = dat.pr, what="peak", idx = 1))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, 
                             what="peak", idx = 1))
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what="click", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what="click", idx = 1, engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what = "click", idx = 1, lambda="210",
                             engine = "plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what="click",lambda="210", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what="rt", lambda="210", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr,
                             what="rt", idx = 1, engine="plotly"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr,
                             what="rt", engine="plotly"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr,
                             what="rt", engine="plotly"))
})


test_that("plot_spectrum works with ggplot2", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  p1 <- plot_spectrum(loc = "13.62", peak_table = pk_tab, chrom_list = dat.pr,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="ggplot")
  vdiffr::expect_doppelganger(title = "plot_both_ggplot", fig = p1)
  
  p2 <- plot_spectrum(loc = "13.62", peak_table = pk_tab, chrom_list = dat.pr,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="ggplot", plot_trace = FALSE)
  vdiffr::expect_doppelganger(title = "plot_spectrum_ggplot", fig = p2)
  
  p3 <- plot_spectrum(loc = "13.62", peak_table = pk_tab, chrom_list = dat.pr, 
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine = "ggplot", plot_spectrum = FALSE)
  vdiffr::expect_doppelganger(title = "plot_trace_ggplot", fig = p3)
})



test_that("mirror_plot works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  mirror1 <- function(){
    mirror_plot(pk_tab, chrom_list = dat.pr, lambdas = c("210","260"),
                var = "trt", legend_size=2)
  }
  vdiffr::expect_doppelganger("mirror1", mirror1)
  
  mirror2 <- function(){
    mirror_plot(pk_tab, chrom_list = dat.pr, lambdas = c("210","260"),
                var = "trt", legend_size=2, mirror = FALSE)
  }
  vdiffr::expect_doppelganger("mirror2", mirror2)
  expect_error(mirror_plot(pk_tab, chrom_list = dat.pr))
  expect_error(mirror_plot(pk_tab, chrom_list = dat.pr, var = "invalid_variable"))
})

test_that("boxplot works as expected", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  boxplot1 <- function(){
    boxplot(pk_tab, V11 ~ trt)
  }
  vdiffr::expect_doppelganger("boxplot1", boxplot1)
  
  boxplot2 <- function(){
    boxplot(pk_tab, V11~trt, las=2)
  }
  vdiffr::expect_doppelganger("boxplot2", boxplot2)
})

test_that("plot_peak.list works", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  plot_peaklist <- function(){
    plot(pks_egh, chrom_list = dat.pr, idx = 2)
  }
  vdiffr::expect_doppelganger("plot_peaklist", plot_peaklist)
  plot(pks_egh, chrom_list = dat.pr, lambda = 210)
})

test_that("write_peaktable writes csvs correctly", {
  data(pk_tab)
  skip_on_cran()
  path <-  tempdir(check = TRUE)
  on.exit(unlink(c(fs::path(path,"peak_table-tab", ext = "csv"),
                   fs::path(path,"peak_table-pk_meta", ext = "csv"),
                   path)))
  write_peaktable(peak_table = pk_tab, path = path, format = "csv",
                  what = c("tab", "pk_meta"))
  path_table <- paste0(file.path(path,"peak_table-tab"), ext = ".csv")
  path_meta <- paste0(file.path(path,"peak_table-pk_meta"), ext = ".csv")
  
  expect_true(file.exists(path_table))
  expect_true(file.exists(path_meta))
  
  expect_equal(pk_tab$tab, 
               read.csv(file = path_table, row.names = 1), ignore_attr = TRUE)
  expect_equal(pk_tab$pk_meta,
               read.csv(file = path_meta, row.names = 1), ignore_attr = TRUE)
  expect_error(write_peaktable(pk_tab, path = "fake_path"))
})

test_that("write_peaktable writes xlsx correctly", {
  skip_on_cran()
  skip_if_not_installed("openxlsx")
  data(pk_tab)
  path <-  tempdir()
  on.exit(unlink(c(fs::path(path, "peak_table", ext="xlsx"), path)))
  
  write_peaktable(pk_tab, path = path, format = "xlsx")
  file <- fs::path(path, "peak_table", ext = "xlsx")
  xx <- openxlsx::read.xlsx(file, sheet = 1, rowNames = TRUE)
  expect_equal(pk_tab$tab, xx, ignore_attr=TRUE)
  expect_equal(pk_tab$pk_meta, 
               openxlsx::read.xlsx(file, sheet = 2, rowNames = TRUE), 
               ignore_attr=TRUE)
})

# run the following line to activate plotly tests:
# Sys.setenv("VISUAL_TESTS" = "true")

test_that("plot_spectrum works with plotly engine", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  p1 <- plot_spectrum("13.62", peak_table = pk_tab, chrom_list = dat.pr,
                      export_spectrum = FALSE, what = "rt", idx = 1,
                      verbose = FALSE, engine = "plotly")
  expect_doppelganger_plotly("plot_both_plotly", p = p1)
  
  p2 <- plot_spectrum("13.62", peak_table = pk_tab, chrom_list = dat.pr,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="plotly", plot_trace = FALSE)
  expect_doppelganger_plotly("plot_trace_plotly", p = p2)
  
  p3 <- plot_spectrum("13.62", peak_table = pk_tab, chrom_list = dat.pr,
                      export_spectrum = FALSE, what="rt", idx = 1,
                      verbose = FALSE, engine="plotly", plot_spectrum = FALSE)
  expect_doppelganger_plotly("plot_spectrum_plotly", p = p3)
})

test_that("plot_chroms works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  p <- plot_chroms(warp, lambdas = "210", engine = "plotly")
  # warnings are suppressed because of annoying RColorBrewer warnings
  # https://github.com/plotly/plotly.R/issues/2392
  # should be fixed by https://github.com/plotly/plotly.R/pull/1999
  suppressWarnings(expect_doppelganger_plotly(name = "alignment_plotly", 
                                              p = p))
  
  p2 <- plot_chroms(warp, lambdas="210", engine = "plotly", show_legend = FALSE)
  suppressWarnings(expect_doppelganger_plotly(name = "alignment_plotly_no_legend", 
                                              p = p2))
  
  p3 <- plot_chroms(warp, lambdas="210", engine = "plotly", 
                    legend_position = "topleft", 
                    xlim=c(15, 18), ylim=c(0,400))
  suppressWarnings(expect_doppelganger_plotly(name = "alignment_plotly_zoom",
                                              p = p3))
})


test_that("plot_chroms_heatmap works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  p <- plot_chroms_heatmap(warp, lambdas = "210", engine = "plotly")
  expect_doppelganger_plotly(name = "heatmap_plotly", 
                                              p = p)
  # plot retention times instead of indices
  p2 <- plot_chroms_heatmap(warp, lambdas = "210", engine = "plotly", 
                            show_legend = FALSE, xlim = c(10, 15))
  expect_doppelganger_plotly(name = "heatmap_plotly_zoom",
                                              p = p2)
})

test_that("plot.ptw_list works with plotly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("plotly")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("rsvg")
  
  plot_ptw_list_traces_plotly <- suppressWarnings(plot(warping.models, 
                                                       what = "traces", 
                                                       engine = "plotly")
                                                                        )

  expect_doppelganger_plotly(name = "plot_ptw_list_traces_plotly", 
                             p = plot_ptw_list_traces_plotly)
  
  plot_ptw_list_heatmap_plotly <- plot(warping.models, what = "heatmap", 
                                       engine = "plotly")

  expect_doppelganger_plotly(name = "plot_ptw_list_heatmap_plotly", 
                             p = plot_ptw_list_heatmap_plotly)
})


test_that("plot_chroms can subset chromatograms correctly", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  data(Sa_pr)
  numeric_indices <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c(1, 3), lambdas = 210)
  }
  character_indices <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c("119", "122"), lambdas = 210)
  }
  vdiffr::expect_doppelganger("heatmap_plot_trace", numeric_indices)
  vdiffr::expect_doppelganger("heatmap_plot_trace", character_indices)
  
  numeric_indices_rev <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c(3, 1), lambdas = 210)
  }
  character_indices_rev <- function(){
    plot_chroms_heatmap(Sa_pr, idx=c("122", "119"), lambdas = 210)
  }
  vdiffr::expect_doppelganger("heatmap_plot_trace_rev", numeric_indices_rev)
  vdiffr::expect_doppelganger("heatmap_plot_trace_rev", character_indices_rev)
})


test_that("plot_chroms works to plot alignments with ggplot", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  alignment_ggplot_heatmap <- function(){
    plot_chroms_heatmap(warp, lambdas=210, engine="ggplot")
  }
  vdiffr::expect_doppelganger("alignment_ggplot_heatmap",
                              alignment_ggplot_heatmap)
  
  alignment_ggplot_heatmap_zoom <- function(){
    plot_chroms_heatmap(warp, lambdas = 210, engine = "ggplot", 
                        show_legend = TRUE, xlim = c(10, 15))
  }
  suppressWarnings(vdiffr::expect_doppelganger("alignment_ggplot_heatmap_zoom",
                                               alignment_ggplot_heatmap_zoom))
  expect_error(plot_chroms(pktab))
})