### test preprocess function ###
path <- "testdata/DAD1.CSV"
x <- read.csv(path, fileEncoding = "utf-16")
x1 <- suppressWarnings(load_chroms(path, format.in = "csv", find_files = FALSE))
folder <- "testdata"
x2 <- suppressWarnings(load_chroms(folder, format.in = "csv", find_files = TRUE))

test_that("load_chroms works", {
  expect_equal(x[,2], x1[[1]][, "220.00000"], ignore_attr = TRUE)
  expect_equal(x1, x2, ignore_attr = TRUE)
  expect_equal(length(x1), length(path))
})

data("Sa")
new.ts <- seq(10,18.66,by=.01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths

out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works on matrix", {
  out <- preprocess_matrix(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas, maxI = 1)
  expect_equal(class(out)[1],"matrix")
  expect_equal(rownames(out), as.character(new.ts))
  expect_equal(colnames(out), as.character(new.lambdas))
})


test_that("preprocess returns correct errors", {
  expect_error(preprocess(X=as.data.frame(Sa[[1]])))
})


test_that("preprocess works without interpolation", {
  dat.pr <- preprocess(X=Sa[[1]], interpolate_cols = F, interpolate_rows=F)
  expect_equal(dim(dat.pr),dim(Sa[[1]]))
  expect_equal(rownames(dat.pr), rownames(Sa[[1]]))
  expect_equal(colnames(dat.pr), colnames(Sa[[1]]))
})

dat.pr <- preprocess(X = Sa[1:2], dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works on a list", {
  expect_equal(length(dat.pr),length(Sa[1:2]))
  expect_equal(rownames(dat.pr[[1]]), as.character(new.ts))
  expect_equal(colnames(dat.pr[[1]]), as.character(new.lambdas))
  expect_equal(names(dat.pr), names(Sa[1:2]))
})

### test correct_rt ###
warping.models <- correct_rt(dat.pr, lambdas = c('210','260','318'), what = "models", progress_bar = FALSE)
warp <- correct_rt(chrom_list = dat.pr, models = warping.models, what = "corrected.values", progress_bar = FALSE)

test_that("correct_rt works", {
  equals(length(warping.models), length(warp), length(Sa[1:2]))
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_equal(t(warping.models[[1]]$warped.sample)[,1], as.numeric(warp[[1]][,"210"]))
  expect_error(correct_rt(dat.pr))
  expect_error(correct_rt(dat.pr, what="x"))
  expect_error(correct_rt(dat.pr, lambdas = "x"))
  expect_error(correct_rt(dat.pr, lambdas = "210", alg="x"))
  expect_error(correct_rt(dat.pr, lambdas = "210", models = "warping.models", alg = "vpdtw"))
})

test_that("correct_rt works with vpdtw", {
  skip_if_not_installed("VPdtw")
  warp <- correct_rt(dat.pr, lambdas = "210", alg="vpdtw")
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_error(correct_rt(dat.pr, lambdas = c("210","260"), alg="vpdtw",  progress_bar = FALSE))
  expect_error(correct_rt(dat.pr, lambdas = c("x"), alg="vpdtw"))
})

test_that("plot_chroms works to plot alignments", {
  skip_if_not_installed("vdiffr")
  alignment <- function(){
    par(mfrow=c(2,1))
    plot_chroms(warp, lambdas="210")
    plot_chroms(dat.pr, lambdas="210")
  }
  vdiffr::expect_doppelganger("alignment", alignment)
})


### test get_peaks ###
lam <- c("210","318")
pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh", smooth_type = "none", progress_bar = FALSE)
pks_gaussian <- get_peaks(dat.pr, lambdas = lam, fit = "gaussian", smooth_type = "none", progress_bar = FALSE)

test_that("get_peaks works", {
  expect_equal(names(pks_egh), names(dat.pr))
  expect_equal(names(pks_gaussian), names(dat.pr))
  expect_equal(names(pks_egh[[1]]), lam)
  expect_equal(names(pks_gaussian[[1]]), lam)
  expect_equal(class(pks_egh), "peak_list")
  expect_error(get_peaks(dat.pr)) # lambdas must be provided
  expect_error(get_peaks(dat.pr, lambdas="210", fit="nonsense"))
})

pks_s <- filter_peaks(pks_egh, min_height = 10, min_area = 10, min_sd = .07)

test_that("filter_peaks works", {
  expect_equal(names(pks_egh), names(pks_s))
  expect_equal(names(pks_s), names(dat.pr))
  expect_equal(names(pks_s[[1]]), lam)
  expect_equal(class(pks_s), "peak_list")
  expect_lt(nrow(pks_s[[1]][[1]]), nrow(pks_egh[[1]][[1]]))
  expect_lt(nrow(pks_s[[2]][[1]]), nrow(pks_egh[[2]][[1]]))
  expect_warning(filter_peaks(pks_egh))
  # expect_warning(filter_peaks(pks_egh, min_sd=0))
  expect_warning(filter_peaks(pks_egh, min_height=0))
  expect_warning(filter_peaks(pks_egh, min_area=0))
  expect_warning(filter_peaks(pks_egh, max_sd=Inf))
})

### test get_peaktable ###

pk_tab <- get_peaktable(pks_egh, dat.pr)
pk_tab_sp <- get_peaktable(pks_egh, dat.pr, clust = "sp.rt")

test_that("get_peaktable works", {
  expect_equal(rownames(pk_tab$tab), names(dat.pr))
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$pk_meta))
  expect_equal(class(pk_tab), "peak_table")
  expect_equal(class(pk_tab$tab), "data.frame")
  expect_equal(class(pk_tab$pk_meta), "data.frame")
  expect_equal(class(pk_tab$args), "list")
})


### test metadata attachment ###

path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(pk_tab, metadata = meta, column="vial")

test_that("attach_metadata works", {
  expect_equal(rownames(pk_tab$tab), pk_tab$sample_meta$vial)
  expect_equal(colnames(pk_tab$sample_meta), colnames(meta))
  ### check errors & warnings ###
  expect_error(attach_metadata(pk_tab))
  expect_error(attach_metadata(pk_tab, metadata=meta, column ="x"))
  expect_error(attach_metadata(pk_tab$tab, metadata=meta, column ="vial"))
  expect_error(attach_metadata(pk_tab, metadata=rbind(meta,meta), column ="vial"))
  expect_warning(attach_metadata(pk_tab, metadata=meta[-1,], column ="vial"))
})

pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.cor")
test_that("attach_ref_spectra works", {
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  expect_equal(pk_tab$args[["reference_spectra"]], "max.cor")
  pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "max.int")
  expect_equal(pk_tab$args[["reference_spectra"]], "max.int")
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  expect_error(attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "x"))
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

test_that("merge_peaks works", {
  pk_tab_m <- merge_peaks(pk_tab, peaks=c("V10","V11"))
  expect_equal(pk_tab_m$tab$V11, pmax(pk_tab$tab$V10, pk_tab$tab$V11))
  expect_equal(ncol(pk_tab_m), ncol(pk_tab)-1)
  expect_equal(ncol(pk_tab_m$pk_meta), ncol(pk_tab$pk_meta)-1)
  expect_equal(rownames(pk_tab_m$tab), rownames(pk_tab$tab))
})

test_that("normalize_data works", {
  pk_tab_norm <- normalize_data(pk_tab, chrom_list = dat.pr, column = "mass")
  expect_equal(rownames(pk_tab_norm$tab), rownames(pk_tab$tab))
  expect_equal(class(pk_tab_norm), class(pk_tab))
  expect_equal(colnames(pk_tab_norm$tab), colnames(pk_tab$tab))
  chroms <- normalize_data(pk_tab, chrom_list=dat.pr, column = "mass", what = "chrom_list")
  expect_equal(dim(chroms), dim(dat.pr))
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "x"))
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "mass", what="x"))
  expect_equal(pk_tab_norm$args[["normalized"]], "TRUE")
  expect_equal(pk_tab_norm$args[["normalization_by"]], "mass")
})

### cluster

# suppressWarnings(cl <- cluster_spectra(pk_tab, chrom_list = dat.pr, nboot = 10,
#                                        parallel = FALSE, verbose = FALSE,
#                                        save = FALSE, output = "both"))
# test_that("cluster_spectra works", {
#   expect_equal(class(cl[[1]]), "pvclust")
#   expect_equal(class(cl[[2]]), "list")
# })

### test plotting functions

test_that("plot.peak_list works", {
  skip_if_not_installed("vdiffr")
  plot_peaks <- function() plot(pks_egh, chrom_list = dat.pr, points = TRUE,
                                ticks=TRUE)
  vdiffr::expect_doppelganger("plot.peak_list", plot_peaks)
})

test_that("plot.peak_table works", {
  skip_if_not_installed("vdiffr")
  plot_peak_table <- function(){
    par(mfrow=c(3,1))
    plot(pk_tab, loc = "V13", chrom_list = dat.pr, box_plot = TRUE,
         vars = "trt", verbose = FALSE, spectrum_labels = TRUE)
  }
  vdiffr::expect_doppelganger("plot.peak_table", plot_peak_table)
})


test_that("plot_all_spectra works", {
  skip_if_not_installed("vdiffr")
  plot_spectra <- function(){
    plot_all_spectra("V13", peak_table=pk_tab, chrom_list = dat.pr, export=TRUE,
                     overlapping=TRUE)
  }
  x <- plot_all_spectra("V13", peak_table=pk_tab, chrom_list = dat.pr, export=TRUE,
                        overlapping=TRUE)
  vdiffr::expect_doppelganger("plot_all_spectra", plot_spectra)
  expect_equal(class(x), "data.frame")
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(colnames(x), rownames(pk_tab$tab))
})

test_that("plot_spectrum works", {
  skip_if_not_installed("vdiffr")
  plot_spec <- function(){
    par(mfrow=c(2,1))
    plot_spectrum("13.62", peak_table=pk_tab, chrom_list = dat.pr, export=TRUE,
                  what="rt", chr=1, verbose = FALSE)
  }
  vdiffr::expect_doppelganger("plot_spectrum", plot_spec)
  x <- plot_spectrum("V13", peak_table=pk_tab, chrom_list = dat.pr, export=TRUE,
                     what="peak", chr=1, verbose = FALSE)
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
})

test_that("plot_spectrum works", {
  x <- plot_spectrum("V13", peak_table=pk_tab, chrom_list = dat.pr, export=TRUE,
                     what="peak", chr=1, verbose = FALSE)
  expect_equal(rownames(x), as.character(new.lambdas))
  expect_equal(class(x), "data.frame")
  expect_equal(ncol(x), 1)
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click", chr=1))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click",lambda="210"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="rt", lambda="210"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="rt", chr=1))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, what="rt"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, what="rt"))
  
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click", chr=1, engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click", chr=1, lambda="210",
                             engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="click",lambda="210", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="rt", lambda="210", engine="plotly"))
  expect_error(plot_spectrum(peak_table = pk_tab, chrom_list = dat.pr, what="rt", chr=1, engine="plotly"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, what="rt", engine="plotly"))
  expect_error(plot_spectrum(loc=12, peak_table = pk_tab, chrom_list = dat.pr, what="rt", engine="plotly"))
})

test_that("mirror_plot works", {
  skip_if_not_installed("vdiffr")
  mirror1 <- function(){
    mirror_plot(pk_tab, chrom_list = dat.pr, lambdas = c("210","260"),
                var = "trt", legend_size=2)
  }
  vdiffr::expect_doppelganger("mirror1", mirror1)
})

test_that("plot_peak.list works", {
  skip_if_not_installed("vdiffr")
  plot_peaklist <- function(){
    plot(pks_egh, chrom_list = dat.pr, index=2)
  }
  vdiffr::expect_doppelganger("plot_peaklist", plot_peaklist)
})

test_that("cluster_spectra works", {
  skip_on_cran()
  cl <- cluster_spectra(pk_tab, chrom_list = dat.pr, nboot = 10,
                          parallel = FALSE, verbose = FALSE, save = FALSE,
                        output = "both", plot_dend = FALSE, plot_spectra = FALSE)
  expect_equal(class(cl[[1]]), "pvclust")
  expect_equal(class(cl[[2]]), "list")
})
