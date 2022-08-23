### test preprocess ###
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

# test_that("load_chroms can read chemstation UV data", {
#   path_csv <- system.file("testdata/DAD1.CSV", package="chromConverter")
#   path_uv <- system.file("testdata/DAD1.uv", package="chromConverter")
#   x <- read.csv(path_csv, fileEncoding = "utf-16")
#   x1 <- load_chroms(path_uv, format.in = "chemstation", find_files = FALSE)
#   folder <- system.file("testdata", package="chromConverter")
#   x2 <- load_chroms(folder, format.in = "chemstation", find_files = TRUE)
#   expect_equal(x[,2], x1[[1]][, "220.0"], ignore_attr = TRUE)
#   expect_equal(x1, x2, ignore_attr = TRUE)
# })

data("Sa")
new.ts <- seq(10,18.66,by=.01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths

out <- preprocess(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)

test_that("preprocess works on matrix", {
  out <- preprocess_matrix(X = Sa[[1]], dim1 = new.ts, dim2 = new.lambdas,maxI = 1)
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

test_that("correct_rt works", {
  warping.models <- correct_rt(dat.pr, lambdas=c('210','260','318'), what = "models")
  warp <- correct_rt(chrom_list=dat.pr, models=warping.models, what = "corrected.values")
  equals(length(warping.models), length(warp), length(Sa[1:2]))
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_equal(t(warping.models[[1]]$warped.sample)[,1], as.numeric(warp[[1]][,"210"]))
  expect_error(correct_rt(dat.pr))
  expect_error(correct_rt(dat.pr, what="x"))
  expect_error(correct_rt(dat.pr, lambdas = "x"))
  expect_error(correct_rt(dat.pr, lambdas = "210", alg="x"))
})

test_that("correct_rt works with vpdtw", {
  skip_if_not_installed("VPdtw")
  warp <- correct_rt(dat.pr, lambdas = "210", alg="vpdtw")
  expect_equal(names(warp), names(dat.pr[1:2]))
  expect_equal(colnames(warp[[1]]), colnames(dat.pr[[1]]), ignore_attr=TRUE)
  expect_error(correct_rt(dat.pr, lambdas = c("210","260"), alg="vpdtw"))
  expect_error(correct_rt(dat.pr, lambdas = c("x"), alg="vpdtw"))
})


### test get_peaks ###
# lam <- c("260") # tests fail if I use multiple wavelengths here
lam <- c("210","318")
pks_egh <- get_peaks(dat.pr, lambdas = lam, fit = "egh")
pks_gaussian <- get_peaks(dat.pr, lambdas = lam, fit = "gaussian")

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
  expect_equal(class(pk_tab$args), "character")
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

pk_tab <- attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref="max.cor")
test_that("attach_ref_spectra works", {
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  pk_tab <- attach_ref_spectra(pk_tab, chrom_list=dat.pr, ref = "max.int")
  expect_equal(colnames(pk_tab$tab), colnames(pk_tab$ref_spectra))
  expect_error(attach_ref_spectra(pk_tab, chrom_list = dat.pr, ref = "x"))
})
# 
# test filter_peaktable

pktab_s <- filter_peaktable(pk_tab, min_rt=6, max_rt=15)

test_that("filter_peaktable works", {
  expect_equal(rownames(pk_tab$tab), rownames(pktab_s$tab))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$pk_meta))
  expect_equal(nrow(pktab_s$tab), nrow(pktab_s$sample_meta))
  expect_equal(colnames(pktab_s$tab), colnames(pktab_s$ref_spectra))
  expect_warning(filter_peaktable(pk_tab))
})

test_that("combine_peaks works", {
  suppressWarnings(pk_tab_c <- combine_peaks(pk_tab))
  expect_lt(ncol(pk_tab_c), ncol(pk_tab))
  expect_lt(ncol(pk_tab_c$pk_meta), ncol(pk_tab$pk_meta))
  expect_equal(rownames(pk_tab_c$tab), rownames(pk_tab$tab))
})

test_that("normalize_data works", {
  pk_tab_norm <- normalize_data(pk_tab, chrom_list = dat.pr, column = "mass")
  expect_equal(rownames(pk_tab_norm$tab), rownames(pk_tab$tab))
  expect_equal(class(pk_tab_norm), class(pk_tab))
  expect_equal(colnames(pk_tab_norm$tab), colnames(pk_tab$tab))
  chroms <- normalize_data(pk_tab, chrom_list=dat.pr, column = "mass",what = "chrom_list")
  expect_equal(dim(chroms), dim(dat.pr))
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "x"))
  expect_error(normalize_data(pk_tab, chrom_list = dat.pr, column = "mass", what="x"))
})

### cluster

suppressWarnings(cl <- cluster_spectra(pk_tab, chrom_list = dat.pr, nboot = 10,
                                       parallel = FALSE, verbose = FALSE,
                                       save = FALSE, output = "both"))
test_that("cluster_spectra works", {
  expect_equal(class(cl[[1]]), "pvclust")
  expect_equal(class(cl[[2]]), "list")
})

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
                                     vars = "trt", verbose = FALSE)
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
})

test_that("mirror_plot works", {
  skip_if_not_installed("vdiffr")
  mirror1 <- function(){
    mirror_plot(pk_tab, chrom_list = dat.pr, lambdas = c("210","260"),
                var = "trt", legend_size=2)
  }
  vdiffr::expect_doppelganger("mirror1", mirror1)
})

#test fit_peaks

test_that("fit_peaks works independently", {
  y<-dat.pr[[1]][,'280']
  pos<-find_peaks(y)
  pks<-fit_peaks(y,pos, max.iter=1000)
  expect_equal(ncol(pos),3)
  expect_equal(ncol(pks),9)
})

