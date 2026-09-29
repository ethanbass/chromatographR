## Fitted peak-model columns carried into a peak table's `pk_meta`, shared by
## the `peak_list` and `matched_peak_list` methods so a new model column is
## registered in one place.
PEAK_META_COLS <- c("sd", "width", "tau", "FWHM", "r.squared", "purity",
                    "window_limited", "converged")

#' Convert peak lists into an aligned peak table
#' 
#' Converts a `peak_list` object containing peaks from multiple samples into
#' an aligned peak table by clustering peaks across samples according to
#' retention time and (optionally) spectral similarity. The tree is cut at
#' height `hmax`, whose meaning depends on `clust` (see `hmax`).
#' 
#' By default, clustering is performed on retention times alone (when 
#' `clust = "rt"`). Alternatively, clustering can also incorporate information 
#' about spectral similarity (when `clust = "sp.rt"`) using a similarity function 
#' adapted from Broeckling et al., 2014:
#' \deqn{
#' S_{ij} =
#'   \exp\left(-\frac{(t_i - t_j)^2}{2\sigma_t^2}\right)
#'    \cdot
#'    \exp\left(-\frac{(1 - c_{ij})^2}{2\sigma_r^2}\right)
#' }
#' where \eqn{c_{ij}} is the spectral correlation between spectra \eqn{i} and 
#' \eqn{j}, and \eqn{\sigma_t} and \eqn{\sigma_r} control the relative 
#' contributions of retention time and spectral similarity, respectively (see
#' `sigma_t` and `sigma_r`).
#' 
#' When two peaks from the same sample fall in the same cluster, only the more
#' intense one enters the peak table, with a warning if `verbose = TRUE`. Such
#' warnings are often benign, but they can also mean `hmax` is too large. Lowering `hmax` merges fewer peaks, at the
#' cost of splitting some across multiple clusters; filtering low-intensity
#' peaks before clustering reduces the warnings without that trade-off. 
#' 
#' Once clusters are formed, peak metadata are summarized within each cluster
#' according to `summarize_by`, by default weighting each peak by its
#' intensity, which reduces the influence of small peaks.
#' 
#' @name get_peaktable
#' @aliases get_peaktable
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom fastcluster hclust
#' @importFrom stats dist cutree as.dist aggregate
#' @importFrom lattice panel.stripplot panel.abline stripplot
#' @importFrom grDevices colorRampPalette 
#' @importFrom scales rescale
#' @importFrom graphics par
#' @param peak_list A `peak_list` object created by [`get_peaks`], or a
#' `matched_peak_list` created by [`match_peaks`].
#' @param chrom_list A list of chromatographic matrices. Only required when
#' `clust = "sp.rt"`.
#' @param response Peak intensity measure: either `"area"` (default) or 
#' `"height"`.
#' @param use.cor Logical. Indicates whether to use corrected retention times
#' (`rt.cor` column) or raw retention times (`rt` column). Unless
#' otherwise specified, the `rt.cor` column will be used by default if it 
#' exists in the provided `peak_list`.
#' @param hmax Height at which the clustering tree is cut. With `clust = "rt"`,
#' peaks are clustered by complete linkage on retention time, so `hmax` is the
#' largest retention time difference allowed between any two peaks in the same
#' cluster, in minutes. With `clust = "sp.rt"`, peaks are clustered
#' by average linkage on the dissimilarity \eqn{1 - S_{ij}}, which runs from 0
#' to 1, and `hmax` is the maximum tree height passed to
#' [`cutreeDynamicTree`][dynamicTreeCut::cutreeDynamicTree]; clusters may be
#' split below it (see `deepSplit`). There, the default of `0.2` lets peaks with
#' identical spectra sit about `0.67 * sigma_t` apart, a third of the mean
#' `end - start` at the default `sigma_t`. Ignored
#' for peak lists produced by `get_peaks(windows = )`, since those peaks are
#' grouped by integration window instead of being clustered.
#' @param summarize_by How to summarize peak metadata within each cluster.
#' Options are `"weighted.mean"` (default, weights peaks by their intensity as
#' given by `response`), `"mean"`, `"median"` (which aggregate metadata across 
#' all peaks in the cluster) or `"max"` (which selects the most intense peak in 
#' the cluster and uses its metadata directly).
#' @param plot_it Logical. If `TRUE`, for every component a strip plot will be
#' shown indicating the clustering.
#' @param ask Logical. Ask before showing new plot? Defaults to `plot_it`.
#' @param clust Specify whether to perform hierarchical clustering based on
#' spectral similarity and retention time (`sp.rt`) or retention time alone
#' (`rt`). Defaults to `rt`. Ignored for peak lists produced by
#' `get_peaks(windows = )`, which are grouped by window.
#' @param sigma_t Width of Gaussian kernel for retention time similarity, in
#' minutes. Controls the weight given to retention time when
#' `clust = "sp.rt"`. Defaults to half the mean `end - start` across all peaks.
#' @param sigma_r Width of Gaussian kernel for spectral similarity. Controls
#' the weight given to spectral correlation when `clust = "sp.rt"`. Defaults
#' to `0.5`.
#' @param sigma.t The `sigma.t` argument is deprecated. Please use `sigma_t`,
#' which is in minutes. `sigma.t` is read in the units of the retention times,
#' as before.
#' @param sigma.r The `sigma.r` argument is deprecated. Please use `sigma_r`
#' instead.
#' @param deepSplit Logical. Controls sensitivity to cluster splitting when
#' `clust = "sp.rt"`. If `TRUE`, returns a larger number of smaller clusters.
#' Defaults to `FALSE`. For additional information, see
#' [`cutreeDynamicTree`][dynamicTreeCut::cutreeDynamicTree].
#' @param verbose Logical. Whether to warn when two peaks from the same sample
#' fall in the same cluster. Defaults to `FALSE`.
#' @param out Output format: `data.frame` (default) or `matrix`.
#' @param na.fill Value used to fill cells where no peak was matched in a given
#'   sample. Defaults to `0`. Only used by the `matched_peak_list`
#'   method.
#' @param ... Additional arguments (currently unused).
#' @return An S3 [`peak_table`] object, containing the following elements:
#' * `tab`: The peak table itself — a `data.frame` (or `matrix`, per `out`) of
#' intensities in a sample × peak configuration.
#' * `pk_meta`: A `data.frame` containing peak metadata (e.g., the spectral
#' component, peak number, and average retention time), summarized across all
#' peaks in each cluster according to `summarize_by`.
#' * `sample_meta`: A `data.frame` of sample metadata. `NA` until added with
#' [`attach_metadata`].
#' * `instrument_meta`: A `data.frame` of instrumental metadata transferred
#' automatically from the chromatograms.
#' * `ref_spectra`: A matrix of reference spectra (in a wavelength × peak
#' configuration). `NA` until added with [`attach_ref_spectra`].
#' * `args`: A list of the arguments given to [`get_peaktable`] to generate
#' the peak table.
#' @author Ethan Bass
#' @note Peak lists produced by `get_peaks(windows = )` are grouped by
#' integration window rather than clustered by retention time, including after
#' [`correct_peaks`]. This gives one peak table column per window, in order of
#' retention time, with every sample contributing to each.
#' @note This function is adapted from the `getPeakTable` function in the alsace
#' package by Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R>.
#' @md
#' @references
#' * Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni J.E. 2014.
#' RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching-Based
#' Annotation for Metabolomics Data. *Anal. Chem.* 86:6812-6817. 
#' \doi{10.1021/ac501530d}.
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. 
#' *Metabolomics* 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @examples
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'))
#' get_peaktable(pks, response = "area")
#' @seealso [`peak_table-class`], [`attach_ref_spectra`], [`attach_metadata`]
#' @export
get_peaktable <- function(peak_list, ...) UseMethod("get_peaktable")

#' @rdname get_peaktable
#' @method get_peaktable default
#' @export
get_peaktable.default <- function(peak_list, chrom_list, response = c("area", "height"),
                          use.cor = NULL, hmax = 0.2, 
                          summarize_by = c("weighted.mean", "median", "mean", "max"),
                          plot_it = FALSE, ask = plot_it, 
                          clust = c("rt", "sp.rt"), 
                          sigma_t = NULL, sigma_r = 0.5,
                          deepSplit = FALSE, verbose = FALSE,
                          out = c('data.frame', 'matrix'), ...,
                          sigma.t = NULL, sigma.r = NULL){
  if (...length()) stop("unused argument(s): ", paste(...names(), collapse = ", "))
  sigma_r <- resolve_deprecated(sigma.r, sigma_r)
  rtf <- rt_factor(get_time_unit(peak_list))
  if (!is.null(sigma.t)){
    warning(paste("The `sigma.t` argument is deprecated. Please use `sigma_t`,",
                  "which is in minutes. `sigma.t` is read in the units of the",
                  "retention times, as before."), immediate. = TRUE)
    sigma_t <- sigma.t / rtf
  }
  if (!is.list(peak_list))
    stop("No get_peaktable method for object of class ",
         paste(class(peak_list), collapse = "/"), ".")
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt", "sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  summarize_by <- match.arg(summarize_by, choices = c("weighted.mean", 
                                                      "median", "mean", "max"))
  if (is.null(use.cor)){
    use.cor <- "rt.cor" %in% colnames(peak_list[[1]][[1]])
  }
  rt <- ifelse(use.cor, "rt.cor", "rt")
  start <- ifelse(use.cor, "start.cor", "start")
  end <- ifelse(use.cor, "end.cor", "end")
  check_duplicated_names(names(peak_list))
  windows <- attr(peak_list, "windows")
  if (clust == "sp.rt" && is.null(windows)){
    if (missing(chrom_list)){
      chrom_list <- get_chrom_list(peak_list)
    } else get_chrom_list(peak_list, chrom_list)
  }
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plot_it) {
    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }
  clusterPeaks <- function(comp, pkLst){
    pkLst <- lapply(pkLst, function(x) lapply(x, function(y){
      if (nrow(y) > 0){
      y[!is.na(y[, rt]), , drop = FALSE]
    } else {
      y
    }}))
    xx <- do.call(rbind, sapply(pkLst, function(samp) samp[comp]))
    file.idx <- xx$sample
    pkcenters <- xx[, rt]
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (!is.null(windows)){
      dist_to_window <- outer(xx[["start"]], windows[, 1], function(a, b) abs(a - b)) +
                        outer(xx[["end"]], windows[, 2], function(a, b) abs(a - b))
      win_idx <- apply(dist_to_window, 1, which.min)
      pkcenters.cl <- match(win_idx, sort(unique(win_idx)))
    } else if (clust == 'rt'){
      pkcenters.hcl <- fastcluster::hclust(dist(pkcenters), method = "complete")
      pkcenters.cl <- cutree(pkcenters.hcl, h = hmax * rtf)
    } else if (clust == 'sp.rt'){
        sigma_t_rt <- if (is.null(sigma_t)){
          0.5 * mean(do.call(rbind, unlist(pkLst, recursive = FALSE))$end -
                       do.call(rbind, unlist(pkLst, recursive = FALSE))$start)
        } else sigma_t * rtf
      ts <- as.numeric(rownames(chrom_list[[1]]))
      sp <- sapply(seq_along(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][
          which(elementwise.all.equal(ts, pkcenters[i])),]))
      }, simplify = TRUE)
      cor.spec <- cor(sp, method = "pearson")
      rt.diff <- abs(outer(unlist(pkcenters), unlist(pkcenters), FUN="-"))
      sim.spec <- exp((-(1 - cor.spec)^2)/(2*sigma_r^2))
      sim.time <- exp(-(rt.diff^2)/(2*sigma_t_rt^2))
      S <- sim.spec*sim.time
      D <- 1 - S
      linkage <- "average"
      pkcenters.hcl <- fastcluster::hclust(as.dist(D), method = linkage)
      pkcenters.cl <- dynamicTreeCut::cutreeDynamicTree(pkcenters.hcl, 
                                                        maxTreeHeight = hmax, 
                                                        deepSplit = deepSplit, 
                                                        minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + seq_along(sing)
    }
    vars <- c(rt, start, end, PEAK_META_COLS)
    vars <- vars[vars %in% colnames(xx)]
    vars.idx <- match(vars, colnames(xx))
    if (summarize_by == "weighted.mean"){
      weights <- xx[, response]
      cl.centers <- do.call(rbind,
                            tapply(seq_along(pkcenters),
                                   pkcenters.cl, function(i){
                                     vapply(vars.idx, function(v){
                                       weighted_mean_finite(xx[i, v], w = weights[i])
                                     }, numeric(1)) |>
                                       setNames(colnames(xx[, vars.idx]))
                                   }))
    } else if (summarize_by == "max"){
      cl.centers <- do.call(rbind, 
                            tapply(seq_along(pkcenters), 
                                   pkcenters.cl, function(i){
        xx[i[which.max(xx[i, response])], vars.idx]
      }))
    } else{
      cl.centers <- aggregate(xx[, vars.idx], by = list(pkcenters.cl),
                              FUN = summarize_by,
                              na.action = "na.pass")[, -1, drop = FALSE]
    }
    ncl <- length(cl.centers[, rt])
    pkcenters.cl <- order(order(cl.centers[, rt]))[pkcenters.cl]
    cl.centers <- cl.centers[order(cl.centers[, rt]),]
    metaInfo <- cbind(lambda = rep(suppressWarnings(
      as.numeric(names(peak_list[[1]])[comp])), ncl),
                      peak = 1:ncl, 
                      round(cl.centers, 2)
                      )
    rownames(metaInfo) <- NULL
    if (plot_it){
      mycols <- myPalette(nrow(cl.centers))
      cl.df <- data.frame(peaks = pkcenters, 
                          files = factor(file.idx), 
                          cluster = pkcenters.cl)
      print(stripplot(files ~ peaks, data = cl.df, 
                      col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl %% 14,
                      xlab = "Retention time", ylab = "",
                      main = paste("Component", comp),
                      panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers[,rt], col = mycols)
                      }))
    }
    if (verbose & max(clusCount <- table(file.idx, pkcenters.cl)) > 1){
      warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    }
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, response]))
    allIs[!is.finite(allIs) | allIs < 0] <- 0
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, names(pkLst)))
    for (i in seq_along(allIs)){
      Iinfo[pkcenters.cl[i],  file.idx[i]] <- 
      max(allIs[i], Iinfo[pkcenters.cl[i], file.idx[i]])
    }
    return(list(Iinfo, metaInfo))
  }
  result <- lapply(seq_len(ncomp), clusterPeaks, peak_list)
  pt <- new_peak_table(result, out, peak_list,
                 args = list(peak_list = deparse(substitute(peak_list)),
                        chrom_list = attr(peak_list, "chrom_list"),
                        chrom_list_aligned = attr(peak_list, "chrom_list_aligned") %||% NA,
                        lambdas = list(names(peak_list[[1]])),
                        response = response,
                        use.cor = use.cor,
                        hmax = hmax,
                        clust = clust,
                        sigma_t = sigma_t,
                        sigma_r = sigma_r,
                        deepSplit = deepSplit,
                        reference_spectra = NA,
                        metadata_path = NA,
                        normalized = FALSE,
                        normalization_by = NA,
                        summarized_by = summarize_by,
                        time_unit = get_metadata_attribute(peak_list, "time_unit"),
                        rt_unit = get_time_unit(peak_list),
                        intensity_unit = get_metadata_attribute(peak_list, "intensity_unit")
                        ))
  attr(pt, "warp_models") <- attr(peak_list, "warp_models")
  pt
}

#' @rdname get_peaktable
#' @param pks A `peak_list` object from [`get_peaks`].
#'   Required when `peak_list` is a `matched_peak_list`. Used to
#'   retrieve per-peak intensities and metadata.
#' @param rt_tol Numeric. Maximum allowed deviation (in retention-time units)
#'   when joining matched peak RTs back to the corresponding row in `pks`.
#'   Defaults to `0.01`.
#' @note When `peak_list` is a `matched_peak_list`, the columns of
#'   the returned [`peak_table`] are defined by the reference
#'   sample's peaks at the wavelength used by [`match_peaks`]; other wavelengths
#'   in `pks` are ignored. Peaks absent from the reference sample will not
#'   appear as columns regardless of their prevalence in other samples. Choose a
#'   reference with a rich peak profile. To discover peaks across all samples,
#'   call `get_peaktable` on the `peak_list` itself, which clusters peaks
#'   across every sample.
#' @method get_peaktable matched_peak_list
#' @export
get_peaktable.matched_peak_list <- function(peak_list, pks,
                                            response = c("area", "height"),
                                            rt_tol = 0.01,
                                            na.fill = 0,
                                            out = c("data.frame", "matrix"),
                                            verbose = FALSE, ...) {
  if (...length()) stop("unused argument(s): ", paste(...names(), collapse = ", "))
  matched  <- peak_list
  response <- match.arg(response, c("area", "height"))
  out      <- match.arg(out, c("data.frame", "matrix"))

  if (missing(pks))
    stop("Please supply the original peak_list from get_peaks() as 'pks'.")
  if (!inherits(pks, "peak_list"))
    stop("'pks' must be a 'peak_list' object.")
  if (!identical(names(matched), names(pks)))
    stop("Sample names of 'peak_list' and 'pks' must match.")

  ref_idx      <- reference_position(matched)
  ref_name     <- names(pks)[ref_idx]
  lambdas      <- attr(matched, "lambda")
  sample_names <- names(pks)

  process_lambda <- function(lam) {
    ref_pks <- pks[[ref_idx]][[lam]]
    ref_pks <- ref_pks[!is.na(ref_pks$rt), , drop = FALSE]
    if (nrow(ref_pks) == 0) return(NULL)

    n_ref      <- nrow(ref_pks)
    rt_ref_all <- ref_pks$rt
    Imat <- matrix(na.fill, nrow = length(sample_names), ncol = n_ref,
                   dimnames = list(sample_names, NULL))

    ref_vals <- ref_pks[, response]
    ref_vals[!is.finite(ref_vals) | ref_vals < 0] <- na.fill
    Imat[ref_idx, ] <- ref_vals

    for (samp_i in seq_along(sample_names)[-ref_idx]) {
      samp_name <- sample_names[samp_i]
      m <- matched[[samp_name]]
      if (is.null(m) || nrow(m$matches) == 0) next

      samp_pks <- pks[[samp_name]][[lam]]
      samp_pks <- samp_pks[!is.na(samp_pks$rt), , drop = FALSE]
      if (nrow(samp_pks) == 0) next

      col_d <- abs(outer(m$matches[, "rt_ref"], rt_ref_all, "-"))
      row_d <- abs(outer(m$matches[, "rt_samp"], samp_pks$rt, "-"))
      col_idx <- max.col(-col_d, "first")
      row_idx <- max.col(-row_d, "first")
      k <- seq_len(nrow(m$matches))
      vals <- samp_pks[[response]][row_idx]
      ok <- col_d[cbind(k, col_idx)] <= rt_tol &
        row_d[cbind(k, row_idx)] <= rt_tol &
        is.finite(vals) & vals >= 0

      Imat[cbind(rep(samp_i, sum(ok)), col_idx[ok])] <- vals[ok]
      n_assigned <- sum(ok)

      if (verbose) {
        n_drop <- n_ref - n_assigned
        message(sprintf("  %s [%s]: %d/%d reference peaks matched%s",
                        samp_name, lam, n_assigned, n_ref,
                        if (n_drop > 0) sprintf(" (%d dropped)", n_drop) else ""))
      }
    }

    present <- c("rt", "start", "end", PEAK_META_COLS)
    present <- present[present %in% colnames(ref_pks)]
    metaInfo <- cbind(
      lambda = rep(suppressWarnings(as.numeric(lam)), n_ref),
      peak   = seq_len(n_ref),
      round(ref_pks[, present, drop = FALSE], 2)
    )
    rownames(metaInfo) <- NULL
    list(t(Imat), metaInfo)
  }

  result <- Filter(Negate(is.null), lapply(lambdas, process_lambda))
  if (length(result) == 0)
    stop("No peaks found in reference sample across any wavelength.")

  new_peak_table(
    result, out, pks,
    args = list(
      peak_list        = deparse(substitute(peak_list)),
      pks              = deparse(substitute(pks)),
      chrom_list       = attr(pks, "chrom_list"),
      lambdas          = list(lambdas),
      reference        = ref_idx,
      reference_sample = ref_name,
      response         = response,
      rt_tol           = rt_tol,
      reference_spectra = NA,
      metadata_path    = NA,
      normalized       = FALSE,
      normalization_by = NA,
      time_unit        = get_metadata_attribute(pks, "time_unit"),
      rt_unit          = get_time_unit(pks),
      intensity_unit   = get_metadata_attribute(pks, "intensity_unit")
    ))
}

#' Assemble a peak table
#'
#' Shared by the `peak_list` and `matched_peak_list` methods so the object
#' layout is defined once. Each method supplies its own `args`, which is
#' recorded verbatim.
#' @noRd
new_peak_table <- function(result, out, source, args){
  as.structure <- switch(out, "data.frame" = as.data.frame, "matrix" = as.matrix)
  pull <- function(k) as.structure(t(do.call("rbind", lapply(result, `[[`, k))))
  structure(
    list(tab = pull(1), pk_meta = pull(2),
         sample_meta = NA,
         instrument_meta = chromConverter::extract_metadata(source),
         ref_spectra = NA,
         args = args),
    class = c("peak_table", "list"),
    pk_args = attr(source, "meta"))
}

#' Reshape peaktable
#' 
#' Reshapes peak table from wide to long format
#' @name reshape_peaktable
#' @importFrom stats reshape
#' @param x A `peak_table` object.
#' @param peaks A character vector specifying the peaks to include. If the
#' character vector is named, the names of the vector elements will be used in
#' place of the original peak names.
#' @param metadata A character vector specifying the metadata fields to include.
#' @param fixed_levels Logical. Whether to fix factor levels of features in the
#' order provided. Defaults to `TRUE`.
#' @return A data.frame containing the information for the specified peaks in
#' long format.
#' @author Ethan Bass
#' @family utility functions
#' @examples
#' data(pk_tab)
#' reshape_peaktable(pk_tab, peaks = c("V10", "V15"))
#' @export
reshape_peaktable <- function(x, peaks, metadata, fixed_levels = TRUE){
  if (!missing(peaks)){
    if (is.numeric(peaks)){
      peaks <- colnames(x$tab)[peaks]
    }
    df <- x$tab[, match(peaks, colnames(x$tab)), drop = FALSE]
    if (!is.null(names(peaks))){
      colnames(df) <- names(peaks)
      peaks <- colnames(df)
    }
  } else {
    df <- x$tab
  }
  if (!missing(metadata)){
    meta_idx <- which(colnames(x$sample_meta) %in% metadata)
    x$sample_meta <- x$sample_meta[, meta_idx, drop = FALSE]
  }
  xx <- reshape(as.data.frame(chr = rownames(df), df), direction = "long",
                varying = list(seq_len(ncol(df))), v.names = x$args[["response"]],
                times = colnames(df), timevar = "peak",
                idvar = "sample", ids = rownames(df))
  rownames(xx) <- NULL
  xx <- merge(xx, data.frame(peak = colnames(x$pk_meta), 
                             t(x$pk_meta[c("lambda", "rt"),])),
              by = "peak", all.x = TRUE)
  xx <- xx[, c(1, 3, 4, 5, 2)]
  if (!is.null(dim(x$sample_meta))){
    xx <- merge(xx, data.frame(sample = row.names(df), x$sample_meta),
                by = "sample", all.x = TRUE)
  }
  if (fixed_levels){
    xx$peak <- factor(xx$peak, levels = peaks)
  }
  xx
}
