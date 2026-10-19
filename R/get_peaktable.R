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
#' @param use_corrected_rt Logical. Indicates whether to use corrected retention times
#' (`rt.cor` column) or raw retention times (`rt` column). Unless
#' otherwise specified, the `rt.cor` column will be used by default if it
#' exists in the provided `peak_list`. Either way, they are stored in the `rt`,
#' `start` and `end` rows of `pk_meta`, and `args$use_corrected_rt` records
#' which were used.
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
#' for peak lists produced by `get_peaks(windows = )` or `get_peaks(targets = )`,
#' since those peaks are grouped by integration window or target instead of
#' being clustered.
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
#' `get_peaks(windows = )` or `get_peaks(targets = )`, which are grouped by window
#' or target.
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
#' @param use.cor The `use.cor` argument is deprecated. Please use
#' `use_corrected_rt` instead.
#' @param deepSplit Logical. Controls sensitivity to cluster splitting when
#' `clust = "sp.rt"`. If `TRUE`, returns a larger number of smaller clusters.
#' Defaults to `FALSE`. For additional information, see
#' [`cutreeDynamicTree`][dynamicTreeCut::cutreeDynamicTree].
#' @param verbose Logical. Whether to warn when two peaks from the same sample
#' fall in the same cluster. Defaults to `FALSE`.
#' @param na.fill Value used to fill cells where no peak was matched in a given
#'   sample. Defaults to `0`. Only used by the `matched_peak_list`
#'   method.
#' @param ... Additional arguments (currently unused).
#' @return An S3 [`peak_table`] object, containing the following elements:
#' * `tab`: The peak table itself — a `data.frame` of
#' intensities in a sample × peak configuration. A peak that was not found in a
#' sample, or whose response is negative, is recorded as `0` (`na.fill` for a
#' `matched_peak_list`). A peak whose model failed, which
#' `get_peaks(targets = )` keeps with a missing `height` and `area`, is
#' recorded as `NA`.
#' * `pk_meta`: A `data.frame` of peak metadata, one column per peak. Its rows
#' are the wavelength (`lambda`), the peak number (`peak`), and whichever of the
#' peak list's `rt`, `start`, `end`, `sd`, `tau`, `FWHM`, `r.squared`, `purity`,
#' `window_limited` and `converged` columns it has, each summarized across the
#' peaks in that column according to `summarize_by`. See the Value section of
#' [`get_peaks`] for their meanings.
#' * `sample_meta`: A `data.frame` of sample metadata: the fields added to the
#' chromatograms with `chromConverter::add_metadata()` (chromConverter 0.11.0
#' or later), and any attached with [`attach_metadata`]. `NA` if there are none.
#' * `instrument_meta`: A `data.frame` of instrumental metadata transferred
#' automatically from the chromatograms.
#' * `ref_spectra`: A matrix of reference spectra (in a wavelength × peak
#' configuration). `NA` until added with [`attach_ref_spectra`].
#' * `args`: A list of the arguments given to [`get_peaktable`] to generate
#' the peak table.
#' @author Ron Wehrens, Ethan Bass
#' @note Peak lists produced by `get_peaks(windows = )` are grouped by
#' integration window rather than clustered by retention time, including after
#' [`correct_peaks`]. This gives one peak table column per window, in order of
#' retention time, with every sample contributing to each.
#' @note Peak lists produced by `get_peaks(targets = )` are grouped by target
#' instead of clustered, each peak going to the target it was selected for.
#' Targets found in no sample get no column. If `targets` has a `name` column, the
#' columns are named after the targets.
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
                          use_corrected_rt = NULL, hmax = 0.2,
                          summarize_by = c("weighted.mean", "median", "mean", "max"),
                          plot_it = FALSE, ask = plot_it,
                          clust = c("rt", "sp.rt"),
                          sigma_t = NULL, sigma_r = 0.5,
                          deepSplit = FALSE, verbose = FALSE, ...,
                          sigma.t = NULL, sigma.r = NULL, use.cor = NULL){
  if (...length()) stop("unused argument(s): ", paste(...names(), collapse = ", "))
  sigma_r <- resolve_deprecated(sigma.r, sigma_r)
  use_corrected_rt <- resolve_deprecated(use.cor, use_corrected_rt)
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
  summarize_by <- match.arg(summarize_by, choices = c("weighted.mean", 
                                                      "median", "mean", "max"))
  if (is.null(use_corrected_rt)){
    use_corrected_rt <- "rt.cor" %in% colnames(peak_list[[1]][[1]])
  }
  rt <- ifelse(use_corrected_rt, "rt.cor", "rt")
  start <- ifelse(use_corrected_rt, "start.cor", "start")
  end <- ifelse(use_corrected_rt, "end.cor", "end")
  check_duplicated_names(names(peak_list))
  windows <- attr(peak_list, "windows")
  targets <- attr(peak_list, "targets")
  if (clust == "sp.rt" && is.null(windows) && is.null(targets)){
    chrom_list <- get_chrom_list(peak_list, if (!missing(chrom_list)) chrom_list)
  }
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plot_it) {
    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }
  pks_found <- lapply(peak_list, function(x) lapply(x, function(y){
    if (nrow(y) > 0){
      y[!is.na(y[, rt]), , drop = FALSE]
    } else {
      y
    }}))
  if (clust == "sp.rt" && is.null(windows) && is.null(targets)){
    sigma_t_rt <- if (is.null(sigma_t)){
      all_pks <- do.call(rbind, unlist(pks_found, recursive = FALSE))
      0.5 * mean(all_pks$end - all_pks$start)
    } else sigma_t * rtf
  }
  clusterPeaks <- function(comp, pkLst){
    xx <- do.call(rbind, sapply(pkLst, function(samp) samp[comp]))
    file.idx <- xx$sample
    pkcenters <- xx[, rt]
    names(pkcenters) <- NULL
    tg <- if (!is.null(targets)){
      intersect(which(is.na(targets$lambda) |
                        targets$lambda == names(peak_list[[1]])[comp]),
                xx$target)
    }
    if (!is.null(tg)){
      if (length(tg) == 0) return(NULL)
      pkcenters.cl <- match(xx$target, tg)
    } else if (length(pkcenters) < 2){
      return(NULL)
    } else if (!is.null(windows)){
      dist_to_window <- outer(xx[["start"]], windows[, 1], function(a, b) abs(a - b)) +
                        outer(xx[["end"]], windows[, 2], function(a, b) abs(a - b))
      win_idx <- apply(dist_to_window, 1, which.min)
      pkcenters.cl <- match(win_idx, sort(unique(win_idx)))
    } else if (clust == 'rt'){
      pkcenters.hcl <- fastcluster::hclust(dist(pkcenters), method = "complete")
      pkcenters.cl <- cutree(pkcenters.hcl, h = hmax * rtf)
    } else if (clust == 'sp.rt'){
      ts <- as.numeric(rownames(chrom_list[[1]]))
      sp <- sapply(seq_along(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][
          get_retention_idx(pkcenters[i], ts, clamp = TRUE),]))
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
        xx[i[c(which.max(xx[i, response]), 1L)[1]], vars.idx]
      }))
    } else{
      cl.centers <- aggregate(xx[, vars.idx], by = list(pkcenters.cl),
                              FUN = summarize_by,
                              na.action = "na.pass")[, -1, drop = FALSE]
    }
    if (is.null(tg)){
      pkcenters.cl <- order(order(cl.centers[, rt]))[pkcenters.cl]
      cl.centers <- cl.centers[order(cl.centers[, rt]), , drop = FALSE]
    }
    ncl <- length(cl.centers[, rt])
    cl_names <- if (!is.null(tg)) target_ids(targets, names(peak_list[[1]]))[tg]
    metaInfo <- cbind(lambda = rep(suppressWarnings(
      as.numeric(names(peak_list[[1]])[comp])), ncl),
                      peak = 1:ncl, 
                      cl.centers
                      )
    if (use_corrected_rt) colnames(metaInfo) <- sub("\\.cor$", "", colnames(metaInfo))
    rownames(metaInfo) <- cl_names
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
    if (verbose && max(clusCount <- table(file.idx, pkcenters.cl)) > 1){
      warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    }
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, response]))
    allIs[which(allIs < 0)] <- 0
    failed <- !is.finite(allIs)
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(cl_names, names(pkLst)))
    Iinfo[cbind(pkcenters.cl, match(file.idx, names(pkLst)))[failed, , drop = FALSE]] <- NA
    for (i in which(!failed)){
      Iinfo[pkcenters.cl[i],  file.idx[i]] <-
      max(allIs[i], Iinfo[pkcenters.cl[i], file.idx[i]], na.rm = TRUE)
    }
    return(list(Iinfo, metaInfo))
  }
  result <- lapply(seq_len(ncomp), clusterPeaks, pks_found)
  if (all(vapply(result, is.null, logical(1)))){
    stop(if (!is.null(targets)) "None of the targets was found in any sample."
         else "Too few peaks to build a peak table.", call. = FALSE)
  }
  pt <- new_peak_table(result, peak_list,
                 args = list(peak_list = deparse(substitute(peak_list)),
                        chrom_list = attr(peak_list, "chrom_list"),
                        chrom_list_aligned = attr(peak_list, "chrom_list_aligned") %||% NA,
                        lambdas = list(names(peak_list[[1]])),
                        response = response,
                        use_corrected_rt = use_corrected_rt,
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
  ids <- target_ids(targets, names(peak_list[[1]]))
  if (!is.null(ids)){
    pt <- pt[, intersect(ids, colnames(pt$tab))]
  }
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
#'   in `pks` are ignored. A `matched_peak_list_set`, from targets at several
#'   wavelengths, gives one table with the columns for each wavelength, each
#'   defined by that wavelength's reference sample. Peaks absent from the reference sample will not
#'   appear as columns regardless of their prevalence in other samples. Choose a
#'   reference with a rich peak profile, or pass a table of `targets` to
#'   [`match_peaks`] to pick one that contains all of them; the
#'   targets' columns are then named after them, and the table can be
#'   restricted to them with `x[, targets$name]`. To
#'   discover peaks across all samples,
#'   call `get_peaktable` on the `peak_list` itself, which clusters peaks
#'   across every sample.
#' @method get_peaktable matched_peak_list
#' @export
get_peaktable.matched_peak_list <- function(peak_list, pks,
                                            response = c("area", "height"),
                                            rt_tol = 0.01,
                                            na.fill = 0,
                                            verbose = FALSE, ...) {
  if (...length()) stop("unused argument(s): ", paste(...names(), collapse = ", "))
  response <- match.arg(response, c("area", "height"))
  check_matched_pks(peak_list, pks)
  use_corrected_rt <- "rt.cor" %in% colnames(pks[[1]][[1]])
  parts <- matched_table_parts(peak_list, pks, response, rt_tol, na.fill, verbose,
                               use_corrected_rt)
  matched_peak_table(parts$result, pks, lambdas = parts$lambdas,
                     reference = parts$ref_idx, reference_sample = parts$ref_name,
                     response = response, rt_tol = rt_tol,
                     use_corrected_rt = use_corrected_rt,
                     labels = list(deparse(substitute(peak_list)), deparse(substitute(pks))))
}

#' @rdname get_peaktable
#' @method get_peaktable matched_peak_list_set
#' @export
get_peaktable.matched_peak_list_set <- function(peak_list, pks,
                                                response = c("area", "height"),
                                                rt_tol = 0.01,
                                                na.fill = 0,
                                                verbose = FALSE, ...) {
  if (...length()) stop("unused argument(s): ", paste(...names(), collapse = ", "))
  response <- match.arg(response, c("area", "height"))
  for (m in peak_list) check_matched_pks(m, pks)
  use_corrected_rt <- "rt.cor" %in% colnames(pks[[1]][[1]])
  parts <- lapply(peak_list, matched_table_parts, pks = pks, response = response,
                  rt_tol = rt_tol, na.fill = na.fill, verbose = verbose,
                  use_corrected_rt = use_corrected_rt)
  result <- unlist(lapply(parts, `[[`, "result"), recursive = FALSE)
  targets <- unlist(lapply(peak_list, function(m) attr(m, "reference_targets")$name))
  cols <- unlist(lapply(result, function(r) rownames(r[[2]])))
  unnamed <- !cols %in% targets
  cols[unnamed] <- paste0("V", seq_len(sum(unnamed)))
  cols <- split(cols, rep(seq_along(result), vapply(result, function(r) nrow(r[[2]]), integer(1))))
  result <- unname(Map(function(r, nm){ rownames(r[[1]]) <- nm; rownames(r[[2]]) <- nm; r },
                       result, cols))
  matched_peak_table(result, pks, lambdas = names(peak_list),
                     reference = vapply(parts, `[[`, integer(1), "ref_idx"),
                     reference_sample = vapply(parts, `[[`, character(1), "ref_name"),
                     response = response, rt_tol = rt_tol,
                     use_corrected_rt = use_corrected_rt,
                     labels = list(deparse(substitute(peak_list)), deparse(substitute(pks))))
}

#' Check that a `peak_list` can supply the values for a `matched_peak_list`
#' @noRd
check_matched_pks <- function(matched, pks){
  if (missing(pks))
    stop("Please supply the original peak_list from get_peaks() as 'pks'.")
  if (!inherits(pks, "peak_list"))
    stop("'pks' must be a 'peak_list' object.")
  if (!identical(names(matched), names(pks)))
    stop("Sample names of 'peak_list' and 'pks' must match.")
}

#' Peak table columns for one `matched_peak_list`, at its wavelength
#' @noRd
matched_table_parts <- function(matched, pks, response, rt_tol, na.fill, verbose,
                                use_corrected_rt = FALSE){
  ref_idx      <- reference_position(matched)
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

      n_assigned <- sum(ok)
      Imat[samp_i, col_idx[ok]] <- vals[ok]

      if (verbose) {
        n_drop <- n_ref - n_assigned
        message(sprintf("  %s [%s]: %d/%d reference peaks matched%s",
                        samp_name, lam, n_assigned, n_ref,
                        if (n_drop > 0) sprintf(" (%d dropped)", n_drop) else ""))
      }
    }

    present <- c(if (use_corrected_rt) c("rt.cor", "start.cor", "end.cor")
                 else c("rt", "start", "end"), PEAK_META_COLS)
    present <- present[present %in% colnames(ref_pks)]
    metaInfo <- cbind(
      lambda = rep(suppressWarnings(as.numeric(lam)), n_ref),
      peak   = seq_len(n_ref),
      ref_pks[, present, drop = FALSE]
    )
    colnames(metaInfo) <- sub("\\.cor$", "", colnames(metaInfo))
    col_names <- target_column_names(attr(matched, "reference_targets"),
                                     rt_ref_all, rt_tol)
    colnames(Imat) <- col_names
    rownames(metaInfo) <- col_names
    list(t(Imat), metaInfo)
  }

  list(result = Filter(Negate(is.null), lapply(lambdas, process_lambda)),
       ref_idx = ref_idx, ref_name = sample_names[ref_idx], lambdas = lambdas)
}

#' Assemble a `peak_table` from matched peak table columns
#' @noRd
matched_peak_table <- function(result, pks, lambdas, reference, reference_sample,
                               response, rt_tol, use_corrected_rt, labels){
  if (length(result) == 0)
    stop("No peaks found in reference sample across any wavelength.")
  pt <- new_peak_table(
    result, pks,
    args = list(
      peak_list        = labels[[1]],
      pks              = labels[[2]],
      chrom_list       = attr(pks, "chrom_list"),
      chrom_list_aligned = attr(pks, "chrom_list_aligned") %||% NA,
      lambdas          = list(lambdas),
      reference        = reference,
      reference_sample = reference_sample,
      response         = response,
      rt_tol           = rt_tol,
      use_corrected_rt = use_corrected_rt,
      reference_spectra = NA,
      metadata_path    = NA,
      normalized       = FALSE,
      normalization_by = NA,
      time_unit        = get_metadata_attribute(pks, "time_unit"),
      rt_unit          = get_time_unit(pks),
      intensity_unit   = get_metadata_attribute(pks, "intensity_unit")
    ))
  attr(pt, "warp_models") <- attr(pks, "warp_models")
  pt
}

#' Column names for a matched peak table, named after any targets
#' @noRd
target_column_names <- function(targets, rt_ref, rt_tol){
  out <- paste0("V", seq_along(rt_ref))
  if (is.null(targets) || all(is.na(targets$name))) return(out)
  k <- vapply(targets$rt, function(r){
    d <- abs(rt_ref - r)
    if (length(d) && min(d) <= rt_tol) which.min(d) else NA_integer_
  }, integer(1))
  if (anyNA(k))
    stop(sprintf(paste("No reference peak found for target%s %s. Build the",
                       "peak table from the same peak list passed to",
                       "`match_peaks`."),
                 plural(sum(is.na(k))),
                 paste(sQuote(targets$name[is.na(k)]), collapse = ", ")),
         call. = FALSE)
  named <- !is.na(targets$name)
  out[k[named]] <- targets$name[named]
  out
}

#' Peak table column names for a targets table
#'
#' `NULL` when unnamed targets are measured at more than one wavelength, since
#' a target then has one column per wavelength.
#' @noRd
target_ids <- function(targets, lambdas){
  if (is.null(targets)) return(NULL)
  if (!anyNA(targets$name)) targets$name
  else if (!anyNA(targets$lambda) || length(lambdas) == 1)
    paste0("V", seq_len(nrow(targets)))
}

#' Assemble a peak table
#'
#' Shared by the `peak_list` and `matched_peak_list` methods so the object
#' layout is defined once. Each method supplies its own `args`, which is
#' recorded verbatim.
#' @noRd
new_peak_table <- function(result, source, args){
  pull <- function(k) as.data.frame(t(do.call("rbind", lapply(result, `[[`, k))))
  added <- unique(unlist(lapply(source, attr, "added_metadata", exact = TRUE)))
  structure(
    list(tab = pull(1), pk_meta = pull(2),
         sample_meta = if (length(added))
           chromConverter::extract_metadata(source, what = added) else NA,
         instrument_meta = chromConverter::extract_metadata(source),
         ref_spectra = NA,
         args = args),
    class = c("peak_table", "list"),
    pk_args = attr(source, "meta"))
}

#' Reshape peaktable
#' 
#' Reshapes a peak table to long format, or to a wide data.frame with the
#' selected peaks alongside sample metadata.
#' @name reshape_peaktable
#' @importFrom stats reshape
#' @param x A `peak_table` object.
#' @param peaks A character vector specifying the peaks to include. If the
#' character vector is named, the names of the vector elements will be used in
#' place of the original peak names.
#' @param metadata A character vector specifying the metadata fields to include.
#' @param fixed_levels Logical. Whether to fix factor levels of features in the
#' order provided. Defaults to `TRUE`. Ignored when `format = "wide"`.
#' @param format Either `"long"` (the default), with one row per sample and
#' peak, or `"wide"`, with one row per sample and one column per peak.
#' @return A data.frame containing the information for the specified peaks. In
#' long format it has columns `sample`, `peak`, `lambda`, `rt`, the response and
#' any metadata; in wide format, `sample`, any metadata, then one column per
#' peak.
#' @author Ethan Bass
#' @family utility functions
#' @examples
#' data(pk_tab)
#' reshape_peaktable(pk_tab, peaks = c("V8", "V13"))
#' reshape_peaktable(pk_tab, peaks = c(a = "V8", b = "V13"), format = "wide")
#' @export
reshape_peaktable <- function(x, peaks, metadata, fixed_levels = TRUE,
                              format = c("long", "wide")){
  format <- match.arg(format)
  if (!missing(peaks)){
    if (is.numeric(peaks)){
      peaks <- colnames(x$tab)[peaks]
    }
    unknown <- setdiff(peaks, colnames(x$tab))
    if (length(unknown))
      stop(sprintf("Peak%s not found in the peak table: %s.", plural(length(unknown)),
                   paste(sQuote(unknown, FALSE), collapse = ", ")), call. = FALSE)
    df <- x$tab[, match(peaks, colnames(x$tab)), drop = FALSE]
    if (!is.null(names(peaks))){
      x$pk_meta <- x$pk_meta[, match(peaks, colnames(x$pk_meta)), drop = FALSE]
      colnames(df) <- colnames(x$pk_meta) <- names(peaks)
      peaks <- colnames(df)
    }
  } else {
    df <- x$tab
  }
  if (!missing(metadata)){
    if (!is_attached(x$sample_meta))
      stop("No sample metadata is attached to `x`, so `metadata` cannot be ",
           "selected. Attach it first with `attach_metadata()`, or omit `metadata`.",
           call. = FALSE)
    absent <- setdiff(metadata, colnames(x$sample_meta))
    if (length(absent))
      stop(sprintf("Metadata column%s not found: %s. Available: %s.",
                   plural(length(absent)), paste(sQuote(absent, FALSE), collapse = ", "),
                   paste(sQuote(colnames(x$sample_meta), FALSE), collapse = ", ")),
           call. = FALSE)
    meta_idx <- which(colnames(x$sample_meta) %in% metadata)
    x$sample_meta <- x$sample_meta[, meta_idx, drop = FALSE]
  }
  if (format == "wide"){
    meta <- if (!is.null(dim(x$sample_meta))) x$sample_meta
            else data.frame(row.names = seq_len(nrow(df)))
    return(data.frame(sample = rownames(df), meta, df, row.names = NULL,
                      check.names = FALSE))
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
