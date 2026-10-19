# Which default should `find_peaks(smooth_window)` have, and how wide should
# the window be relative to the peaks?
#
# Synthetic Gaussian peaks with a true FWHM of 3-40 timepoints, in runs of
# 500-20000 points, under white, 2-scan-smoothed and AR(0.9) noise of unit SD.
# "mixed" traces hold isolated peaks at 3-500 sigma and close pairs 0.75-2
# FWHM apart at 1:1 and 3:1 heights; "crowd" traces have peaks every 2 FWHM.
# Candidates: the old proportional default (0.001), fixed windows of 1-8
# timepoints, and k x the trace's median peak FWHM (`median_peak_fwhm`), each
# at min_prominence_snr = 0, 3 and 4. A true peak is found if a detected apex
# lies within FWHM/2 of it; pair results count only pairs resolvable without
# noise.
#
# Result (2026-10-10): the best window is about a third of the FWHM. One FWHM
# merges close pairs; under a quarter lets the prominence filter drop real
# peaks. No fixed value suits every FWHM; 3 timepoints is best for FWHM 3-10.
#
# Run from the package root:
#   Rscript dev/smooth-window-default.R [n_cores] [noise types]
# e.g. `Rscript dev/smooth-window-default.R 8 white,corr` (about 25 min on 9
# cores) or `... 8 ar0.9`.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
OUT <- tempdir()
args <- commandArgs(trailingOnly = TRUE)
N_CORES <- if (length(args)) as.integer(args[1]) else 8L
NOISES <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("white", "corr")
RES_FILE <- "smooth-window-default.rds"

FWHMS <- c(3, 5, 10, 20, 40)
NS <- c(500, 3000, 20000)
SEPS <- c(0.75, 1, 1.25, 1.5, 2)
SNRS <- c(0, 3, 4)
CANDS <- c("p0.001", "w1", "w2", "w3", "w5", "w8", "ad0.25", "ad0.33", "ad0.5", "ad1")

rlogunif <- function(n, a, b) exp(runif(n, log(a), log(b)))

make_truth <- function(f, N, layout) {
  rows <- list(); x <- 4 * f; k <- 0
  if (layout == "mixed") {
    repeat {
      if (runif(1) < 0.4) {
        if (x > N - 4 * f) break
        k <- k + 1
        rows[[length(rows) + 1]] <- data.frame(c = x, h = rlogunif(1, 3, 500), type = "iso",
                                               feat = k, sep = NA, ratio = NA)
        x <- x + 6 * f
      } else {
        s <- sample(SEPS, 1); r <- sample(c(1, 3), 1)
        if (x + s * f > N - 4 * f) break
        k <- k + 1
        hm <- rlogunif(1, 15, 500)
        h <- if (runif(1) < 0.5) c(hm, hm / r) else c(hm / r, hm)
        rows[[length(rows) + 1]] <- data.frame(c = c(x, x + s * f), h = h, type = "pair",
                                               feat = k, sep = s, ratio = r)
        x <- x + s * f + 6 * f
      }
    }
  } else {
    repeat {
      if (x > N - 4 * f) break
      k <- k + 1
      rows[[length(rows) + 1]] <- data.frame(c = x, h = rlogunif(1, 5, 300), type = "crowd",
                                             feat = k, sep = NA, ratio = NA)
      x <- x + f * (0.6 + rexp(1, 1 / 1.4))
    }
  }
  if (!length(rows)) return(NULL)
  tr <- do.call(rbind, rows)
  tr$c <- tr$c + runif(nrow(tr), -0.5, 0.5)
  tr
}

render <- function(tr, f, N) {
  s <- f / (2 * sqrt(2 * log(2))); t <- seq_len(N); y <- numeric(N)
  for (i in seq_len(nrow(tr))) y <- y + tr$h[i] * exp(-(t - tr$c[i])^2 / (2 * s^2))
  y
}

noise <- function(N, type) {
  if (type == "white") return(rnorm(N))
  if (type == "ar0.9") return(as.numeric(stats::filter(rnorm(N + 200) * sqrt(1 - 0.81), 0.9, "recursive"))[-(1:200)])
  z <- rnorm(N + 1); (z[-1] + z[-(N + 1)]) / sqrt(2)
}

adaptive_fwhm <- function(y) {
  p0 <- find_peaks(y, smooth_window = 1)
  nl <- trace_noise_level(y, seq_along(y), data.frame(start = p0$lower, end = p0$upper),
                          apex = y[p0$pos])
  median_peak_fwhm(y, nl)
}

match_1to1 <- function(det, ctr, tol) {
  found <- rep(NA_integer_, length(ctr))
  if (!length(det)) return(found)
  d <- abs(outer(ctr, det, "-"))
  cand <- which(d <= tol, arr.ind = TRUE)
  if (!nrow(cand)) return(found)
  cand <- cand[order(d[cand]), , drop = FALSE]
  used <- logical(length(det))
  for (r in seq_len(nrow(cand))) {
    i <- cand[r, 1]; j <- cand[r, 2]
    if (is.na(found[i]) && !used[j]) { found[i] <- j; used[j] <- TRUE }
  }
  found
}

run_cell <- function(cell) {
  set.seed(cell$seed)
  f <- cell$f; N <- cell$N
  tr <- make_truth(f, N, cell$layout)
  if (is.null(tr)) return(NULL)
  clean <- render(tr, f, N)
  y <- clean + noise(N, cell$noise)
  tol <- f / 2
  cp <- find_peaks(clean, smooth_type = "none", bounds = FALSE)
  tr$clean_found <- !is.na(match_1to1(cp, tr$c, tol))
  t_ad <- system.time(fw <- adaptive_fwhm(y))[["elapsed"]]
  wins <- c("p0.001" = 0.001, w1 = 1, w2 = 2, w3 = 3, w5 = 5, w8 = 8,
            "ad0.25" = if (is.finite(fw)) max(1, 0.25 * fw) else 1,
            "ad0.33" = if (is.finite(fw)) max(1, fw / 3) else 1,
            "ad0.5" = if (is.finite(fw)) max(1, 0.5 * fw) else 1,
            ad1 = if (is.finite(fw)) max(1, fw) else 1)
  peaks <- list(); traces <- list()
  for (cn in names(wins)) for (snr in SNRS) {
    p <- suppressWarnings(find_peaks(y, smooth_window = wins[[cn]], min_prominence_snr = snr))
    m <- match_1to1(p$pos, tr$c, tol)
    found <- !is.na(m)
    merged <- rep(NA, nrow(tr))
    for (k in unique(tr$feat[tr$type == "pair"])) {
      ii <- which(tr$feat == k)
      lo <- min(tr$c[ii]) - tol; hi <- max(tr$c[ii]) + tol
      nd <- sum(p$pos >= lo & p$pos <= hi)
      merged[ii] <- nd == 1
    }
    peaks[[length(peaks) + 1]] <- data.frame(tr, found = found, merged = merged,
                                             cand = cn, snr = snr)
    traces[[length(traces) + 1]] <- data.frame(cand = cn, snr = snr, window = wins[[cn]],
                                               n_det = nrow(p), spurious = sum(!seq_len(nrow(p)) %in% m),
                                               n_true = nrow(tr))
  }
  meta <- cell[c("f", "N", "noise", "layout", "seed")]
  list(peaks = cbind(meta, do.call(rbind, peaks), row.names = NULL),
       traces = cbind(meta, fw = fw, t_ad = t_ad, do.call(rbind, traces), row.names = NULL))
}

grid <- expand.grid(f = FWHMS, N = NS, noise = NOISES,
                    layout = c("mixed", "crowd"), stringsAsFactors = FALSE)
grid <- grid[grid$N >= 12 * grid$f, ]
cells <- list(); sd0 <- 0
for (g in seq_len(nrow(grid))) {
  f <- grid$f[g]; N <- grid$N[g]
  nfeat <- max(1, (N - 8 * f) / (if (grid$layout[g] == "mixed") 6.8 * f else 2 * f))
  ns <- min(300, max(10, ceiling(150 / (0.6 * nfeat))))
  for (s in seq_len(ns)) { sd0 <- sd0 + 1; cells[[sd0]] <- c(as.list(grid[g, ]), seed = sd0) }
}
cat(length(cells), "traces\n")
t0 <- Sys.time()
res <- parallel::mclapply(cells, run_cell, mc.cores = N_CORES, mc.preschedule = FALSE)
err <- vapply(res, inherits, logical(1), "try-error")
if (any(err)) { print(res[err][1]); res <- res[!err] }
res <- Filter(Negate(is.null), res)
peaks <- do.call(rbind, lapply(res, `[[`, "peaks"))
traces <- do.call(rbind, lapply(res, `[[`, "traces"))
saveRDS(list(peaks = peaks, traces = traces), file.path(OUT, RES_FILE))
print(Sys.time() - t0)

# Summary: one row per FWHM x candidate, one block per min_prominence_snr;
# the 0.001 candidate is split by run length.
files <- RES_FILE
options(width = 220)
pct <- function(x) if (length(x)) round(100 * mean(x)) else NA
for (fl in files) {
  r <- readRDS(file.path(OUT, fl)); pk <- r$peaks; tr <- r$traces
  relab <- function(d) { i <- d$cand == "p0.001"; d$cand[i] <- paste0("p0.001@", d$N[i]); d }
  pk <- relab(pk); tr <- relab(tr)
  lv <- c(paste0("p0.001@", c(500, 3000, 20000)), "w1", "w2", "w3", "w5", "w8", "ad0.25", "ad0.33", "ad0.5", "ad1")
  pr <- aggregate(cbind(both = found, clean = clean_found) ~ seed + feat + f + N + sep + ratio + cand + snr,
                  pk[pk$type == "pair", ], all)
  pr <- pr[pr$clean, ]
  cat("\n\n================", fl, "  (noise:", paste(unique(tr$noise), collapse = ","), ")\n")
  cat("iso3-10/iso10+: isolated-peak recall % by height (sigma); spur/1k: spurious apices per 1000 points (mixed layout);\n",
      "1:1@1.0 etc: % of noiselessly-resolvable pairs with both members found; crowd: recall in crowded layout / noiseless ceiling %\n")
  for (snr in c(0, 3, 4)) {
    rows <- list()
    for (f in sort(unique(pk$f))) for (cn in lv) {
      p <- pk[pk$snr == snr & pk$f == f & pk$cand == cn, ]
      if (!nrow(p)) next
      iso <- p[p$type == "iso", ]; cr <- p[p$type == "crowd", ]
      t <- tr[tr$snr == snr & tr$f == f & tr$cand == cn & tr$layout == "mixed", ]
      q <- pr[pr$snr == snr & pr$f == f & pr$cand == cn, ]
      rows[[length(rows) + 1]] <- data.frame(f = f, cand = cn,
        win = round(median(t$window[t$window >= 1], na.rm = TRUE), 1),
        iso3_10 = pct(iso$found[iso$h < 10]), iso10p = pct(iso$found[iso$h >= 10]),
        spur1k = round(1000 * sum(t$spurious) / sum(t$N), 1),
        `1:1@1.0` = pct(q$both[q$ratio == 1 & q$sep == 1]),
        `1:1@1.25` = pct(q$both[q$ratio == 1 & q$sep == 1.25]),
        `3:1@1.25` = pct(q$both[q$ratio == 3 & q$sep == 1.25]),
        `3:1@1.5` = pct(q$both[q$ratio == 3 & q$sep == 1.5]),
        `3:1@2` = pct(q$both[q$ratio == 3 & q$sep == 2]),
        crowd = pct(cr$found), crowd_ceil = pct(cr$clean_found), check.names = FALSE)
    }
    cat("\n---- min_prominence_snr =", snr, "\n")
    print(do.call(rbind, rows), row.names = FALSE)
  }
  if (!is.null(tr$t_ad)) {
    u <- unique(tr[, c("seed", "N", "t_ad")])
    cat("\nadaptive-width step, seconds per trace (median, max) by N:\n")
    print(aggregate(t_ad ~ N, u, function(z) round(c(median(z), max(z)), 3)))
  }
}
