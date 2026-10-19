devtools::load_all(quiet = TRUE)
set.seed(42)

out_dir <- file.path(tempdir(), "sloped-cap-sim")
if (nzchar(Sys.getenv("CAP_SIM_OUT"))) out_dir <- Sys.getenv("CAP_SIM_OUT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tt <- seq(0, 10, by = 0.005)
emg <- function(x, mu, s, lam) lam * exp(lam * (mu - x) + lam^2 * s^2 / 2) * pnorm((x - mu - lam * s^2) / s)

N <- 200
noise <- as.numeric(Sys.getenv("CAP_SIM_NOISE", "0.01"))
sims <- lapply(seq_len(N), function(i) {
  A <- exp(runif(1, log(0.02), log(2)))
  s <- runif(1, 0.015, 0.05)
  lam <- 1 / runif(1, 0.005, 0.08)
  pk <- A * emg(tt, 5, s, lam)
  has_nb <- runif(1) < 0.5
  nb <- if (has_nb) {
    d <- sample(c(-1, 1), 1) * runif(1, 0.12, 0.35)
    runif(1, 0.2, 1) * A * emg(tt, 5 + d, runif(1, 0.015, 0.05), 1 / runif(1, 0.005, 0.06))
  } else 0
  base <- runif(1, 0, 0.5) + runif(1, -0.05, 0.05) * (tt - 5) + runif(1, -0.01, 0.01) * (tt - 5)^2
  y <- pk + nb + base + rnorm(length(tt), 0, noise)
  list(y = matrix(y, ncol = 1, dimnames = list(tt, "210")), area = A, nb = has_nb)
})
chroms <- setNames(lapply(sims, `[[`, "y"), paste0("s", seq_len(N)))
truth <- vapply(sims, `[[`, numeric(1), "area")
has_nb <- vapply(sims, `[[`, logical(1), "nb")

q <- function(expr) suppressWarnings(suppressMessages(expr))
unit_peak <- matrix(emg(tt, 5, 0.03, 1 / 0.001), ncol = 1, dimnames = list(tt, "210"))
cal <- q(get_peaks(list(a = unit_peak), lambdas = 210, fit = "raw", show_progress = FALSE))[[1]][[1]]
unit <- cal$area[which.min(abs(cal$rt - 5))]
cat("raw area of a unit-area peak:", round(unit, 4), "\n")

orig <- fit_starts
edge_min <- function(...){
  s <- orig(...)
  y <- list(...)[[2]]
  k <- max(3L, ceiling(FLOOR_EDGE_FRAC * length(y)))
  s$upper_edges <- c(floor_s = max(min(head(y, k)), 0),
                     floor_e = max(min(tail(y, k)), 0))
  s
}
uncapped <- function(...){
  s <- orig(...)
  s$upper_edges[] <- Inf
  s
}
variants <- list(uncapped = uncapped, cap_min = edge_min, cap_median = orig)

res <- list()
fits <- strsplit(Sys.getenv("CAP_SIM_FITS", "egh,bemg,gaussian"), ",")[[1]]
for (f in fits) for (v in names(variants)) {
  assignInNamespace("fit_starts", variants[[v]], "chromatographR")
  p <- q(get_peaks(chroms, lambdas = 210, fit = f, baseline = "sloped",
                   show_progress = FALSE, cl = 4))
  a <- vapply(p, function(s){
    x <- s[[1]]
    j <- which(abs(x$rt - 5) < 0.05)
    if (length(j)) x$area[j[which.min(abs(x$rt[j] - 5))]] / unit else NA_real_
  }, numeric(1))
  e <- a / truth - 1
  for (grp in c("all", "isolated", "with neighbour")) {
    k <- switch(grp, all = rep(TRUE, N), isolated = !has_nb, "with neighbour" = has_nb)
    res[[length(res) + 1]] <- data.frame(
      fit = f, variant = v, group = grp, found = sum(!is.na(e[k])),
      bias = median(e[k], na.rm = TRUE),
      med_abs = median(abs(e[k]), na.rm = TRUE),
      p90_abs = quantile(abs(e[k]), 0.9, na.rm = TRUE))
  }
  saveRDS(do.call(rbind, res), file.path(out_dir, paste0("results-", paste(fits, collapse = "-"), "-noise", noise, ".rds")))
  cat("done:", f, v, "\n")
}
assignInNamespace("fit_starts", orig, "chromatographR")

r <- do.call(rbind, res)
r[, 5:7] <- round(r[, 5:7], 3)
options(width = 150)
print(r, row.names = FALSE)
