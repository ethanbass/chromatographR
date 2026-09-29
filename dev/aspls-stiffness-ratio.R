# Backs the sentence in ?baseline_estimators on how the "aspls" stiffness
# relates to the "als" one. For each Sa trace and each "als" stiffness, finds
# the "aspls" stiffness whose baseline is closest (RMSE) to the "als" baseline.
# Run from the package root.
devtools::load_all(quiet = TRUE)
data(Sa)

traces <- list()
for (s in seq_along(Sa)){
  wl <- as.numeric(colnames(Sa[[s]]))
  for (l in c(210, 260, 318))
    traces[[paste(names(Sa)[s], l)]] <- Sa[[s]][, which.min(abs(wl - l))]
}
grid <- 10^seq(5, 12, by = 0.25)

res <- do.call(rbind, lapply(names(traces), function(nm){
  y <- traces[[nm]]
  do.call(rbind, lapply(c(1e5, 1e6, 1e7, 1e8), function(s){
    b_als <- estimate_baseline(y, baseline = "als", stiffness = s)
    rmse <- vapply(grid, function(g){
      sqrt(mean((estimate_baseline(y, baseline = "aspls", stiffness = g) - b_als)^2))
    }, numeric(1))
    k <- which.min(rmse)
    data.frame(trace = nm, als = s, ratio = grid[k] / s,
               rel_rmse = rmse[k] / sd(y))
  }))
}))

aggregate(ratio ~ als, res, function(v) signif(c(median = median(v),
                                                 min = min(v), max = max(v)), 2))
summary(res$rel_rmse)
