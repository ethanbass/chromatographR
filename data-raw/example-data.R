# Regenerates `Sa_pr`, `Sa_warp` and `pk_tab` from `Sa`.

devtools::load_all()
data(Sa)

Sa_pr <- preprocess(Sa, dim1 = seq(10, 18.66, by = 0.02),
                    dim2 = seq(200, 318, by = 2), stiffness = 2e5,
                    zero_floor = FALSE, cl = 1, show_progress = FALSE)

Sa_warp <- correct_rt(Sa_pr, alg = "ptw", lambdas = 210,
                      what = "corrected.values")

pks_egh <- get_peaks(Sa_warp, lambdas = c(210, 260), fit = "egh",
                     smooth_window = 2, cl = 1, show_progress = FALSE)
pk_tab <- get_peaktable(pks_egh, response = "area", hmax = 0.2)

save(Sa_pr, file = "data/Sa_pr.RData", compress = "xz")
save(Sa_warp, file = "data/Sa_warp.RData", compress = "xz")
save(pk_tab, file = "data/pk_tab.RData", compress = "xz")
