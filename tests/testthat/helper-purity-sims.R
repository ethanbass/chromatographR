# Simulated DAD chromatograms for the purity calibration tests: a bilinear
# model, M[t, wl] = profile(t) %o% spectrum(wl), plus white Gaussian noise.
purity_sims <- function(){
  withr::local_seed(42)
  noise_sd <- 5
  wl <- seq(200, 400, by = 5)
  time <- seq(0, 10, length.out = 1000L)

  gaussian_spec <- function(lmax, w) exp(-(wl - lmax)^2 / (2 * w^2))
  gaussian_profile <- function(rt, h, s) h * exp(-(time - rt)^2 / (2 * s^2))
  noise_mat <- function() matrix(rnorm(length(time) * length(wl), sd = noise_sd),
                                 length(time), length(wl))
  label <- function(m){
    dimnames(m) <- list(round(time, 3), wl)
    m
  }

  spec_A <- gaussian_spec(254, 20)
  spec_B <- gaussian_spec(340, 15)
  spec_bg <- gaussian_spec(210, 30)
  apex_t <- which.min(abs(time - 1.0))
  pos <- c(apex_t, apex_t - 25L, apex_t + 25L)

  prof_A <- gaussian_profile(1.0, 500, 0.1)
  mat_pure <- label(outer(prof_A, spec_A) + noise_mat())
  # height 50 keeps the threshold below `max_threshold`
  mat_lo <- label(outer(gaussian_profile(1.0, 50, 0.1), spec_A) + noise_mat())
  prof_B <- gaussian_profile(1.3, 100, 0.1)
  mat_mixed <- label(outer(prof_A, spec_A) + outer(prof_B, spec_B) + noise_mat())
  profiles_w <- cbind(prof_A, prof_B)[pos[2]:pos[3], , drop = FALSE]
  mat_bg <- label(outer(prof_A, spec_A) + outer(50 * time, spec_bg) + noise_mat())

  list(
    pure = list(mat = mat_pure, pos = pos, noise_sd = noise_sd),
    lo_snr = list(mat = mat_lo, pos = pos, noise_sd = noise_sd),
    mixed = list(mat = mat_mixed, pos = pos,
                 impurity_frac = profiles_w[, 2] / rowSums(profiles_w)),
    bg = list(mat = mat_bg, pos = pos, noise_sd = noise_sd)
  )
}
