## Code to reproduce the chromatographR hex sticker (`man/figures/logo.png`).
##
## Run from the package root with `source("dev/logo.R")`. hexSticker and
## ggplot2 are only needed here, so neither is a package dependency.
##
## The sticker shows six wavelengths of the first preprocessed goldenrod
## chromatogram over a 13-15 minute window, drawn on a gold hexagon.

library(chromatographR)
library(ggplot2)
library(hexSticker)

data(Sa_pr)

lambdas <- seq(200, 300, by = 20)
rt_window <- c(12.99, 15.09)

chrom <- reshape_chroms(Sa_pr, idx = 1, lambdas = lambdas)
chrom <- chrom[chrom$rt >= rt_window[1] & chrom$rt <= rt_window[2], ]

p <- ggplot(chrom, aes(x = rt, y = absorbance, group = lambda)) +
  geom_line(aes(color = factor(lambda)), show.legend = FALSE) +
  # a fixed ceiling (rather than the data range) leaves headroom above the
  # tallest peak, so the trace clears the package name below it
  ylim(0, 600) +
  theme_void() +
  theme_transparent()

sticker(p, package = "chromatographR", p_size = 14, p_y = 0.65,
        p_color = "black", s_x = 1, s_y = 1.35, s_width = 1.6, s_height = 1.2,
        h_color = "black", h_fill = "gold",
        filename = "man/figures/logo.png")
