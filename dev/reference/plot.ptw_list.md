# Plot PTW alignments

Plots [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) alignments.

## Usage

``` r
# S3 method for class 'ptw_list'
plot(
  x,
  what = c("traces", "heatmap"),
  engine = c("base", "ggplot", "plotly"),
  lambdas,
  show_legend = TRUE,
  ...
)
```

## Arguments

- x:

  A `ptw_list` object created by
  [`correct_rt`](https://ethanbass.github.io/chromatographR/dev/reference/correct_rt.md).

- what:

  What type of plot to return. Either `traces` or `heatmap`.

- engine:

  What plotting engine to use. Either `base`, `ggplot` or `plotly`

- lambdas:

  Which lambdas to plot.

- show_legend:

  Logical. Whether to include sample legend.

- ...:

  Additional arguments (placeholder).

## Value

No return value, called for side effects.

## Side effects

Plots PTW alignments at the specified wavelength (`lambda`) either as
individual traces or as a heatmap, according to the value of `what`. The
plot can be produced using either base R graphics, `ggplot2`, or
`plotly`, according to the value of `engine`.

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, what = "models", lambdas = 210)
plot(warp)
}
```
