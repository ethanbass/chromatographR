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
  [`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md).

- what:

  What type of plot to return. Either `"traces"` (default) or
  `"heatmap"`.

- engine:

  What plotting engine to use. One of `"base"` (default), `"ggplot"`, or
  `"plotly"`.

- lambdas:

  Wavelengths to plot, from among those used for warping. Defaults to
  all of them. The heatmap shows only the first.

- show_legend:

  Logical. Whether to include sample legend.

- ...:

  Ignored.

## Value

A `plotly` or `ggplot` object when `engine = "plotly"` or
`engine = "ggplot"`, respectively. No value is returned when
`engine = "base"`.

## Details

Plots raw and warped chromatograms, one above the other, either as
individual traces or as a heatmap, according to the value of `what`. The
plot can be produced using either base R graphics, `ggplot2`, or
`plotly`, according to the value of the `engine` argument.

## Side effects

If `engine == "base"`, plots are rendered to the active graphics device.

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, what = "models", lambdas = 210)
plot(warp)
```
