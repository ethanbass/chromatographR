# Get time resolution

Get the sampling interval of a chromatogram: the median gap between
consecutive retention times.

## Usage

``` r
get_time_resolution(x, idx = 1, digits = 6, unit = NULL)
```

## Arguments

- x:

  A chromatogram, a list of chromatograms, or a `peak_table` object.

- idx:

  Index of the chromatogram from which to derive the interval, if `x` is
  a list or `peak_table`. Defaults to `1`.

- digits:

  Number of significant digits to round to. Defaults to `6`, because a
  regular time axis carries floating-point noise in the last digits of
  `diff`.

- unit:

  Unit for the returned interval: `"min"`, `"s"` or `"ms"`. Defaults to
  `NULL`, giving the interval in the units of the chromatogram's own
  time axis, which is what
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)'s
  `dim1` expects. Supply `"min"` to pair with the `time_resolution`
  arguments of
  [`plot_chroms`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md)
  and
  [`reshape_chroms`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md),
  which are in minutes, regardless of the units of `x`. The units of `x`
  are read from its `time_unit` attribute, falling back to minutes when
  it is absent, as
  [`plot_chroms`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md)
  does.

## Value

A length-one numeric giving the time between consecutive timepoints, in
`unit` if supplied and otherwise in the units of the chromatogram's time
axis.

## See also

Other accessors:
[`get_baseline()`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md)

## Examples

``` r
data(Sa_warp)
get_time_resolution(Sa_warp)
#> [1] 0.02
```
